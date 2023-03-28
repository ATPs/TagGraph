#script called: 
# SetUpTagGraphProject.py (scripts)
#	database/Models.py
#	preprocessors/ParsePEAKS7Results.py if not from LADS
# HUMAN_PROTEOME_TAGGRAPH.py (scripts) : now called: scripts/ProteomeTagGraph.py
# Importer.py (database)
# AddPlausibleModAnnotationsDB.py (database)
# CheckProteomeRun.py (database)
#
#slin 20160921  modify to the take value from a init File:  init.config under project directory
#slin 201707xx  added lib path and others.
#slin 201807    rename script from: RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py to: RunTagGraphProteomeEasy.py
#slin 20170720  rename database/CHECK_HUMAN_PROTEOME_RUN.py to database/CheckProteomeRun.py
#slin 20170720  rename scripts/HUMAN_PROTEOME_TAGGRAPH.py to scripts/ProteomeTagGraph.py
#slin 20170816  LADS
#slin 20181004  update print function
#slin 20190115  new version by pass LADs, peaks Denovo.

import os
import sys
import time
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)
DATABASE_SCRIPT_DIR = os.path.join(PAR_DIR, 'database')
import ArgLib
import DataFile
import glob

if __name__ == '__main__':
    print('dtadir is the directory containing the mzXML files to analyze')
    print('peaks is a dictionary mapping {experiment_name: peaks csv}')
    print('output is the directory to move all files to and set up the project in')
    options = ArgLib.parse(['init', 'dtadir', 'peaks', 'output', 'ppmstd', 'modtolerance', 'unimoddict', 'maxcounts', 'modmaxcounts', 'fmindex', 'model', 'config','experimentname','uselads','ladsresults','simple'], optArgs=[{'opts': ('-x', '--splittaxon'), 'attrs': {'dest': 'splittaxon', 'action': 'store_true', 'default': False, 'help': 'Flag. For searches of metaproteomic databases, split identical context entries by taxon for accurate consideration via EM.'}}])
    t0 = time.time()
    args = DataFile.getArgs(options, ['init', 'dtadir', 'peaks', 'output','experimentname','uselads','ladsresults','simple'])
    DataFile.executeProcess(os.path.dirname(__file__), 'SetUpTagGraphProject.py', args)
    base_args = DataFile.getArgs(options, ['init', 'ppmstd', 'modtolerance', 'model', 'config', 'fmindex', 'unimoddict', 'maxcounts', 'modmaxcounts'])
    import_args = DataFile.getArgs(options, ['init', 'fmindex', 'modtolerance', 'ppmstd', 'maxcounts', 'modmaxcounts']) + ['--type', 'experiment']
    experiment_names = eval(options.peaks).keys()
    if options.uselads:
        experiment_names=[options.experimentname]
    for experiment_name in experiment_names:
        t1 = time.time()
        experiment_dir = os.path.join(options.output, experiment_name)
        sqlitedb_loc = os.path.join(options.output, 'results.db')
        #Get all the fraction numbers
        peaks_dir = os.path.join(experiment_dir, 'de_novo')
        peaks_files = glob.glob(peaks_dir+'/'+'*_parsed_F*.tdv')
        fracs=[]
        for file in peaks_files:
            fracs += [int(file.split('_')[-1][1:-4])]
        fracs = sorted(fracs)
        #Run TAG-GRAPH over all fracs
        t2 = time.time()
        print('Running TAGGRAPH for experiment %s'%experiment_name)
        for frac in fracs:
            args = base_args + ['--fraction', str(frac), '--output', 'TAGGRAPH', '--dtadir', experiment_dir]
            args += ['--splittaxon'] if options.splittaxon else []
            DataFile.executeProcess(os.path.dirname(__file__),'ProteomeTagGraph.py', args)
        print('Total runtime for TAG-GRAPH step in experiment %s: %f'%(experiment_name, time.time() - t2))
        # Import results into database
        t2 = time.time()
        print('Importing Results into DB for experiment %s'%experiment_name)
        args = import_args + ['--sqlitedb', sqlitedb_loc, '--taggraph', os.path.join(experiment_dir, 'taggraph'), '--experimentname', experiment_name, '--fracs', ','.join(['F%i'%frac for frac in fracs])]
        DataFile.executeProcess(DATABASE_SCRIPT_DIR, 'Importer.py', args)
        print('Total runtime for import step in experiment %s: %f'%(experiment_name, time.time() - t2))
        t2 = time.time()
        print('Adding plausible mod candidates for experiment %s'%experiment_name)
        args = DataFile.getArgs(options, ['init', 'model', 'config', 'ppmstd', 'modtolerance']) + ['--dtadir', os.path.join(experiment_dir, 'data'),  '--experimentname', experiment_name, '--fracs', 'all', '--sqlitedb', sqlitedb_loc, '--output', options.output+'/'+'%s_addPlausibleMods'%experiment_name]
        DataFile.executeProcess(os.path.dirname(__file__), 'AddPlausibleModAnnotationsDB.py', args)
        print('Total runtime for plausible mod addition step in experiment %s: %f'%(experiment_name, time.time() - t2))
        # Check to make sure all runs finished properly
        print('Checking to make sure run finished properly')
        args = ['--dtadir', options.output, '--experimentname', experiment_name, '--output', options.output+'/'+'%s_CHECK.txt'%experiment_name]
        DataFile.executeProcess(DATABASE_SCRIPT_DIR, "CheckProteomeRun.py", args)
        # TODO: Clean up (remove) taggraph files + any other files I want to delete?
        print('Total runtime for experiment %s: %f'%(experiment_name, time.time() - t1))
    print('Total runtime for project: ', time.time() - t0)
