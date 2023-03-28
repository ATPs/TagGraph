#slin 201707    added lib path and others.
#slin 201807    rename script from: SetUpTAGGRAPHProject.py to: SetUpTagGraphProject.py
#slin 20181004  update print function
#slin 20190116  add simple parse the input csv to tg required format

import os
import sys
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1,PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1,LIB_DIR)
DATABASE_SCRIPT_DIR = os.path.join(PAR_DIR,'database')
PREPROCESSOR_SCRIPT_DIR = os.path.join(PAR_DIR,'preprocessors')
RESOURCES_DIR = os.path.join(PAR_DIR,'resources')
import PepInput
import shutil
import ArgLib
import DataFile
import glob
import pickle
import csv

if __name__ == '__main__':
    print('dtadir is the directory containing the mzXML files to analyze')
    print('peaks is a dictionary mapping {experiment_name: peaks csv}')
    print('output is the directory to move all files to and set up the project in')
    options = ArgLib.parse(['init', 'dtadir', 'peaks', 'output','experimentname','uselads','ladsresults','simple'])
    #Fails with an OSError if directory already exists
    os.makedirs(options.output)
    args = ['--sqlite', os.path.join(options.output, 'results.db')]
    DataFile.executeProcess(DATABASE_SCRIPT_DIR, 'Models.py', args)
    # Make experiment directories
    # Structure
    # /options.output
    # .../ExperimentName
    # ...../data
    # ...../de_novo
    # ...../taggraph
    #Make the directory and subdirectories
    experiment=options.experimentname
    experiment_dir=os.path.join(options.output,experiment)
    os.makedirs(experiment_dir)
    peaks_dir=os.path.join(experiment_dir,'de_novo')
    os.makedirs(peaks_dir)
    data_dir=os.path.join(experiment_dir,'data')
    os.makedirs(data_dir)
    os.makedirs(os.path.join(experiment_dir,'taggraph'))
    if (options.simple=="True"):    #TODO: hardcode for now.
        print("1. nonPEAKS version -setupTagGraphProject")
        shutil.move(options.dtadir + '/fileFractionMapping.pck', options.output + '/')
        parsedPecks=options.dtadir+'/'+experiment+'_PEAKS_parsed.csv'
        new_peaks_file=peaks_dir+'/'+experiment+'_PEAKS_parsed.tdv'
        shutil.move(parsedPecks, new_peaks_file)
        DataFile.PEAKS7_split_by_fraction(peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment)
        print("copy current mzXML/mzML files")
        print options.dtadir + '/' + '*%s*mz*ML'%experiment
        for mzXML in glob.glob(options.dtadir + '/' + '*%s*mz*ML'%experiment):
            print mzXML
            shutil.copy(mzXML, data_dir+'/'+os.path.basename(mzXML))
        for mgf in glob.glob(options.dtadir + '/' + '*%s*mgf'%experiment):
            shutil.copy(mgf, data_dir+'/'+os.path.basename(mgf))
    elif (options.ladsresults!=""):
        print("2. use LADS Result file  -setupTagGraphProject")
        print("copy current mzXML files")
        for mzXML in glob.glob(options.ladsresults + '/' + '*%s*mz*ML'%experiment):
            shutil.copy(mzXML,data_dir+'/'+os.path.basename(mzXML))
        for mgf in glob.glob(options.dtadir + '/' + '*%s*mgf'%experiment):
            shutil.copy(mgf, data_dir+'/'+os.path.basename(mgf))
        sourceFracFile=options.ladsresults+"/fileFractionMapping.pck"
        destFrctFile=options.output+"/fileFractionMapping.pck"
        shutil.copy(sourceFracFile,destFrctFile)
        for root,dirs,files in os.walk(options.ladsresults,topdown=False):
            for curfile in files:
                if curfile.endswith('.txt.tg.txt'):
                    sourcefile=options.ladsresults+"/"+curfile
                    curfile=curfile.replace('.txt.tg.txt','')
                    #A375mzXML_PEAKS_parsed_F1.tdv
                    fracID=curfile[-2:]
                    newfileName=options.experimentname+"_PEAKS_parsed_F"+str(int(fracID))+".tdv"
                    destinationFile=peaks_dir+"/"+newfileName
                    shutil.copy(sourcefile,destinationFile) 
    else:
        print("3. use peaks file -setupTagGraphProject")
        print("copy current mzXML/mzML files")
        for mzXML in glob.glob(options.dtadir + '/' + '*%s*mz*ML'%experiment):
            shutil.copy(mzXML, data_dir+'/'+os.path.basename(mzXML))
        for mgf in glob.glob(options.dtadir + '/' + '*%s*mgf'%experiment):
            shutil.copy(mgf, data_dir+'/'+os.path.basename(mgf))
        for experiment, peaks_file in eval(options.peaks).items():
            new_peaks_file = peaks_dir + '/' + experiment + '_PEAKS.csv' #Copy the peaks file to the de_novo subdirectory
            shutil.move(options.dtadir + '/fileFractionMapping.pck', options.output + '/') #Move the pickled fileFractionMapping file from the symLinkDir to the output directory
            if peaks_file.upper().endswith('.XML') or peaks_file.upper().endswith('.PEPXML'):
                shutil.copy(peaks_file, peaks_dir + "/")   #Copy the original pepXML to the de_novo output directory
                mappingFile = open(options.output + '/fileFractionMapping.pck', 'rb')
                fileFractionMapping = pickle.load(mappingFile)
                pepInput.convertPepxmlToCSV(peaks_file, new_peaks_file, fileFractionMapping)
            else:
                # Note: This assumes it will have .CSV extension, if it doesn't have .[pep]XML
                # Might be good to add an additional check
                shutil.copy(peaks_file, new_peaks_file)
            #preprocess peaks result file (convert to tab-delimited .TDV file)
            args = ['--init', options.init, '--symbolmap', os.path.join(RESOURCES_DIR, 'PEAKSSymbMap.txt'), '--peaks', new_peaks_file, '--output', peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment]
            DataFile.executeProcess(PREPROCESSOR_SCRIPT_DIR,'ParsePEAKS7Results.py', args)
            #Unpack results into fractions
            DataFile.PEAKS7_split_by_fraction(peaks_dir+'/'+'%s_PEAKS_parsed.tdv'%experiment)
    print "done with SetUpTagGraphProject.py file."