#slin 201807    update code for Metaproteome
#slin 201807    update code for Metaproteome
#slin 20170720  rename scripts/RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py to scripts/RunTagGraphProteomeEasy.py
#slin 20170720  rename scripts/RUN_TAGGRAPH_HUMAN_PROTEOME_EASY_NR.py to scripts/RunTagGraphProteomeEasy_NR.py
#slin 20170720  rename scripts/parseResultsDB.py to scripts/ParseResultsDB.py
#slin 20170720  rename scripts/generatePepXMLDBperFrac.py to scripts/GeneratePepXMLDBperFrac.py
#slin 20181004  update print function
#slin 20181214  declare veriable, file handle, simple fail or pass will write to file so wev can parse the file to update the status. (runPassFaile.txt)
#slin 20190115  new version by pass LADs, peaks Denovo
#slin 20190208  writeErr message to file can be 
#slin 20190724  simple version

import os
import sys
import time
import ConfigParser
import getopt
import shutil
import glob
import csv
import pickle
import time, datetime
import tempfile
CUR_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1,CUR_DIR)
LIB_DIR = CUR_DIR+"/lib"
sys.path.insert(1,LIB_DIR)
SCRIPTS_DIR=CUR_DIR+"/scripts"
LADS_DIR=CUR_DIR+"/LADS"
sys.path.insert(1,LADS_DIR)
import ArgLib
import DataFile
import PepInput
import VerifyEM

#un#when publish for LADs
#import runLADS

MSGBORDER="---------------------------------------------------------------------------------------"
TAGGRAPH_CONFIG_HEADER =  '\n\n***************************************\n'
TAGGRAPH_CONFIG_HEADER += '*** START TAGGRAPH FROM CONFIG FILE ***\n'
TAGGRAPH_CONFIG_HEADER += '***************************************\n'
TAGGRAPH_CONFIG_FOOTER =  '\n\n*************************************\n'
TAGGRAPH_CONFIG_FOOTER += '*** END TAGGRAPH FROM CONFIG FILE ***\n'
TAGGRAPH_CONFIG_FOOTER += '*************************************\n'
symLinkBaseDir='/tmp/'
ZTestPctToGrab=0.2

def exitTG(fh,fhfailed,fLog,fErr,outputFolder,sourceDir):
    if (outputFolder=="") or (not os.path.exists(outputFolder)):
        errB4Run="Found error while process the configuration file."
        writeTo2FilesStdout(fh,fhfailed,errB4Run)
    fh.close()
    fhfailed.close()
    if (os.path.exists(outputFolder)):  #cp file out, or remove the tmp file
        fLogDest=outputFolder+"/runReport.log"
        fErrDest=outputFolder+"/runFailed.txt"
        if (os.path.getsize(fErr) > 0):
            shutil.move(fErr,fErrDest)
        else:
            os.remove(fErr)
        shutil.move(fLog,fLogDest)  
    else: #write to sourceDir
        fLogDest=sourceDir+"/runReport.log"
        fErrDest=sourceDir+"/runFailed.txt"
        if (os.path.getsize(fErr) > 0):
            errMsg="True"
            shutil.move(fErr,fErrDest)
        else:
            os.remove(fErr)
        shutil.move(fLog,fLogDest)
    sys.exit(1)

def write2FileStdout(fh,message):
    fh.write("%s\n"%message)
    print(message)

def writeTo2FilesStdout(fh1,fh2,message):
    fh2.write("%s\n"%message)
    write2FileStdout(fh1,message)

def printUsage():
    print('\nUsage: ')
    print('   This HELP:')
    print('      $ python '+sys.argv[0]+' -h')
    print('      $ python '+sys.argv[0]+' --help')
    print('   Run TagGraph:')
    print('      $ python '+sys.argv[0]+' <TagGraph Config File>')

def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def isNumeric(value):
    try:
        float(value)
        return True
    except:
        return False

def main(argv):
    try:
        opts,args = getopt.getopt(argv, "h", ["help"])
    except getopt.GetoptError:
        sys.exit(2)
    for opt,arg in opts:
        print('opt: %s' % (opt))
        print('arg: %s' % (arg))
        if opt in ("-h", "--help"):
            printUsage()
            sys.exit(1)
    errCnt=0
    if len(args) != 1:  #Now process the arguments (INI file path)
        printUsage()
        sys.exit(1)
    tgMetaProteome=False
    runLADSProcess=False
                   
    outputFolder=""
    dataDirectory = ''
    configFileName = args[0]
    #create a output file/handle:
    tmpFolder=tempfile.gettempdir()
    (tm_year,tm_mon,tm_mday,tm_hour,tm_min,tm_sec,tm_wday,tm_yday,tm_isdst)=time.localtime(time.time())
    timeNow=str(tm_mon)+str(tm_mday)+str(tm_hour)+str(tm_min)+str(tm_sec)
    runCaptureFN=tmpFolder+'/RunTG'+timeNow+'.txt'
    runFailFN=tmpFolder+"/RunTGFailed"+timeNow+".txt"
    fh=open(runCaptureFN,'w')
    fhfailed=open(runFailFN,'w')
    write2FileStdout(fh,'**** start TagGraph process: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,TAGGRAPH_CONFIG_HEADER)
    write2FileStdout(fh,configFileName+"\n")
    if os.path.isfile(configFileName) and os.access(configFileName,os.R_OK):
        write2FileStdout(fh,MSGBORDER)
        write2FileStdout(fh,"Using Configuration File: %s"%configFileName)
        write2FileStdout(fh,MSGBORDER)
    else:
        msg=' ** FAILURE ** Could not read configuration file: \'%s\'' % (configFileName)
        writeTo2FilesStdout(fh,fhfailed,msg)
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
    theConfig=ConfigParser.ConfigParser()
    theConfig.optionxform=str
    theConfig.read(configFileName)
    generalSectionMap = ConfigSectionMap(theConfig, "General")
    
    tagGraphSectionMap = ConfigSectionMap(theConfig, "TagGraph")
    #'init', 'dtadir', 'peaks', 'output', 'ppmstd', 'modtolerance', 'unimoddict', 'maxcounts', 'modmaxcounts', 'fmindex', 'model', 'config'
    ## Define our Required Arguments ##
    fatalError = False
    ## Arguments that must exist, and be numbers ##
    requiredTagGraphNumericArgs = ['ppmstd','modtolerance','maxcounts','modmaxcounts']
    ## Arguments that must exist, and be paths that point to files that exist and are Readable ##
    requiredTagGraphExistingFiles = ['unimoddict','model','config','init']

    ## Arguments that must exist, and be directories that can be created on the filesystem ##
    requiredTagGraphToCreate = ['output']
    ## Special Arguments:
    # ExperimentName must be a string (space in the string will be remove)
    # d must be a directory, with mzML/mzXML files in it that start with ExperimentName 
    # f must be an fmindex name of the form <basepath>.fm, where <basepath> is the basename and the following files should exist: <basename>.fm.1, <basename>.seqnames.1, <basename>.offsets
    ## Arguments that must exist, and be numbers ##
    if 'runMetaProteome' in generalSectionMap:
        tgMetaProteome=generalSectionMap['runMetaProteome']
        write2FileStdout(fh,'* Found TagGraph Parameter runMetaProteome: \'%s\'' % (tgMetaProteome))
    tg_uselads=False
    tg_useOther=False
    tgSimple=False
    if 'runLADSProcess' in generalSectionMap:
        if generalSectionMap['runLADSProcess']=="True":
            if not theConfig.has_section("LADS"): 
                fatalError=True
                msg='** FAILURE ** Required LADS Section not found in config file'
                writeTo2FilesStdout(fh,fhfailed,msg) 
            else:
                LADSMap=ConfigSectionMap(theConfig,"LADS")
                runLADSProcess=generalSectionMap['runLADSProcess']
                tg_uselads=runLADSProcess
                write2FileStdout(fh,'* Found TagGraph Parameter runLADSProcess: \'%s\'' % (runLADSProcess))
    if 'nonPEAKS' in generalSectionMap:
        tgSimple=generalSectionMap['nonPEAKS']
        write2FileStdout(fh,'* Found TagGraph Parameter nonPEAKS: \'%s\'' % (tgSimple))
    if  tgSimple == 'False' or tgSimple is False:
        requiredTagGraphExistingFiles.append('de_novo')
    tg_ladsResultsPath=""
    if (runLADSProcess): #run LADS
        tg_ladsResultsPath=runLADS.runLADS(configFileName)
    for currArg in requiredTagGraphNumericArgs:
        if currArg in tagGraphSectionMap:
            if isNumeric(tagGraphSectionMap[currArg]):
                write2FileStdout(fh,'* Found Required Numeric TagGraph Parameter \'%s\'  : \'%s\'' % (currArg, tagGraphSectionMap[currArg]))
            else:
                fatalError = True
                msg='** FAILURE ** Required TagGraph Parameter \'%s\' must be a numeric value, found value \'%s\'' % (currArg,tagGraphSectionMap[currArg])
                writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg))
    ## Arguments that must exist, and be paths that point to files that exist and are Readable ##
    for currArg in requiredTagGraphExistingFiles:
        if currArg in tagGraphSectionMap:
            if os.path.isfile(tagGraphSectionMap[currArg]) and os.access(tagGraphSectionMap[currArg], os.R_OK):
                write2FileStdout(fh,'* Found Required Readable File for TagGraph Parameter \'%s\' : \'%s\'' % (currArg, tagGraphSectionMap[currArg]))
            else:
                if not os.path.isfile(tagGraphSectionMap[currArg]):
                    fatalError = True
                    msg='** FAILURE ** Could not find file for Required Parameter \'%s\' at \'%s\'' % (currArg, tagGraphSectionMap[currArg])
                    writeTo2FilesStdout(fh,fhfailed,msg) 
                elif not os.access(tagGraphSectionMap[currArg], os.R_OK):
                    fatalError = True
                    msg='** FAILURE ** Could not Read file for Required Parameter \'%s\' at \'%s\' (check permissions)' % (currArg, tagGraphSectionMap[currArg])
                    writeTo2FilesStdout(fh,fhfailed,msg)                    
        else:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg))
    ## Arguments that must exist, and be directories that should not already exist but can be created on the file system ##
    for currArg in requiredTagGraphToCreate:
        if currArg in tagGraphSectionMap:
            dirToCreate = tagGraphSectionMap[currArg]
            outputFolder=dirToCreate
            if not os.path.exists(dirToCreate):
                try:
                    ## Should be able to make the directory, and then remove it ##
                    os.makedirs(dirToCreate)
                    os.rmdir(dirToCreate)
                    write2FileStdout(fh,'* Found Required Createable Directory for TagGraph Parameter \'%s\' : \'%s\'' % (currArg, dirToCreate))
                except OSError:
                    fatalError = True
                    msg='** FAILURE ** Unable to Create Directory for Required Parameter \'%s\' at \'%s\'' % (currArg, dirToCreate)
                    writeTo2FilesStdout(fh,fhfailed,msg)
            else:
                fatalError = True
                msg='** FAILURE ** File/Directory for Required Parameter \'%s\' at \'%s\' already exists! Should be created by TagGraph' % (currArg, dirToCreate)
                writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            fatalError = True
            msg='** FAILURE ** Required TagGraph Parameter \'%s\' not found in config file' % (currArg)
            writeTo2FilesStdout(fh,fhfailed,msg)       
    experimentName = ''  #ExperimentName must be a string
    if not 'ExperimentName' in tagGraphSectionMap:
        fatalError = True
        msg='** FAILURE ** Required TagGraph Parameter \'ExperimentName\' not found in config file'
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        experimentName = tagGraphSectionMap['ExperimentName'].replace(' ', '').replace('"', '')
        write2FileStdout(fh,'* Found Required TagGraph Parameter ExperimentName: \'%s\'' % (experimentName))
    ## New Method: numFractions = 2, fraction01 = <path to file 1>, fraction02 = <path to file 2>
    numFractions = 0
    foundNumFractions = False
    symLinkDir = symLinkBaseDir
    symLinkDir += experimentName + '_' + str(os.getpid()) + '/'
    write2FileStdout(fh,'* Created temporary sym-link Directory for TagGraph file fraction/mz[X]ML files \'%s\'' % (symLinkDir))
    os.makedirs(symLinkDir) #make the directory, and then remove it
    fileFractionMapping=[]
    DeNovoScanSplitter='?'
    DeNovoDataStartAtLine='2'
    DeNovoFractionID='?'
    DeNovoRT='?'
    DeNovoALCinPCT='?'
    DeNovoLC='?'
    DeNovoPath='?'
    DenovoPeptideDelimiter='?'
    if (tgSimple is True or tgSimple == 'True'):
        if  theConfig.has_section("DeNovo"): #no peaks denovo, will placed the processed denovo to tmp folder first
            DeNovoSectionMap=ConfigSectionMap(theConfig,"DeNovo")
            DeNovoDenovoMZML=DeNovoSectionMap['denovoMZML']
            DeNovoFileDelimiter=DeNovoSectionMap['denovoFileDelimiter']
            if 'de_novo' in DeNovoSectionMap:
                DeNovoPath=DeNovoSectionMap['de_novo']
            else:
                DeNovoPath=tagGraphSectionMap['de_novo']
            DeNovoScan=DeNovoSectionMap['scan']
            if 'scanSplitter' in DeNovoSectionMap:
                DeNovoScanSplitter=DeNovoSectionMap['scanSplitter']
            DeNovoCharge=DeNovoSectionMap['charge']
            if 'precursor_mz' in DeNovoSectionMap:
                precursor_mz=DeNovoSectionMap['precursor_mz']
            DeNovoPeptide=DeNovoSectionMap['peptide']
            if 'peptideDelimiter' in DeNovoSectionMap:
                DenovoPeptideDelimiter=DeNovoSectionMap['peptideDelimiter']
            if 'ignoreRows' in DeNovoSectionMap:
                DeNovoDataStartAtLine=str(int(DeNovoSectionMap['ignoreRows'])+1)
            if 'fractionID' in DeNovoSectionMap:
                DeNovoFractionID=DeNovoSectionMap['fractionID']
            if 'rt' in DeNovoSectionMap:
                DeNovoRT=DeNovoSectionMap['rt']
            if 'ALCinPCT' in DeNovoSectionMap:
                DeNovoALCinPCT=DeNovoSectionMap['ALCinPCT']
            if 'LC' in DeNovoSectionMap:
                DeNovoLC=DeNovoSectionMap['LC']
            writeDeNovoArgs=[]
            writeDeNovoArgs.extend(['\"'+symLinkDir+'\"'])
            writeDeNovoArgs.extend(['\"'+experimentName+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoPath+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoDenovoMZML+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoFileDelimiter+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoDataStartAtLine+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoFractionID+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoScan+'\"'])
            writeDeNovoArgs.extend(['\"'+DeNovoCharge+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DeNovoRT)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(precursor_mz)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DeNovoPeptide)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DeNovoALCinPCT)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DeNovoLC)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DeNovoScanSplitter)+'\"'])
            writeDeNovoArgs.extend(['\"'+str(DenovoPeptideDelimiter)+'\"'])
            write2FileStdout(fh,'*** CALLING denovo.py from runTG.py')
            write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
            write2FileStdout(fh,MSGBORDER+"\n")
            DataFile.executeProcess(SCRIPTS_DIR, 'Denovo.py',writeDeNovoArgs)
            write2FileStdout(fh,'*** command executed: python Denovo.py %s'%writeDeNovoArgs)  
            write2FileStdout(fh,"\n"+MSGBORDER)
            write2FileStdout(fh,'*** END CALLING denovo.py from runTG.py')
            write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
        else:
            fatalError = True
            msg='** FAILURE ** Missing [Denovo] section '
            writeTo2FilesStdout(fh,fhfailed,msg)
                  
    if not 'numFractions' in tagGraphSectionMap:
        ## Check for dataDirectory and automatically finding data files from the de novo files
        if not 'dataDirectory' in tagGraphSectionMap:
            fatalError = True
            msg='** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' not found in config file'
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            dataDirectory = tagGraphSectionMap['dataDirectory']
            if not dataDirectory.endswith('/'):
                dataDirectory += '/'
            if ( not (dataDirectory.startswith('/'))):
                levelup=dataDirectory.count('../')
                if (levelup==0):
                    dataDirectory=CUR_DIR+'/'+ dataDirectory
                else:
                    splitDataDir=dataDirectory.split("/")
                    splitCurDir=CUR_DIR.split("/")
                    tmpD=''
                    for i in xrange(0,len(splitCurDir)-levelup):
                        tmpD = tmpD+splitCurDir[i]+"/"
                    for i in xrange(levelup,len(splitDataDir)-1):
                        tmpD = tmpD+splitDataDir[i]+"/"
                    dataDirectory=tmpD
            sourceDir=dataDirectory
            write2FileStdout(fh,"dataDirectory: %s" % dataDirectory)
            if not os.path.exists(dataDirectory):
                fatalError = True
                #write2FileStdout(fh,'** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' does not exist at: \'%s\'' % (dataDirectory))
                msg='** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' does not exist at: \'%s\'' % (dataDirectory)
                writeTo2FilesStdout(fh,fhfailed,msg)
            elif not os.path.isdir(dataDirectory):
                fatalError = True
                msg='** FAILURE ** Required Directory TagGraph Parameter \'dataDirectory\' does not point to a directory at: \'%s\'' % (dataDirectory)
                writeTo2FilesStdout(fh,fhfailed,msg) 
            else:
                if (runLADSProcess): #slin: deal with LADS data.
                    pass
                else:
                    deNovoFile = tagGraphSectionMap['de_novo']
                    if (tgSimple=='False' or tgSimple is False):
                        if (deNovoFile.upper().endswith('.XML') or deNovoFile.upper().endswith('.PEPXML') or deNovoFile.upper().endswith('.CSV')):
                            if deNovoFile.upper().endswith('.XML') or deNovoFile.upper().endswith('.PEPXML'):
                                fileFractionMapping = PepInput.getFileFractionMappingFromPepXML(deNovoFile)
                            else: ## deNovoFile.upper().endswith('.CSV'):
                                fileFractionMapping = PepInput.getFileFractionMappingFromCSV(deNovoFile)
                                write2FileStdout(fh,'fileFractionMapping: %s'%fileFractionMapping)
                        else:
                            fatalError = True
                            msg='** FAILURE ** Required de novo TagGraph Parameter \'de_novo\' must be named .CSV or .XML/.PEPXML, found \'%s\'' % (deNovoFile)
                            writeTo2FilesStdout(fh,fhfailed,msg)
                        ## Lets write out the fileFractionMapping, pickled for easy reading/writing
                        mappingFilename = 'fileFractionMapping.pck'
                        mappingFilePath = os.path.join(symLinkDir, mappingFilename)
                        mappingOutput = open(mappingFilePath, 'wb')
                        pickle.dump(fileFractionMapping, mappingOutput)
                        mappingOutput.close()
                    else: #or load the pickle file
                        pickleFile=symLinkDir+'/fileFractionMapping.pck'
                        mappingFile = open(pickleFile,'rb')
                        fileFractionMapping = pickle.load(mappingFile)    
                    ##Create a symbolic link pointing to source named link_name.
                    for currFilledFractionNumber, currFilename, currFilledFractionNumberFromFile  in fileFractionMapping:
                        currFilePath = dataDirectory + currFilename  #Check if source file exists
                        if ((not os.path.exists(currFilePath+".mzML")) and (not os.path.exists(currFilePath+".mzXML")) and (not os.path.exists(currFilePath+".mgf"))):
                            fatalError = True
                            msg='** FAILURE ** Data File \'%s\' referenced in de novo file does not exist in dataDirectory \'%s\'' % (currFilename,dataDirectory)
                            writeTo2FilesStdout(fh,fhfailed,msg)
                        elif ((not os.access(currFilePath+".mzML", os.R_OK)) and (not os.access(currFilePath+".mzXML", os.R_OK)) and (not os.access(currFilePath+".mgf", os.R_OK))):
                            fatalError = True
                            msg='** FAILURE ** Data file \'%s\' Not Readable' % (currFilePath)
                            writeTo2FilesStdout(fh,fhfailed,msg)
                        else:
                            if os.path.exists(currFilePath+".mzML"):
                                currFractionFile = currFilePath+".mzML"
                                dataFileSuffix = 'mzML'
                            elif os.path.exists(currFilePath+".mzXML"):
                                currFractionFile = currFilePath+".mzXML"
                                dataFileSuffix = 'mzXML'
                            elif os.path.exists(currFilePath+".mgf"):
                                currFractionFile = currFilePath+".mgf"
                                dataFileSuffix = 'mgf'
                            else:
                                fatalError = True
                                dataFileSuffix = ''
                                msg='** FAILURE ** Data file \'%s\' must end in .mzML, .mzXML, or .mgf!' % (currFractionFile)
                                writeTo2FilesStdout(fh,fhfailed,msg)
                            if not dataFileSuffix == '':
                                symLinkFile = symLinkDir + experimentName + '_f' + currFilledFractionNumber + '.' + dataFileSuffix
                                os.symlink(currFractionFile, symLinkFile)
                                write2FileStdout(fh,'   * Created symLink \'%s\' to data file \'%s\'' % (symLinkFile, currFractionFile))
    else:
        numFractions = tagGraphSectionMap['numFractions']
        if isNumeric(numFractions):
            if float(numFractions).is_integer():
                foundNumFractions = True
                write2FileStdout(fh,'* Found Required integer TagGraph Parameter \'numFractions\'  : \'%s\'' % (numFractions))
                numFractions = int(numFractions)
            else:
                fatalError = True
                msg='** FAILURE ** Required TagGraph Parameter \'numFractions\' must be an integer value, found value \'%s\'' % (numFractions)
                writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            fatalError = True
            #write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'numFractions\' must be a numeric value, found value \'%s\'' % (numFractions))
            msg='** FAILURE ** Required TagGraph Parameter \'numFractions\' must be a numeric value, found value \'%s\'' % (numFractions)
            writeTo2FilesStdout(fh,fhfailed,msg)
    ## If we found numFractions, lets get the paths to the data files and make sym-links to them in a new directory
    ## sym-links will be named <ExperimentName>_f01.mz[X]ML, etc.
    if foundNumFractions:
        #symLinkDir += experimentName + '_' + str(os.getpid()) + '/'
        dataFileSuffix = "mzML"
        for currFraction in xrange(1,numFractions+1):
            filledFractionNumber = str(currFraction).zfill(2)
            if not str('fraction'+filledFractionNumber) in tagGraphSectionMap:
                fatalError = True
                msg='** FAILURE ** Required TagGraph Parameter \'fraction%s\' not found in config file' % (filledFractionNumber)
                writeTo2FilesStdout(fh,fhfailed,msg)
            currFractionFile = tagGraphSectionMap['fraction'+filledFractionNumber]
            if currFractionFile.endswith('mzML'):
                dataFileSuffix = 'mzML'
            elif currFractionFile.endswith('mzXML'):
                dataFileSuffix = 'mzXML'
            elif currFractionFile.endswith('mgf'):
                dataFileSuffix='mgf'
            else:
                fatalError = True
                msg='** FAILURE ** Data file \'%s\' must end in mzML, mzXML or mgf!' % (currFractionFile)
                writeTo2FilesStdout(fh,fhfailed,msg)
            symLinkFile = symLinkDir + experimentName + '_f' + filledFractionNumber + '.' + dataFileSuffix
            os.symlink(currFractionFile, symLinkFile)
            write2FileStdout(fh,'   * Created symLink \'%s\' to data file \'%s\'' % (symLinkFile, currFractionFile))
    fmindexBase = ''
    if not 'fmindex' in tagGraphSectionMap:
        fatalError = True
        msg='** FAILURE ** Required TagGraph Parameter \'fmindex\' (should be the basename of the fmindex files, ending in \'.fm\') not found in config file'
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        fmParam = tagGraphSectionMap['fmindex']
        write2FileStdout(fh,'* Found Required fmindex TagGraph Parameter \'%s\'' % (fmParam))
        if fmParam.endswith('.fm'):
            fmindexBase = fmParam[:-3]
        else:
            fmindexBase = fmParam
        # Now lets check for 3 fmIndex files ending in: .fm.1, .offsets, and .seqnames.1
        fmFile = fmindexBase + ".fm.1"
        fmOffsetFile = fmindexBase + ".offsets"
        fmSeqnamesFile = fmindexBase + ".seqnames.1"
        if not os.path.isfile(fmFile):
            fatalError = True
            #write2FileStdout(fh,'    ** FAILURE ** Could not find required fmindex file at \'%s\'' % (fmFile))
            msg='    ** FAILURE ** Could not find required fmindex file at \'%s\'' % (fmFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        elif not os.access(fmFile, os.R_OK):
            fatalError = True
            #write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex file \'%s\' (check permissions)' % (fmFile))
            msg='    ** FAILURE ** Could not Read required fmindex file \'%s\' (check permissions)' % (fmFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex file at \'%s\'' % (fmFile))
        if not os.path.isfile(fmOffsetFile):
            fatalError = True
            #write2FileStdout(fh,'    ** FAILURE ** Could not find required fmindex Offset file at \'%s\'' % (fmOffsetFile))
            msg='    ** FAILURE ** Could not find required fmindex Offset file at \'%s\'' % (fmOffsetFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        elif not os.access(fmOffsetFile, os.R_OK):
            fatalError = True
            #write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex Offset file \'%s\' (check permissions)' % (fmOffsetFile))
            msg='    ** FAILURE ** Could not Read required fmindex Offset file \'%s\' (check permissions)' % (fmOffsetFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex Offset file at \'%s\'' % (fmOffsetFile))
        if not os.path.isfile(fmSeqnamesFile):
            fatalError = True
            msg='    ** FAILURE ** Could not find required fmindex Seqnames file at \'%s\'' % (fmSeqnamesFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        elif not os.access(fmSeqnamesFile, os.R_OK):
            fatalError = True
            #write2FileStdout(fh,'    ** FAILURE ** Could not Read required fmindex Seqnames file \'%s\' (check permissions)' % (fmSeqnamesFile))
            msg='    ** FAILURE ** Could not Read required fmindex Seqnames file \'%s\' (check permissions)' % (fmSeqnamesFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'   * Found Required readable fmindex Seqnames file at \'%s\'' % (fmSeqnamesFile))
    ### Now lets Check the EM step parameters that can be checked before TG runs ###
    expectationMaximizationSectionMap = ConfigSectionMap(theConfig, "EM")
    ## Arguments that must exist, and be numbers
    # Special Case: EMFractions must be 'all' or a number. Note: EMFractions is now assumed to always be 'all'
    requiredEMNumericArgs = ['maxIterations','initIterations']#,'EMFractions']
    ## Special Arguments:
    ## -o must be a string, the file prefix for the EM Output files (often 'EM_Results')
    ## Arguments that must exist, and be numbers ('EMFractions' is special, as a number or 'all')
    for currArg in requiredEMNumericArgs:
        if currArg in expectationMaximizationSectionMap:
            if isNumeric(expectationMaximizationSectionMap[currArg]):
                write2FileStdout(fh,'* Found Required EM Numeric Parameter \'%s\'  : \'%s\'' % (currArg, expectationMaximizationSectionMap[currArg]))
            else:
                fatalError = True
                #write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'%s\' must be a numeric value, found value \'%s\'' % (currArg, expectationMaximizationSectionMap[currArg]))
                msg='** FAILURE ** Required EM Parameter \'%s\' must be a numeric value, found value \'%s\'' % (currArg, expectationMaximizationSectionMap[currArg])
                writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            fatalError = True
            #write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'%s\' not found in config file' % (currArg))
            msg='** FAILURE ** Required EM Parameter \'%s\' not found in config file' % (currArg)
            writeTo2FilesStdout(fh,fhfailed,msg)    
    ## Now Lets Handle the Special Cases
    # resultsPrefix (Output Prefix) must be a string
    emResultsPrefix = ''
    if not 'resultsPrefix' in expectationMaximizationSectionMap:
        fatalError = True
        #write2FileStdout(fh,'** FAILURE ** Required EM Parameter \'resultsPrefix\' not found in config file')
        msg='** FAILURE ** Required EM Parameter \'resultsPrefix\' not found in config file'
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        emResultsPrefix = expectationMaximizationSectionMap['resultsPrefix']
        write2FileStdout(fh,'* Found Required EM Parameter \'resultsPrefix\': \'%s\'' % (emResultsPrefix))
    #options = ArgLib.parse(['init', 'dtadir', 'peaks', 'output', 'ppmstd', 'modtolerance', 'unimoddict', 'maxcounts', 'modmaxcounts', 'fmindex', 'model', 'config'], optArgs=[{'opts': ('-x', '--splittaxon'), 'attrs': {'dest': 'splittaxon', 'action': 'store_true', 'default': False, 'help': 'Flag. For searches of metaproteomic databases, split identical context entries by taxon for accurate consideration via EM.'}}])
    tg_cutoff=10
    tg_taxonfastadir=""
    tg_alpha=1
    tg_taxonomy=""
    if tgMetaProteome:
        if not 'taxonfastadir' in tagGraphSectionMap:
            fatalError = True
            write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'taxonfastadir\' not found in config file')
        else:
            tg_taxonfastadir = tagGraphSectionMap['taxonfastadir']
            if tg_taxonfastadir=="":
                fatalError=True
                #write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'taxonfastadir\' is empty')
                msg='** FAILURE ** Required TagGraph Parameter \'taxonfastadir\' is empty'
                writeTo2FilesStdout(fh,fhfailed,msg) 
            else:           
                write2FileStdout(fh,'* Found Required TagGraph Parameter taxonfastadir: \'%s\'' % (tg_taxonfastadir))    
        if not 'cutoff' in tagGraphSectionMap:
            fatalError=True
            #write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'cutoff\' not found in config file')
            msg='** FAILURE ** Required TagGraph Parameter \'cutoff\' not found in config file'
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            tg_cutoff = tagGraphSectionMap['cutoff']
            write2FileStdout(fh,'* Found Required TagGraph Parameter cutoff: \'%s\'' % (tg_cutoff))
        if not 'alpha' in tagGraphSectionMap:
            fatalError=True
            #write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'alpha\' not found in config file')
            msg='** FAILURE ** Required TagGraph Parameter \'alpha\' not found in config file'
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            tg_alpha=tagGraphSectionMap['alpha']
            write2FileStdout(fh,'* Found Required TagGraph Parameter alpha: \'%s\'' % (tg_alpha))
        if not 'taxonomy' in tagGraphSectionMap:
            fatalError=True
            #write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'taxonomy\' not found in config file')
            msg='** FAILURE ** Required TagGraph Parameter \'taxonomy\' not found in config file'
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            tg_taxonomy=tagGraphSectionMap['taxonomy']
            if tg_taxonomy=="":
                fatalError=True
                write2FileStdout(fh,'** FAILURE ** Required TagGraph Parameter \'taxonomy\' is empty')
                msg='** FAILURE ** Required TagGraph Parameter \'taxonomy\' is empty'
                writeTo2FilesStdout(fh,fhfailed,msg)
            else:
                write2FileStdout(fh,'* Found Required TagGraph Parameter taxonomy: \'%s\'' % (tg_taxonomy))
    tg_fasta =""
    if not tgMetaProteome:
        if not 'fasta' in tagGraphSectionMap:
            tg_fasta=fmindexBase+".fasta"
            if not os.path.exists(tg_fasta):
                fatalError = True
                msg='** FAILURE ** Required TagGraph fasta file %s does not exist!!' % (tg_fasta)
                writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            tg_fasta=tagGraphSectionMap['fasta']
    tgArgs = []
    ### If a fatal error was thrown, do not proceed ###
    if fatalError == True:
        write2FileStdout(fh,'*****  HALTING DUE TO FATAL ERROR IN TAGGRAPH OR EM PARAMETERS, SEE OUTPUT ABOVE!!! ')
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
    ## Lets set up the args properly for RUN_TAGGRAPH_HUMAN_PROTEOME_EASY.py ##
    tg_ppmstd = str(tagGraphSectionMap['ppmstd'])
    tg_modtolerance = str(tagGraphSectionMap['modtolerance'])
    tg_maxcounts = str(tagGraphSectionMap['maxcounts'])
    tg_modmaxcounts = str(tagGraphSectionMap['modmaxcounts'])
    tg_config = tagGraphSectionMap['config']
    tg_init = tagGraphSectionMap['init']
    tg_dtadir = symLinkDir # tagGraphSectionMap['d']
    tg_model = tagGraphSectionMap['model']
    tg_output = tagGraphSectionMap['output'].replace(' ','').replace('"','')
    tg_unimoddict = tagGraphSectionMap['unimoddict']
    tg_fmindex = tagGraphSectionMap['fmindex']
    tg_experimentname=tagGraphSectionMap['ExperimentName'].replace(' ','').replace('"','')
    tg_peaks = '{\'' + tg_experimentname + '\': \'' + tagGraphSectionMap['de_novo'] + '\'}' # K = "{'e009133': '/lab/samba/shared/Users/Sam/20160630_Pulldown_dcas9_in_gel_digest_test_DENOVO_5/de_novo_peptides.csv'}"
    ### tg_output directory will now end with a slash
    if not tg_output.endswith('/'):
        tg_output += '/'
    if tgMetaProteome:
        tg_fasta = tg_output+tg_experimentname+"/fmindex/"+tg_experimentname+"_targeted.fasta"
    tgArgs = []
    tgArgs.extend(['-p','\"' + tg_ppmstd + '\"'])
    tgArgs.extend(['-l','\"' + tg_modtolerance + '\"'])
    tgArgs.extend(['-M','\"' + tg_maxcounts + '\"'])
    tgArgs.extend(['-C','\"' + tg_modmaxcounts + '\"'])
    tgArgs.extend(['-c','\"' + tg_config + '\"'])
    tgArgs.extend(['-i','\"' + tg_init + '\"'])
    tgArgs.extend(['-d','\"' + tg_dtadir + '\"'])
    tgArgs.extend(['-m','\"' + tg_model + '\"'])
    tgArgs.extend(['-o','\"' + tg_output + '\"'])
    tgArgs.extend(['-Q','\"' + tg_unimoddict + '\"'])
    tgArgs.extend(['-f','\"' + tg_fmindex + '\"'])
    tgArgs.extend(['-K','\"' + tg_peaks + '\"'])
    tgArgs.extend(['-E','\"' + str(tg_experimentname) + '\"'])
    tgArgs.extend(['-Z','\"' + str(tg_uselads) + '\"'])
    tgArgs.extend(['-j','\"' + str(tgSimple) + '\"'])
    tgArgs.extend(['-z','\"' + tg_ladsResultsPath + '\"'])
    if tgMetaProteome:
        tgArgs.extend(['-H','\"'+tg_taxonfastadir+'\"'])
        tgArgs.extend(['-b','\"'+str(tg_cutoff)+'\"'])
        tgArgs.extend(['-r','\"'+str(tg_alpha)+'\"'])
        tgArgs.extend(['-y','\"'+tg_taxonomy+'\"'])
    write2FileStdout(fh,'\nTG ARGS: %s\n\n'%tgArgs)
    write2FileStdout(fh,MSGBORDER)
    scriptName = ""
    if tgMetaProteome:
        scriptName="RunTagGraphProteomeEasy_NR.py"
    else:
        scriptName="RunTagGraphProteomeEasy.py"
    write2FileStdout(fh,'*** CALLING %s from runTG.py'%scriptName)
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR,'%s'%scriptName, tgArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING %s from runTG.py'%scriptName)
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER)
    ### VERIFY TG RUN ### 
    '''
    Now lets check the TG output to make sure it ran correctly. We'll check for:
    * <output_dir>/results.db should exist and have size > 0 (do actual db check?)
    * The files <output_dir>/<experiment_name>_addPlausibleMods_poss_[combo/single]_mods.tdv both exist and have reasonable sizes
    * Check that output_dir/<experiment_name>/data/ contains directories of DTA files named <experiment_name>_f01/ etc
    * Check that output_dir/<experiment_name>/de_novo/<experiment_name>_PEAKS.csv/PEAKS_parsed.tdv/PEAKS_parsed_F1.tdv etc exist
    * Check that output_dir/<experiment_name>/taggraph/<experiment_name>_PEAKS_parsed_F1_TAGGRAPH.tdv etc exist
    * output_dir/<experiment_name>_CHECK.txt.<numFractions> contains count numbers for each fraction:
    '''
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** VERIFYING TAGGRAPH OUTPUTS in runTG.py ')
    write2FileStdout(fh,MSGBORDER)
    minDBFileSize = 1000000 ## 1Megabyte minimum db size after TG runs?
    minAddPlausibleModsFileSize = 1000 ## 10kBytes min size for <experiment_name>_addPlausibleMods_[combo/single]_mods.tdv files #sarahl 20181213 change from 2000 to 1000
    ## <output_dir>/results.db should exist and have size > 0 (do actual db check?)
    dbFile = tg_output + 'results.db'
    if not os.path.exists(dbFile):
        fatalError = True
        msg='** FAILURE ** Required SQLITE DB File %s does not exist!!' % (dbFile)
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        dbFileSize = os.path.getsize(dbFile)
        if dbFileSize < minDBFileSize:
            fatalError = True
            msg='** FAILURE ** Required SQLITE DB File %s is too small: %d Bytes!!' % (dbFile, dbFileSize)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required SQLITE DB File \'%s\', size %d Bytes OK' % (dbFile, dbFileSize))
    ## The files <output_dir>/<experiment_name>_addPlausibleMods_poss_[combo/single]_mods.tdv both exist
    singleModsFile = tg_output + experimentName + '_addPlausibleMods_poss_single_mods.tdv'
    comboModsFile  = tg_output + experimentName + '_addPlausibleMods_poss_combo_mods.tdv'
    if not os.path.exists(singleModsFile):
        fatalError = True
        msg='** FAILURE ** Required Single Mods File %s does not exist!!' % (singleModsFile)
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        singleModsFileSize = os.path.getsize(singleModsFile)
        if singleModsFileSize < minAddPlausibleModsFileSize:
            fatalError = True
            msg='** FAILURE ** Required Single Mods File %s is too small: %d Bytes!!' % (singleModsFile,singleModsFileSize)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required Single Mods File %s, size %d Bytes OK' % (singleModsFile, singleModsFileSize))
    if not os.path.exists(comboModsFile):
        fatalError = True
        msg='** FAILURE ** Required Combo Mods File %s does not exist!!' % (comboModsFile)
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        comboModsFileSize = os.path.getsize(comboModsFile)
        if comboModsFileSize < minAddPlausibleModsFileSize:
            fatalError = True
            msg='** FAILURE ** Required Combo Mods File \'%s\' is too small: %d Bytes!!' % (comboModsFile, comboModsFileSize)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required Combo Mods File \'%s\', size %d Bytes OK' % (comboModsFile, comboModsFileSize))
    ## Check that output_dir/<experiment_name>/data/ contains directories of DTA files named <experiment_name>_f01/ etc
    dataDir = tg_output + experimentName + '/data/'
    for currFraction in xrange(1,numFractions+1):
        filledFractionNumber = str(currFraction).zfill(2)
        currDtaDirName = dataDir + experimentName + '_f' + filledFractionNumber
        if not os.path.exists(currDtaDirName):
            fatalError = True
            msg='** FAILURE ** Missing directory of DTA files at: \'%s\'' % (currDtaDirName)
            writeTo2FilesStdout(fh,fhfailed,msg)
        elif not os.path.isdir(currDtaDirName):
            fatalError = True
            msg='** FAILURE ** \'%s\' exists but is not a Directory!' % (currDtaDirName)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found DTA directory: \'%s\'' % (currDtaDirName))
    ## Check that output_dir/<experiment_name>/de_novo/<experiment_name>_PEAKS.csv/PEAKS_parsed.tdv/PEAKS_parsed_F1.tdv etc exist
    deNovoDir = tg_output + experimentName + '/de_novo/'
    deNovoCSV = deNovoDir + experimentName + '_PEAKS.csv'
    peaksParsed = deNovoDir + experimentName + '_PEAKS_parsed.tdv'
    fractionsParsedBase = deNovoDir + experimentName + '_PEAKS_parsed_F'
    if (not tg_uselads) and (not tgSimple):
        if not os.path.exists(deNovoCSV):
            fatalError = True
            msg='** FAILURE ** Missing de novo CSV File \'%s\' !!' % (deNovoCSV)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required de novo CSV File \'%s\'' % (deNovoCSV))
        if not os.path.exists(peaksParsed):
            fatalError = True
            msg='** FAILURE ** Missing Parsed de novo File \'%s\' !!' % (peaksParsed)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required Parsed de novo File \'%s\'' % (peaksParsed))
    for currFraction in xrange(1,numFractions+1):
        currParsedFractionFile = fractionsParsedBase + str(currFraction) + '.tdv'
        if not os.path.exists(currParsedFractionFile):
            fatalError = True
            msg='** FAILURE ** Missing Parsed de novo Fraction File \'%s\' !!' % (currParsedFractionFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required Parsed de novo Fraction File \'%s\'' % (currParsedFractionFile))
    ## Check that output_dir/<experiment_name>/taggraph/<experiment_name>_PEAKS_parsed_F1_TAGGRAPH.tdv etc exist
    taggraphDir = tg_output + experimentName + '/taggraph/'
    taggraphParsedBase = taggraphDir + experimentName + '_PEAKS_parsed_F'
    taggraphParsedSuffix = '_TAGGRAPH.tdv'
    for currFraction in xrange(1,numFractions+1):
        currTaggraphFractionFile = taggraphParsedBase + str(currFraction) + taggraphParsedSuffix
        if not os.path.exists(currTaggraphFractionFile):
            fatalError = True
            curErrMsg='** FAILURE ** Missing Parsed TagGraph Fraction File \'%s\' !!'%(currTaggraphFractionFile)
            writeTo2FilesStdout(fh,fhfailed,msg)
        else:
            write2FileStdout(fh,'* Found Required Parsed TagGraph Fraction File \'%s\''%(currTaggraphFractionFile))
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END VERIFYING TAGGRAPH OUTPUTS in runTG.py')
    write2FileStdout(fh,MSGBORDER)
    ### END VERIFY TG RUN ###
    ### If a fatal error was thrown, do not proceed ###
    if fatalError == True:
        msg='*****  HALTING DUE TO FATAL ERROR IN VERIFYING TAGGRAPH RUN, SEE OUTPUT ABOVE!!'
        writeTo2FilesStdout(fh,fhfailed,msg)
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
        sys.exit(1)
    ## Copy configuration file to output tree for safe keeping ##
    configFileBaseName = os.path.basename(configFileName)
    checkConfigDestination = tg_output
    if os.path.exists(checkConfigDestination + configFileBaseName):
        write2FileStdout(fh,'** WARNING ** config file \'%s\' already exists in output directory \'%s\'' % (configFileBaseName,checkConfigDestination))
    else:
        shutil.copy(configFileName, checkConfigDestination)
        write2FileStdout(fh,'* Successfully copied Configuration File \'%s\' to Output Directory \'%s\'' % (configFileName,checkConfigDestination))
    ## Lets set up the args properly for ComputeEMProbabilitiesFromDB.py ##
    '''
    -i: same as TG -i parameter
    -F all
    -M 100
    -C 20
    -B = <-o parameter from TG>/results.db [checked after TG runs]
    -E: Same as TG ExperimentName parameter.
    -o: Output Prefix, will create files with the prefix <EM -o parameter> in the directory specified by the <TG -o parameter>
    '''
    #eg: python /home/slin/TagGraph/taggraphLADS/scripts/ComputeEMProbabilitiesFromDB_NR.py -i ../resources/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini -E Codexis_glycosylation_ITCID_LP6301 -M 100 -C 20 -F all -o /home/slin/TagGraph/sampleInputFiles/TG_Codexis_glycosylation_DENOVO_2_ITCID_LP6301/err/EM_Results -B /home/slin/TagGraph/sampleInputFiles/TG_Codexis_glycosylation_DENOVO_2_ITCID_LP6301/err/results.db -T /home/slin/TagGraph/sampleInputFiles/TG_Codexis_glycosylation_DENOVO_2_ITCID_LP6301/err/Codexis_glycosylation_ITCID_LP6301/fmindex/Codexis_glycosylation_ITCID_LP6301_GetTaxonsFMMATCH.tdv -y /home/slin/TagGraph/taxonomy/genus_info.pck
    em_init = tg_init
    em_fractions = 'all' ## EMFractions is always 'all' now! ## = str(expectationMaximizationSectionMap['EMFractions'])
    em_maxIterations = str(expectationMaximizationSectionMap['maxIterations'])
    em_initIterations = str(expectationMaximizationSectionMap['initIterations'])
    em_dbLocation = tg_output + 'results.db'
    em_experimentName = tg_experimentname
    em_output = tg_output
    if not em_output.endswith('/'):
        em_output += '/'
    em_output += emResultsPrefix
    emTaxonsFMMATCH=tg_output+"/"+em_experimentName+"/fmindex/"+em_experimentName+"_GetTaxonsFMMATCH.tdv"
    em_all='all'
    emArgs = []
    emArgs.extend(['-i', '\"' + em_init + '\"'])
    emArgs.extend(['-F', '\"' + em_fractions + '\"'])
    emArgs.extend(['-M', '\"' + em_maxIterations + '\"'])
    emArgs.extend(['-C', '\"' + em_initIterations + '\"'])
    emArgs.extend(['-B', '\"' + em_dbLocation + '\"']) 
    emArgs.extend(['-E', '\"' + em_experimentName + '\"'])
    emArgs.extend(['-o', '\"' + em_output + '\"']) 
    if tgMetaProteome:
        emArgs.extend(['-T','\"'+emTaxonsFMMATCH+'\"'])
        emArgs.extend(['-y','\"'+tg_taxonomy+'\"'])
    write2FileStdout(fh,'EM ARGS: %s\n'% emArgs)
    write2FileStdout(fh,MSGBORDER)
    if tgMetaProteome:
        scriptName="ComputeEMProbabilitiesFromDB_NR.py"
    else:
        scriptName="ComputeEMProbabilitiesFromDB.py"
    write2FileStdout(fh,'*** CALLING %s from runTG.py'%scriptName)
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR,'%s'%scriptName, emArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING %s from runTG.py'%scriptName)
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER)
    EMProbs_TOPONLY = tg_output + 'EM_Results_EMProbs_END_TOPONLY.tdv'
    if not os.path.exists(EMProbs_TOPONLY):
        fatalError = True
        msg='** FAILURE ** Missing EMProbs END TOPONLY file \'%s\'.' % (EMProbs_TOPONLY)
        writeTo2FilesStdout(fh,fhfailed,msg)
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
        sys.exit(1)
    else:
        write2FileStdout(fh,'* Found EMProbs END TOPONLY file \'%s\'' % (EMProbs_TOPONLY))
    #run CalculateProteinTaxonAbundance.py when it's meta Proteome version
    topResultsFile = tg_output + experimentName + '_TopResults.tdv'
    resultsDBFile = tg_output +"results.db"
    #map back the correct I/L peptide. hardcode path now
    #fasta file path from above.
    #eg: python /home/slin/TagGraph/taggraphLADS/scripts/CorrectILPeptide.py "/home/slin/TagGraph/sampleInputFiles/FMIndices/human_uniprot_12092014_crap.fasta" "/exp/exp_TopResults.tdv" "/exp/results.db"
    ILPeptidesArgs=[]
    ILPeptidesArgs.extend(['\"' + tg_fasta + '\"'])
    ILPeptidesArgs.extend(['\"' + topResultsFile + '\"'])
    ILPeptidesArgs.extend(['\"' + resultsDBFile + '\"'])
    #python scripts/CorrectILPeptide.py 
    write2FileStdout(fh,'*** CALLING CorrectILPeptide.py from runTG.py')
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR,'CorrectILPeptide.py',ILPeptidesArgs)
    #write2FileStdout(fh,'*** command executed: python CorrectILPeptide.py %s'%ILPeptidesArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    if tgMetaProteome:
        #eg: python /home/slin/tgLatest/taggraph/scripts/CalculateProteinTaxonAbundance.py -y /home/slin/taggraph/taxonomy/genus_info.pck -m '' -r 0.99 -T /home/slin/sampleInputFiles/TG/fusion1Frac/fusion1Frac_TopResults.tdv -o /home/slin/sampleInputFiles/TG/fusion1Frac/fusion1Frac_TopResults_v2.txt
        taOutputPrefix=tg_output
        if not taOutputPrefix.endswith('/'):
            taOutputPrefix += '/'
        taOutputPrefix+=tg_experimentname
        taOutput=taOutputPrefix+"_TaxonAbundance"
        taInput=taOutputPrefix+"_TopResults.tdv.new"
        taArgs=[]
        taArgs.extend(['-y','\"'+tg_taxonomy+'\"'])
        taArgs.extend(['-m','\" \"']) 
        taArgs.extend(['-r','\"0.99\"'])
        taArgs.extend(['-T','\"'+taInput+'\"'])
        taArgs.extend(['-o','\"'+taOutput+'\"']) 
        write2FileStdout(fh,'Taxon Abundance ARGS: %s\n'% taArgs)
        write2FileStdout(fh,MSGBORDER)
        scriptName="CalculateProteinTaxonAbundance.py"
        write2FileStdout(fh,'*** CALLING %s from runTG.py'% scriptName)
        write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
        write2FileStdout(fh,MSGBORDER+"\n")
        DataFile.executeProcess(SCRIPTS_DIR,'%s'%scriptName, taArgs)
        write2FileStdout(fh,"\n"+MSGBORDER)
        write2FileStdout(fh,'*** END CALLING %s from runTG.py'%scriptName)
        write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
        write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,"\n\n"+MSGBORDER)
    write2FileStdout(fh,'*** CALLING verify EM result tests from runTG.py')
    write2FileStdout(fh,"\ntime now: @ %s"%datetime.datetime.now())
    result=VerifyEM.verifyEM(tg_output,tgMetaProteome)
    write2FileStdout(fh,result)
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,"\ntime now: @ %s"%datetime.datetime.now())
    write2FileStdout(fh,'*** END CALLING verify EM result tests from runTG.py')
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER)
    topResultsFile=tg_output+experimentName+'_TopResults.tdv'
    if not os.path.exists(topResultsFile):
        fatalError = True
        msg8u='** FAILURE ** Missing TopResult file \'%s\'.' % (topResultsFile)
        writeTo2FilesStdout(fh,fhfailed,msg)
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
        sys.exit(1)
    else:
        write2FileStdout(fh,'* Found TopResult file \'%s\'' % (topResultsFile))
    outputPerFraction="No"
    #eg: python /home/slin/TagGraph/taggraphLADS/scripts/ParseResultsDB.py "PATH of EXPOUTPUT_FOLDER" "/home/slin/TagGraph/sampleInputFiles/TG/TAG_GRAPH_Tryp_CysCarbam_MetOx.ini" "Yes" "0.01" "2" "5" "0"
    write2FileStdout(fh,'**** start ParseResultsDB process: %s'%(datetime.datetime.now()))
    FDRCutoff=0.01
    logEMCutoff=100
    DisplayProteinNum=5
    if "outputPerFraction" in generalSectionMap:
        if True==theConfig.getboolean('General','outputPerFraction'):
            outputPerFraction="Yes"
    if "FDRCutoff" in generalSectionMap:
        if isNumeric(generalSectionMap["FDRCutoff"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("FDRCutoff", generalSectionMap["FDRCutoff"]))
            FDRCutoff=generalSectionMap['FDRCutoff']
    if "logEMCutoff" in generalSectionMap:
        if isNumeric(generalSectionMap["logEMCutoff"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("logEMCutoff", generalSectionMap["logEMCutoff"]))
            logEMCutoff=generalSectionMap['logEMCutoff']
    if "DisplayProteinNum" in generalSectionMap:
        if isNumeric(generalSectionMap["DisplayProteinNum"]):
            write2FileStdout(fh,'* Found  Numeric TagGraph Parameter \'%s\'  : \'%s\'' % ("DisplayProteinNum", generalSectionMap["DisplayProteinNum"]))
            DisplayProteinNum=generalSectionMap['DisplayProteinNum']
    if tgMetaProteome:  #for metaProteome
        metaProteome=1
    else:
        metaProteome=0
    writeTopArgs=[]
    writeTopArgs.extend(['\"'+tg_output+'\"'])
    writeTopArgs.extend(['\"'+tg_init+'\"'])
    writeTopArgs.extend(['\"'+outputPerFraction+'\"'])
    writeTopArgs.extend(['\"'+str(FDRCutoff)+'\"'])
    writeTopArgs.extend(['\"'+str(logEMCutoff)+'\"'])
    writeTopArgs.extend(['\"'+str(DisplayProteinNum)+'\"'])
    writeTopArgs.extend(['\"'+str(metaProteome)+'\"'])
    ## Now lets parse the original TG tab-delimted format ##
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'*** CALLING ParseResultsDB.py from runTG.py')
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER+"\n")
    DataFile.executeProcess(SCRIPTS_DIR, 'ParseResultsDB.py',writeTopArgs)
    write2FileStdout(fh,'*** command executed: python ParseResultsDB.py %s'%writeTopArgs)
    write2FileStdout(fh,"\n"+MSGBORDER)
    write2FileStdout(fh,'*** END CALLING ParseResultsDB.py from runTG.py')
    write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
    write2FileStdout(fh,MSGBORDER)
    topResultsFinalFile = tg_output+experimentName+'_TopResults*.txt'
    foundFile=0
    if len(glob.glob(topResultsFinalFile))>0:
        foundFile=1
    if foundFile==0:
        fatalError=True
        msg='** FAILURE ** Missing result file \'%s\' from ParseResultsDB.py process. Please check.' % (topResultsFinalFile)
        writeTo2FilesStdout(fh,fhfailed,msg)
        exitTG(fh,fhfailed,runCaptureFN,runFailFN,outputFolder,sourceDir)
        sys.exit(1)
    if 'generatePepXML' in generalSectionMap:
        if theConfig.getboolean('General','generatePepXML'):  ## Now lets generate the output in PepXML format
            pepArgs = []
            pepArgs.extend(['\"'+tg_init+'\"'])
            pepArgs.extend(['\"'+tg_ppmstd+'\"'])
            pepArgs.extend(['\"'+tg_modtolerance+'\"'])
            pepArgs.extend(['\"'+tg_maxcounts+'\"'])
            pepArgs.extend(['\"'+tg_modmaxcounts+'\"'])
            pepArgs.extend(['\"'+tg_fmindex+'\"'])  # tagGraphSectionMap['fmindex']
            pepArgs.extend(['\"'+tg_model+'\"'])  # tagGraphSectionMap['model']
            #pepArgs.extend(['\"'+tg_config+'\"'])  # tagGraphSectionMap['config']
            pepArgs.extend(['\"'+tg_unimoddict+'\"'])  # tagGraphSectionMap['unimoddict']
            pepArgs.extend(['\"'+tg_output+'\"'])  # tagGraphSectionMap['output']
            pepArgs.extend(['\"'+tg_dtadir+'\"'])  # symLinkDir
            pepArgs.extend(['\"'+str(FDRCutoff)+'\"'])
            pepArgs.extend(['\"'+str(logEMCutoff)+'\"'])
            write2FileStdout(fh,MSGBORDER)
            write2FileStdout(fh,'*** CALLING GeneratePepXMLDBperFrac.py from runTG.py')
            write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
            write2FileStdout(fh,MSGBORDER+"\n")
            DataFile.executeProcess(SCRIPTS_DIR,'generatePepXMLDBperFrac.py',pepArgs)
            write2FileStdout(fh,'*** command: python generatePepXMLDBperFrac.py %s'%pepArgs)
            write2FileStdout(fh,"\n"+MSGBORDER)
            write2FileStdout(fh,'*** END CALLING generatePepXMLDBperFrac.py from runTG.py')
            write2FileStdout(fh,'*** time now: %s'%(datetime.datetime.now()))
            write2FileStdout(fh,MSGBORDER)
            '''
            Now lets clean up our temporary items and copied data files as configured! ###
            We need to:
            * Remove the sym-link directory in /tmp/ (symLinkDir)
            * If cleanMzDataFilesFromOutput is True, clean the dataDir (<output_dir>/<experiment_name>/data/)
            directory of mz[X]ML files and the DTA directories of the same name
            '''
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'***    CLEANING UP')
    write2FileStdout(fh,MSGBORDER)
    shutil.rmtree(symLinkDir) #Remove the sym-link directory in /tmp/ (symLinkDir)
    if os.path.exists(symLinkDir):
        #write2FileStdout(fh,'** FAILURE ** Failed to removed temporary symbolic link directory \'%s\'' % (symLinkDir))
        msg='** FAILURE ** Failed to removed temporary symbolic link directory \'%s\'' % (symLinkDir)
        writeTo2FilesStdout(fh,fhfailed,msg)
    else:
        write2FileStdout(fh,'* Successfully removed temporary symbolic link directory \'%s\'' % (symLinkDir))
    if 'cleanInputDataFilesFromOutput' in generalSectionMap:
#        if True == theConfig.getboolean('General','cleanInputDataFilesFromOutput'):
        if theConfig.getboolean('General','cleanInputDataFilesFromOutput'):
            shutil.rmtree(dataDir)
            #os.makedirs(dataDir)
            write2FileStdout(fh,'* Removed mz[X]ML and DTA files from data directory \'%s\' (cleanInputDataFilesFromOuput is True)' % (dataDir))
        else:
            write2FileStdout(fh,'* Leaving mz[X]ML and DTA files in data directory \'%s\' (cleanInputDataFilesFromOuput is False)' % (dataDir))
    if 'cleanIntermediateFiles' in generalSectionMap:
        denovoOutputDir=tg_output+'/'+experimentName+'/de_novo/'
        taggraphOutputDir=tg_output+'/'+experimentName+'/taggraph/'
        experimentOutputDir=tg_output+'/'+experimentName
        if theConfig.getboolean('General','cleanIntermediateFiles'):
            shutil.rmtree(denovoOutputDir)
            shutil.rmtree(taggraphOutputDir)
            if os.path.exists(dataDir):
                shutil.rmtree(dataDir)
            shutil.rmtree(experimentOutputDir)
            files = os.listdir(tg_output)
            for file in files:
                if (file.endswith(".tdv") or file.endswith(".tdv.new") or (file.find("_CHECK.txt.") > 0) or file.endswith(".db") or file.endswith(".log")):
                    if (os.path.exists(os.path.join(tg_output,file))):
                        write2FileStdout(fh,"remove %s" % os.path.join(tg_output,file))
                        os.remove(os.path.join(tg_output,file))
                    else:
                        write2FileStdout(fh,"keeper %s" % os.path.join(tg_output,file))
            write2FileStdout(fh,'* Removed mz[X]ML and Intermediate files from output directory \'%s\' (cleanIntermediateFiles is True)' % (dataDir))
        else:
            write2FileStdout(fh,'* Leaving mz[X]ML and Intermediate files in output directory \'%s\' (cleanIntermediateFiles is False)' % (dataDir))
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'***  END CLEANING UP')
    write2FileStdout(fh,MSGBORDER)
    write2FileStdout(fh,'%s'%TAGGRAPH_CONFIG_FOOTER)
    write2FileStdout(fh,'**** end TagGraph process: %s'%(datetime.datetime.now()))
    fh.close()
    #move file back to output folder:
    toDest=tg_output+"runReport.log"
    shutil.move(runCaptureFN,toDest)
    if (os.path.getsize(runFailFN) ==0): #slin since no error, it should always==0
        os.remove(runFailFN)
    sys.exit(0)

if __name__ == '__main__':
    main(sys.argv[1:])
