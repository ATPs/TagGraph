#slin 201807    update code for Metaproteome
#slin 20170720  rename pepInput.py to PepInput.py

import csv
import os
import re
import xml.dom
import xml.etree.ElementTree as ET
from xml.etree import ElementTree
from pyteomics import pepxml, auxiliary

#peaksCSV = "/lab/samba/shared/Users/Sam/20160630_Pulldown_dcas9_in_gel_digest_test_DENOVO_5/de_novo_peptides.csv"
#peaksPepXML = "/lab/samba/shared/Users/Sam/20160630_Pulldown_dcas9_in_gel_digest_test_DENOVO_5/de_novo_peptides.xml"
#rankFile = "c:\pep\A375_RPLC_PEAKS_parsed_F1_TAGGRAPH.30.tdv"
#tagGraphpepXML = "TagGraph.pepXML"
#newmzXMLFile = "newA375_F1_5000_f01.tmp.mzXML"

def writeRootNode(tag,dict):
    rootNode = ET.Element(tag)
    for item in dict:
        rootNode.set(item, dict[item])
    return rootNode

def writeInternalNode(node,tag,dict):
    internalNode = ET.SubElement(node, tag)
    for item in dict:
        internalNode.set(item, dict[item])
    return internalNode

def writeLeafNode(node,tag,dict):
    leafNode = ET.SubElement(node, tag)
    for item in dict:
        leafNode.set(item, dict[item])

def readCSV(infile,delimiter=','):
    rowCnt = 0
    header,data = None,dict()
    with open(infile) as f:
        for line in f:
            line = line.strip()
            splitLine = line.split(delimiter)
            if header is None:
                header = splitLine
                continue
            for i in xrange(len(splitLine)):
                data[rowCnt,header[i]] = splitLine[i]
            rowCnt = rowCnt+1
    return rowCnt, data

def writeNodeToXML(rootNode,outFile):
    output_file = open(outFile,'w' )
    output_file.write( '<?xml version="1.0"?>\n' )
    rootNode = indent(rootNode)
    str = xml.etree.ElementTree.tostring(rootNode)
    output_file.write( str )
    output_file.close()

#obtain online 
def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem

def convertPepxmlToCSV(pepXMLFilepath, csvOutputFilepath, fractionMapping=None):
    protonMass = 1.007276466771
    with open(csvOutputFilepath, "w") as fout:
        with pepxml.read(pepXMLFilepath) as reader:
            columnNames = []
            columnNames.append('Fraction')
            columnNames.append('Scan')
            columnNames.append('Source File')
            columnNames.append('Peptide')
            columnNames.append('Tag Length')
            columnNames.append('ALC (%)')
            columnNames.append('length')
            columnNames.append('m/z')
            columnNames.append('z')
            columnNames.append('RT')
            columnNames.append('Area')
            columnNames.append('Mass')
            columnNames.append('ppm')
            columnNames.append('PTM')
            columnNames.append('local confidence (%)')
            columnNames.append('tag (>=0%)')
            columnNames.append('mode')
            fout.write(','.join(columnNames) + '\n')
            for currSpectrum in reader:
                #Ok,lets gather the fields we need:
                fileName = currSpectrum['spectrum']
                fractionNum = 0
                if not fractionMapping is None:
                    for fractionFilePair in fractionMapping:
                        if fractionFilePair[1] == fileName:
                            fractionNum = int(fractionFilePair[0])
                searchHit = currSpectrum['search_hit']
                searchScore = searchHit[0]['search_score']
                scan = currSpectrum['start_scan']
                PTM = ''
                peptide = searchHit[0]['peptide']
                if searchHit[0]['modified_peptide'] != peptide:
                    peptide = searchHit[0]['modified_peptide']
                    PTM = 'Carbamidomethylation'
                #Remove anything between parentheses, to get rid of modifications from the peptide string
                #e.g. LLEGEEC(+57.02)R --> LLEGEECR
                strippedPeptide = re.sub(r'\([^)]*\)', '', peptide)
                tagLength = len( strippedPeptide )
                deNovoScore = int(100.0 * float(searchScore['PeaksDenovoScore']))
                z = int(currSpectrum['assumed_charge'])
                precursorMass = float(currSpectrum['precursor_neutral_mass'])
                calcMass = float(searchHit[0]['calc_neutral_pep_mass'])
                mOverZ = (precursorMass / float(z)) + protonMass
                ppm = searchHit[0]['massdiff']
                RT = '0' # '73.49'
                Area = '6.92E5'
                #Confidence must be converted from "0.96,0.98,0.99,0.99,0.99,0.99,0.99" in pepXML to "96 98 99 99 99 99 99" for CSV
                localConfidenceString = ''
                localConfidenceList = []
                positionalConfidenceString = searchScore['positional_conf']
                positionalConfidenceList = positionalConfidenceString.split(',')
                for currConfidenceDecimalScore in positionalConfidenceList:
                    localConfidenceList.append( str(int(100.0 * float(currConfidenceDecimalScore))) )
                localConfidenceString = ' '.join(localConfidenceList)
                #Now lets spit out the entry
                fout.write(str(fractionNum) + ',')      # Fraction
                fout.write(str(scan) + ',')             # Scan
                fout.write(fileName + ',')              # Source File
                fout.write(peptide + ',')               # Peptide
                fout.write(str(tagLength) + ',')        # Tag Length
                fout.write(str(deNovoScore) + ',')      # ALC (%)
                fout.write(str(tagLength) + ',')        # length
                fout.write(str(mOverZ) + ',')           # m/z
                fout.write(str(z) + ',')                # z
                fout.write(RT + ',')                    # RT is currently set to a bogus constant 0, as the information is not in the pepXML
                fout.write('' + ',')                    # Area is blank for now
                fout.write(str(calcMass) + ',')         # Mass
                fout.write(str(ppm) + ',')              # ppm
                fout.write('' + ',')                    # PTM name is blank for now
                fout.write(localConfidenceString + ',') # local confidence (%)
                fout.write(peptide + ',')               # tag (>=0%)
                fout.write('' + '\n')                   # mode is blank for now

### These two methods both return an ordered list of Pairs (sorted by filename) of the form:
### [ ('01', 'FirstDataFile.mzML'), ('02', 'SecondDataFile.mzML'), ... ]
def getFileFractionMappingFromPepXML(pepXMLFilepath):
        fileList = []
        fileFractionMapping = []
        with pepxml.read(pepXMLFilepath, read_schema=False) as reader:
                #auxiliary.print_tree(next(reader))
                for currSpectrum in reader:
                        if not currSpectrum['spectrum'] in fileList:
                                fileList.append(currSpectrum['spectrum'])
        sortedFileList = sorted(fileList)
        currFractionNum = 1
        for currDataFile in sortedFileList:
                fileFractionMapping.append( tuple([str(currFractionNum).zfill(2), currDataFile]) )
                currFractionNum += 1
        return fileFractionMapping

def getFileFractionMappingFromCSV(csvFilepath):
    fileList = []
    fileFractionMapping = []
    curCounter=1
    with open(csvFilepath) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            currFileName = row['Source File'].split(".")[0]
            if not currFileName in fileList:
                fileList.append(currFileName)
                fileFractionMapping.append( tuple([str(curCounter).zfill(2), currFileName, str(row['Fraction'])]) )
                curCounter=curCounter+1
        return sorted(fileFractionMapping, key=lambda x: x[1])

if __name__ == '__main__':
    pepFileLocation = '/lab/samba/shared/Users/Sam/de_novo_peptides_small.xml'
    getDataFileNames(pepFileLocation)