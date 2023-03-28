import os
import sys
import argparse
import pickle
import shutil
import csv

H=1.00727638

def createFiles(symLinkDir,experiment,DeNovoPath,DeNovoDenovoMZML,DeNovoDelimiter,dataStartAtLine,DeNovofractionID,DeNovoScan,DeNovoCharge,DeNovoRT,precursor_mz,DeNovoPeptide,DeNovoALCinPCT,DeNovoLC,DeNovoScanSplitter,DeNovoPeptideRemove):
    if DeNovoDelimiter=="":
        DeNovoDelimiter='	'
    fileFractionMapping=[]
    currFractionNum=1
    DeNovoDenovoMZMLPairs=DeNovoDenovoMZML.split(";")
    mappingFilename = 'fileFractionMapping.pck'
    mappingFilePath = os.path.join(symLinkDir, mappingFilename)
    outputfile=symLinkDir+'/'+'%s_PEAKS_parsed.csv'%experiment
    fout=open(outputfile,'w')
    fout.write("ScanF\tCharge\tRT\tObs M+H\tPeptide\tALC (%)\tLC\n")
    for filePair in DeNovoDenovoMZMLPairs:
        (denovo,mzml)=filePair.split("|")
        deNovo(fout,DeNovoPath,DeNovoDelimiter,denovo,dataStartAtLine,currFractionNum,DeNovoScan,DeNovoCharge,DeNovoRT,precursor_mz,DeNovoPeptide,DeNovoALCinPCT,DeNovoLC,DeNovoScanSplitter,DeNovoPeptideRemove)
        fileFractionMapping.append(tuple([str(currFractionNum).zfill(2),mzml,str(currFractionNum)]))
        currFractionNum=currFractionNum+1
    mappingOutput = open(mappingFilePath,'wb')
    pickle.dump(fileFractionMapping, mappingOutput)

def deNovo(fout,DeNovoPath,DeNovoDelimiter,denovo,dataStartAtLine,currFractionNum,DeNovoScan,DeNovoCharge,DeNovoRT,precursor_mz,DeNovoPeptide,DeNovoALCinPCT,DeNovoLC,DeNovoScanSplitter,DeNovoPeptideRemove):
    scanPos=int(DeNovoScan)-1
    chargePos=int(DeNovoCharge)-1
    precursor_mzPos=int(precursor_mz)-1
    PeptidePos=int(DeNovoPeptide)-1
    ALCPos=int(DeNovoALCinPCT)-1
    RTPos=0
    LCPos=0
    if DeNovoLC=="?":
        LC='0'
    else:
        LCPos=int(DeNovoLC)-1
    if DeNovoRT=="?":
        RT='0'
    else:
        RTPos=int(DeNovoRT)-1
    rowID=0
    with open(DeNovoPath+"/"+denovo) as denovofile:
        for line in denovofile:
            row=line.strip().split(DeNovoDelimiter)
            rowID=rowID+1
            if rowID<int(dataStartAtLine):
                pass
            else:
                if DeNovoPeptideRemove=="":
                    peptide=row[PeptidePos].strip()
                else:
                    peptide=row[PeptidePos].replace(DeNovoPeptideRemove,"").strip()
                peptide=peptide.replace("(Carbamidomethylation)","").replace("(Oxidation)","#")
                peptide=peptide.replace("(+15.99)","#").replace("(+57.02)","")
                peptide=peptide.replace("L","I")
                if peptide!="" and len(peptide)>=7:
                    scans=row[scanPos]
                    charge=row[chargePos]
                    precursor_mz=row[precursor_mzPos]
                    obsMH=float(precursor_mz) * int(float(charge)) - ((int(float(charge)) - 1) * H)
                    if DeNovoPeptideRemove=="":
                        ALC=row[ALCPos].strip()
                    else:
                        ALC=row[ALCPos].strip().replace(DeNovoPeptideRemove," ")        
                    if DeNovoLC!="?":
                        LC=row[LCPos]
                    if DeNovoRT!="?":
                        RT=row[RTPos]
                    if DeNovoScanSplitter=="":
                        scanF="F%s:%s"%(currFractionNum,scans.strip())
                        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(scanF,charge,RT,obsMH,peptide,LC,ALC))
                    else:
                        scanInfo=scans.split(DeNovoScanSplitter)
                        for j in range(len(scanInfo)):
                            scanF="F%s:%s"%(currFractionNum,scanInfo[j].strip())
                            fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(scanF,charge,RT,obsMH,peptide,LC,ALC))

if __name__ == '__main__':
    createFiles(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13],sys.argv[14],sys.argv[15],sys.argv[16])
