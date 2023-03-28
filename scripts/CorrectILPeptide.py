#20180501: slin take care the decoy case
#20180524: slin add to return the full accession as the key.
import sys
import subprocess
from itertools import islice
from decimal import Decimal
import os.path
import re
import string
import os.path
import sqlite3 as lite

MARKER=">"
PADDING=0
DELIMITER="\t"

def LIChange(sequence):
    sequence=sequence.replace("L","I")
    return sequence

def proteinAccessionSequenceFromfastaLI(fasta):
    listAnnoationOrder=[]
    dictProteinSequence={}
    dictProteinLISequence={}
    proteinSequence=""
    accession=""
    rowCnt=0
    with open(fasta) as lines:
        for line in lines:
            rowCnt=rowCnt+1
            line=line.strip()
            if (rowCnt>0) and (MARKER not in line):
                proteinSequence=proteinSequence+line
            if MARKER in line:  #only search when ">" exist
                if (rowCnt>1):
                    listAnnoationOrder.append(accession)
                    dictProteinSequence[accession]=proteinSequence
                    dictProteinLISequence[accession]=LIChange(proteinSequence)
                accession=line[1:].replace("","")
                proteinSequence=""
            rowCnt=rowCnt+1    
        listAnnoationOrder.append(accession)
        dictProteinSequence[accession]=proteinSequence   
    return listAnnoationOrder,dictProteinSequence,dictProteinLISequence

def readTopResultFile(inputFile):
    scanInfoList=[]
    peptideList=[]
    proteinsList=[]
    emProbabilityList=[]
    orgProteinsList=[]
    rowNum=0
    scanDataPos=0
    emProbabilityPos=6
    peptidePos=14
    proteinsPos=19
    with open(inputFile,'r') as lines:
        for line in lines:
            line=line.strip()
            lineInfo=line.split("	")
            if rowNum==0:  #header
                rowNum=rowNum+1
            else:
                scanData=lineInfo[scanDataPos].replace("F","").replace(":",",")
                peptide=lineInfo[peptidePos].strip()[2:-2]
                proteins=lineInfo[proteinsPos].strip()[2:-2].replace("\"","'")
                nTerm=lineInfo[peptidePos].strip()[0]
                cTerm=lineInfo[peptidePos].strip()[-1]
                emProbability=lineInfo[emProbabilityPos]
                emProbabilityList.append(emProbability)
                peptideList.append(peptide)
                proteinsList.append(proteins.replace("\\x01",""))  #clean up the protein so it will match to the fasta file
                orgProteinsList.append(proteins)  #clean up the protein so it will match to the fasta file
                scanData=scanData+","+nTerm+","+cTerm
                scanInfoList.append(scanData)
            rowNum=rowNum+1
    return (scanInfoList,emProbabilityList,peptideList,proteinsList,orgProteinsList)

def findCorrectPeptide(scanInfoList,emProbabilityList,listAnnoationOrder,dictProteinSequence,dictProteinLISequence,peptideList,proteinsList,orgProteinsList):
    uniqueCorrectPeptideList=[]
    correctPeptidesList=[]
    sqlList=[]
    fileList=[]
    for i in range(len(peptideList)):
        peptide=peptideList[i]
        lenPeptide=len(peptide)
        proteins=proteinsList[i]
        orgProteinsInfo=orgProteinsList[i].split("', '")
        proteinInfo=proteins.replace("\\x01","")
        proteinInfo=proteins.split("', '")
        fractionID=scanInfoList[i].split(",")[0]
        scanID=scanInfoList[i].split(",")[1]
        nTerm=scanInfoList[i].split(",")[2]
        cTerm=scanInfoList[i].split(",")[3]
        emProbability=emProbabilityList[i]
        peptides=[]
        for j in range(len(proteinInfo)):
            try:
                curPos=dictProteinLISequence[proteinInfo[j]].find(peptide)
                realPeptide=dictProteinSequence[proteinInfo[j]][curPos:curPos+lenPeptide]
            except:
                realPeptide=peptide
            peptides.append(realPeptide)
            if j==0:
                uniqueCorrectPeptideList.append(realPeptide)
                usePeptide=nTerm+"."+realPeptide+"."+cTerm
            correctPeptidesList.append(peptides)
            if orgProteinsInfo[j].find("'")>0:
                orgProtein='["'+orgProteinsInfo[j]+'"]'
            else:
                orgProtein="['"+orgProteinsInfo[j]+"']"
            fileList.append((scanID,emProbability,usePeptide,orgProtein))
            sqlList.append((scanID,fractionID,emProbability,usePeptide,proteinInfo[j],"%"+proteinInfo[j]+"%"))
    return (uniqueCorrectPeptideList,correctPeptidesList,sqlList,fileList)

def sqlProcesses(conn,sqlList):
    conn.execute("DROP TABLE if exists correctedILPeptide")
    conn.commit()
    conn.execute("Create table correctedILPeptide (id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,scan int,fractionID int,EMProbability decimal(18,2), contextILCorrected varchar(255),proteins varchar(2550), likeProteins varchar(2550))")
    conn.commit()
    conn.executemany("insert into correctedILPeptide(scan,fractionID,EMProbability,contextILCorrected,proteins,likeProteins) values (?,?,?,?,?,?)",sqlList)
    conn.commit()
    conn.execute("CREATE INDEX IF NOT EXISTS correctedILPeptide_IL1 ON correctedILPeptide(scan,fractionID,likeProteins)")
    conn.execute("CREATE INDEX IF NOT EXISTS correctedILPeptide_IL2 ON correctedILPeptide(scan,fractionID,contextILCorrected)")
    conn.execute("CREATE INDEX IF NOT EXISTS result_IL1 ON result(top_result,scan,fraction_ID,proteins)")
    conn.commit()
    conn.execute("DROP VIEW IF EXISTS v_CreateCorrectPeptidesWithID")
    conn.commit()
    conn.execute("CREATE VIEW v_CreateCorrectPeptidesWithID AS select il.*, r.id as resultID from correctedILPeptide il, result r where r.top_result=1 and il.scan=r.scan and il.fractionid=r.fraction_id and r.proteins like il.likeProteins")
    conn.commit()
    conn.execute("DROP VIEW IF EXISTS v_CreateCorrectPeptides")
    conn.commit()
    conn.execute("CREATE VIEW v_CreateCorrectPeptides AS select scan,fractionID,EMProbability,contextILCorrected,group_concat(proteins,' ||| ') as groupedProteins,resultID  from v_CreateCorrectPeptidesWithID group by scan,fractionID,contextILCorrected,resultID;")
    conn.commit()
    print("done insert")

def createNewFile(newResultFile,fileList):
    #SCORE_KEY = 'Spectrum Probability Score'
    fout=open(newResultFile,'w')
    fout.write("ScanF\tEM Probability\tContext\tProteins\n")
    for j in range(len(fileList)):
        fout.write("%s\t%s\t%s\t%s\n"%(fileList[j][0],fileList[j][1],fileList[j][2],fileList[j][3]))
    fout.close()
    
if __name__ == '__main__':
    #eg: python CorrectILPeptide.py
    fastaFile=sys.argv[1]
    inputFile=sys.argv[2]
    resultsDBFile=sys.argv[3]
    newResultFile=inputFile+".new"
    #fastaFile="/home/slin/TagGraph/sampleInputFiles/FMIndices/human_uniprot_12092014_crap.fasta"
    #inputFile="/home/slin/TagGraph/sampleInputFiles/TG/mzML_output/A375MZML_TopResults.tdv"
    #resultsDBFile="/home/slin/TagGraph/sampleInputFiles/TG/mzML_output/"+"results.db"
    (listAnnoationOrder,dictProteinSequence,dictProteinLISequence)=proteinAccessionSequenceFromfastaLI(fastaFile)
    (scanInfoList,emProbabilityList,peptideList,proteinsList,orgProteinsList)=readTopResultFile(inputFile)
    (uniqueCorrectPeptideList,correctPeptideList,sqlList,fileList)=findCorrectPeptide(scanInfoList,emProbabilityList,listAnnoationOrder,dictProteinSequence,dictProteinLISequence,peptideList,proteinsList,orgProteinsList)
    conn=lite.connect(resultsDBFile)
    conn.execute("PRAGMA max_page_count=max_page;")
    conn.execute("PRAGMA temp_store=2;")
    sqlProcesses(conn,sqlList)
    conn.commit()
    createNewFile(newResultFile,fileList)