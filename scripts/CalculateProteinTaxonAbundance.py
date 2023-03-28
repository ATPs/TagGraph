#slin 201707    added lib path and others.
#slin 201806    update code for Metaproteome
#slin 201809    add data to db, enable pick the unique topresult.
#slin 20181004  update print function

import os
import sys
import pickle
import glob
import anydbm
import numpy as np
from collections import defaultdict
PAR_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
sys.path.insert(1, PAR_DIR)
LIB_DIR = PAR_DIR+"/lib"
sys.path.insert(1, LIB_DIR)
import ArgLib
import DataFile
import sqlite3 as lite
exclude_taxons = set([32630])

def perform_protein_taxon_abundance_iteration(prot_abundances, taxon_abundances, pept_abundances, peptide_prot_map, prot_taxon_map):
    prev_prot_abundances = prot_abundances
    prev_taxon_abundances = taxon_abundances
    prot_abundances = defaultdict(float)
    taxon_abundances = defaultdict(float)
    # Get nonzero contributors to each set to identify protein and taxon groups
    nonzero_contrib_p = defaultdict(set)
    for peptide in peptide_prot_map:
        sum_abundance = sum([prev_prot_abundances[p] * prev_taxon_abundances[prot_taxon_map[p]] for p in peptide_prot_map[peptide]])
        #Happens when all matching prots have abundance of zero
        #Split spectral counts proportionally based on taxon abundances or, if all taxon abundances are zero, split equally among all proteins
        if sum_abundance == 0:
            taxon_abundance = sum([prev_taxon_abundances[prot_taxon_map[p]] for p in peptide_prot_map[peptide]])
            if taxon_abundance == 0:
                num_prots = len( peptide_prot_map[peptide] )
                for p in peptide_prot_map[peptide]:
                    prot_abundances[p] += (1.0/num_prots) * pept_abundances[peptide]
                    nonzero_contrib_p[p].add(peptide)
            else:
                for p in peptide_prot_map[peptide]:
                    prot_abundances[p] += (prev_taxon_abundances[prot_taxon_map[p]]/taxon_abundance) * pept_abundances[peptide]
                    if prev_taxon_abundances[prot_taxon_map[p]] > 0:
                        nonzero_contrib_p[p].add(peptide)
        else:
            for p in peptide_prot_map[peptide]:
                prot_abundances[p] += (prev_prot_abundances[p] * prev_taxon_abundances[prot_taxon_map[p]] / sum_abundance) * pept_abundances[peptide]
                if prev_prot_abundances[p] * prev_taxon_abundances[prot_taxon_map[p]] > 0:
                    nonzero_contrib_p[p].add(peptide)
    #for p in prot_abundances:
        #prot_abundances[p] = prot_abundances[p] / int( prot_lengths[p] )
    sum_nsaf = sum( prot_abundances.values())
    #print('Sum NSAF',sum_nsaf)
    for p in prot_abundances:
        prot_abundances[p] = prot_abundances[p] / sum_nsaf
    for t in taxon_prot_map:
        for p in taxon_prot_map[t]:
            taxon_abundances[t] += prot_abundances[p]
    dist = 0
    for p in prot_abundances:
        dist += (prot_abundances[p] - prev_prot_abundances[p])**2
    for t in taxon_abundances:
        dist += (taxon_abundances[t] - prev_taxon_abundances[t])**2
    dist = np.sqrt(dist)
    return prot_abundances, taxon_abundances, nonzero_contrib_p, dist

def calculate_initial_maps(tg_files, score_cutoff):
    SCAN_KEY = 'ScanF'
    PEPTIDE_KEY = 'Context'
    SEQ_DELIM = 2
    PROTEIN_KEY = 'Proteins'
    SCORE_KEY = 'EM Probability'
    #SCORE_KEY = 'Spectrum Probability Score'
    pept_scan_map = defaultdict(set)
    prot_peptide_map = defaultdict(set)
    taxon_peptide_map = defaultdict(set)
    peptide_taxon_map = defaultdict(set)
    peptide_prot_map = defaultdict(set)
    prot_taxon_map = {}
    taxon_prot_map = defaultdict(set)
    for tg_file in tg_files:
        for item in DataFile.getScanInfoIterator(tg_file, delimiter='\t'):
            prots = eval(item[PROTEIN_KEY])
            #taxon = int(item['Taxon'])
            # If using topresults file
            taxon = int( prots[0].split('TAXON=')[1] )
            # Don't record taxon data if it is in the set of excluded taxons
            if taxon in exclude_taxons or float(item[SCORE_KEY]) < score_cutoff:
                continue
            peptide = item[PEPTIDE_KEY][SEQ_DELIM:-SEQ_DELIM]
            pept_scan_map[peptide].add( (tg_file, item[SCAN_KEY]) )
            taxon_peptide_map[taxon].add( peptide )
            peptide_taxon_map[peptide].add( taxon )
            for prot in prots:
                prot_peptide_map[prot].add( peptide )
                peptide_prot_map[peptide].add( prot )
                prot_taxon_map[prot] = taxon
                taxon_prot_map[taxon].add( prot )
    return pept_scan_map, prot_peptide_map, peptide_prot_map, taxon_peptide_map, peptide_taxon_map, prot_taxon_map, taxon_prot_map

def initialize_abundances(peptide_prot_map, peptide_taxon_map, pept_scan_map):
    # Initialize protein and taxon abundances
    prot_abundances = defaultdict(float)
    taxon_abundances = defaultdict(float)
    pept_abundances = defaultdict(float)
    num_uniques_t, num_uniques_p = defaultdict(int), defaultdict(int)
    # Maybe change this to initialize based on total peptides matching count, not just peptides which uniquely match. Test it out
    for peptide in peptide_prot_map:
        pept_abundance = len(pept_scan_map[peptide])
        if len(peptide_taxon_map[peptide])==1:
            # sets don't support indexing
            for t in peptide_taxon_map[peptide]:
                taxon_abundances[t] += pept_abundance
                num_uniques_t[t] += 1
        if len(peptide_prot_map[peptide])==1:
            for p in peptide_prot_map[peptide]:
                prot_abundances[p] += pept_abundance
                num_uniques_p[p] += 1
        pept_abundances[peptide] = pept_abundance
    return prot_abundances, taxon_abundances, pept_abundances, num_uniques_t, num_uniques_p

def createTables(conn):
    conn.execute("DROP TABLE if exists taxonsGroup")
    conn.execute("DROP TABLE if exists taxonsProtein")
    conn.execute("DROP TABLE if exists taxonsInfo")
    conn.commit()
    conn.execute("CREATE TABLE taxonsGroup(taxonsGroupID int AUTO_INCREMENT PRIMARY KEY,taxonsProteinID int,contributingPeptide varchar(255))")
    conn.execute("CREATE TABLE taxonsProtein(taxonsProteinID int AUTO_INCREMENT PRIMARY KEY,proteins varchar(2550),taxonInfo varchar(2550),coveringSetSize int,combinedNSAF decimal(14,10))")
    conn.execute("CREATE TABLE taxonsInfo(id int AUTO_INCREMENT PRIMARY KEY,taxon varchar(255),name varchar(255),uniqueMatchingPeptides int,numMatchingPeptides int,abundance decimal(14,10))")
    conn.commit()

def doSQL(conn):
    conn.execute("DROP INDEX IF EXISTS taxonsGroup_id")
    conn.execute("DROP INDEX IF EXISTS taxonsProtein_id")
    conn.execute("DROP INDEX IF EXISTS results_multi1")
    conn.execute("DROP INDEX IF EXISTS results_multi2")
    conn.commit()
    conn.execute("CREATE INDEX taxonsGroup_id ON taxonsGroup(taxonsProteinID)")
    conn.execute("CREATE INDEX taxonsProtein_id ON taxonsProtein(taxonsProteinID)")
    conn.execute("CREATE INDEX results_multi1 ON result(top_result,proteins,context)")
    conn.execute("CREATE INDEX results_multi2 ON result(top_result,scan,fraction_id)")
    conn.commit()
    conn.execute("DROP VIEW IF EXISTS v_proteinsTaxionsPeptides;")
    conn.execute("DROP VIEW IF EXISTS v_topResultNSAF;")
    conn.execute("DROP VIEW IF EXISTS v_uniqueScan;")
    conn.execute("DROP VIEW IF EXISTS v_duplicateScans;")
    conn.execute("DROP VIEW IF EXISTS v_duplicateScansCoveringSetSizeMax")
    conn.execute("DROP VIEW IF EXISTS v_duplicateScansUnique")
    conn.execute("DROP VIEW IF EXISTS v_duplicateScansTie")
    conn.execute("DROP VIEW IF EXISTS v_results;")
    conn.execute("DROP VIEW IF EXISTS v_minIDTopResult;")
    conn.execute("DROP VIEW IF EXISTS v_scanFractionID;")
    conn.execute("DROP VIEW IF EXISTS v_topResults;") 
    conn.commit()
    conn.execute("CREATE VIEW v_proteinsTaxionsPeptides AS select tp.proteins,\"%\"||tp.taxonInfo||\"%\" as taxonInfoLike,tp.taxonInfo,tp.combinedNSAF,tg.contributingPeptide as peptide, replace(tg.contributingPeptide,\"L\",\"I\") as LIPeptide, tp.coveringSetSize as coveringSetSize from taxonsGroup tg, taxonsProtein tp where tp.taxonsProteinID=tg.taxonsProteinID;")
    conn.execute("CREATE VIEW v_topResultNSAF AS select r.*,vr.coveringSetSize,vr.combinedNSAF from result r left outer join v_proteinsTaxionsPeptides vr on (substr(r.context,3,length(r.context)-4)=LIPeptide) and (r.proteins like taxonInfoLike) where r.top_result=1;")
    conn.execute("CREATE VIEW v_uniqueScan AS select * from v_topResultNSAF where top_result=1 group by scan,fraction_id having count(*)=1;")
    conn.execute("CREATE VIEW v_duplicateScans AS select v.* from v_topResultNSAF v, (select scan, fraction_id from v_topResultNSAF group by fraction_id,scan having count(*)>1) q  where v.scan=q.scan and v.fraction_id=q.fraction_id")
    conn.execute("CREATE VIEW v_duplicateScansCoveringSetSizeMax AS select scan,fraction_id,max(coveringSetSize) as coveringSetSize  from v_duplicateScans group by scan,fraction_id")
    conn.execute("CREATE VIEW v_duplicateScansUnique AS select v.scan,v.fraction_id,v.coveringSetSize from v_duplicateScans v, v_duplicateScansCoveringSetSizeMax d where d.scan=v.scan and d.fraction_id=v.fraction_id and v.coveringSetSize=d.coveringSetSize group by v.scan,v.fraction_id having count(*)=1")
    conn.execute("CREATE VIEW v_duplicateScansTie AS select v.scan,v.fraction_id,v.coveringSetSize from v_duplicateScans v, v_duplicateScansCoveringSetSizeMax d where d.scan=v.scan and d.fraction_id=v.fraction_id and v.coveringSetSize=d.coveringSetSize group by v.scan,v.fraction_id having count(*)>1")
    conn.execute("CREATE VIEW v_results AS select * from v_uniqueScan union select v.* from v_duplicateScans v, v_duplicateScansUnique u where v.scan=u.scan and v.fraction_id=u.fraction_id and v.coveringSetSize=u.coveringSetSize union select v.* from v_duplicateScans v, (select v.scan, v.fraction_id,v.coveringSetSize,max(v.combinedNSAF) as combinedNSAF from v_duplicateScans v, v_duplicateScansTie u where v.scan=u.scan and v.fraction_id=u.fraction_id and v.coveringSetSize=u.coveringSetSize group by  v.scan, v.fraction_id,v.coveringSetSize) u where v.scan=u.scan and v.fraction_id=u.fraction_id and v.coveringSetSize=u.coveringSetSize and v.combinedNSAF=u.combinedNSAF")
    conn.execute("CREATE VIEW v_minIDTopResult as select min(id) as minID from v_results group by scan,fraction_id;")
    conn.execute("CREATE VIEW v_scanFractionID as select scan,fraction_id from v_results group by scan,fraction_id;")
    conn.execute("CREATE VIEW v_topResults AS select * from v_results where id in (select min(id) from v_results group by scan,fraction_id) union select v.* from v_topResultNSAF v left join (select scan, fraction_id from v_results group by scan, fraction_id) p on v.scan =p.scan and v.fraction_id=p.fraction_id where p.scan is null and p.fraction_id is null")
    conn.commit()

if __name__ == '__main__':
    print('this program calculates all protein and taxon abundances based on the output of taggraph')
    print('taggraph argument is a list of files to include in the protein abundance calculation, separated by commas')
    print('alpha parameter is cutoff on score (EM Probability) to consider scan in abundance calculation')
    print('model parameter is dbm file containing protein lengths')
    print('taxonomy parameter is genus counts file')
    #python scripts/CalculateProteinTaxonAbundance.py -y "/home/slin/TagGraph/taxonomy/genus_info.pck" -m " " -r "0.99" -T "/home/slin/TagGraph/sampleInputFiles/TG/OutMetaProteome_f1/MetaProteome_f1_TopResults.tdv" -o "/home/slin/TagGraph/sampleInputFiles/TG/OutMetaProteome_f1/MetaProteome_f1_TaxonAbundance"
    options = ArgLib.parse(['taggraph', 'output', 'alpha', 'model', 'taxonomy'])
    #prot_lengths = anydbm.open( options.model )
    CUTOFF = 1e-4
    pept_scan_map, prot_peptide_map, peptide_prot_map, taxon_peptide_map, peptide_taxon_map, prot_taxon_map, taxon_prot_map = calculate_initial_maps( options.taggraph.split(','), options.alpha )
    prot_abundances, taxon_abundances, pept_abundances, num_uniques_t, num_uniques_p = initialize_abundances(peptide_prot_map, peptide_taxon_map, pept_scan_map)
    prot_lengths = {}
    dist = 10
    iter = 0
    while dist > CUTOFF:
        prot_abundances, taxon_abundances, nonzero_contrib_p, dist = perform_protein_taxon_abundance_iteration(prot_abundances, taxon_abundances, pept_abundances, peptide_prot_map, prot_taxon_map)
        iter += 1
        print('Iteration %i, Distance %f'%(iter, dist))
    outputFolder=options.output[:options.output.rfind("/")]
    resultsDBFile=outputFolder+"/"+"results.db"
    conn=lite.connect(resultsDBFile)
    conn.execute("PRAGMA max_page_count=max_page;")
    conn.execute("PRAGMA temp_store=2;")
    createTables(conn)
    taxonsGroupList=[]
    taxonsProteinList=[]
    taxonsProteinID=1
    taxonsGroupID=1
    cols = ['Protein', 'Taxon', 'Unique Matching Peptides', 'Num Matching Peptides', 'NSAF']
    out_file = open(options.output + '.proteins', 'w')
    out_file.write('\t'.join(cols) + '\n')
    for prot in sorted(prot_abundances, key = lambda k: -prot_abundances[k]):
        if prot_abundances[prot] > 0:
            out_file.write('\t'.join([prot, str(prot_taxon_map[prot]), str(num_uniques_p[prot]), str(len(prot_peptide_map[prot])), str(prot_abundances[prot])]) + '\n')
    out_file.close()
    #Combine peptides together if they have the same abundance and set of nonzero contributors
    protein_groups = defaultdict(set)
    for p in nonzero_contrib_p:
        protein_groups[tuple(sorted(nonzero_contrib_p[p])), prot_abundances[p]].add(p)
    cols = ['Proteins', 'Covering Set Size', 'Combined NSAF', 'Contributing Peptides']
    out_file = open(options.output + '.proteins.groups', 'w')
    out_file.write('\t'.join(cols) + '\n')
    for key in sorted(protein_groups, key = lambda k: len(k[0]) * -k[1]):
        out_file.write('\t'.join([str(sorted(tuple(protein_groups[key]))),str(len(key[0])),str(len(protein_groups[key])*key[1]),str(key[0])])+'\n')
        for i in range(0,len(key[0])):
            taxonsGroupList.append((taxonsGroupID,taxonsProteinID,key[0][i]))
            taxonsGroupID=taxonsGroupID+1
        proteins=str(sorted(tuple(protein_groups[key])))
        taxonInfo=""
        taxonInfo=str(proteins[proteins.rfind("TAXON=",1)+1:len(proteins)-2])
        taxonsProteinList.append((taxonsProteinID,proteins,taxonInfo,str(len(key[0])),str(len(protein_groups[key])*key[1])))
        taxonsProteinID=taxonsProteinID+1
    out_file.close()
    conn.executemany("insert into taxonsProtein(taxonsProteinID,proteins,taxonInfo,coveringSetSize,combinedNSAF) values(?,?,?,?,?)",taxonsProteinList)
    conn.executemany("insert into taxonsGroup(taxonsGroupID,taxonsProteinID,contributingPeptide) values(?,?,?)",taxonsGroupList)
    conn.commit()
    with open(options.taxonomy) as fin:
        genus_counts = pickle.load(fin)
    cols = ['Taxon', 'Name', 'Unique Matching Peptides', 'Num Matching Peptides', 'Abundance']
    out_file = open(options.output + '.taxons', 'w')
    out_file.write('\t'.join(cols) + '\n')
    counter=1
    taxonsInfoList=[]
    for taxon in sorted(taxon_abundances, key = lambda k: -taxon_abundances[k]):
        if taxon_abundances[taxon] > 0:
            out_file.write('\t'.join([str(taxon), genus_counts[taxon]['name'], str(num_uniques_t[taxon]), str(len(taxon_peptide_map[taxon])), str(taxon_abundances[taxon])]) + '\n')
            taxonsInfoList.append((counter,str(taxon),genus_counts[taxon]['name'],str(num_uniques_t[taxon]),str(len(taxon_peptide_map[taxon])),str(taxon_abundances[taxon])))
            counter=counter+1
    conn.executemany("insert into taxonsInfo(id,taxon,name,uniqueMatchingPeptides,numMatchingPeptides,abundance) values(?,?,?,?,?,?)",taxonsInfoList)
    conn.commit()
    out_file.close()
    doSQL(conn)
    # TODO: Write code to output which taxons are represented by the same set of peptides?
