# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 12:14:03 2016
@author: ellisrj2

Every now and then using Bioconductor or GEO2R you'll have probes with no corresponding genes. This script finds the genes corresponding to blank probes with the use of probe reference files
"""



def GEO_find_blank_probes(GEO_file, probe_file):
    import pandas
    #probe_file should equal 'EBI_Illumina_probe_data.csv'
    indices = []
    probes_to_find = []
    
#    ebi_genes = []
#    ebi_indices = []
    
    ebi_genes = []
    ilmn_genes = []
    agilent_genes = []
    affy_genes = []

    #open GEO file
    df = pandas.read_excel(GEO_file, skipinitialspace=True)
    #df = df[df['P.Value'] < 0.05] #retain entries where p<0.05 (CAN ALSO USE 'adj.P.Val')
    
    probe_ids = list(df['probe']) #list of probe codes
    gene_symbols = list(df['gene.symbols72517']) #list of gene symbols
    
    #gather list of indices of probes with no corresponding gene symbols
    for i in range(len(gene_symbols)):
        if type(gene_symbols[i]) == float:
            indices.append(i)
            probes_to_find.append(probe_ids[i])
    
    #turn list elements into strings       
    probes_to_find = [str(i) for i in probes_to_find]
    
    #for probes that are all numbers
    if 'Illumina' in probe_file:

                    
        df_EBI2 = pandas.read_csv('EBI_Illumina_probe_data2.csv', 
                                  skipinitialspace=True)
        #list of probes       
        ebi_probe_list2 = list(df_EBI2['Array Design Name'])
        
        #list of genes
        ebi_gene_list2 = list(df_EBI2['Unnamed: 5'])
        
        #gather list of genes for corrresponding probes with unknown genes in GEO file
        for probe in probes_to_find:
            if probe in ebi_probe_list2:
                probe_index = ebi_probe_list2.index(probe)
                ebi_genes.append(ebi_gene_list2[probe_index])
        
        for i in range(len(ebi_genes)):
            df['Gene.symbol'][indices[i]] = ebi_genes[i]
        
        df.to_excel(GEO_file + '_addedgenes.xlsx', index=False)
        
    #for probes that begin with "ILMN"
    elif "ILMN" in probe_file:
        df_ILMN = pandas.read_csv(probe_file, skipinitialspace=True)
        
        #list of probes
        ilmn_probe_list = list(df_ILMN['Probe Set ID'])
        
        #list of genes
        ilmn_gene_list = list(df_ILMN['Gene Symbol'])
        
        #gather list of genes for corrresponding probes with unknown genes in GEO file
        for probe in probes_to_find:
            if probe in ilmn_probe_list:
                probe_index = ilmn_probe_list.index(probe)
                ilmn_genes.append(ilmn_gene_list[probe_index])
                
        for i in range(len(ilmn_genes)):
            df['Gene.symbol'][indices[i]] = ilmn_genes[i]
        
        df.to_excel(GEO_file + '_addedgenes.xlsx', index=False)    
       
    elif "GPL2872" in probe_file:
        df_agilent = pandas.read_csv(probe_file, skipinitialspace=True)
        
        #list of probes
        agilent_probes = list(df_agilent['ID'])
        
        #list of genes
        agilent_gene_list = list(df_agilent['GENE_SYMBOL'])
        
        #gather list of genes for corrresponding probes with unknown genes in GEO file
        for probe in probes_to_find:
            if probe in agilent_probes:
                probe_index = agilent_probes.index(probe)
                agilent_genes.append(agilent_gene_list[probe_index])
                
        for i in range(len(agilent_genes)):
            df['Gene.symbol'][indices[i]] = agilent_genes[i]
        
        df.to_excel(GEO_file + '_addedgenes.xlsx', index=False) 
    elif "GPL1261" in probe_file:
        df_gpl1261 = pandas.read_excel(probe_file)
        
        #list of probes
        affy_probes = list(df_gpl1261['#ID = Affymetrix Probe Set ID'])
        
        #list of genes
        affy_gene_list = list(df_gpl1261['Unnamed: 10'])
        
        #gather list of genes for corrresponding probes with unknown genes in GEO file
        for probe in probes_to_find:
            if probe in affy_probes:
                probe_index = affy_probes.index(probe)
                affy_genes.append(affy_gene_list[probe_index])
                
        for i in range(len(affy_genes)):
            df['gene.symbols72517'][indices[i]] = affy_genes[i]
        
        df.to_excel(GEO_file + '_addedgenes1.xlsx', index=False)