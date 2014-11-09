# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Libraries

# <codecell>

import pandas as pd
import mygene
import numpy as np
import os
import myvariant
import sys
import time
import gc

# <headingcell level=1>

# MyVariant Search

# <codecell>

reload(myvariant)
mv = myvariant.MyVariantInfo()
#%timeit a = mv.query_variant("chr2:179392080-179669337", size=2000)

# <headingcell level=1>

# Gene List - CDS START/END

# <codecell>

mg = mygene.MyGeneInfo()

# <codecell>

#mg.getgenes([7273],'name,symbol,refseq.rna', as_dataframe=True)
#mg.getgene(1017)

# <codecell>

#gene_ids = [mg.query(q='type_of_gene:protein-coding', fields='entrezgene', species='human', skip=i, size=i+1000, as_dataframe=True) for i in range(0,21000,1000)]     
# gene_ids_df = pd.concat(gene_ids)
# gene_ids_df.index = range(0,len(gene_ids_df))
# gene_ids_df.drop_duplicates(inplace=True)
#gene_ids_df = pd.read_table('gene_ids.txt',sep='\t')
#print len(gene_ids_df)

# <headingcell level=1>

# Gene Specific Functional Impact Score Function

# <headingcell level=3>

# Get CDS exon range

# <codecell>


def get_cds_coords(exons, cds_start, cds_end):
    '''
    This function returns coding exon ranges when given a 
    exon coordinate list  [[start1, end1], [start2, end2]],
    the cds_start position and cds_end position 
    '''
    
    def start_change(exon_start_end, cds_start):
        exon_start = int(exon_start_end[0])
        exon_end = int(exon_start_end[1])
        if cds_start > exon_end:  #discard 5' UTR exon
            return 0
        if cds_start >= exon_start and cds_start <= exon_end: #cds start within exon
            return cds_start
        return exon_start  #keep exon start
    
    
    def end_change(exon_start_end, cds_end):
        exon_start = int(exon_start_end[0])
        exon_end = int(exon_start_end[1])
        if cds_end <= exon_start: #discard 3' UTR exon
            return 0
        if cds_end >= exon_start and cds_end <= exon_end:  #cds end within exon
            return cds_end
        return exon_end  #keep exon end
    
    exon_coords = pd.DataFrame.from_records(exons)
    exon_coords[0] = exon_coords.apply(start_change, args=[cds_start], axis=1)
    exon_coords[1] = exon_coords.apply(end_change, args=[cds_end], axis=1)
    exon_coords = exon_coords[(exon_coords[0]!=0) &(exon_coords[1]!=0)]
    return exon_coords

# <headingcell level=3>

# Get VCF Exon df

# <codecell>

class Vcf_metadata(object):
    
    def __init__(self, filename):
        if filename.endswith('.gz'):
            self.compression = 'gzip'
            if filename+'.tbi' in os.listdir(os.path.split(filename)[0]):
                header_lines = os.popen('tabix -H ' + filename).readlines()
                self.header = [l.replace('#CHROM','CHROM') for l in header_lines if l.startswith('#')]
            os.system('tabix -p vcf ' + filename)
            header_lines = os.popen('tabix -H ' + filename).readlines()
            self.header = [l.replace('#CHROM','CHROM') for l in header_lines if l.startswith('#')]
        
        else:
            self.compression = ''
            header_lines = os.popen('head -5000 ' + filename).readlines()
            self.header = [l.replace('#CHROM','CHROM') for l in header_lines if l.startswith('#')]

# <codecell>

class vcf_tabix(object):
    '''
    Loads in a vcf file, aware of gzipped files.
    '''
    
    def __init__(self, filename, chrom, start, end):
        
        #Header
        header_parsed = Vcf_metadata(filename)
        self.header = self.get_header_df(header_parsed.header)  #header parsed into key/values dataframe
    
        #Vcf df indexed on CHROM, POS, REF, ALT
        
        #self.df = pd.read_table(filename, sep="\t", compression=header_parsed.compression, skiprows=(len(self.header)-1))
        self.df = os.popen('tabix ' + vcf_path + ' '+str(chrom)+':' + str(start) + '-' + str(end)).readlines()
        if len(self.df) == 0: return None
        
        self.df = [i.rstrip('\n').split('\t') for i in self.df]
        self.df = pd.DataFrame.from_records(self.df)
        #print header_parsed.header[-1].rstrip('\n').split("\t")
        self.df.columns = header_parsed.header[-1].rstrip('\n').split("\t")[:-1:]
        self.df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True, drop=False)
        self.df_bytes = self.df.values.nbytes + self.df.index.nbytes + self.df.columns.nbytes
        
        #Sample Ids
        self.samples = self.df.columns[8:]
        
        #Format values
        self.FORMAT = self.df[self.df.columns[8]]
        
        #Create Samples df indexed on CHROM, POS, REF, ALT
        split_cols_dict = {}#{'AD':2}#, 'PL':3}
        self.get_vcf_samples_df(split_cols_dict)
        self.sample_df_bytes = self.sample_df.values.nbytes + self.sample_df.index.nbytes + self.sample_df.columns.nbytes
        
    
    
    def get_header_df(self, header_txt):
        '''
        Parses header into pandas DataFrame
        '''
        key_value_header = [i.replace('##','').replace('\n','').split('=',1) for i in header_txt if '##' in i]
        key_value_header.append(['SampleIDs',header_txt[-1].rstrip('\n').split('\t')[9:]])
        header_df =  pd.DataFrame.from_records(key_value_header)
        header_df.set_index(0,inplace=True)
        header_df.index.name = 'header_keys'
        header_df.columns = ['header_values']
        return header_df
    
    
    
    def get_vcf_samples_df(self, split_cols):
        '''
        This function creates a samples_df containing sample
        level information for all non-missing variant calls
        indexed on CHROM, POS, REF, ALT
        
        vcf_df
        
        split_cols
        
        '''
        genotypes = self.df[self.samples].groupby(by=self.FORMAT)  #genotypes grouped by FORMAT variant annotations
        
        #Iterate through genotype groups, dropping missing calls
        master_df = []
        for name,group in genotypes:
            temp_group = genotypes.get_group(name)  #group of interest
            del temp_group['FORMAT']  #remove the format column
            temp_group.replace(to_replace='.', value=nan, inplace=True)  #replace . with none, allows stack to remove null columns, space savings
            
            temp_group = temp_group.stack()  #stacking samples for each variant
        
            #creating sample dataframe
            temp_group_data = pd.DataFrame(temp_group.str.split(':').tolist())
            temp_group_data.index = temp_group.index
            temp_group_data.columns = name.split(':')
            temp_group_data.replace(to_replace='.', value=nan, inplace=True)
            
            master_df.append(temp_group_data)
        
        #Concatenating all genotype groups
        self.sample_df = pd.concat(master_df)
        self.sample_df.index.names = ['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID']
        
        #spliting user-defined columns
        for col in split_cols:
            for i in range(0, split_cols[col]):
                self.sample_df[col + '_' + str(i)] = self.sample_df[col].str.split(',').str[i]
            del self.sample_df[col]
                
        return 0

# <codecell>

def _get_allele(self, line, gt_col):
    '''
    Returns allelic base, handles multi-allelic variants
    '''
    alleles = [line['REF']]
    alleles.extend(line['ALT'].split(","))
    a1 = "."
    try:
        a1 = alleles[int(line[gt_col])]  #returns missing if gt_int_call is "."
    except:
        a1 = "."
    return a1
        
        
def _get_GT_multisample_vcf(self, line, sample_col, gt_index):
    '''
    Slow parser for multisample vcf
    '''
    return int( line[sample_col].split(line['phase'])[int(gt_index)])

# <codecell>

def _get_allele_bases(df, sample_col, single_sample_vcf=True):
    '''
    Adds phase, GT1, GT2, a1, a2 to self._variants_df dataframe
    
    phase : { /, | }, unphased or phased call
    GT1: int, first allele call in the numeric sample genotype column
    GT2: int, second allele call in the numeric sample genotype column
    a1: {A, T, G, C, AA, etc}, nucleotide base representation for GT1
    a2: {A, T, G, C, AA, etc}, nucleotide base representation for GT2
    
    '''
    
    
    if single_sample_vcf:
        df['phase'] = df[sample_col].str[1]
        df['GT1'] = df[sample_col].str[0]
        df['GT1'] = df['GT1'].astype(int)
        df['GT2'] = df[sample_col].str[2]
        df['GT2'] = df['GT2'].astype(int)
    
    
    if not single_sample_vcf:
        df['phase'] = df.apply(get_phase, args=['GT'], axis=1)  #get phase
        df = df[df.phase != "-"]  #likley occurs at sex chromosome sites
        df['GT1'] = df.apply(_get_GT_multisample, args=[sample_col, 0], axis=1)
        df['GT2'] = df.apply(_get_GT_multisample, args=[sample_col, 1], axis=1)
        
        
    
    #SLOW PROCESS MULTIPLE ALLELE GENOTYPES
    df_multi = df[(df.GT1.astype(int)>1) | (df.GT2.astype(int)>1)] #select all multi-alleleic variants
    if len(df_multi) > 0:
        df_multi['a1'] = df_multi.apply(_get_allele, args=['GT1'], axis=1)  #
        df_multi['a2'] = df_multi.apply(_get_allele, args=['GT2'], axis=1)
    
    
    #FAST PROCESS SIMPLE ALLELE GENOTYPES
    df_simple = df[~df.index.isin(df_multi.index)][['REF', 'ALT', 'GT1', 'GT2']]  #dropping multiallele variants, minimize memory usage
    
    df_gt1_ref = df_simple[df_simple.GT1==0][['REF']]  #get a1 ref alleles
    df_gt1_ref.columns = ['a1']
    df_gt2_ref = df_simple[df_simple.GT2==0][['REF']]  #get a2 ref alleles
    df_gt2_ref.columns = ['a2']
    
    
    df_gt1_alt = df_simple[df_simple.GT1==1][['ALT']]  #get a1 alt alleles
    df_gt1_alt.columns = ['a1']
    df_gt2_alt = df_simple[df_simple.GT2==1][['ALT']]  #get a2 alt alleles
    df_gt2_alt.columns = ['a2']
    
    
    gt1_alleles = pd.concat([df_gt1_ref,df_gt1_alt])  #merging GT1 allele bases into a single df

    gt2_alleles = pd.concat([df_gt2_ref,df_gt2_alt])  #merging GT2 allele bases into a single df

    gt1_2_allele_df = gt1_alleles.join(gt2_alleles, how='outer')  #Joining the GT1 and GT2 simple allele bases 
    
    
    df = df.join(gt1_2_allele_df, how='inner')  #Adding simle allele a1 and a2 columns to original df
    df = df.append(df_multi)  #Adding multi-alleleic bases to original df
    
    return df

# <codecell>

#print chrom
#exon_coords.head()

# <codecell>

#exon_coords.tail()

# <codecell>

def get_fi_score(chrom, start, end, score_id1, score_id2):
    '''
    Input mv.query_variant['hits']['hits']
    '''
    assert 'chr' in chrom, 'chr not in chrom'
    start = str(int(start))
    end = str(int(end))
    
    my_var = mv.query_variant( chrom + ":" + start + "-" + end , size=2000)  #  "chr1:100000-1000000")
    hits = my_var['hits']['hits']
    

    fi_scores = []
    counter = 0
    for i in hits:
        try:
            fi_scores.append([i['_source']['dbnsfp'][score_id1][score_id2], i['_source']['_id']])
        except:
            counter +=1
            continue
    
    #print counter
    if len(fi_scores) == 0: return []
    fi_scores = pd.DataFrame.from_records(fi_scores)
    #fi_scores = pd.DataFrame.from_records( [[i['_source']['dbnsfp'][score_id1][score_id2], i['_source']['_id']]  for i in hits] )
    fi_scores.columns = ['score', 'id'] 
    fi_scores['POS'] = fi_scores['id'].str.split('.').str[-1].str.split('>').str[0].str[:-1]
    fi_scores['REF'] = fi_scores['id'].str.split('.').str[-1].str.split('>').str[0].str[-1]
    fi_scores['ALT'] = fi_scores['id'].str.split('.').str[-1].str.split('>').str[-1]
    fi_scores['POS_ALT'] = fi_scores['POS'] + '_' + fi_scores['ALT']
    fi_scores.set_index('POS_ALT', inplace=True)
    fi_scores.sort(columns=['POS'], inplace=True)
    return fi_scores
    

    

# <codecell>



def get_gene_df(gene_id):

    print gene_id
    
    def arrayify(gene_df, gene_id):
        
        gene_df = gene_df['a1_a2_score_sum'].unstack(level=4)
        gene_df['gene_id'] = gene_id
        return gene_df
    
    
    vcf_path = '/Users/erickscott/datasets_raw/snyderome/TB0001907.all.ILLUMNIA.bwa.CEU.high_coverage.20101118.snp.raw.filtered.vcf.gz'
    

    json_txt = mg.query(gene_id, fields='exons_hg19',  as_dataframe=True)
    
    try:
        temp_gene = json_txt.exons_hg19
    except AttributeError:
        return 0  #no hg19 exon information
    
    gene_df = []
    
    gene_coords_tracker = set()
    for tx in json_txt.exons_hg19.values[0].keys():
        
       
        if type(json_txt.exons_hg19.values[0][tx]) == list:
            gene_chrom = json_txt.exons_hg19.values[0][tx][0]['chr']
            if 'chr' not in gene_chrom:
                gene_chrom = 'chr' + gene_chrom
            gene_start = json_txt.exons_hg19.values[0][tx][0]['cdsstart']
            gene_end = json_txt.exons_hg19.values[0][tx][0]['cdsend']
            gene_exons = json_txt.exons_hg19.values[0][tx][0]['exons']
            
        else:
            gene_chrom = json_txt.exons_hg19.values[0][tx]['chr']
            if 'chr' not in gene_chrom:
                gene_chrom = 'chr' + gene_chrom
            gene_start = json_txt.exons_hg19.values[0][tx]['cdsstart']
            gene_end = json_txt.exons_hg19.values[0][tx]['cdsend']
            gene_exons = json_txt.exons_hg19.values[0][tx]['exons']
      
    
        
        exon_coords = get_cds_coords(gene_exons, gene_start, gene_end)
        
        
        for exon_idx in exon_coords.index:  #iterate through gene exons
            
            exon_start = exon_coords.ix[exon_idx][0]
            exon_end = exon_coords.ix[exon_idx][1]
            
            if (exon_start,exon_end) in gene_coords_tracker: continue
            
            #print exon_start, exon_end
            
            #try:
            #vcf_path = '/Users/erickscott/git/fiSSEA/ipynb/SWGR_titin.vcf.gz'
            vcf_df_tabix = vcf_tabix(vcf_path, gene_chrom, exon_start, exon_end)
            if len(vcf_df_tabix.df) == 0: 
                gene_coords_tracker.add((exon_start,exon_end))
                continue  #occurs if the samples have no variants in that exon
            
            
            #FORMAT VCF FOR TRANSLATION
            vcf_df_tabix_sample_df = vcf_df_tabix.sample_df.copy()
            #print len(vcf_df_tabix_sample_df)
            del vcf_df_tabix
            gc.collect()
            
            
            vcf_df_tabix_sample_df.reset_index(inplace=True)
            vcf_df_tabix_sample_df = vcf_df_tabix_sample_df[vcf_df_tabix_sample_df.GT != './.']
            vcf_df_tabix_sample_df = vcf_df_tabix_sample_df[ (vcf_df_tabix_sample_df.ALT.isin(['A', 'T', 'C', 'G'])) & (vcf_df_tabix_sample_df.REF.isin(['A', 'T', 'C', 'G']))]     
            vcf_df_tabix_sample_df.fillna(value='Null', inplace=True)
            vcf_df_tabix_sample_df = vcf_df_tabix_sample_df[vcf_df_tabix_sample_df['GT']!= 'Null']
            vcf_df_tabix_sample_df.set_index(['CHROM', 'POS', 'REF', 'ALT', 'SAMPLE_ID'],inplace=True, drop=False)
            
            
            
            vcf_df_tabix_sample_df = _get_allele_bases(vcf_df_tabix_sample_df, 'GT')
            
             
            vcf_df_tabix_sample_df['pos_a1'] = vcf_df_tabix_sample_df['POS'] + '_' + vcf_df_tabix_sample_df['a1']
            vcf_df_tabix_sample_df['pos_a2'] = vcf_df_tabix_sample_df['POS'] + '_' + vcf_df_tabix_sample_df['a2']
            
            try:
                chrom = vcf_df_tabix_sample_df.index.get_level_values(0).unique()[0]
                start = vcf_df_tabix_sample_df.head(1).index.get_level_values(1).unique()[0]
                end = vcf_df_tabix_sample_df.tail(1).index.get_level_values(1).unique()[0]
                
                cadd = get_fi_score('chr'+chrom, start, end, 'cadd', 'raw')
                if type(cadd) == list: 
                    gene_coords_tracker.add((exon_start,exon_end))
                    continue  #No cadd scores for these variants
                
                cadd['pos_a1'] = cadd['score']
                cadd['pos_a2'] = cadd['score']
                cadd_score_dict = cadd[['pos_a1', 'pos_a2']].to_dict()
                #cadd_score_dict_a2 = cadd[['pos_a2']].to_dict()
                
                vcf_df_tabix_sample_df = vcf_df_tabix_sample_df[['pos_a1','pos_a2']]
                
                vcf_df_tabix_sample_df = vcf_df_tabix_sample_df[ (vcf_df_tabix_sample_df['pos_a1'].isin(cadd_score_dict['pos_a1'].keys())) \
                                                                  | (vcf_df_tabix_sample_df['pos_a2'].isin(cadd_score_dict['pos_a2'].keys()))  ]
                
                vcf_df_tabix_sample_df = vcf_df_tabix_sample_df.replace(cadd_score_dict)
                
                vcf_df_tabix_sample_df['pos_a1'] = vcf_df_tabix_sample_df['pos_a1'].map(lambda x: 0 if '_' in str(x) else x)
                vcf_df_tabix_sample_df['pos_a2'] = vcf_df_tabix_sample_df['pos_a2'].map(lambda x: 0 if '_' in str(x) else x)
                
                if len(vcf_df_tabix_sample_df) == 0: 
                    gene_coords_tracker.add((exon_start,exon_end))
                    continue
            
                gene_df.append(vcf_df_tabix_sample_df) 
            
                print 121
                gene_coords_tracker.add((exon_start,exon_end))
            except: 
                gene_coords_tracker.add((exon_start,exon_end))
                #print 'except cadd'
                continue
    
    #return arrayify( pd.concat(gene_df) )
    if len(gene_df) >0:
        gene_df = pd.concat(gene_df)
        gene_df = gene_df.drop_duplicates()
        gene_df['a1_a2_score_sum'] = gene_df.pos_a1 + gene_df.pos_a2
        
        array_gene_df = arrayify(gene_df , gene_id)
        array_gene_df.to_csv('/Users/erickscott/datasets_raw/snyderome/fiSSEA.txt',sep="\t",mode='a', header=False)
        return 1
    return 0
    

# <codecell>

# vcf_path = '/Users/erickscott/datasets_raw/snyderome/TB0001907.all.ILLUMNIA.bwa.CEU.high_coverage.20101118.snp.raw.filtered.vcf.gz'

# g_df = get_gene_df(entrez_genes[903])

# <codecell>

# import multiprocessing as mp
# pool = mp.Pool(6)
# tasks = entrez_genes[:400]
# results =[]
# r = pool.map_async(get_gene_df, tasks)
# r.wait()
# pool.terminate()
# pool.join()

# <codecell>

entrez_genes = sys.argv[1]  #list of gene_ids to run
entrez_genes = map(int, entrez_genes.split('\t'))

gene_counter = 0
#for g in gene_ids_df['_id']:
for g in entrez_genes:
    gene_counter += 1
    print 'gene_counter', gene_counter
    g_df = get_gene_df( g)
    if g_df > 0:
        print 'NS_found', g
        #g_df.to_csv('/Users/erickscott/datasets_raw/snyderome/fiSSEA.txt',sep="\t",mode='a', header=False)
        

# <headingcell level=2>

# OLD CODE

# <codecell>

#vcf_df_tabix_sample_df.index.get_level_values(0).unique()[0]

# <codecell>

# chrom = vcf_df_tabix_sample_df.index.get_level_values(0).unique()[0]
# start = vcf_df_tabix_sample_df.head(1).index.get_level_values(1).unique()[0]
# end = vcf_df_tabix_sample_df.tail(1).index.get_level_values(1).unique()[0]

# print chrom, start, end

# <codecell>

#vcf_df_tabix_sample_df['hgvs_id'] = vcf_df_tabix_sample_df['CHROM'] + ":g." + vcf_df_tabix_sample_df['POS'] + vcf_df_tabix_sample_df['REF'] + '>' +  vcf_df_tabix_sample_df['ALT']

# <codecell>



# json_txt = mg.query(7273, fields='exons_hg19',  as_dataframe=True)
# temp_gene = pd.DataFrame.from_dict(json_txt.exons_hg19.values[0])
# #json_txt

# gene_chrom = json_txt.exons_hg19.values[0]['NM_133432']['chr']
# gene_start = json_txt.exons_hg19.values[0]['NM_133432']['cdsstart']
# gene_end = json_txt.exons_hg19.values[0]['NM_133432']['cdsend']
# gene_exons = json_txt.exons_hg19.values[0]['NM_133432']['exons']





# exon_coords = get_cds_coords(gene_exons, gene_start, gene_end)



# <codecell>

# chrom = vcf_df_tabix_sample_df.index.get_level_values(0).unique()[0]
# start = vcf_df_tabix_sample_df.head(1).index.get_level_values(1).unique()[0]
# end = vcf_df_tabix_sample_df.tail(1).index.get_level_values(1).unique()[0]

# cadd = get_fi_score('chr'+chrom, start, end, 'cadd', 'raw')

# cadd['pos_a1'] = cadd['score']
# cadd_score_dict = cadd[['pos_a1']].to_dict()

# vcf_df_tabix_sample_df1 = vcf_df_tabix_sample_df.copy()
# vcf_df_tabix_sample_df1[vcf_df_tabix_sample_df1['pos_a1'].isin(cadd_score_dict['pos_a1'].keys())]
# vcf_df_tabix_sample_df.replace(cadd_score_dict, inplace=True)

# <codecell>

# c = vcf_df_tabix_sample_df1[vcf_df_tabix_sample_df1['pos_a1'].isin(cadd_score_dict['pos_a1'].keys())]

# <codecell>

# d = c.replace(cadd_score_dict_a1)

# <codecell>

# def get_fi_score(gene, scaling=False, build=hg19, fi_model=cadd):
	
# gene_df = ??????  

# 	for exon in gene_exons:

# 		vcf_df = tabix[chr, exon.cds_start, exon.cds_end]

# 		myvariant_df = my_variant( [chr, exon.cds_start, exon.cds_end], cadd )
			
# 		intersect(vcf_df, myvariant_df)  #convert genotypes into cadd score

# 		#scale cadd score by allele frequency in 1000g (eur_af)

# 		if gene_df == ??????:
# 			gene_df = intersected_cadd_scores
# 		if gene_df != ??????:
# 			gene_df.append(intersected_cadd_scores)

# 	gene_series = get max cadd score for gene for each person

# 	return gene_series  (index=gene, values=max cadd score)

