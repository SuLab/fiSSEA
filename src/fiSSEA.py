

'''
November 21, 2014
Functional Impact SNP Set Enrichment Analysis (fiSSEA)
Authors: Erick R. Scott, Adam Mark, Ali Torkamani, Chunlei Wu, Andrew Su
The Scripps Research Institute

'''
#Libraries
import os, sys, time, gc, argparse
import json, httplib2
import pandas as pd
from pandas.io.json import json_normalize
import numpy as np

#fiSSEA Libraries
sys.path.append('../pandasVCF/src/single_sample/')
from pdVCFsingle import *

sys.path.append('../myvariant/src/')
import myvariant







def myvariant_post(hgvs_list):
    '''
    Query and Parser for myvariant.info
    Parses Raw Elastic Search results into
    pandas dataframe
    
    Parameters
    -------------
    hgvs_list: list, required
        
        
    Output
    -------------
    pandas df:
        normalized json of myvariant results
    '''
    
    def normalize(con):
        '''
        This function uses the json.loads and json_normalize
        modules to flatten the myvariant.info results
        '''
        temp_df = json_normalize(json.loads(con), 'docs')
        source_parsed = []
        for source_json in temp_df['_source'].values:
            try:
                source_parsed.append(json_normalize(source_json))
            except TypeError:  #occurs with empty results
                pass
        return pd.concat(source_parsed)
    
    
    
    if type(hgvs_list) == list:
        hgvs_list = ','.join(hgvs_list)
        
    assert type(hgvs_list) == str
    
    h = httplib2.Http()
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    params = 'ids=' + hgvs_list +",size=1000"  #can pass a dictionary
    #print params
    res, con = h.request('http://myvariant.info:8001/api/variant/', 'POST', params, headers=headers)
    
    
    mv_df = normalize(con)
    mv_df.index = mv_df['_id']

    return mv_df
    
    


def get_myvariant_snp_annot(vcf_df_mv):
    '''
    Returns myvariant.info annotations for all snps for
    input vcf pandas dataframe. Submits myvariant.info
    POST request in snp batches of 1000.
    
    Parameters
    -------------
    vcf_df_mv: pandas df; required
    
        vcf_df_mv must contain the following columns:
            vartype1, vartype2, CHROM, REF, POS, a1, a2
            *please see pandasVCF docs for column definitions
            

    Output
    -------------
    myvariant.info results: pandas df
    
    '''
    
    def chunks(l, n):
        """ 
        Yield successive n-sized chunks from l.
        """
        l = list(l)
        chunk_list = []
        for i in xrange(0, len(l), n):
            chunk_list.append( l[i:i+n] )
        return chunk_list
    
    
 
    if 'CHROM' not in vcf_df_mv.columns:
        vcf_df_mv.reset_index(inplace=True)
        
    
    vcf_df_mv = vcf_df_mv[(vcf_df_mv['vartype1'] == 'snp') | (vcf_df_mv['vartype2'] == 'snp')] #restrict to snps

    non_ref_a1 = vcf_df_mv[vcf_df_mv['REF'] != vcf_df_mv['a1']][['CHROM', 'POS', 'REF', 'a1']]  #restrict a1 to non-reference alleles
    hgvs_a1 = 'chr' + non_ref_a1['CHROM'].astype(str) + ":g." + non_ref_a1['POS'].astype(str)  \
                                 + non_ref_a1['REF'] + '>' +  non_ref_a1['a1']

    non_ref_a2 = vcf_df_mv[vcf_df_mv['REF'] != vcf_df_mv['a2']][['CHROM', 'POS', 'REF', 'a2']]  #restrict a2 to non-reference alleles
    hgvs_a2 = 'chr' + non_ref_a2['CHROM'].astype(str) + ":g." + non_ref_a2['POS'].astype(str) \
                                 + non_ref_a2['REF'] + '>' +  non_ref_a2['a2']


    hgvs_ids = set(hgvs_a1.values) | set(hgvs_a2.values)  #set of non-reference a1 and a2 snps

    hgvs_ids_chunks = chunks(hgvs_ids, 1000)  #batch queries of hgvs a1 and a2 variants, n=1000

    myvariant_results = map(myvariant_post, hgvs_ids_chunks)
        
    return myvariant_results
    





class fiSSEA(object):

    
    def __init__(self, vcf_file_path, sample_id='', chunksize=100000, fiSSEA_score='dbnsfp.cadd.phred', gene_stat=np.mean):
        
        
        #setting fiSSEA_score
        self.fiSSEA_score = fiSSEA_score
        
        #setting gene statistic
        self.gene_stat = gene_stat
        
        #open Vcf object
        vcf_df_obj = Vcf(vcf_file_path, sample_id, ['#CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'FORMAT'], chunksize=chunksize)

        #read in vcf chunk
        vcf_df_obj.get_vcf_df_chunk()

        #drop duplicate vcf lines, if necessary
        vcf_df_obj.df = vcf_df_obj.df.drop_duplicates(subset=['CHROM', 'POS'])

        #remove empty variant lines
        vcf_df_obj.df = vcf_df_obj.df[vcf_df_obj.df.ALT!='.']

        #replace chr if necessary, depends on format of annotations
        vcf_df_obj.df['CHROM'] = vcf_df_obj.df['CHROM'].str.replace('chr', '')

        #parse vcf lines
        vcf_df_obj.add_variant_annotations( inplace=True, verbose=True)
        
        #settting vcf_df as class object
        self.vcf_df = vcf_df_obj.df

        #get myvariant.info annotations, create fiSSEA dataframe
        self.fiSSEA_status = self.get_fi_scores()
        

    
    def get_multigene_var_records(self, df_line):
        '''
        Splits multigene variants into separate
        rows
        '''
        genes = df_line['dbnsfp.genename']
        gene_series = [] #list of gene variants
        for gene in genes: #iterate through all genes intersecting this variant
            df_line_copy = df_line.copy()  #avoids aliasing
            del df_line_copy['dbnsfp.genename']  #removing genename to allow gene insertion
            df_line_copy['dbnsfp.genename'] = gene
            gene_series.append(df_line_copy)

        return gene_series
    
    
    
    def get_fi_scores(self):
        '''
        Creates fi_scores series, index on gene symbol, value fi_score.
        '''
        
        #Retrieving myvariant.info annotations for all snps
        myvariant_results = get_myvariant_snp_annot(self.vcf_df) #query and parse myvariant.info using all snps in self.vcf_df
        self.mv_annot = pd.concat(myvariant_results)  #setting dataframe of parsed myvariant.info results

        #reducing size of dataframe for fiSSEA
        fiSSEA_cols = ['dbnsfp.genename', self.fiSSEA_score]
        mv_annot_fiSSEA = self.mv_annot[fiSSEA_cols] 
        
        #dropping variants without a genename & filling fi_score NaN with 0
        mv_annot_fiSSEA.dropna(subset=['dbnsfp.genename'], inplace=True)  
        mv_annot_fiSSEA[self.fiSSEA_score].fillna(value=0)  

        
        #Splitting variants intersecting multiple genes into separate rows
        mv_annot_fiSSEA['dbnsfp.genename_type'] = mv_annot_fiSSEA['dbnsfp.genename'].map(type)
        mv_annot_fiSSEA_multigene = mv_annot_fiSSEA[mv_annot_fiSSEA['dbnsfp.genename_type']==list]
        mv_annot_fiSSEA = mv_annot_fiSSEA[~mv_annot_fiSSEA.index.isin(mv_annot_fiSSEA_multigene.index)]
        
        if len(mv_annot_fiSSEA_multigene) > 0:
            mv_annot_fiSSEA_multigene['index'] = mv_annot_fiSSEA_multigene.index
            mv_annot_fiSSEA_multigene.drop_duplicates(subset='index', inplace=True)
            mv_annot_fiSSEA_multigene = mv_annot_fiSSEA_multigene.apply(self.get_multigene_var_records, axis=1)
            mv_annot_fiSSEA_multigene = [n for i in mv_annot_fiSSEA_multigene for n in i]
            mv_annot_fiSSEA = mv_annot_fiSSEA.append(pd.DataFrame(mv_annot_fiSSEA_multigene))
            del mv_annot_fiSSEA_multigene
        
        gc.collect()

        
        #Setting gene specific functional impact score
        self.fi_scores = mv_annot_fiSSEA.groupby('dbnsfp.genename')[self.fiSSEA_score].aggregate(self.gene_stat)
        
        return 0








parser = argparse.ArgumentParser()
parser.add_argument('-i','--vcf_input', type=str,help='Specifies the input coding snp vcf file, /path/to/cds.vcf')
parser.add_argument('-o','--output_path', type=str,help='Specifies the output file path, e.g. /path/to/output/counts_df.txt')

parser.add_argument('-s','--sample_id', type=str,help='Specifies the sample identifier in the vcf file, e.g. NA12878')
parser.add_argument('-c','--chunk', type=str,help='Specifies the chunksize for vcf iteration') 
parser.add_argument('-f','--fi_score', type=str,help='Specifies the name of the functional impact score, e.g. dbnsfp.cadd.phred') 
parser.add_argument('-g','--gene_stat', type=str,help='Specifies the gene statistic, e.g. np.mean')


opts = parser.parse_known_args()
vcf_path, output_path, sample_id = opts[0].vcf_input, opts[0].output_path, opts[0].sample_id 
fi_score, gene_statistic, chunk = opts[0].fi_score, opts[0].gene_stat, opts[0].chunk


mv = myvariant.MyVariantInfo()


if __name__ == '__main__':
    
    
    if not fi_score:
        fi_score = 'dbnsfp.cadd.phred'
        
    if not gene_statistic:
        gene_statistic = np.mean
    
    if not chunk:
        chunk = 100000
    


    fiSSEA_results = fiSSEA(vcf_path, sample_id=sample_id, chunksize=chunk, fiSSEA_score=fi_score, gene_stat=gene_statistic)
    
    
    if not output_path:
        fiSSEA_results.fi_scores.to_csv(sys.stdout)
    else:
        fiSSEA_results.fi_scores.to_csv(output_path, sep='\t')
    
    

    