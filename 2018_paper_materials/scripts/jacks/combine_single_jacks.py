import io, os, sys
import pandas as pd

def mergeFiles(files, outfile):
    merged_data = pd.read_csv(files[0],sep='\t')
    for filename in files[1:]:
        try:
            data = pd.read_csv(filename, sep='\t')
        except pd.errors.EmptyDataError:
            print('No data in', filename)
            continue
        merged_data = pd.merge(merged_data, data, on='Gene')
    merged_data.to_csv(outfile, sep='\t', index=False)
    print(merged_data.shape)

dirname = sys.argv[1]
outdirname = dirname + '_combined'
if not os.path.isdir(outdirname): os.makedirs(outdirname)

all_files = os.listdir(dirname)

file_endings = ['gene', 'gene_pval', 'gene_std','pseudo_combined_gene', 'pseudo_combined_gene_pval','pseudo_combined_gene_std','pseudo_noness_gene','pseudo_noness_gene_std']

for fending in file_endings:

    pseudo = ('pseudo' in fending)
    indiv_files = [dirname + '/' + x for x in all_files if x[-18-len(fending):-18] == fending and ((not pseudo and 'pseudo' not in x) or pseudo)] 
    print(pseudo, fending, len(indiv_files))
    mergeFiles(indiv_files, outdirname + '/' + dirname.split('/')[-1] + '_%s_JACKS_results.txt' % fending)

