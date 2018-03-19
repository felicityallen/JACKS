import sys
import csv
import os
import io

if len(sys.argv) != 2:
    print 'Usage: combine_single_jacks.py single_jacks_dir'
else:

    single_jacks_dir = sys.argv[1]
    all_files = [x for x in os.listdir(single_jacks_dir) if 'gene' in x]
    num_dir_toks = len(single_jacks_dir.split('_'))

    single_results = {}
    for filename in all_files:
        toks = filename.split('_')
        cell_line = '_'.join(toks[3:-3])
        f = io.open(single_jacks_dir + '/' + filename)
        f.readline()
        single_results[cell_line] = {toks[0]:toks[1] for toks in csv.reader(f, delimiter='\t')}
        f.close()

    cell_lines = single_results.keys()
    fout = io.open('%s_gene_JACKS_results.txt' % single_jacks_dir, 'w')
    fout.write(u'Gene\t%s\n' % '\t'.join(cell_lines))
    genes = single_results[cell_lines[0]].keys()
    for gene in genes:
        fout.write(u'%s\t%s\n' % (gene, '\t'.join([single_results[cell_line][gene] for cell_line in cell_lines])))
    fout.close()
