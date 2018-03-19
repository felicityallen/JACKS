
import os
import io
import sys
import random

def readBagelScores(filename):
    f = io.open(filename)
    scores = {row['GENE']: eval(row['BF']) for row in csv.DictReader(f,delimiter='\t')}
    f.close()
    return scores

if len(sys.argv) != 2:
    print 'Usage: collate_bagel_results.py results_dir'
else:

    results_dir = sys.argv[1]
    
    files = os.listdir(results_dir)
    all_results = {}
    for filename in results_dir:
        cell_line = filename.split('_')[-1][:-4]
        all_results[cell_line] = readBagelScores(results_dir + '/' + filename)
    
    fout = io.open(results_dir + '.txt','w')
    cell_lines = all_results.keys()
    fout.write(u'Gene\t%s\n' % '\t'.join(cell_lines))
    genes = all_results[cell_lines[0]].keys()
    for gene in genes:
        fout.write(u'%s\t%s\n' % (gene, '\t'.join([all_results[cellline][gene] for cellline in cell_lines])))
    fout.close()