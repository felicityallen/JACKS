import sys, os, random, csv, io
from jacks.plot_infer import plot_inference_result
from jacks.io_preprocess import loadJacksFullResultsFromPickle
import pickle

if len(sys.argv) != 4:
    print('Usage: gen_heatmap.py picklefile gene(or "random" to randomly pick one) outfile')
else:

    picklefile = sys.argv[1]
    gene = sys.argv[2]
    outfile = sys.argv[3]
    dirname = os.path.dirname(outfile)
    if '/' in outfile and not os.path.exists(dirname): os.makedirs(dirname)
    jacks_results, cell_lines, gene_grnas = loadJacksFullResultsFromPickle(picklefile)

    if gene == 'random':
        gene = random.choice([x for x in jacks_results.keys()])
        print('Selected gene: %s' % gene)
    elif gene not in jacks_results:
        raise Exception('No results for gene %s in JACKS results in %s' % (gene, picklefile))
     
    y, tau, x1, x2, w1, w2 = jacks_results[gene]
    plot_inference_result(y,w1,w2,x1,x2,tau,cell_lines = cell_lines, title=gene, figname=outfile)
    
