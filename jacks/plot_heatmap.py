import sys, os, random, csv, io
from jacks.code.plot_infer import plot_inference_result
from jacks.code.io_preprocess import loadJacksFullResultsFromPickle
import pickle

if len(sys.argv) != 3 and len(sys.argv) != 4 and len(sys.argv) != 5:
    print 'Usage: gen_heatmap.py picklefile gene(or "random" to randomly pick one) (opt)mutationfile del_gene_pickle'
else:

    plot_dir = 'heatmaps'
    if not os.path.isdir(plot_dir):
        os.mkdir(plot_dir)

    picklefile = sys.argv[1]
    gene = sys.argv[2]
    mut_mapping = {}
    if len(sys.argv) >= 4:
        mutfile = sys.argv[3]
        if mutfile != 'None':
    
            f = io.open(mutfile)
            delim = ',' if mutfile[-4:] == '.csv' else '\t'
            mut_mapping = {row['Cell Line']: row for row in csv.DictReader(f, delimiter=delim)}
            f.close()
        
    del_gene_pickle = False
    if len(sys.argv) >= 5:
        del_gene_pickle = eval(sys.argv[4])
    
    genepickle = 'gen_heatmap_%s.pickle' % gene
    if os.path.isfile(genepickle) and not del_gene_pickle:
        y, tau, x1, x2, w1, w2, cell_lines,muts = pickle.load(io.open(genepickle,'rb'))
    else:
    
        jacks_results, cell_lines, gene_grnas = loadJacksFullResultsFromPickle(picklefile)

        if gene == 'random':
            gene = random.choice(jacks_results.keys())
        elif gene not in jacks_results:
            raise Exception('No results for gene %s in JACKS results in %s' % (gene, picklefile))
         
        if len(mut_mapping) > 0:         
            muts = [mut_mapping[x][gene] if gene in mut_mapping[x] else 'Unknown' for x in cell_lines ]
        else: muts = []
         
        y, tau, x1, x2, w1, w2 = jacks_results[gene]
        pickle.dump((y, tau, x1, x2, w1, w2,cell_lines,muts), io.open(genepickle,'wb'))
   
    plot_inference_result(y,w1,w2,x1,x2,tau,cell_lines = cell_lines, title=gene, muts=muts, figname='%s/%s.svg' % (plot_dir,gene))
    
    import pdb; pdb.set_trace()
    
        
    
