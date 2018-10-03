import argparse
import os
import random

from jacks.jacks_io import loadJacksFullResultsFromPickle
from jacks.plot_infer import plot_inference_result


def plot_heatmap(picklefile, gene, outfile):
    dirname = os.path.dirname(outfile)
    if '/' in outfile and not os.path.exists(dirname): os.makedirs(dirname)
    jacks_results, cell_lines, gene_grnas = loadJacksFullResultsFromPickle(picklefile)

    if gene == 'random':
        gene = random.choice([x for x in jacks_results.keys()])
        print('Selected gene: %s' % gene)
    elif gene not in jacks_results:
        raise Exception('No results for gene %s in JACKS results in %s' % (gene, picklefile))

    y, tau, x1, x2, w1, w2 = jacks_results[gene]
    return plot_inference_result(y, w1, w2, x1, x2, tau, cell_lines=cell_lines, title=gene, figname=outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("picklefile",
                        help="Pickle file with JACKS results")
    parser.add_argument("gene",
                        help='Gene name; use "random" to randomly pick a gene')
    parser.add_argument("outfile",
                        help="Name of the output file")
    args = parser.parse_args()
    plot_heatmap(args.picklefile, args.gene, args.outfile)
