
import os
import io
import sys
import random

if len(sys.argv) != 5:
    print 'Usage: run_bagel.py infile cell_line num_rep outprefix'
else:

    infile = sys.argv[1]
    cell_line = sys.argv[2]
    
    f = io.open(infile)
    full_hdrs = [x for x in f.readline()[:-1].split('\t')[2:]]
    hdrs = [x.split('_')[0] for x in full_hdrs]
    f.close()

    test_col_idxs = [idx+1 for idx,x in enumerate(hdrs) if x == cell_line]
    num_rep = eval(sys.argv[3])
    outprefix = sys.argv[4]
    
    #Run Bagel on data from this cell line
    if num_rep < 0: selected_test_idxs = test_col_idxs
    else: selected_test_idxs = random.sample(test_col_idxs, num_rep)
    cmd = 'python BAGEL.py -i %s -o %s -e ../../data/training_essentials.txt -n ../../data/training_nonessential.txt -c %s' % (infile, outprefix, ','.join([str(x) for x in selected_test_idxs]))
    print cmd
    os.system(cmd)
    
    
