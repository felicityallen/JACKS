import os
import io
import sys
import random

if len(sys.argv) != 5:
    print 'Usage: run_mageck.py infile test_line num_rep outprefix'
else:

    infile = sys.argv[1]
    ctrl_col_hdrs = sys.argv[2]
    test_col_hdrs = sys.argv[3].split(',')
    num_rep = eval(sys.argv[4])
    outprefix = sys.argv[5]
    
    #Fetch the cell lines
    f = io.open(infile)
    full_hdrs = [x for x in f.readline()[:-1].split('\t')[2:]]
    hdrs = [x.split('_')[0] for x in full_hdrs]
    f.close()
    ctrl_str = ','.join([x for x in full_hdrs if 'CTRL' in x])
    
    test_cols = [x for x in full_hdrs if test_line in x]
    if num_rep < 0: selected_cols = test_cols
    else: selected_cols = random.sample(test_cols, num_rep)       
    test_str = ','.join(test_cols)
           
    cmd = 'mageck test -k %s -t %s -c %s -n %s' % (infile, test_str, ctrl_str, outprefix)
    print cmd
    os.system(cmd)

