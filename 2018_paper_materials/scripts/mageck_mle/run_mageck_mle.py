import io, os, csv, random
import sys

def prepareFile(filename, hdr):
    #Count any lines before the headers (should be skipped)
    f = io.open(filename)
    skip_lines, line = 0, f.readline()  
    while hdr not in line and skip_lines < 100: skip_lines += 1; line = f.readline()
    f.close()
    
    if skip_lines >= 100:
        raise Exception('Could not find line with header ' + hdr + ' in ' + filename)
    
    #Check for comma vs tab delimited
    delim = ',' if (filename.split('.')[-1] == 'csv') else '\t'
    
    #Reopen the file and skip to the start of the data
    f = io.open(filename); [f.readline() for i in range(skip_lines)]  
    return f, delim

def convertChars(rep):
    for char in [' ','.',';',',']:
        rep = rep.replace(char,'_')
    return rep

def createDesignMatrixFile(repfile, rep_hdr, sample_hdr, ctrl_sample, outfile):
    f, delim = prepareFile(repfile, rep_hdr)
    rdr = csv.DictReader(f, delimiter=delim)
    rep_map = {row[rep_hdr]: row[sample_hdr] for row in rdr}
    cell_lines = [x for x in set([rep_map[x] for x in rep_map]) if x != ctrl_sample]
    f.close()
    fout = io.open(outfile, 'w')
    fout.write(u'Samples\tbaseline\t%s\n' % '\t'.join(cell_lines))
    for rep in rep_map:
        fout.write(u'%s\t1\t%s\n' % (convertChars(rep), '\t'.join(['%d' % int(x==rep_map[rep]) for x in cell_lines])))
    fout.close()
    
def formatCountFile( countfile, guidemappingfile, sgrna_hdr, gene_hdr, outprefix ):
    delim = ',' if countfile[-4:]=='.csv' else '\t'
    f = io.open(countfile); rdr = csv.DictReader(f, delimiter=delim)
    guide_col = rdr.fieldnames[0]
    has_gene_col = (rdr.fieldnames[1].upper() == 'GENE')
    if has_gene_col and delim=='\t' and sum([convertChars(x) != x for x in rdr.fieldnames[1:]])==0: return countfile
    if not has_gene_col:
        fg = io.open(guidemappingfile)
        fg_delim = ',' if guidemappingfile[-4:]=='.csv' else '\t'
        guide_map = {row[sgrna_hdr]: row[gene_hdr] for row in csv.DictReader(fg, delimiter=fg_delim)}
        fg.close()
    newcountfile = outprefix + '_' + countfile.split('/')[-1][:-4] + '.txt'
    fout = io.open(newcountfile, 'w')
    fout.write(u'sgRNA\tgene\t%s\n' % '\t'.join([convertChars(x) for x in rdr.fieldnames[1+has_gene_col:]]))
    for row in rdr:
        if has_gene_col: gene = row[rdr.fieldnames[1]]
        elif row[guide_col] not in guide_map:
            continue
        else: gene = guide_map[row[guide_col]]
        try:
            [eval(row[hdr]) for hdr in rdr.fieldnames[1+has_gene_col:]]
        except:
            print 'EXCEPTION IN ROW:', row
        fout.write(u'%s\t%s\t%s\n' % (row[guide_col], gene, '\t'.join(row[hdr] for hdr in rdr.fieldnames[1+has_gene_col:])))
    fout.close()
    return newcountfile


if len(sys.argv) != 6:
    print 'Usage: run_mageck_mle.py countfile replicatefile:rep_hdr:sample_hdr:ctrl_sample genemapfile:guide_hdr:gene_hdr outputprefix threads'
else:

    #Parse arguments
    countfile = sys.argv[1] #sgrna label is always in the first column, gene in the second
    rep_toks = sys.argv[2].split(':')
    if len(rep_toks) != 4:
        raise Exception('Incorrect replicate file input: expecting "replicatefile:rep_hdr:sample_hdr:ctrl_sample" where replicatefile is a csv or tab delimited file mapping replicates to samples, rep_hdr and sample_hdr specify the column headers for the columns containing the replicate labels and sample labels respectively, and ctrl_sample specifies the name of the control sample')
    replicatefile, rep_hdr, sample_hdr, ctrl_sample = rep_toks
    
    guide_toks = sys.argv[3].split(':')
    if len(guide_toks) != 3:
        raise Exception('Incorrect guidemappingfile input: expecting "guidemappingfile:sgrna_hdr:gene_hdr" where guidemappingfile is a csv or tab delimited file mapping guides to genes, sgrna_hdr and gene_hdr specify the column headers for the columns containing the guide labels and gene labels respectively.')
    guidemappingfile, sgrna_hdr, gene_hdr = guide_toks   
    
    outprefix = sys.argv[4]
    outdir = '/'.join(outprefix.split('/')[:-1])
    if outdir != '' and not os.path.isdir(outdir):
        os.mkdir(outdir)
    threads = eval(sys.argv[5])

    #Create the design matrix file
    designfile = outprefix + '_designmat.txt'
    createDesignMatrixFile(replicatefile, rep_hdr, sample_hdr, ctrl_sample, designfile)
 
    #Adjust the input file (if needed)
    countfile = formatCountFile( countfile, guidemappingfile, sgrna_hdr, gene_hdr, outprefix ) 
    
    #Run Mageck-MLE
    cmd = "mageck mle -k %s -d %s -n %s --threads %d" % (countfile, designfile, outprefix, threads)
    print cmd 
    os.system(cmd)
   

