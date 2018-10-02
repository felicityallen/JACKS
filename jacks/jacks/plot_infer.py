import matplotlib
import numpy as np
import pylab as PL
import scipy as SP
import scipy.stats as ST
from matplotlib import gridspec


def plot_inference_result(y,w1,w2,x1,x2,tau,cell_lines = [], muts=[], title='', figname='test.png'):
    matplotlib.rcParams.update({'font.size': 12}) 
    fig = PL.figure(figsize=(9,6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[len(w1),len(x1)])   

    cell_lines = ['LNCaP' if ('LNCaP' in x) else x for x in cell_lines]
    mut_status = ['(M)' if mut=="True" else ('' if mut=="False" else '(U)') for mut in muts]
    
    #Signal
    ax = PL.subplot(gs[0,0])    
    im = PL.imshow(y, aspect=1.15, interpolation='none', cmap=PL.get_cmap("coolwarm"), vmin=-3, vmax=3)
    ax = PL.gca()
    ax.set_xticks([])      
    ax.set_yticks(range(len(x1)))
    ax.set_yticklabels(['gRNA %d' % (grnano +1) for grnano in range(len(x1))])         
    if len(cell_lines) > 0:
        ax.set_xticks(range(len(cell_lines))) 
        ax.set_xticklabels(cell_lines, rotation='vertical') 
        ax.xaxis.tick_top()  
        for t in ax.xaxis.get_ticklines(): t.set_visible(False)
        for t in ax.yaxis.get_ticklines(): t.set_visible(False)
        for t,mt in zip(ax.xaxis.get_ticklabels(),mut_status): 
            t.set_fontsize(10)
            if mt == '(M)':
                t.set_fontweight('bold')
    
    #x
    PL.subplot(gs[0,1])
    PL.plot([1,1],[-1,len(x1)+1], 'k--')
    PL.plot([0,0],[-1,len(x1)+1], 'k--')
    vdata = [ST.norm.rvs(x1[i], (x2[i]-x1[i]**2)**0.5, size=5000) for i in range(len(x1))]
    vpos = (SP.arange(len(x1))+1)[::-1]
    clrs = ['#FFFFFF','#BBCCEE']
    for i in range(len(x1)):
        PL.axhspan(i+0.5,i+1.5,facecolor=clrs[i%2],alpha=0.1)
    vplot = PL.violinplot(vdata, vpos, widths=0.5*SP.ones(len(x1)), vert=False, showextrema=True, showmeans=True)
    for patch, val in zip(vplot['bodies'], x1): 
        col_val = int(0.8*min(max(256-val*128,0),255))
        patch.set_color('#%02x%02x%02x' % (col_val, col_val, col_val))
    vplot['cmeans'].set_color('darkblue')
    vplot['cmeans'].set_linewidth(2)
    vplot['cmins'].set_color('#444444')
    vplot['cmaxes'].set_color('#444444')
    vplot['cbars'].set_color('#444444')       
    vplot['cbars'].set_visible(False)
    PL.ylim(0.5,len(x1)+0.5); 
    PL.xlim(-0.5,2)
    ax = PL.gca()
    PL.xticks([0,1])
    PL.yticks([])
    PL.xlabel("gRNA efficacy")
    PL.title(title)


    #w
    PL.subplot(gs[1,0])
    clrs = ['#FFFFFF','#BBCCEE']
    for i in range(len(w1)):
        PL.axvspan(i-1,i,facecolor=clrs[i%2],alpha=0.1)  
    vplot = PL.violinplot([ST.norm.rvs(w1[i], (w2[i]-w1[i]**2)**0.5, size=5000) for i in range(len(w1))], SP.arange(len(w1))-0.5, widths=0.9*SP.ones(len(w1)), showmeans=True, showextrema=True)
    for patch, val in zip(vplot['bodies'], w1): 
        col_val = int(0.95*min(max(0,256+val*128),200))
        patch.set_alpha(1.0)
        clr = im.cmap(im.norm(val))
        patch.set_color(clr)
    vplot['cmeans'].set_color('darkblue')
    vplot['cmeans'].set_linewidth(2)
    vplot['cmins'].set_color('#444444')
    vplot['cmaxes'].set_color('#444444')
    vplot['cbars'].set_visible(False)
    PL.plot([-1,len(w1)], [0.0,0.0], 'k--')
    PL.xlim(-1, len(w1)-1)
    PL.ylim(-3.5,1)
    mean_y = np.nanmean(y,axis=0)
    PL.ylabel("Gene essentiality")
    PL.xlabel("Cell lines")
    pws = [1.0-ST.norm.cdf((w1[i])/np.sqrt(w2[i]-w1[i]*w1[i])) for i in range(len(w1))]
    ax = PL.gca()

    if len(cell_lines) > 0:
        ax.set_xticks([])
        for t in ax.xaxis.get_ticklines(): t.set_visible(False)

    PL.subplots_adjust(left=0.08,right=0.94,top=0.82, bottom=0.11, wspace=0.0, hspace=0.0)
    PL.rcParams['svg.fonttype'] = 'none'
    PL.savefig(figname, bbox_inches='tight')
    PL.show(block=False)
    return fig

