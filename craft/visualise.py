import matplotlib
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy
import math

def ld_block(array, names,
             mid = None,
             left = None,
             right = None,
             cmap = 'Reds',
             colorbar = True):
    """Create and return an linkage-disequilibrium block chart.
    `array` is a square numpy array containing LD values. Only the upper
    triangle of the array (above the diagonal) is used.

    `names` is an iterable of names, which should be the same length as the
    number of rows in `array`.

    `mid` is a label to go in the centre of a bar across the top of
    the chart, e.g. 'chr17'.

    `left` is a label to go at the left of that top bar (e.g. 12398013).

    `right` is a label to go at the right of that top bar (e.g. 18290324).

    If there are no labels then no label bar is drawn.

    `cmap` is the name of a Matplotlib colormap.

    `colorbar` controls whether to add a colorbar. The default is True.

    """

    labels = not ((mid is None) and (left is None) and (right is None))
    if labels:
        fig, (chr_ax,ax) = plt.subplots(2, gridspec_kw=dict(height_ratios=(1,10), hspace=0))
    else:
        fig, ax = plt.subplots()

    # how many items
    lds = array.shape[0]

    # remove lower triangle and diagonal of the array
    mask =  numpy.tri(lds)
    array = numpy.ma.array(array, mask=mask)

    # set the colormap
    cmap = matplotlib.cm.get_cmap(cmap)
    cmap.set_bad('w') # so lower triangle and diagonal are white

    # draw the actual block grid and rotate it as needed.
    im = ax.imshow(array, cmap=cmap)
    trans_data = ax.transData
    transform = transforms.Affine2D().rotate_deg(-45) + trans_data
    im.set_transform(transform)
    ax.set_transform(trans_data)

    # add labels for each item.
    s2 = math.sqrt(2)
    for i,n in enumerate(names):
        ax.text(i*s2, 0, n,
                rotation='vertical',
                horizontalalignment='center',
                verticalalignment='bottom')

    # don't draw axes, 
    ax.set_axis_off()

    # set the limits so we only see the necessary part of the grid.
    ax.set_xlim(0, s2*lds)
    ax.set_ylim(-s2*lds/2, s2)

    # draw the (optional) color bar.
    if colorbar:
        fig.colorbar(im, ax=ax, shrink=0.5)

    # draw the label line and labels if required.
    if labels:
        chr_ax.spines['right'].set_visible(False)
        chr_ax.spines['left'].set_visible(False)
        chr_ax.spines['top'].set_visible(False)
        chr_ax.tick_params(which='both', bottom=False, left=False, labelbottom=False, labelleft=False)
        chr_ax.set_xlim(0,1)
        chr_ax.set_ylim(0,1)
        if left is not None:
            chr_ax.text(0,0, str(left), horizontalalignment='left', verticalalignment='bottom')
        if right is not None:
            chr_ax.text(1,0, str(right), horizontalalignment='right', verticalalignment='bottom')
        if mid is not None:
            chr_ax.text(0.5,0, str(mid), horizontalalignment='center', verticalalignment='bottom')

    return fig

def manhattan(df, chr, alpha=5e-8):
    df['minuslog10pvalue'] = -numpy.log10(df.pvalue)
    df['position'] /= 1e6 # megabases, to suppress the annoying useOffset nonsense in the MPL tick formatter.
    fig, ax = plt.subplots()
    ax.scatter(df['position'],df['minuslog10pvalue'], s=1)
    ax.axhline(-numpy.log10(alpha), linestyle='--', linewidth=0.5, color='0.5')

    index_df = df[df.pvalue < alpha]
    ax.scatter(index_df['position'], index_df['minuslog10pvalue'], s=1.5, color='g',marker='D')
    for i, row in index_df.iterrows():
        t = ax.text(row['position'], row['minuslog10pvalue']+0.1, row['rsid'],
                    rotation='vertical', horizontalalignment="center", verticalalignment='bottom')
    ax.set_ylabel(r'$-\log_{10}\ (P)$')
    ax.set_xlabel(f'Chromosome {chr} (Mbp)')
    ax.ticklabel_format(useOffset=False, style='plain')

    return fig

def locus(df, threshold=0.8):
    df['position'] /= 1e6
    fig, (posax, trackax, geneax) = plt.subplots(3, gridspec_kw=dict(height_ratios=(8,2,2), hspace=0.1), sharex=True)
    posax.scatter(df['position'],df['posterior'], s=2)
    posax.axhline(threshold, linestyle='--', linewidth=0.5, color='0.5')

    cred_df = df[df.posterior > threshold]
    posax.scatter(cred_df['position'], cred_df['posterior'], s=1.5, color='r',marker='D')
    for i, row in cred_df.iterrows():
        t = posax.text(row['position'], row['posterior']+0.01, row['rsid'],
                       rotation='vertical', horizontalalignment="center", verticalalignment='bottom')
    posax.set_ylabel('Posterior probability')

    trackax.set_ylabel('tracks')
    trackax.spines['top'].set_visible(False)
    trackax.spines['left'].set_visible(False)
    trackax.spines['right'].set_visible(False)
    trackax.tick_params(which='both', left=False, labelleft=False)

    geneax.set_ylabel('genes')
    geneax.spines['top'].set_visible(False)
    geneax.spines['left'].set_visible(False)
    geneax.spines['right'].set_visible(False)
    geneax.tick_params(which='both', left=False, labelleft=False)

    geneax.set_xlabel(f'position (Mbp)')

    return fig
