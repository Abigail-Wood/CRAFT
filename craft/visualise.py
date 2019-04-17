import math

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import transforms
import numpy
import pandas

def fit_text(ax, *args, **kwargs):
    """Write text into `ax` by passing `args` and `kwargs` to ax.text,
    then change the view limits of `ax` if necessary to include the text.

    Trickier than it sounds, and the results are a little approximate.
    """

    t = ax.text(*args, **kwargs)
    tbox = t.get_tightbbox(ax.get_figure().canvas.get_renderer())
    tbox_in_data = tbox.inverse_transformed(t.get_transform())

    # We should be able to just do:
    #
    # ax.update_datalim_bounds(tbox_in_data)
    # ax.autoscale_view()
    #

    # But the axes data limits are not necessarily what we want. For
    # instance, after doing an ax.imshow with a transform (as we do
    # for ld_block), the ax.dataLim is in the array's data coordinates
    # (whereas the ax.viewLim is in the axes' data coordinates).  That
    # might be a bug but in any case we want really to work with the
    # viewLim, in case the viewLim has previously been deliberately
    # moved. So we do it in a longer way:

    # existing viewLim
    l, b, w, h = ax.viewLim.bounds
    r = l+w
    t = b+h

    # text box bounds in same coordinates
    tl, tb, tw, th = tbox_in_data.bounds
    tr = tl+tw
    tt = tb+th

    # new view limits
    # apply margins where necessary.
    xmargin, ymargin = ax.margins()
    nl = l if tl > l else tl - (r-tl) * xmargin
    nr = r if tr < r else tr + (tr-l) * xmargin
    nb = b if tb > b else tb - (t-tb) * ymargin
    nt = t if tt < t else tt + (tt-b) * ymargin

    ax.set_xlim(nl,nr)
    ax.set_ylim(nb,nt)

def ld_block(array,
             indexes = None,
             names = None,
             labels = None,
             figsize = (8, 5),
             cmap = 'Reds',
             colorbar = True):
    """Create and return an linkage-disequilibrium block chart.  `array`
    is a square numpy array containing LD values (Pearson's
    correlation coefficient r). Only the upper triangle of the array
    (above the diagonal) is used. The values actually plotted are r^2.

    `indexes` is an iterable of index values into the rows and columns
    of `array`, giving the order and identity of the SNPs to
    display. If None, the whole array is displayed in array order.

    `names` is an iterable of names, which should be the same length
    as `indexes`, or the number of rows in `array`, used to display
    row/column names above the LD block. If None, no names are
    displayed.

    `labels`, if present, is a dictionary of labels with available
    keys "mid", "left", and "right", for example
    dict(mid='chr17, left=12398013, right=18290324).

    If there are no labels then no label bar is drawn.

    `figsize` is the size of the figure in inches (width, height).

    `cmap` is the name of a Matplotlib colormap.

    `colorbar` controls whether to add a colorbar. The default is True.

    """

    assert array.ndim == 2

    if labels:
        fig, (lab_ax,ax) = plt.subplots(2,
                                        figsize = figsize,
                                        gridspec_kw = dict(
                                            height_ratios=(1,10),
                                            hspace=0))
    else:
        fig, ax = plt.subplots(figsize=figsize)

    if indexes is not None:
        l = list(indexes) # for consumable iterables
        assert 0 <= min(l)
        assert max(l) < array.shape[0]
        array = array[..., l][l, ...]

    # array contains Pearson's "r" coefficients. We plot r^2.
    # Note: point-wise multiplication, not matrix multiplication!
    array = array*array

    # how many items
    lds = array.shape[0]

    # mask out lower triangle and diagonal of the array
    mask =  numpy.tri(lds)
    array = numpy.ma.array(array, mask=mask)

    # set the colormap
    cmap = matplotlib.cm.get_cmap(cmap)
    cmap.set_bad('w') # so masked values are white

    # draw the actual block grid and rotate it as needed.
    # force colormap range to 0-1.
    im = ax.imshow(array,
                   cmap = cmap,
                   vmin = 0, vmax = 1,
                   transform = (transforms.Affine2D().rotate_deg(-45)
                                + ax.transData))

    # set view limits so we see the necessary part of the grid.
    s2 = math.sqrt(2)
    ax.set_xlim(0, s2*(lds-1))
    ax.set_ylim(-s2*lds/2, 1)

    # add labels for each item, changing the view limits to fit.
    if names is not None:
        for i,n in enumerate(names):
            fit_text(ax, i*s2, 0, n,
                     rotation='vertical',
                     horizontalalignment='center',
                     verticalalignment='bottom')

    # don't draw axes,
    ax.set_axis_off()

    # draw the (optional) color bar.
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, shrink=0.5)
        cbar.set_label(label='$R^2$', rotation=0)

    # draw the label line and labels if required.
    if labels:
        lab_ax.spines['right'].set_visible(False)
        lab_ax.spines['left'].set_visible(False)
        lab_ax.spines['top'].set_visible(False)
        lab_ax.tick_params(which='both',
                           bottom=False,
                           left=False,
                           labelbottom=False,
                           labelleft=False)
        lab_ax.set_xlim(0, 1)
        lab_ax.set_ylim(0, 1)
        for k, x, ha, va in [('left' , 0.0, 'left'  , 'bottom'),
                             ('mid'  , 0.5, 'center', 'bottom'),
                             ('right', 1.0, 'right' , 'bottom'),
                            ]:
            if k in labels:
                lab_ax.text(x, 0,
                            str(labels[k]),
                            horizontalalignment=ha,
                            verticalalignment=va)
    return fig

def manhattan(df, x_label,
              index_df = None,
              alpha = 5e-8,
              figsize = (8, 5),
              size = 1,
              color = 'b',
              marker = '.',
              alpha_line_width = 0.5,
              alpha_line_color = '0.5',
              alpha_line_style = '--',
              good_size = 2,
              good_color = 'g',
              good_marker = 'D',
              good_label_column = 'rsid',
              good_label_rotation = 'vertical'):
    """Draw and return a "Manhattan plot", with values above some
    critical threshold marked distinctly and labelled.

    `df` is a Pandas dataframe with two required columns:

       "pvalue": the P value of each SNP;
       "position": the position of the SNP (in base-pairs);

    If automatic labelling of the most significant SNPs is required,
    (controlled by `index_df` and `alpha` parameters), the DataFrame
    must also have a column labelled by the `good_label_column`
    parameter.

    `index_df`, if present, is a Pandas dataframe like `df`, with the
    same required columns, identifying the SNPs to be labelled.

    `x_label` is a label for the x axis (such as "Chromosome 17").

    `alpha` is a threshold for distinguishing "good" SNPs.

    If `index_df` is present then the SNPs in it will be
    distinguished. If `index_df` is absent, and `alpha` is set, then
    SNPs with pvalue less than `alpha` are distinguished.

    Distinguished SNPs are drawn differently (see the `good_`
    parameters below). If `index_df` is present, or if both `alpha`
    and `good_label_column` are set

    `figsize` is the figure size (width, height) in inches.

    `size` is the point area for the main scatter plot, in square points.

    `color` is the point color for the main scatter plot.

    `marker` is the marker style for the main scatter plot.

    `alpha_line_width` is the width of the horizontal line to be drawn
    at the alpha level.

    `alpha_line_color` is the color specifier for the alpha line.

    `alpha_line_style` is the line style for the alpha line.

    `good_size` is the point area for "good" SNPs, in square points.

    `good_color` is the color for "good" SNPs.

    `good_marker` is the marker style for "good" SNPs.

    `good_label_column` is the column name in `df` or `index_df` for
    labels to be drawn for "good" SNPs, or None for no labels.

    `good_label_rotation` is the rotation for labels for "good" SNPs.

    """

    # calculate -log_10(P) to plot

    df['minuslog10pvalue'] = -numpy.log10(df.pvalue)

     # x axis in megabases, to suppress annoying useOffset nonsense
     # in the MPL tick formatter.
    df['positionMb'] = df.position / 1e6

    fig, ax = plt.subplots(figsize=figsize)

    # plot the main scatter
    ax.scatter(df.positionMb, df.minuslog10pvalue,
               s=size, color=color, marker=marker)

    # if we have a threshold (alpha), draw the line
    if alpha:
        ax.axhline(-numpy.log10(alpha),
                   linestyle=alpha_line_style,
                   linewidth=alpha_line_width,
                   color=alpha_line_color)


    # if we need to distinguish points:
    if alpha or (index_df is not None):
        # which points to distinguish:
        if index_df is None:
            index_df = df[df.pvalue < alpha]
        else:
            index_df['positionMb'] = index_df.position / 1e6
            index_df['minuslog10pvalue'] = -numpy.log10(index_df.pvalue)


        # draw them differently
        ax.scatter(index_df.positionMb, index_df.minuslog10pvalue,
                   s=good_size,
                   color=good_color,
                   marker=good_marker)

        # label them
        if good_label_column:
            for i, row in index_df.iterrows():
                fit_text(ax, row.positionMb, row.minuslog10pvalue +0.1,
                         row[good_label_column],
                         rotation=good_label_rotation,
                         horizontalalignment='left',
                         verticalalignment='bottom')

    # axes and ticks
    ax.set_ylabel(r'$-\log_{10}\ (P)$')
    ax.set_xlabel(f'{x_label} (Mbp)')
    ax.ticklabel_format(useOffset=False, style='plain')

    return fig

def locus(df,
          threshold = 0.8,
          figsize = (5, 8),
          size = 1,
          color = 'b',
          marker = '.',
          alpha_line_width = 0.5,
          alpha_line_color = '0.5',
          alpha_line_style = '--',
          good_size = 2,
          good_color = 'g',
          good_marker = 'D',
          good_label_column = 'rsid',
          good_label_rotation = 'vertical',
          tracks=None,
          track_height = 0.5,
          track_column = 'tracks',
          track_colors = ['r','g','b','c','m','y','k'],
          track_alpha = 0.3,
          track_linelength = 0.7,
          track_good_linelength = 1.0,
          track_lines = False,
          track_line_width=0.5,
          track_line_color='k',
          pos_top = False,
          # future work:
          genes = None,
          gene_height = 0.5,
          ):
    df['positionMb'] = df.position / 1e6
    total_height = 1 + track_height + gene_height
    if tracks and (genes is not None):
        height_ratios = [track_height, gene_height]
        if pos_top:
            height_ratios.push(1)
        else:
            height_ratios.append(1)
        fig, axes = (
            plt.subplots(3, figsize=figsize,
                         gridspec_kw=dict(height_ratios=height_ratios, hspace=0),
                         sharex=True))
        if pos_top:
            posax, trackax, geneax = axes
        else:
            geneax, trackax, posax = axes
        bottomax = axes[-1]

    elif tracks or (genes is not None):
        pos_index = 0 if pos_top else 1
        other_index = 1 - pos_index
        height2 = track_height if tracks else gene_height
        heights = [1,1]
        heights[other_index] = height2
        fig, axes = (
            plt.subplots(2, figsize=figsize,
                         gridspec_kw=dict(height_ratios=heights, hspace=0),
                         sharex=True))
        bottomax = axes[-1]
        posax = axes[pos_index]
        ax2 = axes[other_index]
        if tracks:
            trackax = ax2
        else:
            geneax = ax2

    else: # neither tracks nor genes
        fig, posax = plt.subplots(figsize=figsize)
        bottomax = posax

    posax.scatter(df.positionMb, df.pp,
                  s=size, color=color, marker=marker)

    chromosome = df.chromosome.unique()[0]
    posax.set_ylabel('Posterior probability')
    posax.set_ylim(0,1)
    bottomax.set_xlabel(f'Chromosome {chromosome}; position (Mbp)')

    if threshold:
        posax.axhline(threshold,
                      linestyle=alpha_line_style,
                      linewidth=alpha_line_width,
                      color=alpha_line_color)

        cred_df = df[df.pp > threshold]
        posax.scatter(cred_df.positionMb, cred_df.pp,
                      s=good_size, color=good_color, marker=good_marker)
        if good_label_column:
            for i, row in cred_df.iterrows():
                fit_text(posax,row.positionMb, row.pp+0.01,
                         row[good_label_column],
                         rotation=good_label_rotation,
                         horizontalalignment="left",
                         verticalalignment='bottom')
    print(posax.get_yticks())
    posax.set_yticks([y for y in posax.get_yticks() if y >= 0 and y <= 1])
    print(posax.get_yticks())

    if tracks:
        track_colors = track_colors[:len(tracks)]
        alpha = track_alpha if threshold else 1.0
        colors = [matplotlib.colors.to_rgba(c, alpha) for c in track_colors]
        trackax.set_ylabel('tracks')
        tracks_array = [df[df[track_column].str.contains(t)].positionMb
                        for t in tracks]
        trackax.eventplot(tracks_array, linelengths=track_linelength, colors=colors)
        if threshold:
            if track_lines:
                for p in cred_df.positionMb:
                    trackax.axvline(p, zorder=-1,
                                    linewidth=track_line_width, color=track_line_color)
            tracks_array = [cred_df[cred_df[track_column].str.contains(t)].positionMb
                            for t in tracks]
            trackax.eventplot(tracks_array, linelengths=track_good_linelength, colors=track_colors)

        trackax.set_yticks(range(len(tracks)))
        trackax.set_yticklabels(tracks)
        trackax.set_ylabel('')
        half_linelength = max(track_good_linelength, track_linelength)/2
        trackax.set_ylim(-half_linelength-0.1, len(tracks)-1+half_linelength+0.1)

    if genes is not None:
        geneax.set_ylabel('genes')
        # geneax.spines['top'].set_visible(False)
        # geneax.tick_params(which='both', left=False, labelleft=False)
        xlims = geneax.get_xlim()
        y = 0.25
        for start, end, name, strand in genes:
            geneax.plot((start/1e6, end/1e6), (y,y))
            geneax.text((start + end)/2e6, y+0.05, name)
            y = 1-y
        geneax.set_xlim(xlims)
        geneax.set_ylim(0,1)

    fig.tight_layout()

    return fig

def test():
    fig, ax = plt.subplots()
    ax.eventplot(((1,5,7,8),(3,4,1,6),(5,2,8)))
    fit_text(ax, 3, 2, "Spong")
    fit_text(ax, 2, -3, "Wibble", rotation=-45)
    fit_text(ax, 1, 5, "Foobar", verticalalignment='bottom', rotation='vertical')

    return fig
