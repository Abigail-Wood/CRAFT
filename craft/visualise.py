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

def ld_block(array, names,
             labels = None,
             figsize = (8, 5),
             cmap = 'Reds',
             colorbar = True):
    """Create and return an linkage-disequilibrium block chart.
    `array` is a square numpy array containing LD values. Only the upper
    triangle of the array (above the diagonal) is used.

    `names` is an iterable of names, which should be the same length as the
    number of rows in `array`.

    `labels`, if present, is a dictionary of labels with available
    keys "mid", "left", and "right", for example
    dict(mid='chr17, left=12398013, right=18290324).

    If there are no labels then no label bar is drawn.

    `figsize` is the size of the figure in inches (width, height).

    `cmap` is the name of a Matplotlib colormap.

    `colorbar` controls whether to add a colorbar. The default is True.
    """

    if labels:
        fig, (lab_ax,ax) = plt.subplots(2,
                                        figsize = figsize,
                                        gridspec_kw = dict(
                                            height_ratios=(1,10),
                                            hspace=0))
    else:
        fig, ax = plt.subplots(figsize=figsize)

    # how many items
    lds = array.shape[0]

    # mask out lower triangle and diagonal of the array
    mask =  numpy.tri(lds)
    array = numpy.ma.array(array, mask=mask)

    # set the colormap
    cmap = matplotlib.cm.get_cmap(cmap)
    cmap.set_bad('w') # so masked values are white

    # draw the actual block grid and rotate it as needed.
    im = ax.imshow(array,
                   cmap=cmap,
                   transform = (transforms.Affine2D().rotate_deg(-45)
                                + ax.transData))

    # set view limits so we see the necessary part of the grid.
    s2 = math.sqrt(2)
    ax.set_xlim(0, s2*(lds-1))
    ax.set_ylim(-s2*lds/2, 1)

    # add labels for each item, changing the view limits to fit.
    for i,n in enumerate(names):
        fit_text(ax, i*s2, 0, n,
                 rotation='vertical',
                 horizontalalignment='center',
                 verticalalignment='bottom')

    # don't draw axes,
    ax.set_axis_off()

    # draw the (optional) color bar.
    if colorbar:
        fig.colorbar(im, ax=ax, shrink=0.5)

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

    If labelling of the most significant SNPs is required, (controlled
    by `alpha` and `good_label_column` parameters), the DataFrame must
    also have a column labelled by the `good_label_column` parameter.

    `x_label` is a label for the x axis (such as "Chromosome 17").

    `alpha` is the threshold for distinguishing "good" SNPs (those
    with pvalue less than alpha). If None, no threshold is marked and
    "good" SNPs are not distinguished.

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

    `good_label_column` is the column name in `df` for labels to be
    drawn for "good" SNPs, or None for no labels.

    `good_label_rotation` is the rotation for labels for "good" SNPs.
    """

    # calculate -log_10(P) to plot

    df['minuslog10pvalue'] = -numpy.log10(df.pvalue)

     # x axis in megabases, to suppress annoying useOffset nonsense
     # in the MPL tick formatter.
    df['positionMb'] = df.position / 1e6

    fig, ax = plt.subplots(figsize=figsize)

    # plot the main scatter
    ax.scatter(df['positionMb'],df['minuslog10pvalue'],
               s=size, color=color, marker=marker)

    # if we have a threshold (alpha), draw a horizontal line and
    # redraw the points above it.
    if alpha:
        ax.axhline(-numpy.log10(alpha),
                   linestyle=alpha_line_style,
                   linewidth=alpha_line_width,
                   color=alpha_line_color)

        # identify and redraw the points above the line
        index_df = df[df.pvalue < alpha]
        ax.scatter(index_df['positionMb'], index_df['minuslog10pvalue'],
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
          pos_top = False,
          # future work:
          genes = None,
          gene_height = 0.5,
          ):
    df['positionMb'] = df.position / 1e6
    total_height = 1 + track_height + gene_height
    if tracks and genes:
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

    elif tracks or genes:
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

    posax.scatter(df.positionMb, df.posterior,
                  s=size, color=color, marker=marker)

    posax.set_ylabel('Posterior probability')
    bottomax.set_xlabel(f'position (Mbp)')

    if threshold:
        posax.axhline(threshold,
                      linestyle=alpha_line_style,
                      linewidth=alpha_line_width,
                      color=alpha_line_color)

        cred_df = df[df.posterior > threshold]
        posax.scatter(cred_df.positionMb, cred_df.posterior,
                      s=good_size, color=good_color, marker=good_marker)
        if good_label_column:
            for i, row in cred_df.iterrows():
                fit_text(posax,row.positionMb, row.posterior+0.01,
                         row[good_label_column],
                         rotation=good_label_rotation,
                         horizontalalignment="left",
                         verticalalignment='bottom')

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
                    trackax.axvline(p, zorder=-1, linewidth=0.5, color='k')
            tracks_array = [cred_df[cred_df[track_column].str.contains(t)].positionMb
                            for t in tracks]
            trackax.eventplot(tracks_array, linelengths=track_good_linelength, colors=track_colors)

        trackax.set_yticks(range(len(tracks)))
        trackax.set_yticklabels(tracks)
        trackax.set_ylabel('')
        trackax.set_ylim(-0.5, len(tracks)-0.5)

    if genes:
        geneax.set_ylabel('genes')
        geneax.spines['top'].set_visible(False)
        geneax.spines['left'].set_visible(False)
        geneax.spines['right'].set_visible(False)
        geneax.tick_params(which='both', left=False, labelleft=False)

    return fig

def test():
    fig, ax = plt.subplots()
    ax.eventplot(((1,5,7,8),(3,4,1,6),(5,2,8)))
    fit_text(ax, 3, 2, "Spong")
    fit_text(ax, 2, -3, "Wibble", rotation=-45)
    fit_text(ax, 1, 5, "Foobar", verticalalignment='bottom', rotation='vertical')

    return fig
