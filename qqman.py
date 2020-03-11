import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def qqplot(pvector, filename=None, size=1, pointcolor='black', title=None, linecolor='red'):
    pvector = pvector[~np.isnan(pvector)]
    pvector = pvector[~((pvector >= 1) | (pvector <= 0))]
    pvector.sort()
    o = -(np.log10(pvector))
    e = -(np.log10((np.array(range(1, (len(o) + 1))) - .5) / len(o)))
    fig, ax = plt.subplots(1, 1, figsize=(8 * int(size), 8 * int(size)), dpi=300*(int(size)**2), facecolor='w')
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    maxi = max(np.amax(e), np.amax(o))
    ax.set_xlim([0, maxi + 0.1])
    ax.set_ylim([0, maxi + 0.1])
    ax.set_ylabel('Observed -log10(' + r'$p$' + ')')
    ax.set_xlabel('Expected -log10(' + r'$p$' + ')')
    ax.scatter(e, o, c=pointcolor)
    ax.plot((0, maxi), (0, maxi), linecolor)
    if title is not None:
        ax.set_title(title)
    if filename is not None:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        return fig, ax


def manhattan(df, filename="", sigp=5e-8, sigcolor='black', sugp=1e-5, sugcolor='black',
              pointcolor=['midnightblue', 'goldenrod'], size=1, highlight=[],
              highlightcolor=['orange'], title=None, rainbow=False):
    df.columns = map(str.lower, df.columns)
    if rainbow:
        pointcolor = ['#FF0000', '#FF4000', '#FF8000', '#FFBF00', '#FFFF00', '#BFFF00', '#80FF00', '#40FF00', '#00FF00',
                      '#00FF40', '#00FF80', '#00FFBF', '#00FFFF', '#00BFFF', '#0080FF', '#0040FF', '#0000FF', '#4000FF',
                      '#8000FF', '#BF00FF', '#FF00FF', '#FF00BF']
        highlightcolor = list(reversed(pointcolor))
    df['logp'] = -(np.log10(df['p']))
    df['chromosome'] = df.chr.astype('category')
    df.chromosome = df.chromosome.cat.set_categories(['ch-%i' % i for i in range(23)], ordered=True)
    df = df.sort_values('chr')
    df['bpadd'] = 0
    add = 0
    for c in range(2, 24):
        add += df.loc[df['chr'] == (c - 1), 'bp'].max()
        df.loc[df['chr'] == c, 'bpadd'] = add
    df['ind'] = df['bp'] + df['bpadd']
    if len(highlight) > 0:
        dfhighlight = df[df['rsid'].isin(highlight)]
        dfhighlight_grouped = dfhighlight.groupby('chr')
    df_grouped = df.groupby('chr')
    fig, ax = plt.subplots(1, 1, figsize=(12 * float(size), 6 * float(size)), dpi=300*(float(size)**2), facecolor='w')
    colors = pointcolor
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='logp', color=colors[num % len(colors)],
                   edgecolor=colors[num % len(colors)], ax=ax)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].mean()))
    if len(highlight) > 0:
        for num, (name, group) in enumerate(dfhighlight_grouped):
            group.plot(kind='scatter', x='ind', y='logp', color=highlightcolor[num % len(highlightcolor)],
                       edgecolor=highlightcolor[num % len(highlightcolor)], ax=ax, marker='^', s=45 * size)
    if sigp > 0:
        ax.plot((0, (df['ind'].max() * 1.005)), (-(np.log10(float(sigp))), -(np.log10(float(sigp)))), ls='--', lw=1.2, color=sigcolor)
    if sugp > 0:
        ax.plot((0, (df['ind'].max() * 1.005)), (-(np.log10(float(sugp))), -(np.log10(float(sugp)))), ls='--', lw=0.5, color=sugcolor)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, (df['ind'].max() * 1.005)])
    ax.set_ylim([0, (np.amax(df['logp']) * 1.3)])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(' + r'$p$' + ')')
    if title is not None:
        ax.set_title(title)
    if filename is not None:
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    else:
        return fig, ax


if __name__ == '__main__':
    import argparse
    import json

    parser = argparse.ArgumentParser('QQ and Manhattan plotter. Plots both by default. All column names are case insensitive.')
    parser.add_argument('--csv', help='csv file with summary statistics (supports gzipped files if extension is .gz)', default=None)
    parser.add_argument('--tsv', help='tab separated file with summary statistics (supports gzipped files if extension is .gz)', default=None)
    parser.add_argument('--man', help='Only plot Manhattan plot', default=False, action='store_true')
    parser.add_argument('--qq', help='Only plot QQ-plot', default=False, action='store_true')
    parser.add_argument('--highlight', help='SNPs to highlight. Either comma seperated RSIDs, a txt-file with 1 SNP per line, \'sig\' for significant SNPs (default), or \'none\'', default='sig')
    parser.add_argument('--qq_name', help='Filename for qqplot (default qqplot.png)', default='qqplot.png')
    parser.add_argument('--man_name', help='Filename for manhattan plot (default manhattan.png)', default='manhattan.png')
    parser.add_argument('--col_chr', help='Column name of Chromosome (default \'chr\')', default='chr')
    parser.add_argument('--col_bp', help='Column name of basepair position (default \'bp\')', default='bp')
    parser.add_argument('--col_p', help='Column name of p-value (default \'p\')', default='p')
    parser.add_argument('--col_rsid', help='Column name of RSID (default \'rsid\')', default='rsid')
    parser.add_argument('--options', help='JSON formatted other options (use --adv_help to see options)', default=None)
    parser.add_argument('--adv_help', help='See advanced options', default=False, action='store_true')
    args = parser.parse_args()

    options = {'qq': {'size': 1, 'pointcolor': 'black', 'linecolor': 'red', 'title': None},
               'man': {'sigp': 5e-8, 'sigcolor': 'black', 'sugp': 1e-5, 'sugcolor': 'black',
                       'pointcolor': ['midnightblue', 'goldenrod'], 'size': 1, 'highlightcolor': ['orange'],
                       'title': None, 'rainbow': False}}
    options_help = {'qq': {'size': 'Relative figure size.', 'title': 'Main in-figure title.',
                           'pointcolor': 'Color of plotted points.', 'linecolor': 'Color of plotted lines.'},
                    'man': {'sigp': 'P-value cutoff for significance. (set to negative value to remove)',
                            'sigcolor': 'Color of significance cutoff line.',
                            'sugp': 'P-value cutoff for suggested. (set to negative value to remove)',
                            'sugcolor': 'Color of suggested line.',
                            'pointcolor': 'List of colors (length >=1) to use for points, will loop over list per chromosome.',
                            'size': 'Relative figure size.',
                            'highlightcolor': 'List of colors (length >=1) to use for highlighting, will loop over list per chromosome.',
                            'title': 'Main figure title.',
                            'rainbow': 'Override existing colors and use rainbow instead.'
    }}
    if args.adv_help:
        print('Advanced help:')
        for p in ['qq', 'man']:
            print(p+':')
            for k, v in options_help[p].items():
                print('  {}:'.format(k))
                print('        {} (default: {}) '.format(v, options[p][k]))
            print('')
        print('')
        print('Example usage:')
        print('1: --options "{"qq":{"size": "2", "title":"My QQ-plot title"}, "man": {"title":"My Rainbowplot", "rainbow": true}}"')
        print('2: --options "{"qq":{}, "man": {"title":"My red&blue plot", "pointcolor": ["red", "blue"]}}"')
        print('Note both qq and man should be the first level keys, even though you may not specify custom options for one of them.')
        exit()

    if args.options is not None:
        custom_options = json.loads(args.options)
        for p in ['qq', 'man']:
            for k, v in custom_options[p].items():
                options[p][k] = v

    if ((args.csv is None) and (args.tsv is None)) or ((args.csv is not None) and (args.tsv is not None)):
        raise argparse.ArgumentError("Please supply either csv or tsv file.")

    if args.csv is not None:
        data = pd.read_csv(args.csv)
    else:
        data = pd.read_csv(args.tsv, sep='\t')

    data.columns = [x.lower() for x in data.columns]
    qq = args.qq
    man = args.man
    if (not qq) and (not man):
        qq, man = True, True

    data.rename(columns={args.col_p.lower(): 'p'}, inplace=True)
    if qq:
        qqplot(data['p']. values, filename=args.qq_name, **options['qq'])

    if man:
        if args.highlight.lower() == 'sig':
            data.rename(columns={args.col_rsid.lower(): 'rsid'}, inplace=True)
            highlights = data.loc[data['p'] < options['man']['sigp'], "rsid"].tolist()
        elif args.highlight.lower() == 'none':
            highlights = []
        elif '.' in args.highlight:
            with open(args.highlight, 'r') as f:
                highlights = f.readlines()
        else:
            data.rename(columns={args.col_rsid.lower(): 'rsid'}, inplace=True)
            highlights = args.highlight.split(',')

        data.rename(columns={args.col_chr.lower(): 'chr', args.col_bp.lower(): 'bp'}, inplace=True)
        manhattan(data, highlight=highlights, filename=args.man_name, **options['man'])



data = pd.DataFrame()
data.to_csv('test_{}.csv'.format(len(data)))
for m in [.49, .45, .4, .35, .3, .25, .2, .15, .1, .05, .01, .001]:
    print(m)
    data2 = data.loc[data['MAF'] >= m, :]
    data2.to_csv('test_{}.csv'.format(len(data2)))