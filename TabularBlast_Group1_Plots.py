#!/usr/bin/env python

'''Plots from folder of filtered tabular blasts.

Builds boxplots and hexplots.

File naming scheme and Sample names need to be modified to fit the
current project. 

Modify lines: 75, 116-, 140-, 171, 178-, 193

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 04 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from os import listdir
from os.path import isfile, join
import matplotlib
from matplotlib.patches import Patch
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def func1(fasta, blasts, outpre):
    "reads files builds plots"

    # Get contig lengths from reference fasta
    lengths = {}
    data = {
            'Sample': [], 'Location': [], 'Depth': [],
            'pident': [], 'position': []
            }

    with open(fasta, 'r') as file:
        contig_length = 0
        for name, seq in read_fasta(file):
            lengths[name[1:]] = contig_length
            contig_length += len(seq)

    # Grab list of files
    file_list = [f for f in listdir(blasts) if isfile(join(blasts, f))]
    # remove .DS_Store file for stupid MAC OS
    if '.DS_Store' in file_list: file_list.remove('.DS_Store')
    # Keep a sample list
    Sample_List = []
    depths = {}
    # Read through files and populate the data dict
    for file in file_list:
        # Retreive sample name from file name
        X = file.split('_')
        if X[0] == 'Pacific': X[0] = 'Pac'
        Sample = '-'.join(X[:2]) # adjust to filename
        Sample_List.append(Sample)
        loca = X[0]
        depth = X[1]
        depths[Sample] = depth #int(depth.split('m')[0])

        with open(f'{blasts}/{file}', 'r') as f:
            head = f.readline()
            for line in f:
                X = line.rstrip().split('\t')
                sid = X[1]
                pid = float(X[2])
                pos = int(X[8])
                
                len_adjust = lengths[sid]
                adjusted_pos = pos + len_adjust

                data['Sample'].append(Sample)
                data['Depth'].append(depth)
                data['Location'].append(loca)
                data['pident'].append(pid)
                data['position'].append(adjusted_pos)

    # Save data as tsv
    df = pd.DataFrame(data)
    df.to_csv(f'{outpre}_data.tsv', index=False, sep='\t')


    # Build Plots

    # Fancy Boxen Plot
    '''
    colors = [
                '#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6',
                '#4292c6','#2171b5','#08519c','#08306b', '#023858'
                ]
    '''
    #colors = [ '#6baed6', '#2171b5']

    Cgom = '#c6dbef'
    Cpac = '#2171b5'

    colors = {
                    "GoM-150m": Cgom,
                    "Pac-200m": Cpac,
                    "GoM-300m": Cgom,
                    "Pac-500m": Cpac,
                    "GoM-600m": Cgom,
                    "Pac-1000m": Cpac,
                    "GoM-1000m": Cgom,
                    "GoM-1470m": Cgom,
                    "GoM-2107m": Cgom,
                    "Pac-4580m": Cpac,
                    "Pac-5100m": Cpac,
                    "Pac-5601m": Cpac
                    }

    #sns.set_palette(sns.color_palette(colors))

    plotsize = (5,3)

    # GOM BOX PLOT
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=plotsize)
    sns.boxenplot(
            x='pident', y='Sample', data=df,
            order=[
                    "GoM-150m",
                    "GoM-300m",
                    "GoM-600m",
                    "GoM-1000m",
                    "GoM-1470m",
                    "GoM-2107m"
                    ],
            width=.8,
            #hue='Location',
            palette=colors,
            ax=ax
            )

    #plt.subplots(figsize=(8,5), dpi=300)
    #g.fig.set_figwidth(8)
    #g.fig.set_figheight(5)
    '''
    GoM = Patch(facecolor='#6baed6', label='GoM')
    Pac = Patch(facecolor='#2171b5', label='Pac')

    plt.legend(
                bbox_to_anchor=(0.5, 1.2),
                loc='upper center',
                ncol=2,
                title='Location',
                handles=[GoM, Pac]
                )
    '''

    plt.tight_layout()
    plt.savefig(f'{outpre}-GOM_boxen_plot.png', dpi=300)
    plt.close() 

    # PACIFIC BOX PLOT
    fig, ax = plt.subplots(figsize=plotsize)
    sns.boxenplot(
            x='pident', y='Sample', data=df,
            order=[
                    "Pac-200m",
                    "Pac-500m",
                    "Pac-1000m",
                    "Pac-4580m",
                    "Pac-5100m",
                    "Pac-5601m"
                    ],
            width=.8,
            #hue='Location',
            palette=colors,
            ax=ax
            )

    plt.tight_layout()
    plt.savefig(f'{outpre}-PAC_boxen_plot.png', dpi=300)
    plt.close() 


    # Hex style rec plot
    ANIrlist = open(f'{outpre}_Sample_ANIr.tsv', 'w')
    ANIrlist.write('Sample\tDepth\tANIr\n')
    sns.set_style("white")
    for samp in Sample_List:
        subdf = df[df['Sample'] == samp]
        depth = depths[samp]
        anir = subdf.pident.mean()
        ANIrlist.write(f'{samp}\t{depth}\t{anir}\n')

        xmax = subdf.position.max()

        h = sns.jointplot(
                data=subdf, x="position", y="pident",
                ylim=(80,100), xlim=(0,xmax),
                kind="hex", color="#252525", gridsize=(20,20)
                )

        h.ax_joint.axhline(y=anir, ls='--', color='#d7191c', lw=3)

        h.ax_joint.text(
                xmax,
                80.2,
                #anir - 3,
                f'{anir:.2f}%',
                fontsize=22, color='#d7191c', horizontalalignment='right'
                )

        #plt.subplots(figsize=(2,4), dpi=300)
        h.fig.set_figwidth(2)
        h.fig.set_figheight(4)
        plt.subplots_adjust(left=0.25, bottom=0.1, right=0.9995, top=0.995)
        plt.savefig(f'{outpre}_{samp}_hex_plot.png', dpi=300)
        plt.close() 



def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f', '--input_fasta_file',
        help='Please specify the reference fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-b', '--input_blast_dir',
        help='Please specify the blast directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_pre',
        help='What prefix do you want to use for the output files?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    func1(args['input_fasta_file'], args['input_blast_dir'], args['out_pre'])


if __name__ == "__main__":
    main()
