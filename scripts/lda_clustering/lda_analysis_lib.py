#Author: Timothy Durham (c) 2018

from collections import Counter
import dill
import glob
import gzip
import igraph as ig
import itertools
import leidenalg
#import magic
import matplotlib
from matplotlib import pyplot
import numba
import numpy
import os
import pickle
from plumbum import local
import random
import re
import scipy
from scipy.cluster import hierarchy
import scipy.sparse as sps
from scipy.spatial import distance
import scipy.stats as stats
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD
from sklearn import neighbors
from sklearn import metrics
import sys
import umap

def find_nearest_genes(peak_files, out_subdir, refseq_exon_bed):
    #get unix utilities
    bedtools, sort, cut, uniq, awk = local['bedtools'], local['sort'], local['cut'], local['uniq'], local['awk']

    #process the peak files to find nearest genes
    nearest_genes = []
    for path in sorted(peak_files):
        out_path = os.path.join(out_subdir, os.path.basename(path).replace('.bed', '.nearest_genes.txt'))
        cmd = (bedtools['closest', '-D', 'b', '-io', '-id', '-a', path, '-b', refseq_exon_bed] |
         cut['-f1,2,3,5,9,12'] | #fields are chrom, start, stop, peak sum, gene name, distance
         awk['BEGIN{OFS="\t"}{if($6 > -1200){print($1, $2, $3, $6, $5, $4);}}'] |
         sort['-k5,5', '-k6,6nr'] |
         cut['-f5,6'])()
        with open(out_path, 'w') as out:
            prev_gene = None
            for idx, line in enumerate(str(cmd).strip().split('\n')):
                if prev_gene is None or not line.startswith(prev_gene):
#                    print(line)
                    line_split = line.strip().split()
                    prev_gene = line_split[0]
                    out.write(line + '\n')
        nearest_genes.append(out_path)
    return nearest_genes

def load_expr_db(db_path):
    if os.path.basename(db_path) == 'RepAvgGeneTPM.csv':
        with open(db_path) as lines_in:
            db_headers = lines_in.readline().strip().split(',')[1:]
        db_vals = numpy.loadtxt(db_path, delimiter=',', skiprows=1, dtype=object)[:,1:]
    elif os.path.basename(db_path) == 'gexplore_l2_tissue_expr.jonathan_updated.txt.gz':
        db_headers, db_vals = load_expr_db2(db_path)
    else:
        with open(db_path) as lines_in:
            db_headers = lines_in.readline().strip().split('\t')
        db_vals = numpy.loadtxt(db_path, delimiter='\t', skiprows=1, dtype=object)
    print('Loaded DB shape: {!s}'.format(db_vals.shape))
    return (db_headers, db_vals)

def load_expr_db2(db_path):
    db_vals = numpy.loadtxt(db_path, delimiter='\t', skiprows=1, dtype=object)
    gene_set = sorted(set(db_vals[:,0]))
    tissue_set = sorted(set(db_vals[:,2]))
    db_data = numpy.zeros((len(gene_set), len(tissue_set)))
    gene_idx = None
    gene_name = ''
    gene_vals = []
    tissue_idx = []
    for idx in range(db_vals.shape[0]):
        db_elt = db_vals[idx]
        #rule out any genes with 95% CI lower bound of zero
        if float(db_elt[6]) <= 0:
            continue
        if db_elt[0] != gene_name:
            if gene_idx is not None:
                db_data[gene_idx,tissue_idx] = gene_vals
            gene_name = db_elt[0]
            gene_idx = gene_set.index(gene_name)
            gene_vals = []
            tissue_idx = []
        #use the bootstrap median value
        gene_vals.append(float(db_elt[5]))
        tissue_idx.append(tissue_set.index(db_elt[2]))
    return (['gene_name'] + tissue_set, 
            numpy.hstack([numpy.array(gene_set, dtype=object)[:,None], db_data.astype(object)]))

TOPN=500
def get_gene_data(genes_path, gene_expr_db, topn=TOPN):
    if isinstance(genes_path, list):
        genes_list = genes_path
    else:
        with open(genes_path) as lines_in:
            genes_list = [elt.strip().split()[:2] for elt in lines_in]
    gene_idx = [(numpy.where(gene_expr_db[:,0] == elt[0])[0],elt[1]) for elt in genes_list]
    gene_idx_sorted = sorted(gene_idx, key=lambda x:float(x[1]), reverse=True)
    gene_idx, gene_weights = zip(*[elt for elt in gene_idx_sorted if len(elt[0]) > 0][:topn])
    gene_idx = [elt[0] for elt in gene_idx]
    gene_data = gene_expr_db[:,1:].astype(float)[gene_idx,:]
    denom = numpy.sum(gene_data, axis=1)[:,None] + 1e-8
    gene_norm = gene_data/denom
    return gene_idx, gene_data, gene_norm, len(genes_list), numpy.array(gene_weights, dtype=float)

def sample_db(data_norm, expr_db, data_weights=None, nsamples=1000):
    samples = []
    rs = numpy.random.RandomState(15321)
    random_subset = numpy.arange(expr_db.shape[0])
    num_to_select = data_norm.shape[0]
    for idx in range(nsamples):
        rs.shuffle(random_subset)
        if expr_db.shape[1] - data_norm.shape[1] == 1:
            db_subset = expr_db[random_subset[:num_to_select]][:,1:].astype(float)
        else:
            db_subset = expr_db[random_subset[:num_to_select]].astype(float)
        denom = numpy.sum(db_subset, axis=1)[:None] + 1e-8
        db_subset_norm = numpy.mean((db_subset.T/denom).T, axis=0)
        if data_weights is not None:
            avg_data_norm = numpy.divide(numpy.average(data_norm, axis=0, weights=data_weights) + 1e-5, 
                                         db_subset_norm + 1e-5)
#            avg_data_norm = numpy.divide(numpy.average(data_norm, axis=0, weights=gene_weights), db_subset_norm,
#                                         out=numpy.zeros_like(db_subset_norm), where=(db_subset_norm != 0))
#            samples.append(numpy.log2(numpy.average(data_norm, axis=0, weights=gene_weights)/db_subset_norm))
        else:
            avg_data_norm = numpy.divide(numpy.average(data_norm, axis=0, weights=None) + 1e-5, 
                                         db_subset_norm + 1e-5)
#            avg_data_norm = numpy.divide(numpy.average(data_norm, axis=0, weights=None), db_subset_norm,
#                                         out=numpy.zeros_like(db_subset_norm), where=(db_subset_norm != 0))
#            samples.append(numpy.log2(numpy.average(data_norm, axis=0, weights=None)/db_subset_norm))
        samples.append(numpy.log2(avg_data_norm))
#        samples.append(numpy.log2(avg_data_norm, out=numpy.zeros_like(avg_data_norm), where=(avg_data_norm != 0)))
    samples = numpy.vstack(samples)
    samples_mean = numpy.mean(samples, axis=0)
    samples_sem = stats.sem(samples, axis=0)
    conf_int = numpy.array([stats.t.interval(0.95, samples.shape[0]-1, 
                                             loc=samples_mean[idx], scale=samples_sem[idx])
                            for idx in range(samples.shape[1])]).T
    conf_int[0] = samples_mean - conf_int[0]
    conf_int[1] = conf_int[1] - samples_mean
    return samples_mean, conf_int

def plot_l2_tissues(nearest_genes_glob, refdata, allgenes_path=None, expr_db=None, expr_db_headers=None, ncols=3, 
                    topn=TOPN, weights=False, nsamples=100, savefile=None, display_in_notebook=True, dpi=None,
                    figsize_factor=7):
    if expr_db is None:
        #Get all L2 tissue expression data to normalize the distribution of genes from peaks
        l2_tissue_db_path = os.path.join(refdata,'gexplore_l2_tissue_expr.txt')
        expr_db_headers, expr_db = load_expr_db(l2_tissue_db_path)
    
    gene_lists = glob.glob(nearest_genes_glob)
    if os.path.basename(gene_lists[0]).startswith('peaks'):
        gene_lists.sort(key=lambda x:int(os.path.basename(x).split('.')[0].replace('peaks', '')))
    elif os.path.basename(gene_lists[0]).startswith('topic'):
        gene_lists.sort(key=lambda x:int(os.path.basename(x).split('.')[1].replace('rank', '')))
    else:
        gene_lists.sort(key=lambda x:os.path.basename(x).split('.')[0])
    gene_list_data = [(os.path.basename(path).split('.')[0], get_gene_data(path, expr_db, topn=topn)) for path in gene_lists]
    print('\n'.join(['{!s} nearest genes: found {!s} out of {!s} total'.format(fname, data.shape[0], gene_list_len)
                    for (fname, (data_idx, data, data_norm, gene_list_len, gene_weights)) in gene_list_data]))
    
    if allgenes_path is not None:
        allgenes_norm = get_gene_data(allgenes_path, expr_db, topn=100000)[2]

    l2_tissue_colors = [('Body wall muscle', '#e51a1e'),
                        ('Intestinal/rectal muscle', '#e51a1e'),
                        ('Pharyngeal muscle', '#377db8'),
                        ('Pharyngeal epithelia', '#377db8'),
                        ('Pharyngeal gland', '#377db8'),
                        ('Seam cells', '#4eae4a'),
                        ('Non-seam hypodermis', '#4eae4a'),
                        ('Rectum', '#4eae4a'),
                        ('Ciliated sensory neurons', '#984ea3'),
                        ('Oxygen sensory neurons', '#984ea3'),
                        ('Touch receptor neurons', '#984ea3'),
                        ('Cholinergic neurons', '#984ea3'),
                        ('GABAergic neurons', '#984ea3'),
                        ('Pharyngeal neurons', '#984ea3'),
                        ('flp-1(+) interneurons', '#984ea3'),
                        ('Other interneurons', '#984ea3'),
                        ('Canal associated neurons', '#984ea3'),
                        ('Am/PH sheath cells', '#ff8000'),
                        ('Socket cells', '#ff8000'),
                        ('Excretory cells', '#ff8000'),
                        ('Intestine', '#fcd800'),
                        ('Germline', '#f97fc0'),
                        ('Somatic gonad precursors', '#f97fc0'),
                        ('Distal tip cells', '#f97fc0'),
                        ('Vulval precursors', '#f97fc0'),
                        ('Sex myoblasts', '#f97fc0'),
                        ('Coelomocytes', '#a75629')]
    idx_by_color = {}
    for idx, (name, color) in enumerate(l2_tissue_colors):
        try:
            idx_by_color[color][1].append(idx)
        except KeyError:
            idx_by_color[color] = [name, [idx]]
            
#    rs = numpy.random.RandomState(15321)
#    random_subset = numpy.arange(expr_db.shape[0])
#    rs.shuffle(random_subset)
#    #num_to_select = int(numpy.mean([neuron_data.shape[0], emb_muscle_data.shape[0], l2_muscle_data.shape[0]]))
#    num_to_select = len(random_subset)
#    l2_tissue_db_subset = expr_db[random_subset[:num_to_select]][:,1:].astype(float)
#    denom = numpy.sum(l2_tissue_db_subset, axis=1)[:,None] + 1e-8
#    l2_tissue_db_norm = numpy.mean(l2_tissue_db_subset/denom, axis=0)
    print('Tissue DB norm shape: {!s}'.format(expr_db.shape))

    pyplot.rcParams.update({'xtick.labelsize':14,
                            'ytick.labelsize':14,
                            'xtick.major.pad':8})

    ind = numpy.arange(len(expr_db_headers) - 1)
    width = 0.66
    axis_fontsize = 18
    title_fontsize = 19
    nrows = int(numpy.ceil(len(gene_list_data)/float(ncols)))
    if dpi is not None:
        fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, 
                                    figsize=(figsize_factor * ncols, figsize_factor * nrows), 
                                    sharey=True, dpi=dpi)
    else:
        fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, 
                                    figsize=(figsize_factor * ncols, figsize_factor * nrows), 
                                    sharey=True)
    for idx, (fname, (data_idx, data, data_norm, gene_list_len, gene_weights)) in enumerate(gene_list_data):
        ax_idx = (idx//ncols, idx%ncols) if nrows > 1 else idx
#        to_plot = numpy.log2(numpy.mean(data_norm, axis=0)/l2_tissue_db_norm)
#        import pdb; pdb.set_trace()
        if weights is True:
#            to_plot = numpy.log2(numpy.average(data_norm, axis=0, weights=gene_weights)/l2_tissue_db_norm)
            to_plot, errs = sample_db(data_norm, expr_db if allgenes_path is None else allgenes_norm, 
                                      data_weights=gene_weights, nsamples=nsamples)
        else:
#            to_plot = numpy.log2(numpy.average(data_norm, axis=0, weights=None)/l2_tissue_db_norm)
            to_plot, errs = sample_db(data_norm, expr_db if allgenes_path is None else allgenes_norm, 
                                      data_weights=None, nsamples=nsamples)
        for idx2, (name, color) in enumerate(l2_tissue_colors):
            axes[ax_idx].bar(ind[idx2], to_plot[idx2], width, yerr=errs[:,idx2][:,None], color=color, label=name)
        axes[ax_idx].axhline(0, color='k')
        axes[ax_idx].set_xlim((-1, len(expr_db_headers)))
        axes[ax_idx].set_title('{!s}\n({!s} genes)\n'.format(fname, data.shape[0]), fontsize=title_fontsize)
        axes[ax_idx].set_ylabel('Log2 ratio of mean expr proportion\n(ATAC targets:Random genes)', fontsize=axis_fontsize)
        axes[ax_idx].set_xlabel('L2 tissues', fontsize=axis_fontsize)
        axes[ax_idx].set_xticks(ind + width/2)
        axes[ax_idx].set_xticklabels([])
    #axes[0].set_xticklabels(expr_db_headers[1:], rotation=90)
    if nrows > 1:
        axes[0,ncols-1].legend(bbox_to_anchor=[1.0,1.0])
    else:
        axes[-1].legend(bbox_to_anchor=[1.0,1.0])

    if display_in_notebook is True:
        fig.tight_layout()
    if savefile is not None:
        fig.savefig(savefile, bbox_inches='tight')

def plot_stages(nearest_genes_glob, refdata, expr_db=None, expr_db_headers=None, ncols=3, topn=TOPN, weights=False):
    if expr_db is None:
        #Get all stages expression data to normalize the distribution of genes from peaks
        stage_db_path = os.path.join(refdata,'gexplore_stage_expr.txt')
        expr_db_headers, expr_db = load_expr_db(stage_db_path)

    gene_lists = glob.glob(nearest_genes_glob)
    if os.path.basename(gene_lists[0]).startswith('peaks'):
        gene_lists.sort(key=lambda x:int(os.path.basename(x).split('.')[0].replace('peaks', '')))
    elif os.path.basename(gene_lists[0]).startswith('topic'):
        gene_lists.sort(key=lambda x:int(os.path.basename(x).split('.')[1].replace('rank', '')))
    else:
        gene_lists.sort(key=lambda x:os.path.basename(x).split('.')[0])
    gene_list_data = [(os.path.basename(path).split('.')[0], get_gene_data(path, expr_db, topn=topn)) for path in gene_lists]
    print('\n'.join(['{!s} nearest genes: found {!s} out of {!s} total'.format(fname, data.shape[0], gene_list_len)
                    for (fname, (data_idx, data, data_norm, gene_list_len, gene_weights)) in gene_list_data]))
    
    rs = numpy.random.RandomState(15321)
    random_subset = numpy.arange(expr_db.shape[0])
    rs.shuffle(random_subset)
    #num_to_select = int(numpy.mean([neuron_data.shape[0], emb_muscle_data.shape[0], l2_muscle_data.shape[0]]))
    num_to_select = len(random_subset)
    stage_db_subset = expr_db[random_subset[:num_to_select]][:,1:].astype(float)
    denom = numpy.sum(stage_db_subset, axis=1)[:,None] + 1e-8
    stage_db_norm = numpy.mean(stage_db_subset/denom, axis=0)
    print('Stage DB norm shape: {!s}'.format(stage_db_norm.shape))

    emb_idx = [expr_db_headers[1:].index(elt) for elt in expr_db_headers[1:] 
               if elt.endswith('m') or elt == '4-cell']
    larva_idx = [expr_db_headers[1:].index(elt) for elt in expr_db_headers[1:] 
                 if elt.startswith('L')]
    adult_idx = [expr_db_headers[1:].index(elt) for elt in expr_db_headers[1:]
                if 'adult' in elt]
    dauer_idx = [expr_db_headers[1:].index(elt) for elt in expr_db_headers[1:]
                if 'dauer' in elt]
#    rest_idx = [expr_db_headers[1:].index(elt) for elt in expr_db_headers[1:] 
#                if not elt.endswith('m') and not elt.startswith('L') and elt != '4-cell']

    pyplot.rcParams.update({'xtick.labelsize':20,
                            'ytick.labelsize':20,
                            'xtick.major.pad':8})

    ind = numpy.arange(len(expr_db_headers) - 1)
    width = 0.66
    axis_fontsize = 25
    title_fontsize = 27
    nrows = int(numpy.ceil(len(gene_list_data)/float(ncols)))
    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(7 * ncols, 7 * nrows), sharey=True)
    for idx, (fname, (data_idx, data, data_norm, gene_list_len, gene_weights)) in enumerate(gene_list_data):
        ax_idx = (idx//ncols, idx%ncols) if nrows > 1 else idx
#        to_plot = numpy.log2(numpy.mean(data_norm, axis=0)/stage_db_norm)
        if weights is True:
            to_plot = numpy.log2(numpy.average(data_norm, axis=0, weights=gene_weights)/stage_db_norm)
        else:
            to_plot = numpy.log2(numpy.average(data_norm, axis=0, weights=None)/stage_db_norm)
        axes[ax_idx].bar(ind[emb_idx], to_plot[emb_idx], width, color='orange', label='Embryo')
        axes[ax_idx].bar(ind[larva_idx], to_plot[larva_idx], width, color='blue', label='Larva')
        axes[ax_idx].bar(ind[adult_idx], to_plot[adult_idx], width, color='red', label='Adult')
        axes[ax_idx].bar(ind[dauer_idx], to_plot[dauer_idx], width, color='green', label='Dauer')
#        axes[ax_idx].bar(ind[rest_idx], to_plot[rest_idx], width, color='grey', label='Other')
        axes[ax_idx].axhline(0, color='k')
        axes[ax_idx].set_xlim((-1, len(expr_db_headers)))
        axes[ax_idx].set_title('{!s}\n({!s} genes)\n'.format(fname, data.shape[0]), fontsize=title_fontsize)
        axes[ax_idx].set_ylabel('Log2 Ratio of Mean Expr Proportion\n(ATAC Targets:All Genes)', fontsize=axis_fontsize)
        axes[ax_idx].set_xlabel('Developmental Stage', fontsize=axis_fontsize)
        axes[ax_idx].set_xticks(ind + width/2)
        axes[ax_idx].set_xticklabels([])

    fig.tight_layout()

def leiden_clustering(umap_res, resolution_range=(0,1), random_state=2, kdtree_dist='euclidean'):
    tree = neighbors.KDTree(umap_res, metric=kdtree_dist)
    vals, i, j = [], [], []
    for idx in range(umap_res.shape[0]):
        dist, ind = tree.query([umap_res[idx]], k=25)
        vals.extend(list(dist.squeeze()))
        j.extend(list(ind.squeeze()))
        i.extend([idx] * len(ind.squeeze()))
    print(len(vals))
    ginput = sps.csc_matrix((numpy.array(vals), (numpy.array(i),numpy.array(j))), 
                            shape=(umap_res.shape[0], umap_res.shape[0]))
    sources, targets = ginput.nonzero()
    edgelist = zip(sources.tolist(), targets.tolist())
    G = ig.Graph(edges=list(edgelist))
    optimiser = leidenalg.Optimiser()
    optimiser.set_rng_seed(random_state)
    profile = optimiser.resolution_profile(G, leidenalg.CPMVertexPartition, resolution_range=resolution_range, number_iterations=0)
    print([len(elt) for elt in profile])
    return profile

def write_peaks_and_map_to_genes(data_array, row_headers, c_labels, out_dir, refseq_exon_bed, 
                                 uniqueness_threshold=3, num_peaks=1000):
    #write the peaks present in each cluster to bed files
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    else:
        local['rm']('-r', out_dir)
        os.makedirs(out_dir)
    
    #write a file of peaks per cluster in bed format
    peak_files = []
    for idx, cluster_name in enumerate(sorted(set(c_labels))):
        cell_coords = numpy.where(c_labels == cluster_name)
        peak_sums = numpy.mean(data_array[:,cell_coords[0]], axis=1)
        peak_sort = numpy.argsort(peak_sums)
#        sorted_peaks = peak_sums[peak_sort]
#        print('Cluster {!s} -- Present Peaks: {!s}, '
#              'Min Peaks/Cell: {!s}, '
#              'Max Peaks/Cell: {!s}, '
#              'Peaks in {!s}th cell: {!s}'.format(cluster_name, numpy.sum(peak_sums > 0), 
#                                                  sorted_peaks[0], sorted_peaks[-1], 
#                                                  num_peaks, sorted_peaks[-num_peaks]))
        out_tmp = os.path.join(out_dir, 'peaks{!s}.tmp.bed'.format(cluster_name))
        out_path = out_tmp.replace('.tmp', '')
        peak_indices = peak_sort[-num_peaks:]
        with open(out_tmp, 'w') as out:
            out.write('\n'.join('chr'+'\t'.join(elt) if not elt[0].startswith('chr') else '\t'.join(elt) 
                                for elt in numpy.hstack([row_headers[peak_indices],
                                                         peak_sums[peak_indices,None].astype(str)])) + '\n')
        (local['sort']['-k1,1', '-k2,2n', out_tmp] > out_path)()
        os.remove(out_tmp)
        peak_files.append(out_path)

    bedtools, sort, cut, uniq, awk = local['bedtools'], local['sort'], local['cut'], local['uniq'], local['awk']
    out_subdir = os.path.join(out_dir, 'nearest_genes')
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    nearest_genes = []
    for path in sorted(peak_files):
        out_path = os.path.join(out_subdir, os.path.basename(path).replace('.bed', '.nearest_genes.txt'))
        cmd = (bedtools['closest', '-D', 'b', '-io', '-id', '-a', path, '-b', refseq_exon_bed] |
         cut['-f1,2,3,5,9,12'] | #fields are chrom, start, stop, peak sum, gene name, distance
         awk['BEGIN{OFS="\t"}{if($6 > -1200){print($1, $2, $3, $6, $5, $4);}}'] |
         sort['-k5,5', '-k6,6nr'] |
         cut['-f5,6'])()
        with open(out_path, 'w') as out:
            prev_gene = None
            for idx, line in enumerate(str(cmd).strip().split('\n')):
                if prev_gene is None or not line.startswith(prev_gene):
#                    print(line)
                    line_split = line.strip().split()
                    prev_gene = line_split[0]
                    out.write(line + '\n')
        nearest_genes.append(out_path)

    all_genes = []
#    for idx in range(len(nearest_genes)):
#        nearest_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes.txt'.format(idx))
    for nearest_genes_path in nearest_genes:
        with open(nearest_genes_path) as lines_in:
            all_genes.append([elt.strip().split() for elt in lines_in.readlines()])

#    count_dict = Counter([i[0] for i in itertools.chain(*[all_genes[elt] for elt in range(len(nearest_genes))])])
    count_dict = Counter([i[0] for i in itertools.chain(*all_genes)])
    #print unique genes
    for idx, nearest_genes_path in enumerate(nearest_genes):
        unique_genes = [elt for elt in all_genes[idx] if count_dict[elt[0]] < uniqueness_threshold]
        print(idx, len(unique_genes))
#        unique_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes_lt_{!s}.txt'.
#                                         format(idx, uniqueness_threshold))
        unique_genes_path = os.path.splitext(nearest_genes_path)[0] + '_lt_{!s}.txt'.format(uniqueness_threshold)
        with open(unique_genes_path, 'w') as out:
            out.write('\n'.join(['\t'.join(elt) for elt in unique_genes]) + '\n')
    #print shared genes
    shared_genes_by_cluster = []
    all_genes = [dict([(k,float(v)) for k,v in elt]) for elt in all_genes]
    for gene_name in sorted(count_dict.keys()):
        if count_dict[gene_name] < uniqueness_threshold:
            continue
        shared_genes_by_cluster.append([gene_name])
        for cluster_dict in all_genes:
            shared_genes_by_cluster[-1].append(cluster_dict.get(gene_name, 0.0))
    shared_out = os.path.join(out_subdir, 'non-unique_genes_lt_{!s}.txt'.
                              format(uniqueness_threshold))
    numpy.savetxt(shared_out, shared_genes_by_cluster, fmt='%s')
#                  fmt=('%s',)+tuple('%18f' for _ in range(len(all_genes))))

    return

def write_peaks_and_map_to_genes2(data_array, peak_topic_specificity, row_headers, c_labels, out_dir, 
                                  refseq_exon_bed, uniqueness_threshold=3, num_peaks=1000):
#    import pdb; pdb.set_trace()
    #write the peaks present in each cluster to bed files
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    else:
        local['rm']('-r', out_dir)
        os.makedirs(out_dir)
    
    #write a file of peaks per cluster in bed format
    peak_files = []
    for idx, cluster_name in enumerate(sorted(set(c_labels))):
        cell_coords = numpy.where(c_labels == cluster_name)
        peaks_present = numpy.sum(data_array[cell_coords[0],:], axis=0)
        out_tmp = os.path.join(out_dir, 'peaks{!s}.tmp.bed'.format(cluster_name))
        out_path = out_tmp.replace('.tmp', '')
#        peak_indices = peak_sort[-num_peaks:]
        peak_scores = (peak_topic_specificity ** 2) * peaks_present
        sort_idx = numpy.argsort(peak_scores[peaks_present.astype(bool)])
        peak_indices = sort_idx[-num_peaks:]
        with open(out_tmp, 'w') as out:
#            out.write('\n'.join('chr'+'\t'.join(elt) if not elt[0].startswith('chr') else '\t'.join(elt) 
#                                for elt in numpy.hstack([row_headers[peaks_present.astype(bool)][peak_indices],
#                                                         peak_scores[peaks_present.astype(bool)][peak_indices,None].astype(str)])) + '\n')
            out.write('\n'.join('\t'.join(elt) for elt in 
                                numpy.hstack([row_headers[peaks_present.astype(bool)][peak_indices],
                                              peak_scores[peaks_present.astype(bool)][peak_indices,None].astype(str)])) + '\n')
        (local['sort']['-k1,1', '-k2,2n', out_tmp] > out_path)()
        os.remove(out_tmp)
        peak_files.append(out_path)

    bedtools, sort, cut, uniq, awk = local['bedtools'], local['sort'], local['cut'], local['uniq'], local['awk']
    out_subdir = os.path.join(out_dir, 'nearest_genes')
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    nearest_genes = []
    for path in sorted(peak_files):
        out_path = os.path.join(out_subdir, os.path.basename(path).replace('.bed', '.nearest_genes.txt'))
        cmd = (bedtools['closest', '-D', 'b', '-io', '-id', '-a', path, '-b', refseq_exon_bed] |
         cut['-f1,2,3,5,9,12'] | #fields are chrom, start, stop, peak sum, gene name, distance
         awk['BEGIN{OFS="\t"}{if($6 > -1200){print($1, $2, $3, $6, $5, $4);}}'] |
         sort['-k5,5', '-k6,6nr'] |
         cut['-f5,6'])()
        with open(out_path, 'w') as out:
            prev_gene = None
            for idx, line in enumerate(str(cmd).strip().split('\n')):
                if prev_gene is None or not line.startswith(prev_gene):
#                    print(line)
                    line_split = line.strip().split()
                    prev_gene = line_split[0]
                    out.write(line + '\n')
        nearest_genes.append(out_path)

    all_genes = []
#    for idx in range(len(nearest_genes)):
#        nearest_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes.txt'.format(idx))
    for nearest_genes_path in nearest_genes:
        with open(nearest_genes_path) as lines_in:
            all_genes.append([elt.strip().split() for elt in lines_in.readlines()])

#    count_dict = Counter([i[0] for i in itertools.chain(*[all_genes[elt] for elt in range(len(nearest_genes))])])
    count_dict = Counter([i[0] for i in itertools.chain(*all_genes)])
    #print unique genes
    for idx, nearest_genes_path in enumerate(nearest_genes):
        unique_genes = [elt for elt in all_genes[idx] if count_dict[elt[0]] < uniqueness_threshold]
        print(idx, len(unique_genes))
#        unique_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes_lt_{!s}.txt'.
#                                         format(idx, uniqueness_threshold))
        unique_genes_path = os.path.splitext(nearest_genes_path)[0] + '_lt_{!s}.txt'.format(uniqueness_threshold)
        with open(unique_genes_path, 'w') as out:
            out.write('\n'.join(['\t'.join(elt) for elt in unique_genes]) + '\n')
    #print shared genes
    shared_genes_by_cluster = []
    all_genes = [dict([(k,float(v)) for k,v in elt]) for elt in all_genes]
    for gene_name in sorted(count_dict.keys()):
        if count_dict[gene_name] < uniqueness_threshold:
            continue
        shared_genes_by_cluster.append([gene_name])
        for cluster_dict in all_genes:
            shared_genes_by_cluster[-1].append(cluster_dict.get(gene_name, 0.0))
    shared_out = os.path.join(out_subdir, 'non-unique_genes_lt_{!s}.txt'.
                              format(uniqueness_threshold))
    numpy.savetxt(shared_out, shared_genes_by_cluster, fmt='%s')
#                  fmt=('%s',)+tuple('%18f' for _ in range(len(all_genes))))

    return

def write_peaks_and_map_to_genes3(data_array, row_headers, c_labels, out_dir, 
                                  refseq_exon_bed, uniqueness_threshold=3, num_peaks=1000):
#    import pdb; pdb.set_trace()
    #write the peaks present in each cluster to bed files
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    else:
        local['rm']('-r', out_dir)
        os.makedirs(out_dir)

    agg_clusters = numpy.vstack([numpy.sum(data_array[numpy.where(c_labels == cluster_idx)[0]], axis=0)
                                 for cluster_idx in sorted(set(c_labels))])
    tfidf = TfidfTransformer(norm='l2', use_idf=True, smooth_idf=True, sublinear_tf=False)
    agg_clusters_tfidf = tfidf.fit_transform(agg_clusters).toarray()

    #write a file of peaks per cluster in bed format
    peak_files = []
    for idx, cluster_name in enumerate(sorted(set(c_labels))):
        out_tmp = os.path.join(out_dir, 'peaks{!s}.tmp.bed'.format(cluster_name))
        out_path = out_tmp.replace('.tmp', '')
        sort_idx = numpy.argsort(agg_clusters_tfidf[idx])
        peak_indices = sort_idx[-num_peaks:]
        with open(out_tmp, 'w') as out:
#            out.write('\n'.join('chr'+'\t'.join(elt) if not elt[0].startswith('chr') else '\t'.join(elt) 
#                                for elt in numpy.hstack([row_headers[peaks_present.astype(bool)][peak_indices],
#                                                         peak_scores[peaks_present.astype(bool)][peak_indices,None].astype(str)])) + '\n')
            out.write('\n'.join('\t'.join(elt) for elt in 
                                numpy.hstack([row_headers[peak_indices],
                                              agg_clusters_tfidf[idx][peak_indices,None].astype(str)])) + '\n')
        (local['sort']['-k1,1', '-k2,2n', out_tmp] > out_path)()
        os.remove(out_tmp)
        peak_files.append(out_path)

    bedtools, sort, cut, uniq, awk = local['bedtools'], local['sort'], local['cut'], local['uniq'], local['awk']
    out_subdir = os.path.join(out_dir, 'nearest_genes')
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    nearest_genes = []
    for path in sorted(peak_files):
        out_path = os.path.join(out_subdir, os.path.basename(path).replace('.bed', '.nearest_genes.txt'))
        cmd = (bedtools['closest', '-D', 'b', '-io', '-id', '-a', path, '-b', refseq_exon_bed] |
         cut['-f1,2,3,5,9,12'] | #fields are chrom, start, stop, peak sum, gene name, distance
         awk['BEGIN{OFS="\t"}{if($6 > -1200){print($1, $2, $3, $6, $5, $4);}}'] |
         sort['-k5,5', '-k6,6nr'] |
         cut['-f5,6'])()
        with open(out_path, 'w') as out:
            prev_gene = None
            for idx, line in enumerate(str(cmd).strip().split('\n')):
                if prev_gene is None or not line.startswith(prev_gene):
#                    print(line)
                    line_split = line.strip().split()
                    prev_gene = line_split[0]
                    out.write(line + '\n')
        nearest_genes.append(out_path)

    all_genes = []
#    for idx in range(len(nearest_genes)):
#        nearest_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes.txt'.format(idx))
    for nearest_genes_path in nearest_genes:
        with open(nearest_genes_path) as lines_in:
            all_genes.append([elt.strip().split() for elt in lines_in.readlines()])

#    count_dict = Counter([i[0] for i in itertools.chain(*[all_genes[elt] for elt in range(len(nearest_genes))])])
    count_dict = Counter([i[0] for i in itertools.chain(*all_genes)])
    #print unique genes
    for idx, nearest_genes_path in enumerate(nearest_genes):
        unique_genes = [elt for elt in all_genes[idx] if count_dict[elt[0]] < uniqueness_threshold]
        print(idx, len(unique_genes))
#        unique_genes_path = os.path.join(out_subdir, 'peaks{!s}.nearest_genes_lt_{!s}.txt'.
#                                         format(idx, uniqueness_threshold))
        unique_genes_path = os.path.splitext(nearest_genes_path)[0] + '_lt_{!s}.txt'.format(uniqueness_threshold)
        with open(unique_genes_path, 'w') as out:
            out.write('\n'.join(['\t'.join(elt) for elt in unique_genes]) + '\n')
    #print shared genes
    shared_genes_by_cluster = []
    all_genes = [dict([(k,float(v)) for k,v in elt]) for elt in all_genes]
    for gene_name in sorted(count_dict.keys()):
        if count_dict[gene_name] < uniqueness_threshold:
            continue
        shared_genes_by_cluster.append([gene_name])
        for cluster_dict in all_genes:
            shared_genes_by_cluster[-1].append(cluster_dict.get(gene_name, 0.0))
    shared_out = os.path.join(out_subdir, 'non-unique_genes_lt_{!s}.txt'.
                              format(uniqueness_threshold))
    numpy.savetxt(shared_out, shared_genes_by_cluster, fmt='%s')
#                  fmt=('%s',)+tuple('%18f' for _ in range(len(all_genes))))

    return

class DistanceException(Exception):
    pass
class NoPeakException(Exception):
    pass

def get_closest_peaks(gene_name, row_headers, verbose=False, dist_threshold=1200, dist_excpt=False):
    gene_coord = gene_locs[gene_name][0] if gene_locs[gene_name][0][-1] == '+' else gene_locs[gene_name][-1]
    if verbose:
        print(gene_coord)
    if gene_coord[-1] == '+':
        try:
            nearest_peak = numpy.where(numpy.logical_and(row_headers[:,0] == gene_coord[0], 
                                                         row_headers[:,1].astype(int) <= gene_coord[1]))[0][-1]
        except IndexError:
            raise NoPeakException()
        alt_peak = nearest_peak - 1
#        peak_dist = numpy.absolute(gene_coord[1] - row_headers[[nearest_peak, alt_peak],1].astype(int))
        peak_dist = gene_coord[1] - row_headers[[nearest_peak, alt_peak],2].astype(int)
        if verbose:
            print(row_headers[[nearest_peak, alt_peak]])
            print(peak_dist)
    else:
        try:
            nearest_peak = numpy.where(numpy.logical_and(row_headers[:,0] == gene_coord[0], 
                                                         row_headers[:,2].astype(int) >= gene_coord[2]))[0][0]
        except IndexError:
            raise NoPeakException()
        alt_peak = nearest_peak + 1
#        peak_dist = numpy.absolute(gene_coord[2] - row_headers[[nearest_peak, alt_peak],2].astype(int))
        peak_dist = row_headers[[nearest_peak, alt_peak],1].astype(int) - gene_coord[2]
        if verbose:
            print(row_headers[[nearest_peak, alt_peak]])
            print(peak_dist)
    if peak_dist[0] > dist_threshold:
        msg = 'Warning: nearest peak to {!s} is far away! ({!s} bp)'.format(gene_name, peak_dist[0])
        if dist_excpt:
            raise DistanceException(msg)
        else:
            print(msg)
    return nearest_peak, alt_peak

def get_gene_cells(gene_name, row_headers, peak_data_array, **kwargs):
    nearest_peak, alt_peak = get_closest_peaks(gene_name, row_headers, **kwargs)
    cells_idx = peak_data_array[:,nearest_peak].astype(bool)
    return cells_idx

def get_gene_idx(gene_name, row_headers, peaktopic_frac, topic_prob_threshold=0.5, **kwargs):
    nearest_peak, alt_peak = get_closest_peaks(gene_name, row_headers, **kwargs)
    topic_idx = numpy.argsort(peaktopic_frac[nearest_peak])[::-1]
    num_to_get = numpy.where(numpy.cumsum(peaktopic_frac[nearest_peak][topic_idx]) > topic_prob_threshold)[0][0] + 1
    return nearest_peak, topic_idx[:num_to_get]

def get_gene_topn_topics(gene_name, row_headers, peaktopic_frac, ntopics=1, **kwargs):
    nearest_peak, alt_peak = get_closest_peaks(gene_name, row_headers, **kwargs)
    topic_idx = numpy.argsort(peaktopic_frac[nearest_peak])[::-1]
    return nearest_peak, topic_idx[:ntopics]

def find_nearest_genes_per_peak(peak_row_headers_subset, refseq_exon_bed, dist_threshold=1200):
    #get unix utilities
    bedtools, uniq, awk = local['bedtools'], local['uniq'], local['awk']

    peak_str = '\n'.join(['\t'.join(elt) for elt in peak_row_headers_subset]) + '\n'
    cmd = ((bedtools['closest', '-D', 'b', '-id', '-k', '10', '-a', 'stdin', '-b', refseq_exon_bed] << peak_str) 
           | uniq
           | awk['BEGIN{{OFS="\t"}}{{if($11 > -{!s}){{print($0);}}}}'.format(dist_threshold)]
           | awk['!seen[$1,$2,$3,$8]++'])
    res = cmd()
    res = numpy.array([elt.split() for elt in res.strip().split('\n')], dtype=object)
    return res

def get_peak_idx(coord, peak_row_headers):
    idx1 = numpy.where(peak_row_headers[:,0] == coord[0])[0]
    idx2 = numpy.where(peak_row_headers[idx1,1] == coord[1])[0]
    idx3 = numpy.where(peak_row_headers[idx1,2][idx2] == coord[2])[0]
    if len(idx3) == 0:
        raise Exception('Could not find the peak {!s}:{!s}-{!s}.'.format(*coord[:3]))
    return idx1[idx2[idx3[0]]]

def get_gene_cells2(gname, peak_row_headers, peak_data_array, verbose=True, peak_to_gene_mapping=None,
                    dist_threshold=1200,
                    refseq_exon_bed='ATAC_sequencing/2018_worm_atac/ref_data/WS235/c_elegans.WS235.exon_info.bed.gz'):
    if peak_to_gene_mapping is None:
        peak_to_gene_mapping = find_nearest_genes_per_peak(peak_row_headers, refseq_exon_bed,
                                                           dist_threshold=dist_threshold)
    gene_idx = numpy.where(peak_to_gene_mapping[:,7] == gname)[0]
    gene_peaks = peak_to_gene_mapping[gene_idx]
    if verbose is True:
        print(gene_peaks)
    peak_idx = [get_peak_idx(elt, peak_row_headers) for elt in gene_peaks]
    cell_idx = numpy.where(numpy.any(peak_data_array[:,peak_idx], axis=1))[0]
    return cell_idx

REFDATA = 'ATAC_sequencing/2018_worm_atac/ref_data/WS235'
#refseq_exon_bed = os.path.join(REFDATA, 'c_elegans.WS272.canonical_geneset.genes.common_names.sorted.bed.gz')
refseq_exon_bed = os.path.join(REFDATA, 'c_elegans.WS235.exon_info.bed.gz')
#ucsc = True if peak_row_headers[0][0].startswith('chr') else False
ucsc = True
with gzip.open(refseq_exon_bed, 'rb') as lines_in:
    exon_locs = []
    for line in lines_in:
        line = line.decode()[3:].strip().split()
        if ucsc is True:
            line[0] = 'chr{!s}'.format(line[0])
        line[1] = int(line[1])
        line[2] = int(line[2])
        exon_locs.append(line)

gene_locs = {}
for exon in exon_locs:
    gene_locs.setdefault(exon[3], []).append(exon)
for gene, locs in gene_locs.items():
    gene_locs[gene] = sorted(locs, key=lambda x:(x[1],x[2]))