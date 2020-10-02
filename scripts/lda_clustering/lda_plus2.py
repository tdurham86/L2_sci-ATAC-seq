#!/usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import glob
import matplotlib
from matplotlib import pyplot
import numpy
import os
from plumbum import local
import pwd
import scipy.sparse as sps
from scipy import stats
from sklearn import neighbors
import sys
import time
import umap

sys.path.append(os.path.dirname(__file__))
import lda_analysis_lib as lal

def wait_for_jobs(jobname, uge_proj='sage', num_jobs=1):
    qstat_cmd = local['qstat']['-u', pwd.getpwuid(os.getuid())[0], '-j', jobname] | local['grep']['-P', 'project:\s*{!s}'.format(uge_proj)]
    while len([elt for elt in qstat_cmd(retcode=None).strip().split('\n') if elt]) >= num_jobs:
        time.sleep(10)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cuts')
    parser.add_argument('--bow')
    parser.add_argument('--indextable')
    parser.add_argument('--peaks_bed')
    parser.add_argument('--outdir')
    parser.add_argument('--topics', default='2,5,10,15,20,30,40,50,60,70')
    parser.add_argument('--seed', type=int)
    parser.add_argument('--job_threads', type=int)
    parser.add_argument('--job_coremem', type=int)
    parser.add_argument('--job_proj')
    parser.add_argument('--job_name')
    parser.add_argument('--no_peak_prevalence_filter', action='store_true', default=False)
    parser.add_argument('--alpha', type=float, default=3.0)
    parser.add_argument('--beta', type=float, default=2000.0)
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    #remove high- and low-prevalence peaks and low-coverage cells
    try:
        peak_data_sparse = numpy.loadtxt(args.bow, dtype=int, skiprows=3)
    except StopIteration:
        #probably NFS lag; just wait a few seconds and try again
        time.sleep(10)
        peak_data_sparse = numpy.loadtxt(args.bow, dtype=int, skiprows=3)
    peak_data = sps.csr_matrix((peak_data_sparse[:,2], 
                                (peak_data_sparse[:,0] - 1, 
                                 peak_data_sparse[:,1] - 1)))

    try:
        cell_names = numpy.loadtxt(args.indextable, dtype=object)[:,0]
    except StopIteration:
        time.sleep(10)
        cell_names = numpy.loadtxt(args.indextable, dtype=object)[:,0]

    try:
        peak_row_headers = numpy.loadtxt(args.peaks_bed, dtype=object)
    except StopIteration:
        time.sleep(10)
        peak_row_headers = numpy.loadtxt(args.peaks_bed, dtype=object)

    peak_row_headers = numpy.hstack([peak_row_headers, numpy.array(['name'] * peak_row_headers.shape[0])[:,None]])

    peak_data_array = peak_data.toarray().astype(numpy.int8)
    del(peak_data)

    peak_prevalence = numpy.sum(peak_data_array, axis=0)
    if args.no_peak_prevalence_filter is False:
        peak_prev_sort_idx = numpy.argsort(peak_prevalence)
        peak_prev_sort_data = numpy.log10(peak_prevalence[peak_prev_sort_idx]/peak_data_array.shape[0])

        test_step = numpy.convolve(numpy.hstack([numpy.ones(500), 
                                                 0 - numpy.ones(500)]), 
                                   numpy.ones(200)/200, mode='valid')
        conv_dist = numpy.convolve(peak_prev_sort_data - numpy.mean(peak_prev_sort_data), test_step, mode='valid')
        iqr_threshold = max((stats.iqr(conv_dist) * 4) + numpy.mean(conv_dist), 1.0)
        outlier_idx = numpy.where(conv_dist > iqr_threshold)[0] + 400

        try:
            thresh1 = outlier_idx[numpy.where(outlier_idx < peak_prevalence.shape[0]/2)[0][-1]]
        except IndexError:
            thresh1 = 0
        prev_thresh1 = 10**peak_prev_sort_data[thresh1]

        try:
            thresh2 = outlier_idx[numpy.where(outlier_idx > peak_prevalence.shape[0]/2)[0][0]]
        except IndexError:
            thresh2 = -1
        prev_thresh2 = 10**peak_prev_sort_data[thresh2]

        fig, axes = pyplot.subplots(nrows=1, ncols=1, figsize=(6,4))
        axes.plot(numpy.arange(peak_prev_sort_data.shape[0]), 
                  peak_prev_sort_data)
        ax2 = axes.twinx()
        ax2.plot(numpy.arange(peak_prev_sort_data.shape[0] - len(test_step) + 1) + len(test_step)/2, conv_dist, color='orange')
        axes.axvline(thresh1, color='r')
        axes.axvline(thresh2, color='r')
        fig.savefig(os.path.join(args.outdir, 'peak_prev_filter.pdf'),
                    bbox_inches='tight')
    else:
        prev_thresh1 = numpy.min(peak_prevalence) - 1
        prev_thresh2 = numpy.max(peak_prevalence) + 1

    filtered_peaks_idx = numpy.where(numpy.logical_and(peak_prevalence/peak_data_array.shape[0] > prev_thresh1, peak_prevalence/peak_data_array.shape[0] < prev_thresh2))[0]
    filtered_peaks = peak_data_array[:,filtered_peaks_idx]

    cell_coverage = numpy.sum(filtered_peaks, axis=1)
    cell_cov_sort_idx = numpy.argsort(cell_coverage)
    cell_cov_sort_data = numpy.log10(cell_coverage[cell_cov_sort_idx])

    test_step = numpy.convolve(numpy.hstack([numpy.ones(50), 0 - numpy.ones(50)]), numpy.ones(20)/20, mode='valid')
    conv_dist = numpy.convolve(cell_cov_sort_data - numpy.mean(cell_cov_sort_data), test_step, mode='valid')
    iqr_threshold = max((stats.iqr(conv_dist) * 4) + numpy.mean(conv_dist), 1.0)
    outlier_idx = numpy.where(conv_dist > iqr_threshold)[0] + 40

    try:
        thresh1 = outlier_idx[numpy.where(outlier_idx < cell_coverage.shape[0]/2)[0][-1]]
    except IndexError:
        thresh1 = 0
    cov_thresh1 = 10**cell_cov_sort_data[thresh1]

    try:
        thresh2 = outlier_idx[numpy.where(outlier_idx > cell_coverage.shape[0]/2)[0][0]]
    except IndexError:
        thresh2 = -1
    cov_thresh2 = 10**cell_cov_sort_data[thresh2]

    fig, axes = pyplot.subplots(nrows=1, ncols=1, figsize=(6,4))
    axes.plot(numpy.arange(cell_cov_sort_data.shape[0]), 
              cell_cov_sort_data)
    ax2 = axes.twinx()
    ax2.plot(numpy.arange(cell_cov_sort_data.shape[0] - len(test_step) + 1) + len(test_step)/2, conv_dist, color='orange')
    axes.axvline(thresh1, color='r')
    axes.axvline(thresh2, color='r')
    fig.savefig(os.path.join(args.outdir, 'cell_cov_filter.pdf'),
                bbox_inches='tight')

    filtered_low_cells_idx = numpy.where(cell_coverage > cov_thresh1)[0]
    filtered_peaks_low_cells = filtered_peaks[filtered_low_cells_idx]

    filtered_bow = os.path.join(args.outdir, 'filtered_peaks_iqr4.0_low_cells.bow')
    filt_peak_data_array = filtered_peaks_low_cells
    with open(filtered_bow, 'w') as out:
        out.write('{!s}\n{!s}\n1000\n'.format(*filtered_peaks_low_cells.shape))
        for cell_idx, peak_idx in zip(*numpy.where(filtered_peaks_low_cells)):
            out.write('{!s}\t{!s}\t{!s}\n'.format(cell_idx + 1, peak_idx + 1, 
                                                  filtered_peaks_low_cells[cell_idx, peak_idx]))

    filtered_peaks_bed = os.path.join(args.outdir, 'filtered_peaks_iqr4.0_low_cells.bed')
    filt_peak_row_headers = peak_row_headers[filtered_peaks_idx]
    with open(filtered_peaks_bed, 'w') as out:
        for idx in filtered_peaks_idx:
            out.write('\t'.join(peak_row_headers[idx]) + '\n')

    filtered_indextable = os.path.join(args.outdir, 'filtered_peaks_iqr4.0_low_cells.indextable.txt')
    filt_cell_names = cell_names[filtered_low_cells_idx]
    with open(filtered_indextable, 'w') as out:
        out.write('\n'.join(['{!s}\tfilt'.format(elt) for elt in filt_cell_names]) + '\n')

    #split BOW for the hyperparameter search
    lda_jar = 'lda_commandline/dist9/LatentDirichletAllocation.jar'
    hypsearch_out = os.path.join(args.outdir, 'hypsearch_iters')
    if not os.path.isdir(hypsearch_out):
        os.makedirs(hypsearch_out)
    num_parts = 5
    cmd = local['java']['-Xms16G', '-cp', lda_jar, 'org.rhwlab.lda.cache.LDA_CommandLine', '-part', '-p', str(num_parts), '-ib', filtered_bow, '-o', hypsearch_out, '-s', str(args.seed)] > os.path.join(hypsearch_out, 'bow_split.out')
#    cmd()

    #hyperparameter search over topics
    script_path = 'ATAC_sequencing/all_sciATAC/src_me/lda_hyperparameter_search.py'
    cmd = local['python'][script_path, 
                          '--topics={!s}'.format(args.topics), 
                          '--iterations=4000', 
                          '--alphas={!s}'.format(args.alpha), 
                          '--betas={!s}'.format(args.beta), 
                          '--parts={!s}'.format(os.path.join(hypsearch_out, '*part*.bow')), 
                          '--seed={!s}'.format(args.seed), 
                          '--out_root={!s}'.format(hypsearch_out), 
                          '--lda_burnin=1000', 
                          '--threads={!s}'.format(args.job_threads), 
                          '--coremem={!s}'.format(args.job_coremem), 
                          '--lda_thinning=5', 
                          '--lda_cache=0', 
                          '--pe_dist=topic', 
                          '--pe_statistic=mode', 
                          '--chib_iters=1000', 
                          '--chib_burn=200', 
                          '--cleanup', 
                          '--job_name={!s}'.format(args.job_name), 
                          '--uge_proj={!s}'.format(args.job_proj)]
#    cmd()
    wait_for_jobs(args.job_name + '*', args.job_proj, num_jobs=1)

    #calculate a reasonable number of topics to use
    results_paths = glob.glob(os.path.join(hypsearch_out, '*/*.chib'))
    sortkey = lambda x: (int(x.split('_')[-4].replace('topics','')), 
                         int(os.path.basename(os.path.dirname(x)).split('_')[0]))
    all_results = sorted([elt for elt in results_paths], key=sortkey)
    topic_num = len(all_results)//num_parts
    all_results = numpy.array(all_results).reshape(topic_num,num_parts).T.flatten()
    shortnames = ['part{!s}_{!s}'.format(idx//topic_num, elt.split('_')[-4]) 
                  for idx, elt in enumerate(all_results)]

    all_res_list = [[numpy.loadtxt(elt, delimiter=',')[:,1] 
                     for elt in all_results[idx:idx+topic_num]]
                    for idx in range(0, len(shortnames), topic_num)]
    all_res_mean = [numpy.vstack(elt).mean(axis=1) for elt in all_res_list]
    topic_xvals = sorted(set([int(elt.split('_')[1].replace('topics','')) 
                              for elt in shortnames]))
    mean_min_val = numpy.mean([topic_xvals[numpy.argmin(elt)] for elt in all_res_mean])
    topic_num_to_use = max(int(1.5 * mean_min_val), 5)
#    alpha_to_use = 3.0/topic_num_to_use

    fig, axes = pyplot.subplots(nrows=1, ncols=1, figsize=(7,5))
    for idx in range(num_parts):
#        part_data = numpy.vstack(all_res_list[idx])
#        axes.plot(topic_xvals, numpy.mean(part_data, axis=1), 
#                  label='Part {!s}'.format(idx))
        part_data = all_res_mean[idx]
        axes.plot(topic_xvals, part_data, 
                  label='Part {!s}'.format(idx))
    axes.axvline(mean_min_val, color='k')
    axes.axvline(topic_num_to_use, color='k', linestyle='--')
    axes.set_xticks(topic_xvals)
    axes.set_xticklabels(topic_xvals, rotation=90, fontsize=10)
    axes.set_xlabel('Number of topics')
    axes.set_ylabel('Mean test set perplexity')
    axes.legend(loc='best')
    fig.savefig(hypsearch_out + '.pdf', bbox_inches='tight')

    #train detailed LDA model
    cmd = local['python'][script_path, 
                          '--topics={!s}'.format(topic_num_to_use), 
                          '--iterations=4000', 
                          '--alphas={!s}'.format(args.alpha), 
                          '--betas={!s}'.format(args.beta), 
#                          '--alphas={!s}'.format(alpha_to_use), 
#                          '--betas=0.23',
                          '--parts={!s}'.format(filtered_bow), 
                          '--seed={!s}'.format(args.seed), 
                          '--out_root={!s}'.format(args.outdir), 
                          '--lda_burnin=1000', 
                          '--threads={!s}'.format(args.job_threads), 
                          '--coremem={!s}'.format(args.job_coremem), 
                          '--lda_thinning=5', 
                          '--lda_cache=100', 
                          '--pe_dist=topic', 
                          '--pe_statistic=mode', 
                          '--job_name={!s}'.format(args.job_name), 
                          '--uge_proj={!s}'.format(args.job_proj)]
    cmd()
    wait_for_jobs(args.job_name + '*', args.job_proj, num_jobs=1)

    #generate umap plot
    doctopic_path = os.path.join(args.outdir, '0000_topics{!s}_alpha{:1.3f}_beta{:1.3f}/topic_mode.theta'.format(topic_num_to_use, args.alpha, args.beta))
    doctopic_peaks = numpy.loadtxt(doctopic_path, delimiter=',', dtype=float)

    col_means = numpy.mean(doctopic_peaks.T, axis=0)
    doctopic_peaks_norm = doctopic_peaks.T - col_means
    l2_for_norm = (doctopic_peaks_norm ** 2).sum(axis=0).flatten() ** 0.5
    doctopic_peaks_norm /= l2_for_norm
    doctopic_peaks_norm = doctopic_peaks_norm.T
    doctopic_peaks_frac = (doctopic_peaks.T/doctopic_peaks.sum(axis=1).astype(float)).T

    peaktopic_path = os.path.join(args.outdir, '0000_topics{!s}_alpha{:1.3f}_beta{:1.3f}/topic_mode.phi'.format(topic_num_to_use, args.alpha, args.beta))
    peaktopic = numpy.loadtxt(peaktopic_path, delimiter=',', dtype=float).T

    col_means = numpy.mean(peaktopic.T, axis=0)
    peaktopic_norm = peaktopic.T - col_means
    l2_for_norm = (peaktopic_norm ** 2).sum(axis=0).flatten() ** 0.5
    peaktopic_norm /= l2_for_norm
    peaktopic_norm = peaktopic_norm.T
    peaktopic_frac = (peaktopic.T/peaktopic.sum(axis=1).astype(float)).T

    doctopic_peaks_umap2_obj = umap.UMAP(n_components=2, random_state=165)
    doctopic_peaks_umap2_res = doctopic_peaks_umap2_obj.fit_transform(doctopic_peaks_norm)

    fig, axes = pyplot.subplots(nrows=1, ncols=1, figsize=(8,7))
    im = axes.scatter(doctopic_peaks_umap2_res[:,0], 
                      doctopic_peaks_umap2_res[:,1], 
                      s=1)
    axes.set_xlabel('UMAP 1')
    axes.set_ylabel('UMAP 2')
    fig.savefig(os.path.join(args.outdir, 'umap.pdf'), bbox_inches='tight')

    #assign cells to topic clusters
    topic_max_fracs = numpy.max(doctopic_peaks_frac, axis=0)
    topic_mean_sims = []
    for topic in range(doctopic_peaks_frac.shape[1]):
        cellsort = numpy.argsort(doctopic_peaks_frac[:,topic])
        cells_to_use = 50
        celldata = doctopic_peaks_frac[cellsort[-cells_to_use:]]
        meancell = numpy.mean(celldata, axis=0)
        topic_mean_sims.append(numpy.mean([meancell.dot(celldata[idx]) for idx in range(celldata.shape[0])]))
    topic_mean_sims = numpy.array(topic_mean_sims)
    sorted_topics = numpy.argsort(topic_mean_sims)[::-1]

    fig, axes = pyplot.subplots(nrows=1, ncols=1, figsize=(8,4))
    axes.bar(numpy.arange(topic_num_to_use), topic_mean_sims[sorted_topics])
    axes.axhline(0.2, color='r')
    axes.set_xticks(numpy.arange(topic_num_to_use))
    axes.set_xticklabels(numpy.arange(topic_num_to_use)[sorted_topics], 
                         rotation=90, fontsize=6)
    axes.set_xlabel('Topics')
    axes.set_ylabel('Mean sims value')
    fig.savefig(os.path.join(args.outdir, 'topic_cells_mean_similarities.pdf'),
                bbox_inches='tight')

    topics_to_use = numpy.where(topic_mean_sims > 0.2)[0]
    doctopic_frac_subset = doctopic_peaks_frac[:,topics_to_use]
    max_to_use = numpy.max(doctopic_frac_subset, axis=1)
    cell_topic_clusts = numpy.argmax(doctopic_peaks_frac, axis=1)
    cell_topic_clusts[numpy.where(max_to_use <= 0.5)[0]] = -1

    doctopic_peaks_umap10_res = umap.UMAP(n_components=10, random_state=args.seed + 5).fit_transform(doctopic_peaks_norm)
    cell_topic_clusts_expanded = cell_topic_clusts.copy()

    tree = neighbors.KDTree(doctopic_peaks_umap10_res, metric='euclidean')
    min_cluster_size = 150

    test_step = numpy.hstack([numpy.ones(5), 0 - numpy.ones(5)])
    for idx in sorted(set(cell_topic_clusts)):
        topic_cell_idx = list(numpy.where(cell_topic_clusts == idx)[0])
        if idx > -1 and len(topic_cell_idx) < min_cluster_size:
            cells_added = 0
            deficit = min_cluster_size - len(topic_cell_idx)
            while True:
                kde = neighbors.KernelDensity(bandwidth=1.0,
                                              kernel='gaussian',
                                              metric='euclidean').fit(doctopic_peaks_umap10_res[tuple(topic_cell_idx),:])
                meancell = numpy.mean(kde.sample(200, random_state=args.seed + idx), axis=0)
                nn_num = max(len(topic_cell_idx) + int(0.25 * len(topic_cell_idx)), 20)
                orig_dist, orig_ind = tree.query([meancell], k=nn_num)
                dist = orig_dist.squeeze()
                ind = orig_ind.squeeze()
            
                #check dist to avoid jumping to another cluster
                if len(dist) >= 1.5 * len(test_step):
                    conv_dist = numpy.convolve(dist - numpy.mean(dist), test_step, mode='valid')
                    iqr_threshold = max((stats.iqr(conv_dist) * 1.5) + numpy.mean(conv_dist), 1.0)
                    outlier_idx = numpy.where(conv_dist > iqr_threshold)[0]
                    if len(outlier_idx) > 0:
                        outlier_idx = int(len(test_step)/2) + outlier_idx[0]
                        ind = ind[:outlier_idx]
                if -1 in cell_topic_clusts_expanded[ind]:
                    to_add_idx = ind[numpy.where(cell_topic_clusts_expanded[ind] == -1)[0][0]]
                else:
                    break
                topic_cell_idx.append(to_add_idx)
                cell_topic_clusts_expanded[to_add_idx] = idx
                cells_added += 1

            unique_vals, unique_counts = numpy.unique(cell_topic_clusts_expanded[orig_ind], return_counts=True)
            if len(topic_cell_idx) < 50:
                cell_topic_clusts_expanded[topic_cell_idx] = -1

    final_topic_clusters = sorted(set(cell_topic_clusts_expanded))
    print('{!s} topic clusters passed filtering.'.format(len(final_topic_clusters) - 1))
    final_total_cell_count = len(numpy.where(cell_topic_clusts_expanded > -1)[0])
    print('{!s} cells ({!s}%) passed filtering.'.format(final_total_cell_count, 
                                                     round((final_total_cell_count/doctopic_peaks_norm.shape[0]) * 100, 2)))
    for idx in final_topic_clusters:
        print('{!s}: {!s}'.format(idx, len(numpy.where(cell_topic_clusts_expanded == idx)[0])))
    cluster_map_path = os.path.join(args.outdir, 'cell_name_to_topic_cluster_map.txt')
    valid_cell_idx = numpy.where(cell_topic_clusts_expanded > -1)[0]
    numpy.savetxt(cluster_map_path, 
                  numpy.vstack([filt_cell_names.astype(object), 
                                cell_topic_clusts_expanded.astype(object)]).T[valid_cell_idx],
                  delimiter='\t', fmt='%s')

    #generate topic tissue expression distribution
    ncols=3
    nrows=int(numpy.ceil(topic_num_to_use/ncols))
    xlocs = numpy.arange(peaktopic_frac.shape[0])
    peak_num = []
    peak_sum = []
    topic_top_peaks_idx = []
    all_topic_peaks_idx = []

    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(5*ncols,3*nrows))
    for idx, topic in enumerate(numpy.arange(peaktopic_frac.shape[1])):
        row_idx, col_idx = idx//ncols, idx%ncols
        peak_sort_by_topic = numpy.argsort(peaktopic_frac[:,topic])[::-1]
        sorted_data = peaktopic_frac[peak_sort_by_topic, topic]
        axes[row_idx, col_idx].plot(xlocs, sorted_data)
    
        try:
            zero_idx = numpy.where(sorted_data == 0)[0][0]
        except IndexError:
            zero_idx = sorted_data.shape[0]
        zero_sum = numpy.sum(sorted_data[:zero_idx])
        peak_num.append(zero_idx)
        peak_sum.append(zero_sum)

        try:
            info_idx = max(250, numpy.where(sorted_data > 0.5)[0][-1])
        except IndexError:
            info_idx = 250
        info_sum = numpy.sum(sorted_data[:info_idx])

        topic_top_peaks_idx.append(peak_sort_by_topic[:info_idx])
        all_topic_peaks_idx.append(peak_sort_by_topic[:zero_idx])

        axes[row_idx, col_idx].axvline(peak_num[-1], color='r')
        axes[row_idx, col_idx].annotate('{!s} Peaks\nSum: {!s}'.format(peak_num[-1], round(peak_sum[-1])), (zero_idx+100, sorted_data[0] * 0.85), color='r')

        axes[row_idx, col_idx].axvline(info_idx, color='purple')
        axes[row_idx, col_idx].annotate('{!s} Peaks\nSum: {!s}'.format(info_idx, round(info_sum)), (info_idx+100, sorted_data[0] * 0.75), color='purple')

        axes[row_idx, col_idx].set_xlabel('Sorted peaks')
        axes[row_idx, col_idx].set_ylabel('Topic {!s} fraction'.format(topic))
    fig.savefig(os.path.join(args.outdir, 'peaktopic_specificity.pdf'),
                bbox_inches='tight')

    out_dir = os.path.join(args.outdir, 'topic_mode_top_topic_peaks_sorted_topics')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    peak_files = []
    for idx, topic in enumerate(sorted_topics):
        peak_file = os.path.join(out_dir, 'topic{!s}.rank{!s}.bed'.format(topic, idx))
        peak_files.append(peak_file)
        peak_data = numpy.hstack([filt_peak_row_headers[:,:4], peaktopic_frac[:,topic][:,None].astype(str)])
        with open(peak_file, 'w') as out:
            out.write('\n'.join(['\t'.join(list(elt)) 
                                 for elt in peak_data[sorted(topic_top_peaks_idx[topic])]]) + '\n')

    REFDATA = 'ATAC_sequencing/2018_worm_atac/ref_data/WS235'
    refseq_exon_bed = os.path.join(REFDATA, 'c_elegans.WS235.exon_info.bed.gz')
    out_subdir = os.path.join(out_dir, 'nearest_genes')
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    nearest_genes = lal.find_nearest_genes(peak_files, out_subdir, refseq_exon_bed)

    allpeaks_path = filtered_peaks_bed
    allpeaks_extracols = os.path.splitext(allpeaks_path)[0] + '.extra_cols.bed'
    cmd = local['awk']['BEGIN{OFS="\t"}{print($0,"0.0");}', allpeaks_path] > allpeaks_extracols
    cmd()

    allpeaks_nearest_genes = lal.find_nearest_genes([allpeaks_extracols], 
                                                    args.outdir, 
                                                    refseq_exon_bed)[0]

    l2_tissue_db_path = os.path.join(os.path.dirname(REFDATA),'gexplore_l2_tissue_expr.txt')
    expr_db_headers, expr_db = lal.load_expr_db(l2_tissue_db_path)
    genes = lal.get_gene_data(allpeaks_nearest_genes, expr_db, topn=100000)

    nearest_genes_glob = os.path.join(out_subdir, '*.nearest_genes.txt')
    lal.plot_l2_tissues(nearest_genes_glob, os.path.dirname(REFDATA), 
                        expr_db=None, topn=250, nsamples=100,
                        savefile=out_dir + '_peaksnorm.pdf',
                        display_in_notebook=False, 
                        allgenes_path=allpeaks_nearest_genes)

    ncols=3
    nrows = int(numpy.ceil(doctopic_peaks_frac.shape[1]/ncols))
    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols,3*nrows))
    for idx, topic in enumerate(sorted_topics):
        row_idx, col_idx = int(idx/ncols), int(idx%ncols)
        if nrows > 1 and ncols > 1:
            ax = axes[row_idx, col_idx]
        elif nrows > 1 or ncols > 1:
            ax = axes[idx]
        else:
            ax = axes
        im = ax.scatter(doctopic_peaks_umap2_res[:,0], 
                        doctopic_peaks_umap2_res[:,1],
                        cmap='viridis',
                        c=doctopic_peaks_frac[:,topic],
                        s=2)
        ax.set_ylabel('UMAP2')
        ax.set_xlabel('UMAP1')
        ax.set_title('Topic {!s}'.format(topic))
        fig.colorbar(im, ax=ax)
    fig.savefig(os.path.join(args.outdir, 'topic_dist.png'),
                bbox_inches='tight')

    out_dir = os.path.join(args.outdir, 'topic_mode_top_topic_peaks_sorted_filtered_topics')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    peak_files = []
    for idx, topic in enumerate([elt for elt in sorted_topics if elt in final_topic_clusters]):
        peak_file = os.path.join(out_dir, 'topic{!s}.rank{!s}.bed'.format(topic, idx))
        peak_files.append(peak_file)
        peak_data = numpy.hstack([filt_peak_row_headers[:,:4], peaktopic_frac[:,topic][:,None].astype(str)])
        with open(peak_file, 'w') as out:
            out.write('\n'.join(['\t'.join(list(elt)) 
                                 for elt in peak_data[sorted(topic_top_peaks_idx[topic])]]) + '\n')
        
    out_subdir = os.path.join(out_dir, 'nearest_genes')
    if not os.path.isdir(out_subdir):
        os.makedirs(out_subdir)
    nearest_genes = lal.find_nearest_genes(peak_files, out_subdir, refseq_exon_bed)

    nearest_genes_glob = os.path.join(out_subdir, '*.nearest_genes.txt')

    lal.plot_l2_tissues(nearest_genes_glob, os.path.dirname(REFDATA), 
                        expr_db=None, topn=250, nsamples=100,
                        savefile=out_dir + '_peaksnorm.pdf',
                        display_in_notebook=False, 
                        allgenes_path=allpeaks_nearest_genes)

    ncols=3
    nrows = int(numpy.ceil(len(final_topic_clusters)/ncols))
    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, figsize=(3*ncols,3*nrows))
    for idx, topic in enumerate([elt for elt in sorted_topics if elt in final_topic_clusters]):
        row_idx, col_idx = int(idx/ncols), int(idx%ncols)
        if nrows > 1 and ncols > 1:
            ax = axes[row_idx, col_idx]
        elif nrows > 1 or ncols > 1:
            ax = axes[idx]
        else:
            ax = axes
        ax.scatter(doctopic_peaks_umap2_res[:,0], 
                   doctopic_peaks_umap2_res[:,1],
                   cmap='viridis',
                   c=doctopic_peaks_frac[:,topic],
                   s=2)
        ax.set_ylabel('UMAP2')
        ax.set_xlabel('UMAP1')
        ax.set_title('Topic {!s}'.format(topic))
    fig.savefig(os.path.join(args.outdir, 'topic_dist.filtered_topics.png'),
                bbox_inches='tight')

    #do a topic cluster analysis
    topic_analysis_dir = os.path.join(args.outdir, 'topic_analysis')
    if not os.path.isdir(topic_analysis_dir):
        os.makedirs(topic_analysis_dir)
    local['python']('ATAC_sequencing/all_sciATAC/src_me/split_bed_script.py', '--cell_name_to_cluster_map={!s}'.format(cluster_map_path), '--cuts_bed_file={!s}'.format(args.cuts), '--out_dir={!s}'.format(topic_analysis_dir))

    macs_temp = 'macs2 callpeak -t {cuts_bed!s} -c ATAC_sequencing/all_sciATAC/data_reprocess_ce11/neg_ctrl/gdna_neg.blacklist_filtered.cuts.bed --format=BED -n {basename!s} --outdir={outdir!s} -g 9e7 --nomodel --qvalue=0.05 --SPMR --tsize=60 --bdg --keep-dup all --call-summits;\n\nbedtools merge -i {peaks_base!s}.narrowPeak > {peaks_base!s}.merged.bed;\n\nbedtools slop -i {summits_base!s}.bed -g ATAC_sequencing/2018_worm_atac/ref_data/WS235/c_elegans.WS235.chrom.sizes -b 50 | bedtools merge -i stdin > {summits_base!s}.slop50.merged.bed;'

    for cuts_bed in glob.glob(os.path.join(topic_analysis_dir, '*/*.cuts.bed')):
        macs_cmd = macs_temp.format(cuts_bed=cuts_bed,
                                    basename=os.path.basename(os.path.dirname(cuts_bed)),
                                    outdir=os.path.dirname(cuts_bed),
                                    peaks_base=cuts_bed.replace('.cuts.bed', '_peaks'),
                                    summits_base=cuts_bed.replace('.cuts.bed', '_summits'))
        qsub_path = os.path.join(topic_analysis_dir, 'qsub_macs2.sh')
        with open(qsub_path, 'w') as out:
            out.write('#! /usr/bin/env bash\n\n')
            out.write('source activate merged_env;\n\n')
            out.write(macs_cmd + '\n\n')
        os.chmod(qsub_path, 0o774)

        jobname = '{!s}_macs2'.format(args.job_name)
        local['qsub']('-j', 'y',
                      '-N', jobname,
                      '-o', os.path.join(topic_analysis_dir, 'job.out'),
                      '-cwd',
                      '-P', args.job_proj,
                      '-l', 'mfree=16G',
                      qsub_path)
    wait_for_jobs(jobname, args.job_proj, num_jobs=1)
    
    merged_peaks = os.path.join(topic_analysis_dir, 'all_peaks.merged.bed')
    peak_files = glob.glob(os.path.join(topic_analysis_dir, '*/*_peaks.merged.bed'))
    peak_merge_cmd = (local['cat'].__getitem__(peak_files) 
                      | local['sort']['-k1,1', '-k2,2n']
                      | local['bedtools']['merge', '-i', 'stdin'] > merged_peaks)
    peak_merge_cmd()

    to_bow_cmd = 'bedtools intersect -wo -a {cuts_bed!s} -b {peaks_bed!s} | python ATAC_sequencing/all_sciATAC/src_me/peaks_to_bow2.py {peaks_bed!s} {indextable!s} --out_file={outfile!s};'

    qsub_path = os.path.join(topic_analysis_dir, 'qsub_peakbow_allcells.sh')
    out_bow = merged_peaks.replace('.bed', '.allcells.bow')
    with open(qsub_path, 'w') as out:
        out.write('#! /usr/bin/env bash\n\n')
        out.write('source activate merged_env;\n\n')
        out.write(to_bow_cmd.format(cuts_bed=args.cuts, 
                                    peaks_bed=merged_peaks, 
                                    indextable=args.indextable,
                                    outfile=out_bow) + '\n\n')
    os.chmod(qsub_path, 0o774)

    jobname = '{!s}_peakbow_allcells'.format(args.job_name)
    local['qsub']('-j', 'y',
                  '-N', jobname,
                  '-o', os.path.join(topic_analysis_dir, '{!s}.out'.format(jobname)),
                  '-cwd',
                  '-P', args.job_proj,
                  '-l', 'mfree=16G',
                  qsub_path)

    wait_for_jobs(args.job_name + '*bow', args.job_proj, num_jobs=1)
