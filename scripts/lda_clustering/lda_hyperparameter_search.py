#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import glob
import numpy
import os
from plumbum import local
import pwd
import subprocess
import time

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--topics', help='Comma-separated list of topic nums.')
    parser.add_argument('--iterations', help='Comma-separated list of iteration numbers to test.')
    parser.add_argument('--alphas', help='Range from which to pick alphas. (e.g. 0.001-0.1)')
    parser.add_argument('--betas', help='Range from which to pick betas. (e.g. 0.001-0.1)')
    parser.add_argument('--parts', help='Glob of cross validation set bows.')
    parser.add_argument('--search_len', help='Number of hyperparameter combos to test (i.e. the number of models to train). By default, will automatically compute this as the number of topics times the number of bow partitions.', type=int, default=-1)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--out_root')
    parser.add_argument('--start_from', type=int, default=0, help='Start numbering the search iterations from this number. Useful if continuing a search started by a previous round of model training.')
    parser.add_argument('--lda_burnin', type=int, default=750)
    parser.add_argument('--logsearch', action='store_true', default=False)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--coremem', type=int, default=3)
    parser.add_argument('--lda_thinning', type=int, default=4)
    parser.add_argument('--lda_cache', type=int, default=0)
    parser.add_argument('--pe_dist', default='topic')
    parser.add_argument('--pe_statistic', default='mode')
    parser.add_argument('--kde_precision', type=float, default=1.0)
    parser.add_argument('--chib_iters', type=int, default=1000)
    parser.add_argument('--chib_burnin', type=int, default=200)
    parser.add_argument('--cleanup', action='store_true', default=False)
    parser.add_argument('--random_topics', action='store_true', default=False)
    parser.add_argument('--num_jobs', type=int, default=100)
    parser.add_argument('--uge_proj', default='sage')
    parser.add_argument('--job_name', default='ldahyp')

    args = parser.parse_args()

    out_root = args.out_root if args.out_root else os.getcwd()
    if not os.path.isdir(out_root):
        os.makedirs(out_root)

    #seed the random number generator
    numpy.random.seed(args.seed)

    data_parts = numpy.array(sorted(glob.glob(args.parts)))
#    peak_num = int(local['head']('-n', '2', data_parts[0]).strip().split()[1])

    #generate the lists of hyperparameter values to test
    topic_options = [int(elt) for elt in args.topics.split(',')]
    if args.search_len == -1:
        args.search_len = len(topic_options) * data_parts.shape[0]

    if args.random_topics is True:
        topics = numpy.random.choice(topic_options, size=args.search_len, replace=True)
    else:
        num_each = int(numpy.ceil(args.search_len/len(topic_options)))
        topics = numpy.array([topic_options] * num_each).T.flatten()[::-1]

    iter_options = [int(elt) for elt in args.iterations.split(',')]
    iters = numpy.random.choice(iter_options, size=args.search_len, replace=True)

    if '-' in args.alphas:
        alpha_low, alpha_high = numpy.array(args.alphas.split('-'), dtype=float)
        if args.logsearch is True:
            alphas = (numpy.log10(alpha_high) - numpy.log10(alpha_low)) * numpy.random.random_sample(size=args.search_len) + numpy.log10(alpha_low)
            alphas = 10 ** alphas
        else:
            alphas = (alpha_high - alpha_low) * numpy.random.random_sample(size=args.search_len) + alpha_low
    else:
        alphas = [float(args.alphas)] * args.search_len

    if '-' in args.betas:
        beta_low, beta_high = numpy.array(args.betas.split('-'), dtype=float)
        if args.logsearch is True:
            betas = (numpy.log10(beta_high) - numpy.log10(beta_low)) * numpy.random.random_sample(size=args.search_len) + numpy.log10(beta_low)
            betas = 10 ** betas
        else:
            betas = (beta_high - beta_low) * numpy.random.random_sample(size=args.search_len) + beta_low
    else:
        betas = [float(args.betas)] * args.search_len

    if len(data_parts) > 1:
        cmd_template = 'java -Xms{mem!s}G -cp lda_commandline/dist10/LatentDirichletAllocation.jar org.rhwlab.lda.cache.matrix.LDA_CommandLine -lda -a {alpha!s} -b {beta!s} {train!s} -li {lda_iters!s} -o {out_dir!s} -s {seed!s} -t {topic!s} -th {threads!s} -tn {thinning!s} -ch {cache!s} -rid {run_id!s} -d {dist!s} -st {statistic!s} -sk {skip!s} -v 1 -id {phi_dir!s} -pr {kde_precision!s} -chib -ic {test!s} -ip {input_phi!s} -ci {chib_iters!s} -cb {chib_burn!s};'
    else:
        cmd_template = 'java -Xms{mem!s}G -cp lda_commandline/dist10/LatentDirichletAllocation.jar org.rhwlab.lda.cache.matrix.LDA_CommandLine -lda -a {alpha!s} -b {beta!s} {train!s} -li {lda_iters!s} -o {out_dir!s} -s {seed!s} -t {topic!s} -th {threads!s} -tn {thinning!s} -ch {cache!s} -rid {run_id!s} -pe -d {dist!s} -st {statistic!s} -sk {skip!s} -v 1 -id {phi_dir!s} -pr {kde_precision!s};'

#    if len(data_parts) > 1:
#        cmd_template = 'java -Xms{mem!s}G -cp lda_commandline/dist9/LatentDirichletAllocation.jar org.rhwlab.lda.cache.LDA_CommandLine -lda -a {alpha!s} -b {beta!s} {train!s} -li {lda_iters!s} -o {out_dir!s} -s {seed!s} -t {topic!s} -th {threads!s} -tn {thinning!s} -ch {cache!s} -rid {run_id!s} -d {dist!s} -st {statistic!s} -sk {skip!s} -v 1 -id {phi_dir!s} -pr {kde_precision!s} -chib -ic {test!s} -ip {input_phi!s} -ci {chib_iters!s} -cb {chib_burn!s};'
#    else:
#        cmd_template = 'java -Xms{mem!s}G -cp lda_commandline/dist9/LatentDirichletAllocation.jar org.rhwlab.lda.cache.LDA_CommandLine -lda -a {alpha!s} -b {beta!s} {train!s} -li {lda_iters!s} -o {out_dir!s} -s {seed!s} -t {topic!s} -th {threads!s} -tn {thinning!s} -ch {cache!s} -rid {run_id!s} -pe -d {dist!s} -st {statistic!s} -sk {skip!s} -v 1 -id {phi_dir!s} -pr {kde_precision!s};'

    #submit jobs
    javamem = args.coremem * args.threads
    skip = args.lda_burnin//args.lda_thinning
    for idx, (iter_len, topic, alpha, beta) in enumerate(zip(iters, topics, alphas, betas)):
        idx += args.start_from
        test_set_idx = numpy.zeros(len(data_parts)).astype(bool)
        if len(data_parts) > 1:
            test_set_idx[idx % len(data_parts)] = True

#        alpha /= topic
#        beta /= peak_num

        run_id = '{:04d}'.format(idx)
        phi_dir = os.path.join(out_root, '{!s}_topics{!s}_alpha{:1.3f}_beta{:1.3f}'.format(run_id, topic, alpha, beta))
        phi_path = os.path.join(phi_dir, '{!s}_{!s}.phi'.format(args.pe_dist, args.pe_statistic))
#        outdir = os.path.join(out_root, run_id)

        if not os.path.isdir(phi_dir):
            os.makedirs(phi_dir)
        else:
            raise Exception('Improbably, outdir already exists: {!s}'.format(phi_dir))

        params = {'mem':javamem,
                  'alpha':alpha,
                  'beta':beta,
                  'topic':topic,
                  'train':'-ib ' + ' -ib '.join(data_parts[~test_set_idx]),
                  'out_dir':out_root,
                  'test':data_parts[test_set_idx][0] if numpy.any(test_set_idx) else None,
#                  'test_no_bow':os.path.basename(data_parts[test_set_idx][0].replace('.bow', '')),
                  'seed':args.seed + idx,
                  'lda_iters':iter_len,
                  'threads':args.threads,
                  'thinning':args.lda_thinning,
                  'cache':args.lda_cache,
                  'run_id':run_id,
                  'phi_dir':phi_dir,
                  'kde_precision':args.kde_precision,
                  'input_phi':phi_path,
                  'skip':skip,
                  'dist':args.pe_dist,
                  'statistic':args.pe_statistic,
                  'chib_iters':args.chib_iters,
                  'chib_burn':args.chib_burnin}
        cmd = cmd_template.format(**params)
        qsub_path = os.path.join(phi_dir, 'qsub.sh')
        with open(qsub_path, 'w') as out:
            out.write('#! /usr/bin/env bash\n\n')
            out.write('source activate python3;\n\n')
            out.write(cmd + '\n\n')
            if args.cleanup is True:
                out.write('rm {!s};\n\n'.format(os.path.join(phi_dir, '*.Z')))
                out.write('rm {!s};\n\n'.format(os.path.join(phi_dir, '*.xml')))
                out.write('rm {!s};\n\n'.format(os.path.join(phi_dir, '*.random')))
        os.chmod(qsub_path, 0o774)

        qsub_cmd = ['qsub', '-j', 'y', 
                    '-N', '{!s}_{!s}_{:1.3f}_{:1.3f}'.format(args.job_name, topic, alpha, beta), 
                    '-o', os.path.join(phi_dir, 'job.out'),
                    '-cwd',
                    '-P', args.uge_proj,
                    '-l', 'mfree={!s}G'.format(args.coremem), 
                    '-pe', 'serial', str(args.threads), 
                    qsub_path]

#        qstat_cmd = local['qstat']['-u', pwd.getpwuid(os.getuid())[0], '-j', 'ldahyp_*'] | local['grep']['job_number:']
        qstat_cmd = (local['qstat']['-u', pwd.getpwuid(os.getuid())[0], 
                                    '-j', '{!s}_*'.format(args.job_name)] 
                     | local['grep']['-P', 'project:\s*{!s}'.format(args.uge_proj)])
        while len(qstat_cmd(retcode=None).strip().split('\n')) > args.num_jobs - 1:
            time.sleep(10)

        subprocess.check_call(qsub_cmd)
