#!/usr/bin/env python
"""
Deconvolution for antibody pools
joel.bader@jhu.edu
"""

import logging
import os
from gpr import GPR
from dataframe import DataFrame
import numpy as np
import numpy.ma

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)

def get_data_dir():
    """ directory holding the pooled experiment results """
    # data_dir = '/Users/joel/Dropbox/Pooled data and individual retests_12511/Pools'
    data_dir = '../data'
    # data_dir = '/Users/joel/Dropbox/GPR files'
    logger.info('data_dir %s', data_dir)
    return(data_dir)
    
def get_results_dir():
    """ directory to write the results """
    results_dir = '../results'
    results_dir = '/Users/joel/Dropbox/GPR files/results'
    logger.info('results_dir %s', results_dir)
    return(results_dir)

def get_channel():
    """ foreground and background channel """
    (channel_fg, channel_bg) = [ 'F635 Median', 'B635 Median' ]
    return(channel_fg, channel_bg)
    
def get_pool_filename():
    """ file with pool-file map """
    pool_filename = '../data/pool_to_file.txt'
    
    return(pool_filename)

def get_naive(fg, bg):
    """
    calculate ratio as fg / bg
    calculate mean and stdev of ratio
    calculate z-score as (ratio - mean)/stdev
    return ratio and zscore
    """
    ratio = np.array(fg, dtype=float) / np.array(bg, dtype=float)
    mean = ratio.mean()
    stdev = ratio.std(ddof = 1) # subtract 1 ddof for mean to be consistent with excel stdev
    zscore = (ratio - mean)/stdev
    return(ratio, zscore)
    
def get_mean_by_group(group_list, data_list):
    """
    group_list provides the group for each row in data_list
    for each group, calculate the mean of the corresponding rows
    return:
    the means for each key (len = number of distinct keys)
    the corresponding mean for each row
    """
    group_to_sum = dict()
    group_to_cnt = dict()
    for (i, grp) in group_list.iteritems():
        group_to_sum[grp] = group_to_sum.get(grp, 0) + data_list[i]
        group_to_cnt[grp] = group_to_
        for i in range(len(zscore_naive)):
            if not mask[i]:
                this_id = id[i]
                id_to_sum[this_id] = id_to_sum.get(this_id, 0) + zscore_naive[i]
                id_to_cnt[this_id] = id_to_cnt.get(this_id, 0) + zscore_naive[i]
        id_to_mean = dict()
        for this_id in id_to_sum.keys():
            id_to_mean[this_id] = float(id_to_sum[this_id]) / float(id_to_cnt[this_id])
    
        zscore_mean = [ id_to_mean.get(this_id, np.nan) for this_id in id ]

def extract_by_group(row_group, row_data):
    """
    return a dict indexed by unique values in row_group
    dict value = list of values in row_data belonging to this group
    """
    data_by_group = dict()
    for (g, d) in zip(row_group, row_data):
        if g not in data_by_group:
            data_by_group[g] = [ ]
        # data_by_group[g] = data_by_group.get(g, [ ]) + [ d ]
        data_by_group[g].append(d)
    return(data_by_group)
        
def apply_by_group(fn, row_group, row_data):
    """
    apply function fn to groups defined by row_group with data row_data
    fn reduces list to a scalar
    """
    data_by_group = extract_by_group(row_group, row_data)
    fn_by_group = dict()
    for (grp, data_list) in data_by_group.items():
        fn_by_group[grp] = fn(data_list)
    fn_by_row = [ fn_by_group[x] for x in row_group ]
    return(fn_by_group, fn_by_row, data_by_group)

def get_good_ids_rows(id_list, zscore_list, z_threshold = 2.5):
    """
    ad hoc definition: retain rows where mask is good and zscore is above a threshold of 2.5
    then retain the ids corresponding to these rows
    """
    n_row = len(id_list)
    row_num = range(1, n_row + 1)
    id_subset = [ ]
    row_subset = [ ]
    for (r, i, z) in zip(row_num, id_list, zscore_list):
        if z < z_threshold:
            continue
        row_subset.append(r)
        if i not in id_subset:
            id_subset.append(i)
    return(id_subset, row_subset)

def process_gpr_file(input_file, control_cnt):
    """
    open input_file as a gpr
    columns Flags == -100 marks a control
    update control_cnt to count how often each (id, name) pair is/isnot a control
    """
    FLAG_BAD = -100
    logger.info('%s => %s', input_file, output_file)
    gpr = GPR(input_file)
    gpr.print_summary()

    # extract the columns we need: id, name, flags
    (id, name, flags) = gpr.get_columns(['ID', 'Name', 'Flags'])
    n_row_orig = len(flags)
    logger.info('n_row_orig %d', n_row_orig)

    # do all the work
    for (this_id, this_name, this_flag) in zip(id, name, flag):
        status = 'control' if this_flag <= FLAG_BAD else 'exptl'
        key = (this_id, this_name)
        if key not in control_cnt:
            control_cnt[key] = dict()
            control_cnt[key]['control'] = 0
            control_cnt[key]['exptl'] = 0
        control_cnt[key][status] = control_cnt[key][status] + 1
        
def process_gpr_dir(data_dir):
    """
    process each gpr file in the data_dir, writing results to results_dir
    keep track of ids and names that are used as controls
    for each (id, name) pair
        count how often the id is/isnot a control
        count how often name is/isnot a control
    """
    file_list = sorted(os.listdir(data_dir))
    control_cnt = dict()
    for file_name in file_list:
        (base, ext) = os.path.splitext(file_name)
        if (ext == '.gpr') or (ext == '.GPR'):
            logger.info('dir %s file %s base %s ext %s', data_dir, file_name, base, ext)
            input_file = os.path.join(data_dir, file_name)
            logger.info('input %s', input_file)
            process_gpr_file(input_file, control_cnt)
            


def main():
    
    # for each gpr file in the data directory,
    #   analyze the file and generate results for that file
    data_dir = get_data_dir()
    control_list = process_gpr_dir(data_dir)
    control_list.write('controls.txt')


if __name__ == '__main__':
    main()