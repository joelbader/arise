#!/usr/bin/env python
"""
Deconvolution for antibody pools
joel.bader@jhu.edu
"""

import logging
import os
from gpr import GPR
from gpr import DataFrame
import numpy as np
import numpy.ma

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)

# code for pools is possibly dead
# i originally thought that pools would be encoded in the file name
# in fact, there will actually be a mapping table from file to pool (or vice versa)
# begin dead code

def get_expected_pools():
    expected_pools = [ ]
    for char in ('H', 'V'):
        for i in range(12):
            pool_name = char + str(i+1)
            expected_pools.append(pool_name)
    logger.info('expected pools: %s', ' '.join(expected_pools))
    return(expected_pools)

def is_valid_pool_name(name):
    """ a valid pool name is the letter H or V followed by an integer 1 thoguh 12 """
    n = len(name)
    first = name[0]
    last = name[1:n]
    ret = first in ['H', 'V']
    if ret:
        ret = last.isdigit()
    if ret:
        ret = (1 <= int(last)) and (int(last) <= 12)
    return(ret)

def get_pool_to_file(data_dir, expected_pools):
    """ create a map from a pool to the file that has the data """
    pool_to_file = dict()
    file_list = os.listdir(data_dir)
    for file_name in file_list:
        (base, ext) = os.path.splitext(file_name)
        toks = base.split('-')
        # the pool name should be the final token
        pool_name = toks[-1]
        # if is_valid_pool_name(pool_name):
        if pool_name in expected_pools:
            pool_to_file[pool_name] = file_name
            logger.info('pool %s file %s', pool_name, file_name)
        else:
            logger.info('%s is not a pool', file_name)
    return(pool_to_file)
    
def analyze_file(full_path):
    wb = open_workbook(full_path, on_demand=True)
    logger.info('%s has %d sheets', full_path, wb.nsheets)
    sheet = wb.sheet_by_index(0)
    (nr, nc) = (sheet.nrows, sheet.ncols)
    logger.info('sheet name %s nr %d nc %d', sheet.name, nr, nc)
    header = get_header(sheet)
    for row in range(10):
        vals = [ ]
        for col in range(nc):
            vals.append( '-'.join([cellname(row, col), sheet.cell(row, col).value]))
        print " ".join(vals)
    
def analyze_pools(expected_pools, pool_to_file, data_dir):
    for pool_name in expected_pools:
        if pool_name in pool_to_file:
            file_name = pool_to_file[pool_name]
            full_path = os.path.join(data_dir, file_name)
            analyze_file(full_path)

# end dead code            

def get_data_dir():
    """ directory holding the pooled experiment results """
    # data_dir = '/Users/joel/Dropbox/Pooled data and individual retests_12511/Pools'
    data_dir = '../data'
    logger.info('data_dir %s', data_dir)
    return(data_dir)
    
def get_results_dir():
    """ directory to write the results """
    results_dir = '../results'
    logger.info('results_dir %s', results_dir)
    return(results_dir)

def get_channel():
    """ foreground and background channel """
    (channel_fg, channel_bg) = [ 'F635 Median', 'B635 Median' ]
    return(channel_fg, channel_bg)

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
        id_subset.append(i)
    return(id_subset, row_subset)

def process_gpr_file(input_file, output_file, channel_fg, channel_bg):
    """
    open input_file as a gpr
    extract columns corresponding to F635 Median and B635 Median (fore- and back-ground)
    add new column fg/bg ratio
    extract Flags column as a mask
    mask out values with Flags == -100
    calculate mean and standard deviation of the ratio
    calculate z-score for each row
    calculate stouffer's z-score ?or mean z-score? for probes with same ID
    print probes with (mean) z-score >= 2.5
    """
    FLAG_BAD = -100
    logger.info('%s => %s', input_file, output_file)
    gpr = GPR(input_file)
    gpr.print_summary()

    # keep track of which columns we've added
    columns_added = [ ]

    # start by extracting the flags and adding an index for the original row number
    (flags,) = gpr.get_columns(['Flags'])
    n_row_orig = len(flags)
    logger.info('n_row_orig %d', n_row_orig)
    row_number_orig = np.array(range(1, n_row_orig + 1))
    
    gpr.add_columns( ('row_number_orig', row_number_orig))
    columns_added += ['row_number_orig']

    # identify rows with bad flags and delete them
    # follow the semantics of a numpy masked array: delete where mask is True
    mask = flags <= FLAG_BAD
    logger.info('deleting %d rows with flag <= %d', sum(mask), FLAG_BAD)
    gpr.delete_rows(mask)
    
    # re-extract just the good columns
    columns_extracted = [ 'Name', 'ID', channel_fg, channel_bg, 'Flags' ]
    (name, id, fg, bg, flags) = gpr.get_columns(columns_extracted)
    n_row = len(name)
    assert(sum(bg == 0) == 0), 'bg has %d zero values' % sum(bg==0)
    
    # group rows by id
    id_to_name = gpr.get_id_to_name()
    
    (ratio_naive, zscore_naive) = get_naive(fg, bg)    
    (id_to_mean_naive, row_to_mean_naive, id_to_zscores) = apply_by_group(np.mean, id, zscore_naive)
    (id_to_mean_ratio, row_to_mean_ratio, id_to_ratios) = apply_by_group(np.mean, id, ratio_naive)

    gpr.add_columns(('ratio_naive', ratio_naive),
        ('zscore_naive', zscore_naive),
        ('zscore_mean_naive', row_to_mean_naive))
    columns_added += ['ratio_naive', 'zscore_naive', 'zscore_mean_naive' ]
    
    # collect rows where flag is good and either zscore is above a threshold
    (id_subset, row_subset) = get_good_ids_rows(id, zscore_naive)
    
    columns_display = columns_extracted + columns_added
    gpr.write(output_file, rows=row_subset, columns=columns_display)
    
    # gather data for each good id:
    # id, name, zscore_mean, zscores
    name_list = [ id_to_name[i] for i in id_subset ]
    zscore_list = [ id_to_mean_naive[i] for i in id_subset ]
    ratio_list = [ id_to_mean_ratio[i] for i in id_subset ]
    zscores_list = [ ';'.join([ str(x) for x in id_to_zscores[i] ]) for i in id_subset]
    ratios_list = [ ';'.join([ str(x) for x in id_to_ratios[i] ]) for i in id_subset]
    id_data = DataFrame( data=[('ID', id_subset), ('Name', name_list),
        ('zscore', zscore_list), ('ratio', ratio_list),
        ('zscores', zscores_list), ('ratios', ratios_list)] )
        


def process_gpr_dir(data_dir, results_dir, channel_fg, channel_bg):
    """ process each gpr file in the data_dir, writing results to results_dir """
    file_list = sorted(os.listdir(data_dir))
    for file_name in file_list:
        (base, ext) = os.path.splitext(file_name)
        if (ext == '.gpr') or (ext == '.GPR'):
            logger.info('dir %s file %s base %s ext %s', data_dir, file_name, base, ext)
            input_file = os.path.join(data_dir, file_name)
            output_file = os.path.join(results_dir, base + '-top.txt')
            process_gpr_file(input_file, output_file, channel_fg, channel_bg)

def main():
    
    # for each gpr file in the data directory,
    #   analyze the file and generate results for that file
    data_dir = get_data_dir()
    results_dir = get_results_dir()
    (channel_fg, channel_bg) = get_channel()
    process_gpr_dir(data_dir, results_dir, channel_fg, channel_bg)
    
    # pool_to_file = get_pool_to_file(data_dir, expected_pools)
    # analyze_pools(expected_pools, pool_to_file, data_dir)

if __name__ == '__main__':
    main()