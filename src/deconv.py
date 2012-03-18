#!/usr/bin/env python
"""
Deconvolution for antibody pools
joel.bader@jhu.edu
"""

import logging
import os
from gpr import GPR
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
    (name, id, fg, bg, flags) = gpr.get_columns(['Name', 'ID', channel_fg, channel_bg, 'Flags'])
    assert(sum(bg == 0) == 0), 'bg has %d zero values' % sum(bg==0)
    ratio = np.array(1.0 * fg) / np.array(1.0 * bg)

    # mask flagged data when computing mean and stdev
    # numpy masked arrays exclude values where mask is True
    mask = flags <= FLAG_BAD
    logger.info('masking %d rows with flag <= %d', sum(mask), FLAG_BAD)
    ratio_masked = np.ma.masked_array(ratio, mask=mask)
    mean = ratio_masked.mean()
    stdev = ratio_masked.std(ddof = 1) # subtract 1 ddof for mean to be consistent with excel stdev

    zscore = (ratio - mean)/stdev
    # gpr.add_column('ratio', ratio)
    # gpr.add_column('zscore', zscore)
    
    # calculate the mean z-score by id, excluding masked z-scores
    id_to_sum = dict()
    id_to_cnt = dict()
    for i in range(len(zscore)):
        if not mask[i]:
            this_id = id[i]
            id_to_sum[this_id] = id_to_sum.get(this_id, 0) + zscore[i]
            id_to_cnt[this_id] = id_to_cnt.get(this_id, 0) + zscore[i]
    id_to_mean = dict()
    for this_id in id_to_sum.keys():
        id_to_mean[this_id] = float(id_to_sum[this_id]) / float(id_to_cnt[this_id])
    
    zscore_mean = [ id_to_mean.get(this_id, np.nan) for this_id in id ]
    gpr.add_columns(['ratio', 'zscore', 'zscore_mean'], ratio, zscore, zscore_mean)
    
    
    
    
    
    

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