#!/usr/bin/env python
"""
Deconvolution for antibody pools
joel.bader@jhu.edu
"""

import logging
import os
from xlrd import open_workbook, cellname

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)

def get_expected_pools():
    expected_pools = [ ]
    for char in ('H', 'V'):
        for i in range(12):
            pool_name = char + str(i+1)
            expected_pools.append(pool_name)
    logger.info('expected pools: %s', ' '.join(expected_pools))
    return(expected_pools)

def get_data_dir():
    """ directory holding the pooled experiment results """
    data_dir = '/Users/joel/Dropbox/Pooled data and individual retests_12511/Pools'
    logger.info('data_dir %s', data_dir)
    return(data_dir)

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
    
def get_header(sheet):
    
    
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
            

def main():
    expected_pools = get_expected_pools()
    data_dir = get_data_dir()
    pool_to_file = get_pool_to_file(data_dir, expected_pools)
    analyze_pools(expected_pools, pool_to_file, data_dir)

if __name__ == '__main__':
    main()