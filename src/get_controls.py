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
    
def process_gpr_file(input_file, control_cnt):
    """
    open input_file as a gpr
    columns Flags == -100 marks a control
    update control_cnt to count how often each (id, name) pair is/isnot a control
    """
    FLAG_BAD = -100
    logger.info('%s', input_file)
    gpr = GPR(input_file)
    gpr.print_summary()

    # extract the columns we need: id, name, flags
    (ids, names, flags) = gpr.get_columns(['ID', 'Name', 'Flags'])
    n_row_orig = len(flags)
    logger.info('n_row_orig %d', n_row_orig)

    # do all the work
    for (id, name, flag) in zip(ids, names, flags):
        status = 'control' if flag <= FLAG_BAD else 'exptl'
        key = (id, name)
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
            
    # create a dataframe
    keys = sorted(control_cnt.keys())
    id = [ x[0] for x in keys ]
    name = [ x[1] for x in keys ]
    control = [ control_cnt[x]['control'] for x in keys ]
    exptl = [ control_cnt[x]['exptl'] for x in keys ]
    control_df = DataFrame(data= [ ('id', id), ('name', name), ('control', control), ('exptl', exptl) ] )
    return(control_df)
       


def main():
    
    # for each gpr file in the data directory,
    #   analyze the file and generate results for that file
    data_dir = get_data_dir()
    control_list = process_gpr_dir(data_dir)
    control_list.write('control.txt')


if __name__ == '__main__':
    main()