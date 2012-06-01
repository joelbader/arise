#!/usr/bin/env python
"""
list_to_grid
copyright (c) 2012
joel.bader@jhu.edu
"""
import os
import subprocess
import copy
from dataframe import DataFrame

import logging
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)


def make_grid_for_file(results_dir, list_file, grid_file):
    POOL_DIRECTIONS = ['H', 'V']
    POOL_RANGE = range(1, 13)
    logger.info('%s => %s', list_file, grid_file)
    df = DataFrame(filename=os.path.join(results_dir, list_file))
    (pair, id) = df.get_columns('pair', 'id')
    logger.info('%d hits', len(pair))
    col_names = [ 'V' + str(i) for i in POOL_RANGE ]
    row_names = [ 'H' + str(i) for i in POOL_RANGE ]
    data_dict = dict() # will hold a list of the hits for each row, column pair
    for r in row_names:
        for c in col_names:
            data_dict[(r,c)] = [ ]
    for (mypair, myid) in zip(pair, id):
        (horiz, vert) = mypair.split(' x ')
        data_dict[(horiz, vert)] = data_dict[(horiz, vert)] + [ myid ]
    # now build a new data frame as a list of tuples, column name and column list
    data_by_column = [ ]
    # first column is the row names
    data_by_column += [(grid_file, row_names)]
    # subsequent columns are by vertical pool
    for c in col_names:
        col_data = [ ]
        for r in row_names:
            col_data.append(' '.join(sorted(data_dict[(r,c)])))
        data_by_column += [ (c, col_data)]
    grid_dataframe = DataFrame(data=data_by_column)
    grid_dataframe.write(os.path.join(results_dir, grid_file))
    

def match_prefix(mystr, prefix):
    if (len(mystr) < len(prefix)):
        return False
    if mystr[0:len(prefix)] != prefix:
        return False
    return True

def get_files_with_prefix(directory, prefix):
    all_files = sorted(os.listdir(directory))
    match_files = [ f for f in all_files if match_prefix(f, prefix)]
    return match_files

def make_grid_for_dir(results_dir, old_prefix, new_prefix):
    # get the list of files to process
    logger.info('working on directory %s', results_dir)
    match_files = get_files_with_prefix(results_dir, old_prefix)
    for list_file in match_files:
        grid_file = new_prefix + list_file[ len(old_prefix) : ]
        make_grid_for_file(results_dir, list_file, grid_file)

def main(results_dir, old_prefix, new_prefix):
    make_grid_for_dir(results_dir, old_prefix, new_prefix)

if __name__ == '__main__':
    results_dir = '/Users/joel/Dropbox/GPR_files/results'
    old_prefix = 'intersection_hit'
    new_prefix = 'intersection_grid'
    main(results_dir, old_prefix, new_prefix)
