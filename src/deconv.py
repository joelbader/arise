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

import argparse # command line arguments

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)

ap = argparse.ArgumentParser(description ='Deconvolute protein array pooling data.', \
                             epilog = 'copyright (c) 2012 joel.bader@jhu.edu')
ap.add_argument('data_dir', help='directory for reading gpr data files')
ap.add_argument('results_dir', help='directory for writing results')
ap.add_argument('--control_filename', default='~/Dropbox/deconv/src/control_seth_2012_04_26.txt', \
                help='file with controls, header "id name" then one row for each id and name (default: %(default)s)')
ap.add_argument('--create_map', action='store_true', help='create pool-to-file map by parsing gpr filenames (default: %(default)s)')
ap.add_argument('--pool_filename', default = 'pool_to_file.txt', help='RESULTS_DIR/POOL_FILENAME has the pool-to-file map (default: %(default)s)')
ap.add_argument('--signal_fg', default = 'F635 Median', help='gpr signal foreground (default: %(default)s)')
ap.add_argument('--signal_bg', default = 'B635 Median', help='gpr signal background (default: %(default)s)')
ap.add_argument('--norm_fg', default = 'F532 Median', help='gpr normalization foreground (default: %(default)s)')
ap.add_argument('--norm_bg', default = 'B532 Median', help='gpr normalization background (default: %(default)s)')
ap.add_argument('--normalize', action='store_true', help='normalize signal fg/bg by norm fg/bg (default: %(default)s)')
ap.add_argument('--log', action='store_true', help='take the log before calculating z-scores (default: %(default)s)')

def get_data_dir():
    """ directory holding the pooled experiment results """
    # data_dir = '/Users/joel/Dropbox/Pooled data and individual retests_12511/Pools'
    data_dir = '../data'
    data_dir = '/Users/joel/Dropbox/GPR files'
    data_dir = '/Users/joel/Dropbox/GPR files/2012-04-25-IgM'
    data_dir = '/Users/joel/Dropbox/GPR files/2012-05-03 Validations'    
    logger.info('data_dir %s', data_dir)
    return(data_dir)
    
def get_results_dir(data_dir=None):
    """ directory to write the results """
    results_dir = '../results'
    results_dir = '/Users/joel/Dropbox/GPR files/results'
    if (data_dir is not None):
        results_dir = os.path.join(data_dir, 'results')
    logger.info('results_dir %s', results_dir)
    return(results_dir)

def get_channel():
    """ foreground and background channel """
    (channel_fg, channel_bg) = [ 'F635 Median', 'B635 Median' ]
    return(channel_fg, channel_bg)
    
def get_control_filename():
    """
    file with four columns: id, name, control, expt
    counts how many times this (id, name) pair is marked as a control
    """
    filename = 'control_seth_2012_04_26.txt'
    return(filename)

def get_control_from_file(filename, simple=True):
    """
    read the file as a data frame
    for each id, check how many times it occurs as control or experimental
    make a dict with (id, name) as key where pair is often as controls, or name is nd
    """
    logger.info('reading controls from %s', filename)
    control = DataFrame(filename=filename)
    control_dict = dict()
    
    if (simple):
        (ids, names) = control.get_columns('id','name')
        for (id, name) in zip(ids, names):
            control_dict[(id,name)] = True
    else:
        (id, name, control, exptl) = control.get_columns('id', 'name', 'control', 'exptl')
        id_to_name = dict()
        for (i, n, c, e) in zip(id, name, control, exptl):
            isND = n in [ 'ND', 'nd', 'N.D.' ]
            isControl = (i == 'CONTROL')
            isIgg = (n == 'IgG')
            if ((c >= e) or isND or isControl or isIgg):
                control_dict[(i, n)] = True
                
        # insert some special cases
        control_dict[('CONTROL', 'IgG')] = True
    
        for (i, n) in zip(id, name):
            if i not in id_to_name:
                id_to_name[i] = dict()
            id_to_name[i][n] = True
        
        id_to_names = dict()
        for i in id_to_name:
            names = sorted(id_to_name[i].keys())
            cnt = len(names)
            name_str = ','.join(names)
            id_to_names[i] = dict()
            id_to_names[i]['cnt'] = cnt
            id_to_names[i]['names'] = name_str
        ids = sorted(id_to_names.keys())
        cnts = [ id_to_names[x]['cnt'] for x in ids ]
        names = [ id_to_names[x]['names'] for x in ids ]
        df = DataFrame(data=[ ('id', ids), ('cnt', cnts), ('names', names)])
        df.write('id_to_names.txt')
    return(control_dict)

def print_control_dict(control_dict):
    keys = sorted(control_dict.keys())
    ids = [ x[0] for x in keys ]
    names = [ x[1] for x in keys ]
    df = DataFrame(data=[('id', ids), ('name', names)])
    df.write(filename='control_dict.txt')
    

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

def process_gpr_file(input_file, output_file, summary_file, channel_fg, channel_bg, control_dict=None):
    """
    open input_file as a gpr
    extract columns corresponding to F635 Median and B635 Median (fore- and back-ground)
    add new column fg/bg ratio
    extract Flags column as a mask
    
    if control_dict is None:
        mask out values with Flags == -100
    else:
        mask out values based on control_dict
    
    calculate mean and standard deviation of the ratio
    calculate z-score for each row
    calculate stouffer's z-score ?or mean z-score? for probes with same ID
    print probes with (mean) z-score >= 2.5
    """
    FLAG_BAD = -100
    logger.info('%s => %s', input_file, output_file)
    gpr = GPR(input_file)
    # print debug information for a gpr file
    # gpr.print_summary()

    # keep track of which columns we've added
    columns_added = [ ]

    # start by extracting the flags and adding an index for the original row number
    (flags, ids, names, fg, bg) = gpr.get_columns(['Flags', 'ID', 'Name', channel_fg, channel_bg])
    n_row_orig = len(flags)
    logger.info('n_row_orig %d', n_row_orig)
    row_number_orig = np.array(range(1, n_row_orig + 1))
    
    gpr.add_columns( ('row_number_orig', row_number_orig))
    columns_added += ['row_number_orig']

    # identify rows with bad flags and delete them
    # follow the semantics of a numpy masked array: delete where mask is True
    
    # controls from a dictionary
    mask_control = [ False for x in ids ]
    if (control_dict is not None):
        control_ids = dict()
        # for controls, just worry about ID, not name
        for (i, n) in control_dict.keys():
            control_ids[i] = True
        mask_control = [ i in control_ids for i in ids ]
        
    # user interface permits manual flagging of bad data, usually -100
    mask_flag = flags <= FLAG_BAD
    
    # some text values are clearly controls
    mask_text = [ id == 'CONTROL' for id in ids ]
    
    # bad signal
    mask_signal = [ (x[0] <= 0) or (x[1] <= 0) for x in zip(fg, bg) ]
    
    mask = [ x[0] or x[1] or x[2] or x[3] for x in zip(mask_control, mask_flag, mask_text, mask_signal) ]
    logger.info('deleting %d control rows', sum(mask))
    gpr.delete_rows(mask)


    # re-extract just the good columns
    columns_extracted = [ 'Name', 'ID', channel_fg, channel_bg ]
    (name, id, fg, bg) = gpr.get_columns(columns_extracted)
    n_row = len(name)
    assert(sum(bg == 0) == 0), 'bg has %d zero values' % sum(bg==0)
    
    # create a new index, idname, combining id with name
    # this avoids having one id map to multiple names, which could reflect a difference in probes, etc.
    idname = [ '_'.join([i,n]) for (i, n) in zip(id, name) ]
    idname_to_id = dict()
    idname_to_name = dict()
    for (idn, i, n) in zip(idname, id, name):
        idname_to_id[idn] = i
        idname_to_name[idn] = n
    
    gpr.add_columns( ('idname', idname))
    columns_added += ['idname']
    
    (ratio_naive, zscore_naive) = get_naive(fg, bg)    
    (id_to_mean_naive, row_to_mean_naive, id_to_zscores) = apply_by_group(np.mean, idname, zscore_naive)
    (id_to_mean_ratio, row_to_mean_ratio, id_to_ratios) = apply_by_group(np.mean, idname, ratio_naive)

    gpr.add_columns(('ratio_naive', ratio_naive),
        ('zscore_naive', zscore_naive),
        ('zscore_mean_naive', row_to_mean_naive))
    columns_added += ['ratio_naive', 'zscore_naive', 'zscore_mean_naive' ]
    
    # collect rows where flag is good and either zscore is above a threshold
    (id_subset, row_subset) = get_good_ids_rows(idname, zscore_naive)
    
    columns_display = columns_extracted + columns_added
    gpr.write(output_file, rows=row_subset, columns=columns_display)
    
    # gather data for each good id:
    # id, name, zscore_mean, zscores
    id_list = [ idname_to_id[i] for i in id_subset ]
    name_list = [ idname_to_name[i] for i in id_subset ]
    zscore_list = [ id_to_mean_naive[i] for i in id_subset ]
    ratio_list = [ id_to_mean_ratio[i] for i in id_subset ]
    zscores_list = [ ';'.join([ str(x) for x in id_to_zscores[i] ]) for i in id_subset]
    ratios_list = [ ';'.join([ str(x) for x in id_to_ratios[i] ]) for i in id_subset]
    id_data = DataFrame( data=[('IDName', id_subset),
        ('ID', id_list), ('Name', name_list),
        ('zscore', zscore_list), ('ratio', ratio_list),
        ('zscores', zscores_list), ('ratios', ratios_list)] )
    id_data.write(summary_file)
        


def process_gpr_dir(data_dir, results_dir, channel_fg, channel_bg, control_dict):
    """ process each gpr file in the data_dir, writing results to results_dir """
    file_list = sorted(os.listdir(data_dir))
    for file_name in file_list:
        (base, ext) = os.path.splitext(file_name)
        if (ext == '.gpr') or (ext == '.GPR'):
            logger.info('dir %s file %s base %s ext %s', data_dir, file_name, base, ext)
            input_file = os.path.join(data_dir, file_name)
            output_file = os.path.join(results_dir, base + '-top.txt')
            summary_file = os.path.join(results_dir, base + '-summary.txt')
            logger.info('input %s output %s summary %s', input_file, output_file, summary_file)
            process_gpr_file(input_file, output_file, summary_file, channel_fg, channel_bg, control_dict)

POOL_DIRECTIONS = ['H', 'V']
POOL_RANGE = range(1, 13)

def parse_pool_name(p):
    direction = p[0]
    number = int(p[1:])
    return(direction, number)

def is_valid_pool_name(p):
    (direction, number) = parse_pool_name(p)
    ret = True
    if direction not in ['H', 'V']:
        ret = False
    if number not in POOL_RANGE:
        ret = False
    return(ret)
    
def is_available(f):
    ret = os.path.isfile(f)
    return(ret)

def validate_pools(pool_to_file):
    """ check that pool names are valid and that all the summary files exist """
    for (p, f) in zip(pool_to_file.data['pool'], pool_to_file.data['full_path']):
        assert(is_valid_pool_name(p)), 'bad pool name: %s' % p
        assert(is_available(f)), 'results unavailable: %s' % f
        logger.info('pool %s file %s validated', p, f)
    return True

def write_pool_hit(pool_to_file, pool_hit):
    pool_list = [ ]
    file_list = [ ]
    id_list = [ ]
    zscore_list = [ ]
    ratio_list = [ ]
    for (p, f) in zip(pool_to_file.data['pool'], pool_to_file.data['file']):
        for h in pool_hit[p]:
            pool_list.append(p)
            file_list.append(f)
            id_list.append(h)
            zscore_list.append(pool_hit[p][h]['zscore'])
            ratio_list.append(pool_hit[p][h]['ratio'])
    df = DataFrame( data=[('pool', pool_list), ('file', file_list), ('id', id_list), ('zscore', zscore_list), ('ratio', ratio_list)] )
    df.write('pool_hit.txt')

def get_pool_hit(pool_to_file):
    pool_hit = dict()
    for (p, f) in zip(pool_to_file.data['pool'], pool_to_file.data['full_path']):
        pool_data = DataFrame(filename=f)
        for (idname, zscore, ratio) in zip(pool_data.data['IDName'], pool_data.data['zscore'], pool_data.data['ratio']):
            if p not in pool_hit:
                pool_hit[p] = dict()
            assert(idname not in pool_hit[p]), 'pool %s file %f duplicated id %s' % (p, f, idname)
            pool_hit[p][idname] = dict()
            pool_hit[p][idname]['zscore'] = zscore
            pool_hit[p][idname]['ratio'] = ratio
    return pool_hit

def get_intersection_hit(pool_hit, horizontal_pools, vertical_pools):
    intersection_hit = dict()
    pair_order = [ ]
    THRESHOLD = 2.5
    for h in horizontal_pools:
        if h not in pool_hit:
            continue
        
        for v in vertical_pools:
            if v not in pool_hit:
                continue
            pair = h + ' x ' + v
            for id in pool_hit[h].keys():
                if id not in pool_hit[v].keys():
                    continue
                (zscoreh, zscorev) = [pool_hit[x][id]['zscore'] for x in [h, v]]
                (ratioh, ratiov) = [ pool_hit[x][id]['ratio'] for x in [h, v]]
                if (zscoreh < THRESHOLD) or (zscorev < THRESHOLD):
                    continue
                logger.info('%s %s %f %f', pair, id, zscoreh, zscorev)
                if pair not in pair_order:
                    pair_order.append(pair)
                if pair not in intersection_hit:
                    intersection_hit[pair] = dict()
                intersection_hit[pair][id] = dict()
                ptr = intersection_hit[pair][id]
                ptr['zscore_H'] = zscoreh
                ptr['zscore_V'] = zscorev
                ptr['ratio_H'] = ratioh
                ptr['ratio_V'] = ratiov
                        
    # convert to a data frame
    pair_list = [ ]
    id_list = [ ]
    zscoreh_list = [ ]
    zscorev_list = [ ]
    ratioh_list = [ ]
    ratiov_list = [ ]
    for pair in pair_order:
        for id in sorted(intersection_hit[pair].keys()):
            ptr = intersection_hit[pair][id]
            (zscoreh, zscorev, ratioh, ratiov) = (ptr['zscore_H'],
                                                  ptr['zscore_V'],
                                                  ptr['ratio_H'],
                                                  ptr['ratio_V'])
            pair_list.append(pair)
            id_list.append(id)
            zscoreh_list.append(zscoreh)
            zscorev_list.append(zscorev)
            ratioh_list.append(ratioh)
            ratiov_list.append(ratiov)
    intersection_hit_df = DataFrame(data = [ ('pair', pair_list), ('id', id_list),
        ('zscore_h', zscoreh_list), ('zscore_v', zscorev_list),
        ('ratio_h', ratioh_list), ('ratio_v', ratiov_list) ])
    
    return(intersection_hit, intersection_hit_df)
                        
                    
                    
def deconv_pools(results_dir, pool_to_file):
    # create full path to file
    full_path = [ os.path.join(results_dir, str(f) + '-summary.txt') for f in pool_to_file.data['file'] ]
    pool_to_file.add_columns( ('full_path', full_path) )
    
    # check that pool names are correct and that summary files exist
    validate_pools(pool_to_file)
    
    # for each pool, get the hits' zscore and ratio
    pool_hit = get_pool_hit(pool_to_file)
    write_pool_hit(pool_to_file, pool_hit)
    
    HORIZONTAL = [ 'H' + str(i) for i in range(1, 13) ]
    VERTICAL = ['V' + str(i) for i in range(1, 13) ]
    (intersection_hit_dict, intersection_hit_df) = get_intersection_hit(pool_hit, HORIZONTAL, VERTICAL)
    filename = os.path.join(results_dir, 'intersection_hit.txt')
    intersection_hit_df.write(filename=filename)

def main(args):
    

    
    # get a dictionary of controls
    control_dict = get_control_from_file(args.control_filename)
    print_control_dict(control_dict)
    
    # for each gpr file in the data directory,
    #   analyze the file and generate results for that file
    if doFiles:
        process_gpr_dir(args.data_dir, args.results_dir, args.signal_fg, args.signal_bg, \
                        args.norm_fg, args.norm_bg, args.normalize, args.log, \
                        control_dict)

    doDeconv = False
    if (doDeconv):
        pool_fullpath = os.path.join(args.results_dir, args.pool_filename)
        pool_to_file = DataFrame(filename=pool_fullpath)
        deconv_pools(args.results_dir, pool_to_file)

if __name__ == '__main__':
    args = ap.parse_args()
    main(args)