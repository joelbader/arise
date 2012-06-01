#!/usr/bin/env python
"""
deconv driver
copyright (c) 2012
joel.bader@jhu.edu
"""
import os
import subprocess
import copy

import logging
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)


dropbox_base = '/Users/joel/Dropbox'
deconv_cmd = os.path.join(dropbox_base, 'deconv', 'src', 'deconv.py')
gpr_base = os.path.join(dropbox_base, 'GPR_files')
if __name__ == '__main__':
    data_subdirs = ['2012-03-05', '2012-03-05_IgG_FF', '2012-04-25-IgG', \
                   '2012-04-25-IgG_FF', '2012-04-25-IgM', '2012-05-03 Validations']
#    data_subdirs = ['2012-04-25-IgG_FF']
    args_base = ['python', deconv_cmd,'--create_map']
    for subdir in data_subdirs:
        for norm_flag in ['', '--do_norm']:
            for log_flag in ['', '--do_log']:
                print 'norm_flag ' + norm_flag
                print 'log_flag ' + log_flag
                args = copy.deepcopy(args_base)
                print args
                if norm_flag != '':
                    args.append(norm_flag)
                if log_flag != '':
                    args.append(log_flag)
                data_dir = os.path.join(gpr_base, subdir)
                results_dir = os.path.join(gpr_base, 'results', subdir)
                if '--do_norm' in args:
                    results_dir =  results_dir + '_norm'
                if '--do_log' in args:
                    results_dir = results_dir + '_log'
                args = args + [ data_dir, results_dir ]
                cmd_str = ' '.join(args)
                print '***\n' + cmd_str + '\n***\n'
                retval = subprocess.call(args)
                print 'return value: ' + str(retval)        