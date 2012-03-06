#!/usr/bin/env python
"""
Deconvolution for antibody pools
joel.bader@jhu.edu
"""

import logging
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.DEBUG)

def get_data_dir():
    data_dir = '/Users/joel/Dropbox/Pooled data and individual retests_12511/Pools'
    logger.info('data_dir %s', data_dir)
    return(data_dir)

def main():
    data_dir = get_data_dir()

if __name__ == '__main__':
    main()