#!/usr/bin/env python
"""
Utilities for gpr (GenePix Results) files
joel.bader@jhu.edu
"""

import logging
import os
import numpy as np

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger(name='gpr')
logger.setLevel(logging.INFO)

class GPR:
    """
    Utilities for GenePix Results (GPR) files
    """
    def __init__(self, filename):
        logger.info('reading from %s', filename)
        fp = open(filename, 'r')
        
        # line 1 should be 'ATF 1.0'
        line1 = fp.readline()
        toks = line1.strip().split()
        assert(len(toks) == 2), 'gpr line 1 should be ATF 1.0, got %s' % line1
        (self.file_type, self.version_number) = toks
        assert(self.file_type == 'ATF'), 'expecting file_type ATF got %s' % self.file_type
        assert(self.version_number == '1.0'), 'expecting version_number 1.0 got %s' % self.version_number
        
        # line 2 has number of optional header records and data fields (columns)
        line2 = fp.readline()
        toks = line2.strip().split()
        assert(len(toks) == 2), 'gpr line 2 should be <n_header> <n_column>, got %s' % line2
        (self.n_header, self.n_column) = (int(toks[0]), int(toks[1]))
        
        # process the header lines
        # expected format is "<key>=<value>" with double quotes and value possibly empty
        
        self.header_list = list() # to remember the order of the header keys
        self.header_dict = dict()
        for cnt in range(self.n_header):
            header_line = fp.readline()
            toks = header_line.strip().strip('"').split('=')
            assert(len(toks) == 2), 'header expected "<key>=<value>" got %s' % header_line
            (k, v) = toks
            self.header_list.append(k)
            self.header_dict[k] = v
            logger.info('header line %d %s = %s', cnt + 1, k, v)
            
        column_line = fp.readline()
        toks = column_line.strip().split('\t')
        assert(len(toks) == self.n_column), 'expected %d columns got %d: %s' % (self.n_column, len(toks), column_line)
        
        self.column_list = [ x.strip('"') for x in toks ]
        
        self.column_type = dict()
        for c in self.column_list:
            c_type = np.int
            if (c == 'Name') or (c == 'ID'):
                c_type = type('')
            elif (c == 'SNR 525') or (c == 'Log Ratio (1/525)'):
                c_type = np.float
            self.column_type[c] = c_type
        
        for (i, c) in enumerate(self.column_list):
            logger.info('%d\t%s\t%s', i+1, c, str(self.column_type[c]))
        
        # store each column as a separate object in a dict
        # most objects will be numpy int arrays, a few will be numpy float arrays
        self.data = dict()
        rows = fp.readlines()
        self.n_row = len(rows)
        logger.info('%d data rows', self.n_row)
        
        # initialize the data matrix
        for c in self.column_list:
            c_type = self.column_type[c]
            if c_type == type(''):
                self.data[c] = [''] * self.n_row
            elif c_type == np.int:
                self.data[c] = np.zeros(self.n_row, dtype=c_type)
            elif c_type == np.float:
                self.data[c] = np.array([np.nan] * self.n_row, dtype=c_type)
            else:
                assert 1 == 0, 'column %s has bad data type %s' % (c, str(c_type))
                
        # read the data
        for (i, row) in enumerate(rows):
            if (i + 1) % 10000 == 0:
                logger.info('... %d', i+1)
            toks = row.strip().split('\t')
            assert(len(toks) == self.n_column), 'expected %d columns got %d: %s' % (self.n_column, len(toks), row)
            for (j, c) in enumerate(self.column_list):
                c_type = self.column_type[c]
                tok = toks[j]
                if c_type == type(''):
                    value = tok.strip('"')
                elif c_type == np.int:
                    value = 0 if tok == '0.000' else c_type(tok)
                elif c_type == np.float:
                    value = np.nan if tok == 'Error' else c_type(tok)
                else:
                    assert 1 == 0, 'column %s has bad data type %s' % (c, str(c_type))
                self.data[c][i] = value
                
        return None
        
    def write(self, filename):
        fp = open(filename, 'w')
        logger.info('writing %d by %d data matrix to %s', self.n_row, self.n_column, filename)
        fp.write('\t'.join(self.column_list) + '\n')
        for i in range(self.n_row):
            if (i + 1) % 10000 == 0:
                logger.info('... %d', i+1)
            toks = [ str(self.data[c][i]) for c in self.column_list ]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
                
                
        

def main():
    gpr = GPR('../data/tmp.gpr')
    gpr.write('../data/tmp.txt')

if __name__ == '__main__':
    main()
