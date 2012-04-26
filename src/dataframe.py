#!/usr/bin/env python
"""
Utilities for data frames
joel.bader@jhu.edu
"""

import logging
import os
import numpy as np

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger(name='dataframe')
logger.setLevel(logging.INFO)

class DataFrame:
    """
    somewhat like an R data frame
    store a dictionary of numpy arrays
    each array must be the same length
    the column order is stored in the headers list
    """
    def __init__(self, data=None, filename=None, headers=None, sep='\t'):
        """
        initialize either from data or from a file
        if from data, data is a list of tuples (header_name, data_list)
        if from a file, extract the headers as the first line unless headers is not null
        """
        
        def is_type(str_val, typecast):
            try:
                typecast(str_val)
                return True
            except ValueError:
                return False
            
        def try_to_convert(str_val):
            new_val = str_val
            if is_type(str_val, int):
                new_val = int(str_val)
            elif is_type(str_val, float):
                new_val = float(str_val)
            return(new_val)
        
        n_err = 0
        if (data is None) and (filename is None):
            logger.error('data and filename both none')
            ++n_err
        if (data is not None) and (filename is not None):
            logger.error('data and filename both not none')
            ++n_err
        elif (data is not None) and (headers is not None):
            logger.error('data and headers both not none')
            ++n_err
        if (n_err > 0):
            logger.error('stopping after %d errors', n_err)
            assert(1==0)
        
        # initialize from a file by creating the same data format
        if filename is not None:
            logger.info('reading from %s', filename)
            fp = open(filename, 'r')
            if headers is None:
                header_line = fp.readline()
                headers = header_line.strip().split(sep)
            n_column = len(headers)
        
            # check that headers are unique
            header_dict = dict()
            for h in headers:
                assert(h not in header_dict), 'reused header %s' % h
                header_dict[h] = True

            # read the data        
            data_lines = fp.readlines()
            fp.close()
            data_tmp = dict()
            for h in headers:
                data_tmp[h] = [ ]
            for data_line in data_lines:
                toks = data_line.strip().split(sep)
                assert(len(toks) == n_column), 'expected %d columns found %d: %s' % (n_column, len(toks), data_line)
                for (j, h) in enumerate(headers):
                    new_val = try_to_convert(toks[j])
                    data_tmp[h].append(new_val)
            data = [ ]
            for h in headers:
                data.append( (h, data_tmp[h]))
        
        # initialize from data
        self.headers = [ ]
        self.data = dict()
        self.n_row = None
        for (h, data_list) in data:
            assert(h not in self.headers), 'repeated header %s' % h
            self.headers.append(h)
            if (self.n_row is None):
                self.n_row = len(data_list)
            else:
                assert(self.n_row == len(data_list)), 'column %s expected %d rows found %d' % (h, self.n_row, len(data_list))
            self.data[h] = np.array(data_list)
        self.n_column = len(self.headers)

    def get_columns(self, *args):
        """
        each argument is the name of a column
        return a list-of-lists, one list for each column requested
        error if the requested column does not exist
        """
        ret = [ ]
        for arg in args:
            assert(arg in self.headers), 'column name %s does not exist' % arg
            ret = ret + [ self.data[arg] ]
        return(ret)
        
        
    def add_columns(self, *args):
        """
        each argument is a tuple (column name, column data)
        error if a column name already exists
        error if number of data elements is not n_row
        modify
        n_column
        headers
        data
        """
        for (hdr, data_list) in args:
            assert(hdr not in self.headers), 'column name %s already exists' % this_name
            assert(len(data_list) == self.n_row), 'expected %d rows but found %d' % (self.n_row, len(this_data))
            self.n_column += 1
            self.headers.append(hdr)
            self.data[hdr] = np.array(data_list)
        n_new = len(args)
        logger.info('added %d columns: %s', n_new, ' '.join(self.headers[-n_new:]))
        
    def write(self, filename, rows=None, columns=None, sep='\t'):
        """
        rows is a list of row numbers, with the first list element being row 1 (not row 0)
        columns is a list of column headers
        """
        if rows is None:
            rows = range(1, self.n_row + 1)
        if columns is None:
            columns = self.headers
        fp = open(filename, 'w')
        n_col = len(columns)
        logger.info('writing %d by %d data frame to %s', len(rows), n_col, filename)
        fp.write('\t'.join(columns) + '\n')
        for i in rows:
            if i % 10000 == 0:
                logger.info('... %d', i+1)
            toks = [ str(self.data[c][i-1]) for c in columns ]
            fp.write('\t'.join(toks) + '\n')
        fp.close()         

def main():
    return 1

if __name__ == '__main__':
    main()
