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
                    data_tmp[h].append(toks[j])
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

class GPR:
    """
    Utilities for GenePix Results (GPR) files
    """
    def __init__(self, filename):
        """
        components:
        file_type 'ATF'
        version_number '1.0'
        header_list
        header_dict
        n_column
        column_list
        column_type
        n_row
        data
        """
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
            logger.debug('header line %d %s = %s', cnt + 1, k, v)
            
        column_line = fp.readline()
        toks = column_line.strip().split('\t')
        assert(len(toks) == self.n_column), 'expected %d columns got %d: %s' % (self.n_column, len(toks), column_line)
        
        self.column_list = [ x.strip('"') for x in toks ]
        
        self.column_type = dict()
        str_columns = ['Name', 'ID']
        int_columns = ['Block', 'Column', 'Row', 'X', 'Y', 'Dia.', 'F Pixels', 'B Pixels', 'Circularity', 'Flags', 'Normalize', 'Autoflag']
        for c in self.column_list:
            c_type = np.float
            if c in str_columns:
                c_type = type('')
            elif c in int_columns:
                c_type = np.int
            self.column_type[c] = c_type
        
        for (i, c) in enumerate(self.column_list):
            logger.debug('%d\t%s\t%s', i+1, c, str(self.column_type[c]))
        
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
                    value = c_type(tok)
                elif c_type == np.float:
                    value = np.nan if tok == 'Error' else c_type(tok)
                else:
                    assert 1 == 0, 'column %s has bad data type %s' % (c, str(c_type))
                self.data[c][i] = value
                
        return None
    
    def delete_rows(self, mask):
        """ delete rows where mask is true """
        assert(len(mask) == self.n_row), 'expected mask with %d rows but found %d' % (self.n_row, len(mask))
        row_subset = [ i for i in range(self.n_row) if not mask[i] ]
        for c in self.column_list:
            self.data[c] = np.array([ self.data[c][r] for r in row_subset ], dtype=self.column_type[c])
        self.n_row = len(row_subset)
        logger.info('%d rows deleted, new length is %d', sum(mask), self.n_row)        
        
    
    def get_columns(self, request_list):
        ret = [ ]
        for c in request_list:
            assert c in self.data, 'requested column %s missing' % c
            ret.append(self.data[c])
        return(ret)
        
    def add_columns(self, *args):
        """
        each argument is a tuple (column name, column data)
        error if a column name already exists
        error if number of data elements is not n_row
        modify
        n_column
        column_list
        column_type
        data
        """
        for (this_name, this_data) in args:
            assert(this_name not in self.column_list), 'column name %s already exists' % this_name
            assert(len(this_data) == self.n_row), 'expected %d rows but found %d' % (self.n_row, len(this_data))
            self.n_column += 1
            self.column_list.append(this_name)
            self.column_type[this_name] = type(this_data[0])
            self.data[this_name] = np.array(this_data)
        n_new = len(args)
        logger.info('added %d columns: %s', n_new, ' '.join(self.column_list[-n_new:]))
    
    def get_id_to_name(self):
        id_to_name = dict()
        for (i, n) in zip(self.data['ID'], self.data['Name']):
            if i not in id_to_name:
                id_to_name[i] = [ ]
            if n not in id_to_name[i]:
                id_to_name[i].append(n)
        for i in id_to_name:
            id_to_name[i] = ';'.join(id_to_name[i])
        return id_to_name
        
        
    def print_summary(self):
        # count how many ids for each name, how many rows for each name and id
        
        def create_hist(value_list):
            cnt = dict()
            for v in value_list:
                cnt[v] = cnt.get(v, 0) + 1
            return cnt
        
        def print_hist(hist, key_str, value_str):
            print '%s\t%s' % (key_str, value_str)
            for k in sorted(hist.keys()):
                v = hist[k]
                print '%d\t%d' % (k, v)

        id_to_mask = dict()
        id_to_name = dict()

        for masked in (False, True):       
            name_cnt = dict()
            id_cnt = dict()
            name_to_id = dict()
            for i in range(self.n_row):
                (name, id, flag) = (self.data['Name'][i], self.data['ID'][i], self.data['Flags'][i])
                mymask = flag <= -100
                if (mymask == masked):
                    if id not in id_to_mask:
                        id_to_mask[id] = dict()
                    id_to_mask[id][masked] = True
                    id_to_name[id] = name
                    name_cnt[name] = name_cnt.get(name, 0) + 1
                    id_cnt[id] = id_cnt.get(id, 0) + 1
                    if name not in name_to_id:
                        name_to_id[name] = dict()
                    name_to_id[name][id] = True
            name_to_idcnt = dict()
            for name in name_to_id.keys():
                name_to_idcnt[name] = len(name_to_id[name].keys())
            nameid_hist = create_hist(name_to_idcnt.values())
            name_hist = create_hist(name_cnt.values())
            id_hist = create_hist(id_cnt.values())
            print '\nhistogram for mask = ' + str(masked)
            print_hist(nameid_hist, 'ids_per_name', 'number_of_names')
            print_hist(name_hist, 'rows_per_name', 'number_of_names')
            print_hist(id_hist, 'rows_per_id', 'number_of_ids')
            
            print '\nnames with many ids for mask = ' + str(masked)
            for name in sorted(name_to_idcnt.keys()):
                cnt = name_to_idcnt[name]
                if (cnt < 6):
                    continue
                print '%s\t%d' % (name, cnt)
            
            print '\nids with many rows for mask = ' + str(masked)
            for id in sorted(id_cnt.keys()):
                cnt = id_cnt[id]
                if (cnt < 6):
                    continue
                print '%s\t%s\t%d' % (id_to_name[id], id, cnt)
        
        print 'checking for ids with multiple mask values'
        for id in sorted(id_to_mask.keys()):
            nkey = len(id_to_mask[id].keys())
            if (nkey != 1):
                name = id_to_name[id]
                print 'id %s named %s has %d masks' % (id, name, nkey)
                
    
    def write(self, filename, rows=None, columns=None):
        """
        rows is a list of row numbers, with the first list element being row 1 (not row 0)
        columns is a list of column headers
        """
        if rows is None:
            rows = range(1, self.n_row + 1)
        if columns is None:
            columns = self.column_list
        fp = open(filename, 'w')
        n_col = len(columns)
        logger.info('writing %d by %d data matrix to %s', len(rows), n_col, filename)
        fp.write('\t'.join(columns) + '\n')
        for i in rows:
            if i % 10000 == 0:
                logger.info('... %d', i+1)
            toks = [ str(self.data[c][i-1]) for c in columns ]
            fp.write('\t'.join(toks) + '\n')
        fp.close()
        

def main():
    gpr = GPR('../data/tmp.gpr')
    gpr.write('../data/tmp.txt')

if __name__ == '__main__':
    main()
