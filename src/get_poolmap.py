#!/usr/bin/env python

def main():
    listing_file = 'list.txt'
    output_file = 'pool_to_file.txt'
    fp = open(listing_file, 'r')
    gp = open(output_file, 'w')
    gp.write('pool\tfile\n')
    for line in fp:
        toks = line.strip().split()
        assert(len(toks) == 1), 'expected one token: %s' % line
        tok = toks[0]    
        (prefix, suffix) = tok.split('.')
        (base, discard, pool) = prefix.split('-')
        gp.write('\t'.join([pool, base]) + '\n')
    fp.close()
    gp.close()

if __name__ == '__main__':
    main()
