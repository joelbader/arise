#!/usr/bin/env python

import logging

#logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logging.basicConfig(format='%(funcName)s: %(message)s')
logger = logging.getLogger(name='deconv')
logger.setLevel(logging.INFO)

def mysub(*args):
    print str(len(args)) + ' args'
    for a in args:
        print a

def main():
    mysub(1,2,'a','b')
    logger.error('throwing an error')
    print 'after the error'

if __name__ == '__main__':
    main()
    
