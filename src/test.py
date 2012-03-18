#!/usr/bin/env python

def mysub(*args):
    print str(len(args)) + ' args'
    for a in args:
        print a

def main():
    mysub(1,2,'a','b')

if __name__ == '__main__':
    main()
    
