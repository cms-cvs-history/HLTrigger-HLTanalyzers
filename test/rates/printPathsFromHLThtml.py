#!/usr/bin/python

#import pprint
import sys,HTMLTableParser

def main(argv) :
    """
    arguments: infile
    """

    #infname = 'HLT_208_V3.html'
    infname = 'HLT_V9.html'
    alen = len(argv)
    if alen > 0: infname = argv[0]

    f = open(infname,'r')

    p = HTMLTableParser.TableParser()
    p.feed(f.read())
    f.close()
    
    count = 0
    for l in p.doc:
        for row in l:
            count += 1
            if count == 1: continue
            row0 = row[0]
            row1 = row[1].replace('OR','OR ')
            print row0+': "'+row1+'"'
            #print '"'+row1+'"'
            #for entry in row:

    #print p.doc # Get to the data

if __name__ == '__main__' :
    main(sys.argv[1:])
