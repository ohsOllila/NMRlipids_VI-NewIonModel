#!/usr/bin/env python
"""
  Reads in specific analysis output containing average Box dimensions
  and calculates APL (defaulting 64 lipids per leaflet)
"""

import sys, os # math # os.path, string, re

def main():
   from optparse import OptionParser
   # help message is automatically provided
   # type=string, action=store is default
   parser = OptionParser()
   parser.add_option('-f', '--infile', dest='infilename', help='input filename of *.xvg file', default='')
   parser.add_option('-n', '--nlip', dest='nLipids', help='nuber of lipids per leaflet', default='64', type=int)
   options, args = parser.parse_args(sys.argv[1:])  #argument 'sys.argv..' not necessary, it's the default

   nLipids=options.nLipids

   #open the inputfile
   if len(options.infilename)==0 :
      print 'Provide filename \n' ;  sys.exit(1)

   try: infile = open(options.infilename, 'r') # open file for reading 
   except: print '%s:  reading file "%s" failed; does it exist?' % (sys.argv[0],options.infilename) ; sys.exit(1)

   # get the dimensions
   lines = infile.readlines()
   for line in lines:
       if line.startswith("Box-X"):
           xdim=float(line.split()[1])
       if line.startswith("Box-Y"):
           ydim=float(line.split()[1])

   print "apl= (%6.3f*%6.3f)/%3d = %6.3f nm2" % (xdim,ydim,nLipids, xdim*ydim/nLipids)

   sys.exit(0)

if __name__ == '__main__':
    main()
