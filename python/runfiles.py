""" @file runfiles.py
@brief select a list of data files to process

@author: Matthew Kerr <kerm@u.washington.edu>

$Header$

"""

import os, glob

class RunFiles:

   def __init__(self, inpath, runlist):
      self.inpath=inpath
      if not os.path.exists(inpath):
         raise Exception('path "%s" is not valid'%inpath)
      if not os.path.exists(runlist):
         raise Exception('run list "%" not found'%runlist)
      self.runs=[line.strip() for line in open(runlist) if line[0]!='#' and len(line)>4]
      if len(self.runs)==0:
         print 'warning: no runs specified'
      

   def __call__(self, filetype=None):
      if filetype is not None: return self.__get_files__(filetype)
      ev_files=self.__get_files__('ph')
      sc_files=self.__get_files__('pt')
      return [ [ev_files[i],sc_files[i] ] for i in xrange(len(self.runs)) ]

   def __get_files__(self,filetype):
      import glob
      
      files=[0]*len(self.runs)
      for i in xrange(len(self.runs)):
         cands=glob.glob('%s*%s*%s*'%(self.inpath,filetype,self.runs[i]))
         if len(cands)==0:
            print 'warning: no files found with filetype %s '% filetype
            continue
         cands.sort()
         files[i]=cands[-1]
      return files

if __name__=='__main__':
   f=RunFiles(r'f:/glast/downloads/', r'f:/glast/data/first_light/runlist.txt')
   files=f.get_files()
   for x in files:
      print 'Event file: %s\n\tPointing file: %s'%tuple(x)

