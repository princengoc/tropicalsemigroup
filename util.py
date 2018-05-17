#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 00:01:16 2018

@author: ngoc

Utility functions.
"""

import numpy as np
import gzip
import cPickle as cPickle

#--- input/output functions
def _getWord(line):
  """Take a line, return a word w if it encodes a word.
  Otherwise, return None"""
  if len(line) < 3:
    return None
  else:
    w = map(int, list(line[0:-2]))
  return w

def getPairs(filename):
  """get pairs of words from a file. 
  Example file format: input3.txt
"""
  f = open(filename,'r')
  line = f.next()
  wordList = []
  while line is not None:
    w = _getWord(line)
    if w is not None:
      wordList += [w]
    try: 
      line = f.next()
    except:
      line = None
      f.close()
  #regroup in pairs: 
  m = len(wordList)
  wList = [wordList[i] for i in xrange(0,m,2)]
  w2List = [wordList[i] for i in xrange(1,m,2)]
  return (wList,w2List)

def save_zipped_pickle(obj, filename, protocol=-1):
  """Save an object as a zipped pickle file"""
  with gzip.open(filename, 'wb') as f:
    cPickle.dump(obj, f, protocol)

def load_zipped_pickle(filename):
  """Open a zipped pickle file"""
  with gzip.open(filename, 'rb') as f:
    loaded_object = cPickle.load(f)
    return loaded_object

#--- function to load dict
def getDict(na,n=20,prefix = 'pickle-zip/dict'):
  """Load a dictionary of UT2 identities. 
  Default name: pickle-zip/dict(na)(nb).p"""
  nb = n-na
  filename = prefix + str(na) + str(nb) +'.p'
  return load_zipped_pickle(filename)


#---- Small utility functions
def _reverseTest(w,w2):
  """Returns True if w is the reverse of w2 (also called twin). False otherwise"""
  w = np.array(w)
  w2 = np.array(w2)
  return all(np.flip(w,0) == w2) 

def isAdjacent(w,w2):
  """Returns true if w differs from w2 by one adjacent swap."""
  if ell1(w,w2) != 2:
    return False
  #find the first position they differ, and if adjacent, must also differ in the next.
  for i in range(len(w)):
    if w[i] != w2[i]:
      return w[i+1] != w2[i+1]  

#---- ell1 distance
def ell1(w,w2):
  """Returns the ell-1 distance between w and w2"""
  return sum(np.array(w) != np.array(w2))

#--- conversion between path, Sa/Sb and the word
def pathToSa(path,na):
  """Returns Sa from path with na's many a's"""
  S = [(i,max(path[path[:,0]==i,1])) for i in xrange(na)]
  return np.array(S)

def pathToSb(path,nb):
  """Returns Sb from path with nb's many b's"""
  T = [(max(path[path[:,1]==i,0]),i) for i in xrange(nb)]
  return np.array(T)

def SatoW(S,nb):
  """Returns word from S(a) and nb, the number of b's"""
  dS = [S[0,1]] + [S[i,1] - S[i-1,1] for i in range(1,len(S))]
  w = reduce(lambda x,y: x+[0]*y+[1], dS,[])
  w += [0]*(nb-S[len(S)-1,1])
  return w

def SbtoW(S,na):
  """Returns word from S(b) and na, the number of a's"""
  dS = [S[0,0]] + [S[i,0] - S[i-1,0] for i in range(1,len(S))]
  w = reduce(lambda x,y: x+[1]*y+[0], dS,[])
  w += [1]*(na-S[len(S)-1,0])
  return w  

def _vinP(v,P):
  """True iff v lies in the polytope P
  P = ConvexHull(S)
  v in P iff P.equations*(v,1) <= 1e-13"""
  return np.all(np.dot(P.equations,np.array([v[0],v[1],1])) <= 1e-13)

def _vindP(v,P):
  """True if v lies on the boundary of P"""
  b = np.dot(P.equations,np.array([v[0],v[1],1]))
  return np.sum(abs(b) <= 1e-13) > 0

