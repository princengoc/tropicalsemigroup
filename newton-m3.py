#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 17:10:09 2018

@author: ngoc

UT2 for THREE LETTER
Analogue of newton.py for THREE letter alphabet. 

Functions to compute the path, degree-1 signatures of a word
their UT2 equivalence classes. 
***This code is NOT optimized***. Instead, its goal is to verify the counter-example to the Adjacent Swap conjecture for UT2 in three letters.

Words represented as integer strings (a = 0, b = 1, c = 2)

"""

import numpy as np
#import itertools as itt
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from collections import deque  
import util as util

def _path(w):
  """Input: w = string of integers 0/1/2
  Output: path of w"""
  path = np.array([0,0,0])
  s = np.array([0,0,0]) #last point
  for i in range(len(w)):
    increase = np.array([0,0,0])
    increase[w[i]] = 1
    s += increase
    path = np.vstack((path,s))
  return path

def _S(w):
  """Input: path of a word w
  Output: points in path belong to each specific letter.
  Comment: path naturally sorted, so taking last appearance should be ok."""
  gamma = _path(w)
  A = map(lambda x: sum(gamma[:,0] <= x)-1, xrange(max(gamma[:,0])))
  B = map(lambda x: sum(gamma[:,1] <= x)-1, xrange(max(gamma[:,1])))
  C = map(lambda x: sum(gamma[:,2] <= x)-1, xrange(max(gamma[:,2])))
  return((gamma[A],gamma[B],gamma[C]))

def _removeCollinear(Q):
  """Remove collinear points in a qhull object by checking change in volume after removing each point."""
  P = Q.points[Q.vertices]
  m = P.shape[0]
  vertices = []
  for i in range(m):
    Pprime = P[range(0,i)+range(i+1,m),]
    Qprime = ConvexHull(Pprime)
    if(abs(Q.volume-Qprime.volume) > 1e-10):
      vertices += [Q.vertices[i]]
  return(vertices)

def lowdimConvex(P):
  """Compute convex hull of a potentially low dimensional polytope."""
  pmat = np.matrix(P-P[0,])
  dim = np.linalg.matrix_rank(pmat)
  #if full-dim:
  if dim == pmat.shape[1]:
    Q = ConvexHull(P)
    return(Q.vertices)
  else:
    #two-dim: QR decomposition
    (u,s,v) = np.linalg.svd(pmat.T,full_matrices=False)
    #orthonormal basis: u2 = u[:,:2]
    #coefficients
    v2 = v[:dim,:]
    Q = ConvexHull(v2.T)
    #remove fake vertices due to numerical issues. 
    vertices = _removeCollinear(Q)
    #return: index of vertices
    return(vertices)

def signature(w):
  """ Input: w = a word in 0/1/2
  Output: the signature of w as list of vertices"""
  abc = _S(w)
  sig = map(lambda x: map(tuple,x[lowdimConvex(x)]), abc)
  #convert to sets for comparisons
  return(map(set,sig))

def charToWord(s):
  """Convert string s to word"""
  w = map(lambda x: 0 if x == 'a' else 1 if x == 'b' else 2,list(s))
  return(w)

def wordToChar(w):
  """Convert word to string"""
  s =''.join(map(lambda x: 'a' if x == 0 else 'b' if x == 1 else 'c', w))
  return(s)

def equal(w,v):
  """Returns True if two words have equal signature, false otherwise"""
  sigw = signature(w)
  sigv = signature(v)
  return(sigw == sigv)
  
#--- compute the adjacent swap tree. 
def _idloop(w,knownList):
  """
  Internal function for idloop in the module globalid
  Returns a list of w' which are one swap away from w, and which are NOT in knownList.
  Returns: newList.
  idloop is implemented with a global variable in globalid.py
  Note: knownList is a SET"""
  newList = []
  wshift = w[1:]
  wright = w[:-1]
  for i in range(len(wshift)):
    if wshift[i] != wright[i]:
      w2 = np.array(w)
      w2[i] = w[i+1]
      w2[i+1] = w[i]
      if equal(w,w2) and tuple(w2) not in knownList:
        newList += [tuple(w2)]
  return newList

global globeList
globeList = set()

def clear():
  """Clear the globeList."""
  global globeList
  globeList = set()

def addToGlobeList(val):
  global globeList
  globeList = globeList.union(set(val))

def idloop(w):
  """  Return a list of all w' that can be obtained as a sequence of 
    adjacent swaps, each of which is still equivalent to w. """
  newList = _idloop(w,globeList)
  addToGlobeList(newList)
  for w2 in newList:
    idloop(w2)

if False:
  """Counter example to the adjacent swap conjecture for UT2 in three letters."""
  w = charToWord('abccbaabcabcabccbaabc')
  v = charToWord('abccbaabccbaabccbaabc')
  equal(w,v)
  clear()
  addToGlobeList(set([tuple(w)]))  
  idloop(w)
  listw = globeList
  clear()
  addToGlobeList(set([tuple(v)]))  
  idloop(v)
  listv = globeList
  listv.intersection(listw)
  

#--- Next: check for adjacent swaps 

#-- double check with polymake. 
#abc = _S(v)
#A = abc[2]
#ones = np.ones((A.shape[0],1),dtype='int8')
#Ap = np.hstack((ones,A))
#print(Ap)
#COMMENT: has issue with _S(v)[2]

