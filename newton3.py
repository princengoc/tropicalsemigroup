#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 16:48:18 2018

@author: ngoc

Code for UT3.
"""

import numpy as np
import itertools as itt
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import newton as newton

def _Saa(S):
  """ Input: S from newton._S(w,a)
  Output: the list of points for Saa in \R^4"""
  #take the pairs of S, with the constraint 
  #0 <= i < k < n_a
  #this is exactly the same as taking itt.combinations of S
  pairGen = itt.combinations(S,2)
  Saa = map(np.concatenate, pairGen)
  return Saa

def _Sab(Sa,Sb,a=1):
  """Input: Sa,Sb from newton._S -- see _newton3 for example usage.
  Output: the list of points for Sab in \R^4.
  This version uses the transformation of Proposition 2.6
  """    
  if a == 0:
    Sa = np.flip(Sa,1)
    Sb = np.flip(Sb,1) 
  #take combinations, then remove irrelevant points
  genPair = itt.product(Sa,Sb)
  points = map(np.concatenate,genPair)
  keep = [x for x in points if (x[0] < x[2]) & (x[1] <= x[3])]
  return keep

#--- compute from definition: quadratic in word length
def _getCrossSignature(w):
  """Input: word w. Output: points in \gamma_{ab},\gamma_{ba}
  This function uses a raw computation from definition (content between indices)
  """
  w = np.array(w)
  n = len(w)
  idxpair = itt.combinations(xrange(n),2)
  gamab = []
  gamba = []
  for rs in idxpair:
    r,s = rs
    if w[r] == 1 and w[s] == 0:
      br = sum(1-w[:r])
      ar = sum(w[:r])
      bs = sum(1-w[:s])
      ass = sum(w[:s])
      gamab += [(ar,br,ass-ar-1,bs-br)]
    if w[r] == 0 and w[s] == 1:
      br = sum(1-w[:r])
      ar = sum(w[:r])
      bs = sum(1-w[:s])
      ass = sum(w[:s])
      gamba += [(br,ar,bs-br-1,ass-ar)]      
  return (map(np.array,gamab),map(np.array,gamba))
    
#--- create function for testing equality
def _equal3(S,T):
  """Compare two polytopes to see if they are the same. 
  P,Q = convexhull of points S,T.
  Will throw an exception if either P or Q is lower dimensional. 
  In this case, recompute: sequentially drop eachc coordinate
  and compare P and Q again. (This results in m number of comparisons,
  where m is the number of columns of S)
  """
  try:
    P = ConvexHull(S)
    Q = ConvexHull(T)
  except:
    S = np.reshape(S,(len(S),len(S[0])))
    T = np.reshape(T,(len(T),len(T[0])))
    cols = S.shape[1]
    for c in range(cols):
      S2 = S[:,[i for i in range(cols) if i != c]]
      T2 = T[:,[i for i in range(cols) if i != c]]
      if not _equal3(S2,T2):
        return False
    return True
  vp= [S[i] for i in P.vertices]
  vq = [T[i] for i in Q.vertices]
  #sort in increasing x coordinate
  vp = np.sort(vp,axis = -1)
  vq = np.sort(vq, axis = -1)
  if(len(vp) != len(vq)):
    return False
  return all(np.ndarray.flatten(vp == vq))  

def _newton3(w):
  """Returns a tuple of (Paa,Saa) etc for the word w"""
  (Sa,path) = newton._S(w,1)
  (Sb,path) = newton._S(w,0)
  tup = []
  tup += [_Saa(Sa)]
  tup += [_Sab(Sa,Sb,a=1)]
  tup += [_Sab(Sb,Sa,a=0)]
  tup += [_Saa(Sb)]
  return tup
  
def equal3(w,w2):
  """Check if w ~ w2 in UT3"""
  if not newton.equal(w,w2):
    return False
  else:
    tup = _newton3(w)
    tup2 = _newton3(w2)  
  for i in range(len(tup)):
    S = tup[i] 
    T = tup2[i]
    if not _equal3(S,T):
      return False
  return True

def equal3fast(w,w2):
  """Check if w ~ w2 in UT3.
  Faster version: applicable to w ~ w2 in UT2. 
  """
  tup = _newton3(w)
  tup2 = _newton3(w2)  
  for i in range(len(tup)):
    S = tup[i] 
    T = tup2[i]
    if not _equal3(S,T):
      return False
  return True

def convexHull3(S):
  """Modified version of convex hull.
  Returns the vertices of S if S is full-dim. 
  Raise an exception if P is low dimensional.
  """
  vp = []
  try:
    P = ConvexHull(S)
    vp= [S[i] for i in P.vertices]        
  except:
    raise Exception('convexhull3 does not handle low dimensions')
  #sort in increasing x coordinate
  vp = np.sort(vp,axis = -1)    
  return tuple(map(tuple,vp))      
    
def getSignature(w):
  """Returns the signature of w = vertices of (Saa,Sab,Sba,Sbb).
  """
  tup = _newton3(w)
  try:
    vp = map(convexHull3, tup)
  except:
    raise
  return tuple(vp)
 
def compareSignature(vp,vp2):
  """"Compare two signatures."""
  if len(vp) != len(vp2):
    return False
  for i in xrange(len(vp)):
    if vp[i] != vp2[i]:
      return False
  return True

def hasFriend3(w):
  """Return True if w has a UT3 friend that is one swap away. If False, conjectured to be an isoterm."""
  wshift = w[1:]
  wright = w[:-1]
  for i in range(len(wshift)):
    if wshift[i] != wright[i]:
      w2 = np.array(w)
      w2[i] = w[i+1]
      w2[i+1] = w[i]
      if equal3(w,w2):
        return True
  return False

#example run:
#w = [1,0,0,1,0,1,1,0,0,1]
#w2 = [1,0,0,1,1,0,1,0,0,1]
#equal3(w,w2)
