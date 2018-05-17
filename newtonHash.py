#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 18:59:55 2018

@author: ngoc

Generate the dictionary of all UT2 equivalence classes for W(na,nb). 

The list of vertices is the keyword for each word. This allows easy checking if a word is an isoterm for UT2.

Program idea:
  * fix na,nb
  * iterate over all words with na 1's and nb 0's
  * for each word, compute (P,Q) (the vertices of each polygon, ordered in increasing natural coordinate)
  * store the sequence of vertices as two tuples (),(). This is the key for w. 
  * look up a dictionary:
    * if key is in the dictionary, add w to the value associated to this key (as a tuple)
    * if key is NOT in the dictionary, add (key,[w]) in
  * after this computation: pickle the dictionary. 
  
Version 13/04: faster generation using minmax  
  
"""
import newton as newton
import numpy as np
import itertools as itt
from util import *

def _addToDict(d,kw,w):
  if d.has_key(kw):
    d[kw] += [w]
  else:
    d[kw] = [w]  

def addToDict(d,w):
  keyword = newton.getSignature(w)
  _addToDict(d,keyword,w)

def generateDict(na=6,nb=6):
  """Generate the dictionary of all UT2 equivalence classes for W(na,nb).
Dictionary entries are of the form (key,list), where key is the signature of this equivalence class. 
 If na = nb, then only enumerate words that start with an `a' (1)."""
  identityDict = {}
  n = na+nb
  if na != nb:
    agen = itt.combinations(xrange(n),na)
  else:
    agen = itt.combinations(xrange(1,n),na-1)    
  while True:
    try:
      aloc = agen.next()
      w = [0]*n
      for i in aloc:
        w[i] = 1
      if na == nb:
        w[0] = 1
      addToDict(identityDict,w)
    except:
      return identityDict
 
def listIdentities2(na=6,nb=6):
  """Generate the dictionary using the min/max element algorithm
  (ListIdentities2)
  as outlined in the paper. 
  Note: does not do half for the symmetric case. 
  COMMENT: for small values, due to these functions being in Python, 
  not significantly faster."""
  L = {}
  wlist = [[1]*na+[0]*nb]
  while len(wlist) > 0:
    wlist2 = []
    for w in wlist:
      (wmin,wmax) = newton.getMinMax(w)      
      if not L.has_key(wmin):
        (Smin,p) = newton._S(wmin,a=1,details=False)
        (Smax,p) = newton._S(wmax,a=1,details=False)          
        L[wmin] = wmax
        #add in BIGGER neighbors of wmin only
        bigger = newton.swaps(wmin,forward=[1])
        for w2 in bigger:
          (S,p) = newton._S(w2,a=1,details=False)
          if sum(S[:,1] < Smin[:,1]) == 0:      
            wlist2 += [w2]
    wlist = wlist2
  return L            