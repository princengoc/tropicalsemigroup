#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:22:27 2018

@author: ngoc

Examples and documentation guide for the UT2/UT3 package.

"""
import newtonHash as newtonHash
import newtonHash3 as newtonHash3
import newton as newton
import newton3 as newton3
import numpy as np
from util import *
import globalid as gid
import cPickle as cPickle
import matplotlib.pyplot as plt
import itertools as itt



if False:
  """Example: plot the polygons of a single word"""
  w = map(int,list('100010101010101110')) 
  #plot this word against itself
  newton.plotPair(w,w,title =''.join(map(str,w)))
  #--- another word, just with the path and either A or B polygons but not both
  w = map(lambda x: 1 if x == 'a' else 0, list('abbababaab'))
  newton.plotPath(w)
  plt.savefig('./talks/figures/adjan1.pdf')
  newton.plotPath(w,a=1)
  plt.savefig('./talks/figures/adjan2.pdf')
  newton.plotPath(w,a=1)
  newton.plotWord(w,a=1)
  plt.savefig('./talks/figures/adjan3.pdf')
  newton.plotPath(w,a=0)
  plt.savefig('./talks/figures/adjan4.pdf')
  newton.plotPath(w,a=0)
  newton.plotWord(w,a=0)
  plt.savefig('./talks/figures/adjan5.pdf')
  w = map(lambda x: 1 if x == 'a' else 0, list('abbaabbaab'))
  newton.plotPath(w,a=1)
  newton.plotWord(w,a=1)
  plt.savefig('./talks/figures/adjan3-v.pdf')  
  newton.plotPath(w,a=0)
  newton.plotWord(w,a=0)
  plt.savefig('./talks/figures/adjan5-v.pdf')    
  
if False:
  """Example: read in identities from a file and plot them"""
  filename = './input/input3.txt'  
  wList, w2List = getPairs(filename)
  for i in range(len(wList)):
    w = wList[i]
    w2 = w2List[i]
    newton.plotPair(w,w2)
    #check that these are indeed UT3 identities
    print newton3.equal3(w,w2)

if False:
  """Example: convert a pair of words from 'ab' to '10'
  """
  w = map(lambda x: 1 if x == 'a' else 0, list('ababbabaababbabaababba'))
  w2 = map(lambda x: 1 if x == 'a' else 0, list('ababbabaabbababaababba'))
  newton.plotPair(w,w2)
  newton.equal(w,w2)
  print newton3.equal3(w,w2)
  newton.equal(np.flip(w,0),np.flip(w2,0))
  print newton3.equal3(np.flip(w,0),np.flip(w2,0))
  
 
if False:
  """Example: read in non-identities from a file and plot the words individually. 
  input file: inputnon.txt"""
  filename = './input/inputnon.txt'
  wList,w2List = getPairs(filename)  
  for i in range(len(wList)):
    w = wList[i]
    w2 = w2List[i]
    wtext = ''.join(map(str,w))
    w2text = ''.join(map(str,w2))
    #compare the Q polygon of these two words
    newton.plotCustomPair(w,w2,a=0,title = wtext + ',' +w2text)
    #compare the P polygon of these two words
    newton.plotCustomPair(w,w2,a=1,title = wtext + ',' +w2text)
    #if save:
    #plt.savefig('./figures/' + wtext + '.pdf')  

if False:
  """Example: compute the min and max element of a UT2 equivalence class
  in two different ways.
  The first uses the Structural Theorem and Theorem 3.18 to compute min/max.
  The second enumerates the entire equivalence class by a recursion (breadth-first search tree)
  and take the min and max element"""  
  k = 4
  prefix = [1,0,0,1]
  suffix = [1-x for x in prefix]
  between = [1,0]*k
  w = prefix + between + suffix        
  #first method: use Theorem 3.18
  (w1,w2) = newton.getMinMax(w)
  newton.plotPair(w1,w2)
  #second method:
  #compute the equivalence class with a recursion
  gid.clear()
  gid.addToGlobeList(set([tuple(w)]))  
  gid.idloop(w)
  #list out all words equivalent to w
  globeList = list(gid.getGlobeList())
  print len(globeList)
  #plot min/max words
  (minword,maxword) = newton._getMinMaxPair(globeList)
  wtext = ''.join(map(str,minword))
  w2text = ''.join(map(str,maxword))
  newton.plotPair(minword,maxword,title = wtext + ',' +w2text)  
  #check that the two methods give the same answer
  print w1 == minword and w2 == maxword  
 

if False:
  """Example: find a random non-isoterm w in UT2, 
  and check the fraction of its neighbors which also
  form an identity with w in UT3
  """
  na = 100
  nb = 100
  w = [1]*na + [0]*nb
  np.random.shuffle(w)
  while not newton.hasFriend(w):    
    np.random.shuffle(w)
  #compute adjacent friends and check the fraction which are UT3 identities
  friendList = newton.returnFriendFast(w)
  print "number of UT2 neighbors = " + str(len(friendList)) + "\n"
  print "number of UT3 neighbors = " + str(sum([newton3.equal3fast(w,w2) for w2 in friendList]))

if False:
  """Double-check Proposition 2.6 for the case of degree-2 signature. """
  #w = map(int,list('1000101010101011100011111100')) 
  w = map(lambda x: 1 if x == 'a' else 0, list('abbbabababababaaab'))
  (Sa,path) = newton._S(w,a=1)
  (Sb,path) = newton._S(w,a=0)
  #Sab = _Sab(Sa,Sb,a=1)
  (Sab,Sba) = newton3._getCrossSignature(w)
  keepab = newton3._Sab(Sa,Sb,a=1)
  keepba = newton3._Sab(Sb,Sa,a=0)
  #do a transformation
  pi = [1,0,0,0,0,1,0,0,1,0,1,0,0,1,0,1]
  pi = np.reshape(pi,(4,4))
  Sab = np.matrix(Sab)
  Sba = np.matrix(Sba)
  #comment: getCrossSignature counts the intermediate letter, so does not need to add (0,0,1,0)
  pix = np.dot(pi,Sab.T).T + np.matrix((0,0,1,0))
  piy = np.dot(pi,Sba.T).T + np.matrix((0,0,1,0))
  print np.all(pix == keepab)
  print np.all(piy == keepba)  

  

if False:
  """Find an instance to produce Figure 9, Example 5.4"""
  #get a random word. Compute up to recursion depth 2 only, to get a smallest counter example possible
  na = 22
  nb = 22
  w = [1]*na + [0]*nb
  found = False
  while not found:
    np.random.shuffle(w)
    #compute the equivalence class with a recursion depth 2 only
    wmin,wmax = newton.getMinMax(w)
    globeList = set([tuple(w)])
    neighbors = newton._idloopFast(w,wmin,wmax,globeList)
    globeList = globeList.union(set(neighbors))
    print len(globeList)
    #go through globeList: compute ut3 signatures
    d2 = {}
    d2[newton.getSignature(w)] = globeList
    d3 = newtonHash3.allUT3Hash(d2,d2.keys())
    biglist = [x[0] for x in d3.values() if len(x[0]) > 2]
    for wlist in biglist:
      wlist = biglist[0]
      degree = len(wlist)      
      if degree == 3:
        wlist = diamond(wlist)
        if not np.all(map(lambda w: newton3.equal3(w,wlist[0]), wlist)):
          print "found counter example"
          found = True
          break

  map(lambda w: newton3.equal3(w,wlist[0]), wlist)
  wmin,wmax = newton.getMinMax(wlist[0])
  newton.plotPair(wmin,wmax)
  plt.savefig('./figures/fib-minmax.pdf')
