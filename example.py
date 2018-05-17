#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:22:27 2018

@author: ngoc

Examples and documentation guide for the UT2/UT3 package.

"""
import newtonHash as newtonHash
import newton as newton
import newton3 as newton3
import numpy as np
from util import *
import globalid as gid
import cPickle as cPickle
import matplotlib.pyplot as plt

if False:
  """Example: plot the polygons of a single word"""
  w = map(int,list('100010101010101110')) 
  #plot this word against itself
  newton.plotPair(w,w,title =''.join(map(str,w)))

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
    
