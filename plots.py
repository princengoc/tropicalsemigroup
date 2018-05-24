#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 16:17:12 2018

@author: ngoc

Scripts to produce various plots, given input file
consisting of pairs of words w ~ w'.

Need: newton.pyc
"""

import numpy as np
import matplotlib.pyplot as plt
import newton as newton
import newton3 as newton3
import newtonHash as newtonHash
from util import *

#--------- Main plot codes for the paper. 
 
if False:
  """Plot: proof of Ardjan's identity. 
  Code for Figure 2, Example 3.4
  """
  w = map(int,list('1001011001'))
  (w1,w2) = newton.getMinMax(w)
  newton.plotPair(w1,w2)
  plt.savefig('./figure/na5.pdf')
  w = map(int,list('1001010110'))
  (w1,w2) = newton.getMinMax(w)
  newton.plotPair(w1,w2)  
  plt.savefig('./figure/na5-2.pdf')
  w = map(int,list('0110100110'))
  (w1,w2) = newton.getMinMax(w)
  newton.plotPair(w1,w2)    
  #smallest identity with na = 4
  w = [1]+[0]*4+[1]*2+[0]*2+[1,0]
  (w1,w2) = newton.getMinMax(w) 
  newton.plotPair(w1,w2)
  plt.savefig('./figure/na4.pdf')
  
if False:
  """Figure 7, Example 3.20. """
  random = False
  if random:
    na = 14
    nb = 10
    w = [1]*na + [0]*nb
    np.random.shuffle(w)  
  else: #actual word used in the paper:
    w = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1]
  (w1,w2) = newton.getMinMax(w)  
  wtext = ''.join(map(str,w1))
  w2text = ''.join(map(str,w2))
  #do plots of the gray region
  (P,S) = newton.newton(w,1)
  A = newton.getVertices(P,S)
  (Q,T) = newton.newton(w,0)
  B = newton.getVertices(Q,T)  
  na = sum(w)
  nb = len(w)-na    
  #Compute the control points
  points = np.sort(np.array(A+B), axis = 0)
  points= np.append(points, [[na,nb]],axis=0)  
  newton._plotHull(P,S,na,nb,linecol='b--')
  newton._plotHull(Q,T,na,nb,linecol='r--')
  for i in xrange(len(points)-1):
    p = points[i]
    q = points[i+1]
    plt.plot((p[0],q[0]),(p[1],p[1]),'xkcd:dark lavender')
    plt.plot((p[0],p[0]),(p[1],q[1]),'xkcd:dark lavender')
    plt.plot((p[0],q[0]),(q[1],q[1]),'xkcd:dark lavender')
    plt.plot((q[0],q[0]),(p[1],q[1]),'xkcd:dark lavender')  
  plt.savefig('./figure/'+wtext+'-R.pdf') 
  #lower word
  (S1,path1) = newton._S(w1,a=1,details=True)
  newton.plotWord(w1)
  for i in xrange(len(points)-1):
    p = points[i]
    q = points[i+1]
    plt.plot((p[0],q[0]),(p[1],p[1]),'xkcd:dark lavender')
    plt.plot((p[0],p[0]),(p[1],q[1]),'xkcd:dark lavender')
    plt.plot((p[0],q[0]),(q[1],q[1]),'xkcd:dark lavender')
    plt.plot((q[0],q[0]),(p[1],q[1]),'xkcd:dark lavender')   
  plt.plot(path1[:,0],path1[:,1],'black')      
  plt.savefig('./figure/'+wtext+'-low.pdf')
  #upper word
  (S2,path2) = newton._S(w2,a=1,details=True)
  newton.plotWord(w2)
  for i in xrange(len(points)-1):
    p = points[i]
    q = points[i+1]
    plt.plot((p[0],q[0]),(p[1],p[1]),'xkcd:dark lavender')
    plt.plot((p[0],p[0]),(p[1],q[1]),'xkcd:dark lavender')
    plt.plot((p[0],q[0]),(q[1],q[1]),'xkcd:dark lavender')
    plt.plot((q[0],q[0]),(p[1],q[1]),'xkcd:dark lavender')   
  plt.plot(path2[:,0],path2[:,1],'black')      
  plt.savefig('./figure/'+wtext+'-high.pdf')
  newton.plotPair(w1,w2)  
  plt.savefig('./figure/'+wtext+'-both.pdf')
  print ''.join(map(str,['a' if i == 1 else 'b' for i in w]))
  print ''.join(map(str,['a' if i == 1 else 'b' for i in w1]))
  print ''.join(map(str,['a' if i == 1 else 'b' for i in w2]))

if False:
  """Figure 8, Example 3.21. Plot one example from the the Dyck path family"""
  k = 4
  r = 4
  prefix = [1] + [0]*r +[1]
  suffix = [1-i for i in prefix]
  w = prefix + [0,1]*k + suffix
  w2 = prefix + [1,1,1,0,1,0,0,0] + suffix
  newton.equal(w,w2)
  wtext = ''.join(map(str,w))
  w2text = ''.join(map(str,w2))
  newton.plotPair(w,w2,title ='')
  plt.savefig('figures/fibonacci.pdf')  
  
  
if False:
  """Figure 9, Example 5.4. 
  Plot a particular UT3 equivalence class of size 3, where the swap order matters.
  """
  w0 = map(lambda x: 0 if x == 'b' else 1, list('baaaabaaaaaababbbbbabbbaabbabbabbabbaabbabaa'))
  w2 = map(lambda x: 0 if x == 'b' else 1, list('baaaabaaaaaababbbbbabbbababbababbabbaabbabaa'))
  w3 = map(lambda x: 0 if x == 'b' else 1, list('baaaabaaaaaababbbbbabbbabababbabbabbaabbabaa'))  
  w4 = map(lambda x: 0 if x == 'b' else 1, list('baaaabaaaaaababbbbbabbbaabbbababbabbaabbabaa'))
  #here: w0 ~ w2, min and max elements. but not ~ w4.
  wtext = ''.join(map(str,w0))
  newton.plotPair(w0,w2)
  plt.savefig('./figures/fib-ut3.pdf')
  #-- double-check
  map(lambda w: newton3.equal3(w,w0), [w0,w2,w3,w4])

  




