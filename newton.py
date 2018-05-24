#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 14:29:43 2018

@author: ngoc

Input: two words as two binary strings (1 = a, 0 = b)
Output: True if w ~ w', False else. 

Subfunctions:
  * given w, compute (N_a, N_b)
  * given a pair of polytopes P, Q, check if they are equal
  * plot polygons (or pairs of polygons)
  
Next part: given a w, find all w' such that w ~ w'
"""

import numpy as np
#import itertools as itt
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from collections import deque  
import util as util

def _S(w,a=1,details=False):
  """Input: w = a binary string
  Output: list of points S, and if details, the path"""
  path = np.array([0,0])
  s = np.array([0,0]) #last point
  for i in range(len(w)):
    if w[i] == a:
      increase = np.array([1,0])
    else:
      increase = np.array([0,1])
    s += increase
    path = np.vstack((path,s))
  loc = [i for i in range(len(w)) if w[i] == a]
  if(a == 1):
    S = np.array([(i,loc[i]-i) for i in range(len(loc))])
  else:
    S = np.array([(loc[i]-i,i) for i in range(len(loc))])
  if details:
    return (S,path)
  else:
    return (S,None)

def newton(w,a=1,details=False):
  """ Input: w = a binary string. 
  Output: the single polygon P, list of points S, and the path (mostly for debug)"""
  #OK.
  (S,path) = _S(w,a,details)
  try: 
    P = ConvexHull(S)
  except: #degenerate case: S is a line. 
    #In this case, one of the coordinates is increasing, so it is sufficient to keep the (0.) and (last,.) points
    P = None
  if details:
    return (P,S,path)
  else:
    return (P,S)

#plot code

def _plotHull(P,points,na,nb,linecol='k-'):
  #plt.plot(points[:,0], points[:,1], 'o')
  for simplex in P.simplices:
      plt.plot(points[simplex, 0], points[simplex, 1], linecol)
  plt.plot(points[P.vertices,0], points[P.vertices,1],'o')
  plt.grid(b=True, which='major', axis='both')
  plt.xticks(np.arange(0, na+1, 1))
  plt.yticks(np.arange(0, nb+1, 1))  
  

def plotWord(w):
  na = sum(np.array(w))
  nb = len(w) - na
  (P,S) = newton(w,1)
  (Q,T) = newton(w,0)
  _plotHull(P,S,na,nb,linecol='b--')
  _plotHull(Q,T,na,nb,linecol='r--')

def plotPair(w,w2,title=""):
  """Plot an identity with the two paths"""
  na = sum(np.array(w))
  nb = len(w) - na  
  (P,S,path1) = newton(w,1,True)
  (P2,S2,path2) = newton(w2,1,True)
  (Q,T) = newton(w,0)
  fig, ax = plt.subplots()
  ax.set_color_cycle(['black','red','blue','green'])
  plt.plot(path1[:,0], path1[:,1])
  plt.plot(path2[:,0], path2[:,1])
  _plotHull(P,S,na,nb,linecol='b-.')
  _plotHull(Q,T,na,nb,linecol='g-.')
  plt.title(title)

  
def plotCustomPair(w,w2,a=1,title=""):
  """Plot a non-identity, with custom polygons. 
  if a = 1, plot the P pair. 
  if a = 0, plot the Q pair."""
  na = sum(np.array(w))
  nb = len(w) - na    
  (P,S,path1) = newton(w,a,True)
  na2 = sum(np.array(w2))
  nb2 = len(w2) - na2  
  (P2,S2,path2) = newton(w2,a,True)
  fig, ax = plt.subplots()
  ax.set_color_cycle(['black','red','blue','green'])
  if a == 1:    
    plt.plot(path1[:,0], path1[:,1])
    plt.plot(path2[:,0], path2[:,1])
  else:
    plt.plot(path1[:,1], path1[:,0])
    plt.plot(path2[:,1], path2[:,0])    
  _plotHull(P,S,na,nb,linecol='b-.')
  _plotHull(P2,S2,na2,nb2,linecol='g-.')
  plt.title(title)
  
 
def getVertices(P,S):
  """Returns the vertices of P as a tuple ()"""
  if P == None:
    vp = np.vstack((S[0],S[-1]))
  else:
    vp= S[P.vertices]
  #sort in increasing x coordinate
  vp = np.sort(vp,axis=0)
  return tuple(map(tuple,vp))

def _equal(P,S,Q,T):
  """Compare two polygons to see if they are the same.
  P,Q = convexhull of points S,T.
  """
  if P == None:
    vp = np.vstack((S[0],S[-1]))
  else:
    vp= S[P.vertices]
  
  if Q == None:
    vq = np.vstack((T[0],T[-1]))
  else:
    vq = T[Q.vertices]
  #sort in increasing x coordinate
  vp = np.sort(vp,axis=0)
  vq = np.sort(vq,axis=0)
  if(len(vp) != len(vq)):
    return False
  return all(np.ndarray.flatten(vp == vq))

def equal(w,w2):
  """Compare two words as an identity in UT2. 
  Return TRUE if w ~ w', FALSE else"""
  (P,S) = newton(w,1)
  (Q,T) = newton(w,0)
  (P2,S2) = newton(w2,1)
  (Q2,T2) = newton(w2,0)
  return _equal(P,S,P2,S2) and _equal(Q,T,Q2,T2)

def getSignature(w):
  """Returns the UT2 signature of a word. 
  This is a pair of tuple of vertices of P and Q. """
  (P,S) = newton(w,1)
  vp = getVertices(P,S)
  (Q,T) = newton(w,0)
  vq = getVertices(Q,T)
  return (vp,vq)

def hasFriend(w):
  """ Return True if w has a UT2 friend. False if it is an isoterm. 
  Sufficient to check for adjacent swaps."""
  wshift = w[1:]
  wright = w[:-1]
  (P,S) = newton(w,1)
  (Q,T) = newton(w,0)
  for i in range(len(wshift)):
    if wshift[i] != wright[i]:
      w2 = np.array(w)
      w2[i] = w[i+1]
      w2[i+1] = w[i]
      (P2,S2) = newton(w2,1)
      (Q2,T2) = newton(w2,0)
      if _equal(P,S,P2,S2) and _equal(Q,T,Q2,T2):
        return True
  return False

def returnFriend(w):
  """Return a list of all w' ~ w where w' differs from w by one swap"""
  friends = []
  wshift = w[1:]
  wright = w[:-1]
  (P,S) = newton(w,1)
  (Q,T) = newton(w,0)
  for i in range(len(wshift)):
    if wshift[i] != wright[i]:
      w2 = np.array(w)
      w2[i] = w[i+1]
      w2[i+1] = w[i]
      (P2,S2) = newton(w2,1)
      (Q2,T2) = newton(w2,0)
      if _equal(P,S,P2,S2) and _equal(Q,T,Q2,T2):
        friends += [w2]
  return friends      

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

def _idloopFast(w,wmin,wmax,knownList):
  """Same as _idloop but use the min/max pair of w and thus faster"""
  newList = []
  wshift = w[1:]
  wright = w[:-1]
  for i in range(len(wshift)):
    if wshift[i] != wright[i]:
      w2 = np.array(w)
      w2[i] = w[i+1]
      w2[i+1] = w[i]
      if tuple(w2) not in knownList:
        inclass = (isLess(wmin,w2) == 1) and (isLess(w2,wmax) == 1)
        if inclass:
          newList += [tuple(w2)]
  return newList  


def _getMinMaxPair(wList):
  """Input: wList = list of all elements in a UT2 equivalence class.
  Return: the min and max element of a UT2 equivalent class.
  The min element is given by the lower most path (most b's first)
  The max element is given by the upper-most path (most a's first)
  Obsolete by the function getMinMax - used for debugging purposes only. 
  """
  patha = np.matrix(map(lambda w: _S(w,a=1)[0][:,1], wList))
  maxpath = patha.max(axis=0)
  #find the word
  maxword = next(wList[i] for i in range(len(wList)) if np.all(patha[i] == maxpath))
  pathb = np.matrix(map(lambda w: _S(w,a=0)[0][:,0], wList))
  minpath = pathb.max(axis=0)
  minword = next(wList[i] for i in range(len(wList)) if np.all(pathb[i] == minpath))
  return (minword,maxword)      


def _rPath(v,w,upper = True):
  """Compute the r path in the square with opposite vertices
  v and w. Assume v < w.
  Upper = True: max path on this square. 
  Upper = False: min path on this square."""
  if upper:
    #add the vertical points
    up = [(v[0],v[1]+i) for i in range(w[1]-v[1]+1)]
    #add the horizontal points
    right = [(v[0]+i,w[1]) for i in range(1,w[0]-v[0]+1)]
    return up + right
  else:
    up = [(w[0],v[1]+i) for i in range(1,w[1]-v[1]+1)]
    right = [(v[0]+i,v[1]) for i in range(w[0]-v[0]+1)]
    return right+up

def _maxsegment(p,q,labp,labq,P,Q):
  """The Maxsegment procedure in the paper.
  Input: points p, q = p_{i+1}, labels of p,q, and signature polytopes A,B
  Output: point p',label(p').
  label = 1 --> a; = 0 --> b
  """
  if p[0] == q[0] or p[1] == q[1]:
    return (q,labq)
  if labp == 0:
    for k in range(q[1],p[1],-1):
      if util._vinP((p[0],k),P) and util._vinP((p[0],k-1),Q):
        newp = (p[0],k)
        return (newp,1-labp)  
  else:
    for k in range(p[0]+1,q[0]+1):
      if util._vinP((k,p[1]),Q) and util._vinP((k-1,p[1]),P):
        newp = (k,p[1])
        return (newp,1-labp)
  #catch: iterating through an emptyset. 
  #this means the label on p needs to be flipped
  return (p, 1-labp)
  
def getMax(w):
  """Input: a word w.  
  Returns the max element overline{w}. Version May 1st, new algorithm.
  """   
  (P,S) = newton(w,1)
  A = getVertices(P,S)
  (Q,T) = newton(w,0)
  B = getVertices(Q,T)  
  na = sum(w)
  nb = len(w)-na    
  #Compute the control points
  points = np.sort(np.array(A+B), axis = 0)
  points= np.append(points, [[na,nb]],axis=0)
  #compute labels
  labels = map(lambda x: int(tuple(x) in A), points)
  labels[-1] = None
  #build the word directly (instead of building the path as written in the algorithm)
  bigw = []
  for i in xrange(len(points)-1):
    p = points[i]
    labp = labels[i]
    q = points[i+1]
    labq = labels[i+1]
    while tuple(p) != tuple(q):
      (newp, labnewp) = _maxsegment(p,q,labp,labq,P,Q)
      db = newp[1-labp] - p[1-labp]
      bigw += [labp]*db
      p = newp
      labp = labnewp
  return tuple(bigw)

def getMin(w):
  """Returns the minimum element. Version May 1st. """
  dualw = 1-np.array(w)
  return tuple(1-np.array(getMax(dualw)))
  
def getMinMax(w):
  return (getMin(w), getMax(w))

def _getMinMax(w,plot=None):
  """Input: a word. 
  Returns the min and max element in the orbit of w.  Old version: has errors. 
  Currently only used to plot boxes.
  """  
  (P,S) = newton(w,1)
  A = getVertices(P,S)
  (Q,T) = newton(w,0)
  B = getVertices(Q,T)  
  na = sum(w)
  nb = len(w)-na    
  #Compute the control points
  points = np.sort(np.array(A+B), axis = 0)
  points= np.append(points, [[na,nb]],axis=0)
  qu = deque(points)
  p = qu.popleft()
  ur = [tuple(p)]
  lr = [tuple(p)]  
  while(len(qu) > 0):
    if tuple(p) in A:
      v = (p[0]+1,p[1])
    else:
      v = (p[0],p[1]+1)
    u = qu.popleft()
    ur += _rPath(v,u,upper=True)
    lr += _rPath(v,u,upper=False)
    p = u
  ur = np.array(ur)
  lr = np.array(lr)
  #recode as lower and upper bounds for the region
  ls = util.pathToSa(lr,na)
  us = util.pathToSa(ur,na)
  ut = util.pathToSb(lr,nb)
  lt = util.pathToSb(ur,nb)
  #compute upper path
  upath = []
  for i in range(na):
    if np.all(ls[i] == us[i]):
      upath += [tuple(ls[i])]
    else:
      for j in range(us[i,1],ls[i,1]-1,-1):
        v = (i,j)
        if util._vinP(v,P) and (util._vinP(v,Q) or util._vinP((i,j-1),Q)):
          upath += [v]
          break       
  bigw = util.SatoW(np.array(upath),nb)
  #compute the lower path
  lpath = []
  for i in range(nb):
    if np.all(lt[i] == ut[i]):
      lpath += [tuple(lt[i])]
    else:
      for j in range(ut[i,0],lt[i,0]-1,-1):
        v = (j,i)
        if util._vinP(v,Q) and (util._vinP(v,P) or util._vinP((j-1,i),P)):
          lpath += [v]
          break
  smallw = util.SbtoW(np.array(lpath),na)
  #--- various plot types
  if plot is not None:
    _plotHull(P,S,na,nb,linecol='b--')
    _plotHull(Q,T,na,nb,linecol='r--')    
    plt.plot(ur[:,0],ur[:,1],'xkcd:dark lavender')
    plt.plot(lr[:,0],lr[:,1],'xkcd:dark lavender')  
    if plot == "low":
      (smallS,smallPath) = _S(smallw,a=1,details=True)
      plt.plot(smallPath[:,0],smallPath[:,1],'black')
    if plot == "high":
      (bigS,bigPath) = _S(bigw,a=1,details=True)
      plt.plot(bigPath[:,0],bigPath[:,1],'black')      
  return (tuple(smallw),tuple(bigw))

def isLess(w1,w2):
  """Compare two words. Returns: 1 if g1 <= g2, 0 if g2 < g1, None if they are not comparable."""
  if len(w1) == len(w2) and sum(np.array(w1)) == sum(np.array(w2)):
    (S1,p) = _S(w1,a=1,details=False)
    (S2,p) = _S(w2,a=1,details=False)
    count = sum(S1[:,1] <= S2[:,1])
    if count == S1.shape[0]:
      return 1
    if count == 0:
      return 0
    return None

def swaps(w,forward=[1]):
  """Returns a list of neighbors w' such that w' <-> w. 
  forward = list of stuff [1,0], [0] or [1]
  If forward = 1, require w' > w. ie: has more 0's to the front. 
  If forward = 0, require w' < w. """
  wlist = []
  for f in forward:
    (Tw,p) = _S(w,a=1-f,details=False)
    #find the first occurence
    (acount,bcount) = np.unique(Tw[:,1-f],True)
    for i in xrange(len(acount)):
      height = acount[i]
      loc = bcount[i]
      if height > 0:
        w2 = np.array(w)
        w2[loc+height-1] = w[loc+height]
        w2[loc+height] = w[loc+height-1]
        wlist += [w2]        
  return wlist    

def returnFriendFast(w):
  """Faster version of returnFriend, using minmax computations
  AND just use forward and backward swaps (not iterate through the entire word length)"""
  friends = []
  (wmin,wmax) = getMinMax(w)
  (Smin,p) = _S(wmin,a=1,details=False)
  (Smax,p) = _S(wmax,a=1,details=False)  
  wlist = swaps(w,forward=[0,1])
  for w2 in wlist:
      (S,p) = _S(w2,a=1,details=False)
      if sum(S[:,1] > Smax[:,1]) == 0 and sum(S[:,1] < Smin[:,1]) == 0:      
        friends += [w2]
  return friends  

if __name__ == "__main__":
  w = [1]*100 + [0]*100
  np.random.shuffle(w)
  _getMinMax(w) ==  getMinMax(w)
  #correct, though no significant computational improvement.