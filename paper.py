#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:22:27 2018

@author: ngoc

All codes to reproduce experimental results in the paper

"""
import newtonHash as newtonHash
import newtonHash3 as newtonHash3
import newton as newton
import newton3 as newton3
import numpy as np
from util import *
import cPickle as cPickle

if False:
  """Proof of Theorem 1.3.
  Method: Computed UT2 dictionary of equivalence classes for W(na,n-na) for na = 4,5,6,...11; 
  and n = 21 and n = 22.
  Then we load these dictionaries one by one, and hunt for UT3 identities. 
  Store the non-isoterm identities in a dictionary. 
  For n = 21, all dictionaries are empty, meaning that there are no UT3 identities 
  For n = 22, the only non-empty dictionary is the case na = nb = 11, and 
  these are the identities displayed in the paper.
  """   
  #Generate the dictionary of all UT2 equivalence classes for W(na,nb), then save to a zipped pickle file.
  #file name format: dict(na)(nb).p
  if False:
    for na in range(4,12):
      print "do ListAll in UT2 for na = " + str(na) + ", n = 21"
      nb = 21-na
      #compute the dictionary for (na,21)
      identityDict = newtonHash.generateDict(na,nb)
      #zip and save
      filename = './output/dict'+str(na)+str(nb)+'.p'
      save_zipped_pickle(identityDict,filename) 
      #compute the dictionary for (na,22)
      nb = 22-na
      print "do ListAll in UT2 for na = " + str(na) + ", n = 22"
      identityDict = newtonHash.generateDict(na,nb)
      filename = './output/dict'+str(na)+str(nb)+'.p'
      save_zipped_pickle(identityDict,filename)     
  #reload the dictionaries and hunt for UT3 identities amongst these UT2 identities
  #--- list all UT3 identities for n = 21 and n = 22
  #metadict21[na] = a dictionary of all UT3 identities for W(na,21-na)
  #metadict22[na] = a dictionary of all UT3 identities for W(na,22-na)
  metadict21 = {}
  metadict22 = {}
  for na in xrange(11,3,-1):
    print "na = " + str(na) + ",n = 21"
    dictn = getDict(na,21,prefix='./output/dict')
    dickey = dictn.keys()
    metadict21[na] = newtonHash3.allUT3Hash(dictn,dickey,start=0)
    print "na = " + str(na) + ",n = 22"
    dictn = getDict(na,22,prefix='./output/dict')
    dickey = dictn.keys()
    metadict22[na] = newtonHash3.allUT3Hash(dictn,dickey,start=0)    

  #count number of UT3 identities for each setting of na and n
  [len(metadict21[na]) for na in metadict21.keys()]
  [len(metadict22[na]) for na in metadict22.keys()]
  #display all ut3 identities for na = 11, n = 22
  metadict22[11]
  #save for future reference
  filename = './output/metadict-ut3-22.pk'
  f = open(filename, 'wb')
  cPickle.dump(metadict22,f)
  filename = './output/metadict-ut3-21.pk'
  f = open(filename, 'wb')
  cPickle.dump(metadict21,f)  
  #print the identities for n = 22, na = 10
  d10 = metadict22[10].values()
  map(lambda x: map(stringToWord,x[0]),d10)
  #print the identities for n = 22, na = 11
  d11 = metadict22[11].values()
  map(lambda x: len(x[0]), d11)
  #each class has size 2. 
  map(lambda x: map(stringToWord,x[0]),d11)
  #merge class by flips and reversing the role of a's and b's. 
  #for the case na = 10 it is easy to do by inspection that the two classes are equal up to a flip.
  #for the other case, we recompute the keys up to reversals and flips. 
  d = metadict22[11]
  dflip = {}
  dkeys = d.keys()
  seen = set()
  for i in xrange(len(dkeys)):
    key = dkeys[i]
    if key not in seen:      
      seen.add(key)
      w1,w2 = d[key][0]
      f1 = np.flip(w1,0)
      f1sig = newton3.getSignature(f1)
      r1 = map(lambda x: 1 if x == 0 else 0, w1)
      r1sig = newton3.getSignature(r1)
      #need to flip the reverse as well
      f2 = np.flip(r1,0)
      f2sig = newton3.getSignature(f2)      
      r2 = map(lambda x: 1 if x == 0 else 0, f2)
      r2sig = newton3.getSignature(r2)
      dflip[key] = [w1,w2]
      if f1sig in d:
        seen.add(f1sig)
      if f2sig in d:
        seen.add(f2sig)
      if r1sig in d:
        seen.add(r1sig)
      if r2sig in d:
        seen.add(r2sig)
    #print out the words up to flip and reversals
    map(lambda x: map(stringToWord,x),dflip.values())
    len(dflip.values())

if False: 
  """ Generate data for Figure 10, Example 5.4"""
  def randomLong3(na=30, ntrials = 50000):
    """Draw ntrials many words from W(na,na), and count 
    the number of locally isolated UT3 words."""
    nb = na
    w = [1]*na + [0]*nb
    isocount = 0    
    printrate = ntrials/10
    for trialid in xrange(ntrials):
      if trialid % printrate == 0:
        print "done up to trial #" + str(trialid) + " out of " + str(ntrials) + ",na = " + str(na)
      np.random.shuffle(w)
      if not newton3.hasFriend3(w):
        isocount += 1    
    return (isocount,ntrials)  
  #collect statistics
  from multiprocessing import Pool
  pool = Pool()
  ss = pool.map(randomLong3, xrange(15,41))
  pool.close()
  pool.join()
  #reorganize asa dictionary
  statdict= {}
  for na in xrange(15,41):
    statdict[na] = ss[na-15]
  #-- save for future reference. Note that the experiments are random 
  #so recomputed data may not be identical to the figures used the paper
  filename = './output/random-ut3-rate-15to40.pk'
  f = open(filename, 'wb')
  cPickle.dump(statdict,f)
  f.close() 
  #load files
  filename = './output/random-ut3-rate-15to40.pk'
  stat3 = cPickle.load(open(filename,'r'))
 
if False: 
  """Generate data for Conjecture 5.6 (Figure 11). 
  Fraction of long ut2 adjacent identities which are also ut3 identities"""
  def ut3fraction(ell):
    repeat = 100
    count2 = []
    count3 = []
    w = [1]*ell + [0]*ell
    while repeat > 0:
      #progress counter
      if repeat % 10 == 0:
        print "iter = " + str(repeat) + "ell = " + str(ell)
      np.random.shuffle(w)
      wlist = newton.returnFriendFast(w)
      count2 += [len(wlist)]
      wlist3 = map(lambda v: newton3.equal3fast(w,v),wlist)
      count3 += [sum(wlist3)]
      repeat -= 1
    return (count2,count3)
  
  #collect statistics
  from multiprocessing import Pool
  pool = Pool()  
  metadict = {}
  ss = pool.map(ut3fraction, xrange(100,700,100))
  pool.close()
  pool.join()
  
  for i in xrange(1,7):
    metadict[ell*100] = ss[i-1]
  #save data
  filename = './output/longut3-partial.pk'
  f = open(filename, 'wb')
  cPickle.dump(metadict,f)  
  #reload, generate data for plotting
  filename = './output/longut3-partial.pk'
  f = open(filename, 'r')
  metadict = cPickle.load(f)    
  for i in metadict.keys():
    (count2,count3) = metadict[i]
    data = 1.0*np.array(count3)/np.array(count2)
    f = open("./output/longut3-"+str(i)+'.txt','wb')
    np.savetxt(f,data,fmt='%.10f')  
    f.close()  
  
if False:
  """Data for Figure 12 and 13, Conjecture 5.8"""
  def getStats(dictn):
    """Returns a dictionary with the following 4 statistics:
      nidentity = number of identities (= equiv - iso)
      niso = number of isoterms.
      nbin2 = number of bins with 2 words
      nbin2r = number of bins with 2 words, where w' is reverse of w.
      """
    stats = {}
    dicval = dictn.values()
    count = np.array(map(len,dicval))
    stats['nidentity'] = len(count)-sum(count==1)
    stats['niso'] = sum(count == 1)
    stats['nbin2'] = sum(count == 2)
    count2r = 0
    for i in xrange(len(count)):
      if len(dicval[i]) == 2:
        (w,w2) = dicval[i]
        if util._reverseTest(w,w2):
          count2r += 1
    stats['nbin2r'] = count2r
    return stats
  #analyze statistics for UT2 equivalence classes for n = 22
  n = 22
  metastats = {}
  for na in xrange(4,12):
    print "computing na = " + str(na)
    nb = n - na
    dictn = getDict(na,22,prefix='./output/dict')
    metastats[na] = getStats(dictn)
  #--- save this
  filename = 'stats-ut2-22.pk'
  f = open(filename, 'wb')
  cPickle.dump(metastats,f)
  
