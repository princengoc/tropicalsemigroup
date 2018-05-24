#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:53:15 2018

@author: ngoc

Generate the UT2 equivalence class for a given word w. This uses a global variable. 
"""
import newton as newton

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
  newList = newton._idloop(w,globeList)
  addToGlobeList(newList)
  for w2 in newList:
    idloop(w2)

def idloopFast(w,wmin,wmax):
  """Same as idloop, but uses the min/max elements and
  do a search using word comparison instead of signature comparison"""
  newList = newton._idloopFast(w,wmin,wmax,globeList)
  addToGlobeList(newList)
  for w2 in newList:
    idloopFast(w2,wmin,wmax)  

def getGlobeList():
  return globeList