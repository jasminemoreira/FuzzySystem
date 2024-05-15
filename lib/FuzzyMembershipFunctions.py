# -*- coding: utf-8 -*-
"""
Fuzzy membership functions module.

This module is intended to provide the basic fuzzy membership functions 

Created on Sun Oct 15 11:55:02 2023

@author: Jasmine Moreira
"""
import math

class Triangular:
    def __init__(self, a,m,b,h=1):
        self.a = a
        self.m = m
        self.b = b
        self.h = h
        self.type = 'triangular'

    def mu(self, x):
        if x <= self.a: 
            return 0
        if x > self.a and x <= self.m:
            return self.h*(x-self.a)/(self.m-self.a)
        if x > self.m and x <= self.b:
            return self.h*(self.b-x)/(self.b-self.m)
        if x > self.b:
            return 0
    
class Trapezoidal:
    def __init__(self, a,m,n,b,h=1):
        self.a = a
        self.m = m
        self.n = n
        self.b = b
        self.h = h
        self.type = 'trapezoidal'
        
    def mu(self, x):
        if x <= self.a: 
            return 0
        if x > self.a and x <= self.m:
            return self.h*(x-self.a)/(self.m-self.a)
        if x > self.m and x <= self.n:
            return self.h
        if x > self.n and x <= self.b:
            return self.h*(self.b-x)/(self.b-self.n)
        if x > self.b:
            return 0

        
class Gaussian:
    def __init__(self, m,sigma,h=1):
        self.m = m
        self.sigma = sigma
        self.h = h
        self.type = 'gaussian'
        
    def mu(self, x):
        return self.h*math.exp(-self.sigma*math.pow(x-self.m,2))
    
        
class Singleton:
    def __init__(self, m,h=1):
        self.m = m
        self.h = h
        self.type = 'singleton'
        
    def mu(self,x):
        return self.h if x == self.m else 0