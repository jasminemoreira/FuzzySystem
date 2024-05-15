# -*- coding: utf-8 -*-
"""
Fuzzy basic sets module.

This module is intended to provide basic fuzzy sets.

Created on Thu Oct 26 23:39:35 2023

@author: Jasmine Moreira
"""
import numpy as np
from sys import float_info as fi
from FuzzyBase import _base
from FuzzyMembershipFunctions import Gaussian, Singleton,Trapezoidal, Triangular

class FuzzySet(_base):
    """A flexible class to create fuzzy sets.

    Attributes
    ----------
    function : MembershipFunction or function
        A membership function or a user defined function to support mu calculation.
    term : str
        The term used to identify the set.
        
    Methods
    -------
    mu(x)
        Use to calculate the membership degree of x (provided by the user function).
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, function, term='', limiter=1):
        """Use to construct the object.
        
        Parameters
        ----------
        function : function
            A user defined membership function to support mu calculation.
        term : str
            The term used to identify the set.
        limiter:
            The max mu value the set can return
        """
        self.limiter = limiter
        self.term = term
        self.mf = ''
        try:
            self.mf = function
            self.muf = function.mu # MembershipFunctions
        except AttributeError:
            self.muf = function    # Single functions
        
    def mu(self, x, y=0):
        """Use to calculate the membreship degree of x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
        y: float
            An optional number within the Universe to calculate a mu value when 2 arguments are needed.
            
        Returns
        -------
            mu(x) if mu(x) > alpha and 0 otherwise, considering U. Where mu(x) is the membership degree of x to the set.
        """
        mu = self.muf(x)
        return mu if mu <= self.limiter else self.limiter


class TriangularSet(FuzzySet):
    """A triangular fuzzy set class.

    Attributes
    ----------
    a : float
        The set starting point.
    m : float
        The set middle point (triangle vertex).
    b : float
        The set final point.
    h : float
        The set height.
    term : str
        The term used to identify the set.
        
    Methods
    -------
    mu(x)
        Use to calculate the membership degree of x to the set.
        
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, a, m, b, h=1, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        a : float
            The set starting point.
        m : float
            The set middle point (triangle vertex).
        b : float
            The set final point.
        h : float
            The set height.
        term : str
            The term used to identify the set.
        """          
        a = a-fi.epsilon if a==m else a
        b = b+fi.epsilon if b==m else b
        self.a = a
        self.b = b
        self.m = m
        self.h = h
        super().__init__(Triangular(a, m, b, h), term)
        
    def updateParam(self, param, value):
        param = setattr(self, param, value)
        self.mf = Triangular(self.a, self.m, self.b, self.h)
        self.muf = self.mf.mu

class TrapezoidalSet(FuzzySet):
    """A trapezoidal fuzzy set class.

    Attributes
    ----------
    a : float
        The set start point.
    m : float
        The upper base start point.
    n : float
        The upper base end point.        
    b : float
        The set final point.
    h : float
        The set height.
    term : str
        The term used to identify the set.
        
    Methods
    -------
    mu(x)
        Use to calculate the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, a, m, n, b, h=1, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        a : float
            The set start point.
        m : float
            The upper base start point.
        n : float
            The upper base end point.        
        b : float
            The set final point.
        h : float
            The set height.
        term : str
            The term used to identify the set.
        """           
        a = a-fi.epsilon if a==m else a
        b = b+fi.epsilon if b==n else b
        self.a = a
        self.b = b
        self.m = m
        self.n = n
        self.h = h
        super().__init__(Trapezoidal(a, m, n, b, h), term)
        
    def updateParam(self, param, value):
        param = setattr(self, param, value)
        self.mf = Trapezoidal(self.a, self.m, self.n, self.b, self.h)
        self.muf = self.mf.mu
        
class GaussianSet(FuzzySet):
    """A gaussian fuzzy set class.

    Attributes
    ----------
    m : float
        The gaussian mean.
    sigma : float
        The gaussian variance.        
    h : float
        The set height.
    term : str
        The term used to identify the set.
        
    Methods
    -------
    mu(x)
        Use to calculate the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, m, sigma, h=1, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        m : float
            The gaussian mean.
        sigma : float
            The gaussian variance.        
        h : float
            The set height.
        term : str
            The term used to identify the set.
        """
        self.m = m
        self.sigma = sigma
        self.h = h
        super().__init__(Gaussian(m, sigma, h), term)
        
    def updateParam(self, param, value):
        param = setattr(self, param, value)
        self.mf = Triangular(self.m, self.sigma, self.h)
        self.muf = self.mf.mu
        
class SingletonSet(FuzzySet):
    """A singleton fuzzy set class.

    Attributes
    ----------
    m : float
        The singleton definition point.    
    h : float
        The singleton set height.
    term : str
        The term used to identify the set.
        
    Methods
    -------
    mu(x)
        Use to calculate the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, m, h=1, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        m : float
            The singleton definition point.    
        h : float
            The singleton set height.
        term : str
            The term used to identify the set.
        """
        self.m = m
        self.h = h
        super().__init__(Singleton(m, h), term)
        
    def updateParam(self, param, value):
        param = setattr(self, param, value)
        self.mf = Triangular(self.m, self.h)
        self.muf = self.mf.mu
        
        
class DiscreteSet(FuzzySet):
    """A class to create a discrete fuzzy set.

    Attributes
    ----------
    a : [float]
        A discrete list of points.    
    b : [float]
        A discrete list of mu corresponding values.
    term : str
        The term used to identify the set.
    """
    
    def __init__(self, a, b, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        a : [float]
            A discrete list of points.    
        b : [float]
            A discrete list of mu corresponding values.
        term : str
            The term used to identify the set.
        """
        if type(a) is Universe:
            self.a = a.getU()
        else:
            self.a = a
        self.b = b
        sets = []
        for x,y in zip(self.a, self.b):
            sets.append(SingletonSet(x, y))
        
        def muf(x):
            mus = []
            for s in sets:
                mus.append(s.mu(x))
            return max(mus)
            
        super().__init__(muf, term)
        

class Universe(FuzzySet):
    """A class to define a fuzzy universe.

    Attributes
    ----------
    a : float
        The universe starting point.    
    b : float
        The universe final point.
    npoints : str
        The number of points in the universe.
        
    Methods
    -------
    mapSet(A)
        Use to calculate the membership degree of x to the set.
    """
        
    def __init__(self, a, b='', npoints=2000, discrete=False, term=''):
        self.a = a
        self.b = b
        self.npoints = npoints
        self.extra = []
        self.discrete = discrete
        if type(self.a) is list or type(self.a) is np.ndarray:
            super().__init__(lambda x: 1 if x in a else 0, term)
        else:
            super().__init__(lambda x: 1 if x >= a and x <=b else 0, term)
             
    def _appendPoint(self, A, attrib):
        u = self.getU()
        if hasattr(A.mf,attrib):
            mfp = getattr(A.mf, attrib)
            if mfp not in u and mfp > min(u) and mfp < max(u):
                self.addPoints([mfp])
                
    def mapSet(self,A):
        """Use to map the main points of a FuzzyBasicSet object into the universe instance.
        
        Parameters
        ----------
        A : FuzzyBasicSet
            The set to be mapped into the universe instance.   
        """
        if not hasattr(A, 'mf'):
            return
        
        self._appendPoint(A, 'a')
        self._appendPoint(A, 'm')
        self._appendPoint(A, 'n')
        self._appendPoint(A, 'b')     
             
    def addPoints(self,a):
        """Use to map a list of points into the universe instance.
        
        Parameters
        ----------
        a : [float]
            A list of points to be mapped into the universe instance.   
        """ 
        u = self.getU()
        for p in a:
            if p not in u:
                self.extra = np.append(self.extra,[p])
       
    def getU(self):
        """Use to get the points of the universe instance.
        
        Returns
        -------
            A numpy array containing the points of the universe instance.   
        """ 
        if type(self.a) is list:
            u = np.array(self.a, dtype=float)
        elif type(self.a) is np.ndarray:
            u = self.a
        else:
            u = np.linspace(self.a, self.b, self.npoints, dtype=float)
        for p in self.extra:
            if p not in u:
                u = np.append(u,[p])
        u.sort()
        return u

    def getDenseU(self):
        """Use to get finer step points for the universe instance.
        
        This method aims to provide more points in order to increase the precision on calculations.  
        
        Returns
        -------
            A numpy array containing 100.000 points of the universe instance.   
        """
        if self.discrete:
            return self.getU()
        
        if self.b == '' or self.npoints > 100000:
            return self.getU()
        npoints = self.npoints
        self.npoints = 100000
        U = self.getU()
        self.npoints = npoints
        
        #TODO: retornar os pontos mapeados
        
        return U    