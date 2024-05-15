# -*- coding: utf-8 -*-
"""
Fuzzy operations module.

This module is intended to provide basic fuzzy operations.

Created on Sun Oct 22 16:15:58 2023

@author: Jasmine Moreira
"""
import pandas as pd
from sys import float_info as fi
from FuzzyBase import _base
from FuzzyBasicSets import SingletonSet

epsilon = fi.epsilon


#################################################################
#
# ONE ARGUMENT OPERATIONS
#
#################################################################


class AlphaCut(_base):
    """A class used to calculate an alpha-cut for a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    U : Universe
        A Universe instance defining the interval or the points to calculate the alpha-cut.
    alpha : float
        A number representing the desired alpha for the cut, usually between 0 and 1.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Use to check if x belongs to the alpha-cut set, returns the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, U, alpha, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        U : Universe
            A Universe instance defining the interval or the points to calculate the alpha-cut.
        alpha : float
            A number representing the desired alpha for the cut, usually between 0 and 1.
        term : str
            The term used to identify the resulting set.
        """    
        self.A = A
        self.term = term 
        self.alpha = alpha
        self.ucut = []
        U = U.getDenseU()
        for u in U:
            if self.mu(u) >= alpha:
                self.ucut.append(u)
        self.inf = min(self.ucut)
        self.sup = max(self.ucut)        
        
    def mu(self, x):
        """Use to check if x belongs to the alpha-cut set.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            mu(x) if mu(x) > alpha and 0 otherwise, considering U. Where mu(x) is the membership degree of x to the set.
        """
        mu = self.A.mu(x)
        return mu if mu >= self.alpha else 0


class Normalization(_base):
    """A class used to normalize a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    U : Universe
        A Universe instance defining the interval or the points to calculate the normalization.
    h : float
        A number representing the set height, usually between 0 and 1.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the normalized mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, U, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        U : Universe
            A Universe instance defining the interval or the points to calculate the alpha-cut.
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.term = term 
        self.h = max(A.vmu(U.getDenseU()))
        
    def mu(self, x):
        """Use to calculate the normalized mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            Normalized mu(x) considering U points. Where mu(x) is the membership degree of x to the set.
        """
        return self.A.mu(x)/self.h

class Concentration(_base):
    """A class used to concentrate a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    p : float
        The concentration factor must higher than 1 (p>1).
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the concentrated mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, p=2, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        p : float
            A number representing the concentration factor (p>1).
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.term = term
        self.p = p
        
    def mu(self, x):
        """Use to calculate the concentrated mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            Concentrated mu(x). Where mu(x) is the membership degree of x to the set.
        """
        return self.A.mu(x)**self.p

class Dilatation(_base):
    """A class used to dilatate a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    r : float
        The dilatation factor (0<r<1).
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the concentrated mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, r=0.5, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        r : float
            A number representing the dilatation factor (r<1).
        term : str
            The term used to identify the resulting set.
        """         
        self.A = A
        self.term = term
        self.r = r
        
    def mu(self, x):
        """Use to calculate the dilatated mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            Dilatated mu(x). Where mu(x) is the membership degree of x to the set.
        """        
        return self.A.mu(x)**self.r

class Contrast(_base):
    """A class used to contrast a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    p : float
        The constrast factor (p>1).
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the contrasted mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, p=2, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        p : float
            A number representing the contrast factor (p>1).
        term : str
            The term used to identify the resulting set.
        """          
        self.A = A
        self.term = term
        self.p = p
        
    def mu(self, x):
        """Use to Calculate the contrasted mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            Constrasted mu(x). Where mu(x) is the membership degree of x to the set.
        """         
        setmu = self.A.mu(x)
        if setmu <= 0.5:
            ret = (2**(self.p-1))*(setmu**self.p)
        else:
            ret = 1-2**(self.p-1)*(1-setmu)**self.p         
        return ret

class Fuzzyfication(_base):
    """A class used to the fuzzyfication of a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the Fuzzy mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        term : str
            The term used to identify the resulting set.
        """         
        self.A = A
        self.term = term
        
    def mu(self, x):
        """Use to calculate the mu fuzzyfication for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            Fuzzyfication of mu(x). Where mu(x) is the membership degree of x to the set.
        """    
        setmu = self.A.mu(x)
        if setmu <= 0.5:
            ret = (setmu/2)**0.5
        else:
            ret = 1-((1-setmu)/2)**0.5        
        return ret    

class Complement(_base):
    """A class used to calculate the complement of a given set.

    Attributes
    ----------
    A : FuzzySet
        A FuzzySet instance defining the mu function.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the complementary mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """ 
    
    def __init__(self, A, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining the mu function.
        term : str
            The term used to identify the resulting set.
        """           
        self.A = A
        self.term = term 
        
    def mu(self, x):
        """Use to calculate the complementary mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value
        Returns
        -------
            Complementary mu(x). Where mu(x) is the membership degree of x to the set.
        """           
        return 1-self.A.mu(x)

#################################################################
#
# TWO ARGUMENT OPERATIONS (SAME UNIVERSE)
#
#################################################################

class Intersection(_base):
    """A class used to calculate the intersection of many sets.

    Attributes
    ----------
    sets : FuzzySet
        A list of FuzzySet instances, [set1, set2, set3, ...]
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the resulting set mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, SETS, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        SETS : [FuzzySet]
            A list of FuzzySet instances, [set1, set2, set3, ...]
        term : str
            The term used to identify the resulting set.
        """         
        self.SETS = SETS
        self.term = term # union term
        
        
        #self.max_idx = 0
        #self.max_value = 0
        
        
    def mu(self, x):
        """Use to calculate the intersection mu for a given x.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            The resulting set mu(x). Where mu(x) is the membership degree of x to the set.
        """          
        mus = []
        for s in self.SETS:
            mus.append(s.mu(x))
        return min(mus)

class Union(_base):
    """A class used to calculate the union of many sets.

    Attributes
    ----------
    sets : FuzzySet
        A list of FuzzySet instances, [set1, set2, set3, ...]
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the resulting set(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, SETS, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        SETS : [FuzzySet]
            A list of FuzzySet instances, [set1, set2, set3, ...]
        term : str
            The term used to identify the resulting set.
        """            
        self.SETS = SETS
        self.term = term 
        
    def mu(self, x):
        """Use to calculate the union mu(x). Where mu(x) is the membership degree of x to the set.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            The resulting set mu(x). Where mu(x) is the membership degree of x to the set.
        """           
        mus = []
        for s in self.SETS:
            mus.append(s.mu(x))
        return max(mus)
       
# TODO compensation
# TODO distance
# TODO equal

class Inclusion(_base):
    """A class used to calculate the inclusion of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.
    U : Universe
        A Universe instance defining the interval or the points to calculate the degree of inclusion of A in B.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(x)
        Used to calculate the resulting set mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, B, U, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining a mu function.
        B : FuzzySet
            A FuzzySet instance defining a mu function.       
        U : Universe
            A Universe instance defining the interval or the points to calculate the inclusion.
        term : str
            The term used to identify the resulting set.
        value : 
            The inclusion value for A in B in the universe U (read only).
        """        
        self.A = A
        self.B = B
        self.U = U
        self.term = term
        ca = A.cardinality(U)
        MUA = pd.DataFrame(A.vmu(U.getDenseU()))
        MUB = pd.DataFrame(B.vmu(U.getDenseU()))
        self.value = (ca-(MUA-MUB)[0].apply(lambda x: max(0,x)).sum())/ca
        
    def mu(self, x):
        """Use to calculate the inclusion resulting set mu.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            The resulting set mu(x). Where mu(x) is the membership degree of x to the set.
        """        
        return Intersection([self.A,self.B]).mu(x)

class Possibility(_base):
    """A class used to calculate the possibility of two given sets.
    
    The possibility result is available at "value".
    The class mu method provides the intersection between the sets A and B in the universe U.

    Important Notes
    ---------------
    The points related to the possibility region will be added to the universe U.
    
    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.
    U : Universe
        A Universe instance defining the interval or the points to calculate the possibility.
    term : str
        The term used to identify the resulting set.
    value : 
        The possibility value for A and B in the universe U (read only).
        
    Methods
    -------
    mu(x)
        Used to calculate the possibility mu(x). Where mu(x) is the membership degree of x to the set.
            
    Inherited Methods
    -----------------
    vmu(X)
        Use to calculate the membership degree for a list of points X. (FuzzyBase._base)        
    cardinality(U)
        Use to calculate the cardinality of the set in the universe U. (FuzzyBase._base)
    """
    
    def __init__(self, A, B, U, term='', precision_mode=False):
        """Use to construct the object.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet instance defining a mu function.
        B : FuzzySet
            A FuzzySet instance defining a mu function.        
        U : Universe
            A Universe instance defining the interval or the points to calculate the possibility.
        term : str
            The term used to identify the resulting set.
        """          
        self.A = A
        self.B = B
        self.U = U
        self.term = term
        U.mapSet(A)
        U.mapSet(B)

        # mapear intervalo, calcular mu e mapear no universo
        if type(A) is SingletonSet:
            self.value = B.mu(A.m)
            return 
        
        if type(B) is SingletonSet:
            self.value = A.mu(B.m)
            return
        
        self.I = Intersection([self.A, self.B])       
    
        if precision_mode:
            u = U.getDenseU()
        else:
            u = U.getU()

        muvals = self.I.vmu(u)
        self.U.addPoints([u[muvals.argmax()]])
        self.value = max(muvals)

        
        
    def mu(self, x):
        """Use to calculate the possibility  resulting set mu.
        
        Parameters
        ----------
        x: float
            A number within the Universe to calculate a mu value.
            
        Returns
        -------
            The resulting set mu(x). Where mu(x) is the membership degree of x to the set.
        """          
        v = self.I.mu(x)
        return  v if v==self.value else 0
    
