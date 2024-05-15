# -*- coding: utf-8 -*-
"""
Fuzzy relations module.

This module is intended to provide some fuzzy rules.

Created on Thu Nov  2 13:21:06 2023

@author: Jasmine Moreira
"""
import numpy as np
from FuzzyTriangularNorms import T_Norms, S_Norms
from FuzzyBase import _basev

class FuzzyRelation(_basev):
    """A flexible class to create a fuzzy relation.
    
    This class takes a function to calculate the membership degree of a provided (x,y) pair.

    Attributes
    ----------
    function : function
        A user defined relation function to support mu calculation.
    term : str
        The term used to identify the relation.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the S norm mu for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (x,y) from matrices or vectors m1 and m2. (FuzzyBase._basev)
    """
    
    def __init__(self, function, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        function : function
            A user defined relation function to support mu calculation.
        term : str
            The term used to identify the relation.
        """          
        self.term = term
        self.mf = function
        
    def mu(self, a, b):
        """Use to calculate the membership degree for a pair (a,b).
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """        
        return self.mf(a, b)


#################################################################
#
# TWO ARGUMENT OPERATIONS (DIFFERENT UNIVERSES)
#
#################################################################
        
class TNorm(_basev):
    """A class used to calculate the TNorm of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.
    norm : int
        The selected T norm to calculate the mu, 1 to 9.
    pt : float
        Factor used in the triangular norms as follows:
        pt>0 for t1, t3 and t7.
        pt >= -1 for t2.
        pt >=0 for t5.
        pt not used for t4, t6, t8 and t9.      
    term : str
        The term used to identify the resulting set.
    symmetric: bool
        Use to apply the same points to each set in order to get mu.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the T norm mu for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation. 
            
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)
    """    
    
    def __init__(self, A, B, norm, pt=1, term='', symmetric=False):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.       
        pt : float
            Factor used in the triangular norms as follows:
            pt>0, t1, t3 and t7.
            pt >=0, t2 and t5.
            pt not used for t4, t6, t8 and t9.          
        term : str
            The term used to identify the resulting set.
        symmetric: bool
            Use to apply the same points to each set in order to get mu.
        """
        self.A = A
        self.B = B
        self.norm = norm
        self.pt = pt
        self.term = term
        self.symmetric = symmetric
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the TNorm for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """          
        t_norm = T_Norms()
        norm = getattr(t_norm, 't'+str(self.norm))
        
        if type(self.A) is FuzzyRelation and type(self.B) is FuzzyRelation:
            return norm(self.A.mu(a,b), self.B.mu(a,b), self.pt)
    
        if self.symmetric:
            if c=='':
                return norm(self.A.mu(a, b), self.B.mu(a, b), self.pt)
            else:
                return norm(self.A.mu(a, b, c), self.B.mu(a, b, c), self.pt)       
        else:
            if c=='':
                return norm(self.A.mu(a), self.B.mu(b), self.pt)
            elif d=='':
                return norm(self.A.mu(a, b), self.B.mu(c), self.pt)
            else:
                return norm(self.A.mu(a, b), self.B.mu(c, d), self.pt)
        
    
class SNorm(_basev):
    """A class used to calculate the SNorm of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.
    norm : int
        The selected S norm to calculate the mu, 1 to 9.
    pt : float
        Factor used in the triangular norms as follows:
        pt>0, s1, s3 and s7.
        pt >=0, s2 and s5.
        pt not used for s4, s6, s8 and s9.       
    term : str
        The term used to identify the resulting set.
    symmetric: bool
        Use to apply the same points to each set in order to get mu.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the S norm mu for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, norm, pt=1, term='', symmetric=False):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.       
        pt : float
            Factor used in the triangular norms as follows:
            pt>0, s1, s3 and s7.
            pt >=0, s2 and s5.
            pt not used for s4, s6, s8 and s9.        
        term : str
            The term used to identify the resulting set.
        symmetric: bool
            Use to apply the same points to each set in order to get mu.
        """        
        self.A = A
        self.B = B
        self.norm = norm
        self.pt = pt
        self.term = term
        self.symmetric = symmetric
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the SNorm for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """        
        s_norm = S_Norms()
        norm = getattr(s_norm, 's'+str(self.norm))
        
        if type(self.A) is FuzzyRelation and type(self.B) is FuzzyRelation:
            return norm(self.A.mu(a,b), self.B.mu(a,b), self.pt)
        
        if self.symmetric:
            if c=='':
                return norm(self.A.mu(a, b), self.B.mu(a, b), self.pt)
            else:
                return norm(self.A.mu(a, b, c), self.B.mu(a, b, c), self.pt)       
        else:
            if c=='':
                return norm(self.A.mu(a), self.B.mu(b), self.pt)
            elif d=='':
                return norm(self.A.mu(a, b), self.B.mu(c), self.pt)
            else:
                return norm(self.A.mu(a, b), self.B.mu(c, d), self.pt)
        
    

#####################################################################
#
# Fuzzy Conjunctions (cartesian product generalizations using T norm)
#
#####################################################################

class MamdaniConjunction(_basev):
    """A class used to calculate the Mamdani conjunction of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the Mamdani's Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.             
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Mamdani's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """   
        if c=='':
            return T_Norms().t9(self.A.mu(a), self.B.mu(b), pt=1)
        elif d=='':
            return T_Norms().t9(self.A.mu(a, b), self.B.mu(c), pt=1)
        else:
            return T_Norms().t9(self.A.mu(a, b), self.B.mu(b, c), pt=1)
        

class LarsenConjunction(_basev):
    """A class used to calculate the Larsen conjunction of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
    symmetric : bool
        If True, a, b, c points will be applied to A and B in order to get the mu's.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the Larsen's Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.             
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Larsen's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """
        if c=='':
            return T_Norms().t4(self.A.mu(a), self.B.mu(b), pt=1)
        elif d=='':
            return T_Norms().t4(self.A.mu(a, b), self.B.mu(c), pt=1)
        else:
            return T_Norms().t4(self.A.mu(a, b), self.B.mu(b, c), pt=1)


#####################################################################
#
# Fuzzy Disjunctions (cartesian product generalizations using S norm)
#
#####################################################################

class MaxDisjunction(_basev):
    """A class used to calculate the max disjunction of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the max Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.             
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the max disjunction mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """

        if c=='':
            return S_Norms().s9(self.A.mu(a), self.B.mu(b), pt=1)
        elif d=='':
            return S_Norms().s9(self.A.mu(a, b), self.B.mu(c), pt=1)
        else:
            return S_Norms().s9(self.A.mu(a, b), self.B.mu(b, c), pt=1)



#######################################################################
#
# Fuzzy Implications (generalization of multi values Lukasiewicz logic)
#
#######################################################################

# TODO: terminar de implementar para a, b, c, d

class LukasiewiczImplication(_basev):
    """A class used to calculate the Lukasiewicz implication for two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the max Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        term : str
            The term used to identify the resulting set.

        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Lukasiewicz mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """
        if c=='':
            return min(1, (1-self.A.mu(a)) + self.B.mu(b))
        elif d=='':
            return min(1, (1-self.A.mu(a, b)) + self.B.mu(b))
        else:
            return min(1, (1-self.A.mu(a, b)) + self.B.mu(b, c))



class GodelImplication(_basev):
    """A class used to calculate the Godel's implication for two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
    symmetric : bool
        If True, a, b, c points will be applied to A and B in order to get the mu's.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the max Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Godel's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """
        if c=='':
            return 1 if self.A.mu(a) <= self.B.mu(b) else self.B.mu(b)
        elif d=='':
            return 1 if self.A.mu(a, b) <= self.B.mu(c) else self.B.mu(c)
        else:
            return 1 if self.A.mu(a, b) <= self.B.mu(c, d) else self.B.mu(c, d)

    
class KleeneImplication(_basev):
    """A class used to calculate the Kleene's implication for two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the max Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, A, B, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        term : str
            The term used to identify the resulting set.
        """        
        self.A = A
        self.B = B
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Kleene's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """   
        if c=='':
            return max(1-self.A.mu(a), self.B.mu(a))
        elif d=='':
            return max(1-self.A.mu(a, b), self.B.mu(c))
        else:
            return max(1-self.A.mu(a, b), self.B.mu(c, d))


#####################################################################
#
# Fuzzy Comparing (same dimensions comparing)
#
#####################################################################

class Minimum(_basev):
    """A class used to calculate the Minimum of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
    symmetric : bool
        If True, a, b, c points will be applied to A and B in order to get the mu's.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the Mamdani's Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, SETS, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.             
        term : str
            The term used to identify the resulting set.
        """        
        self.SETS = SETS
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Mamdani's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """
        if c=='':
            ret = self.SETS[0].mu(a, b)
            for cnt in range(1, len(self.SETS)):
                ret = min(ret, self.SETS[cnt].mu(a, b))
            return ret
        else:
            ret = self.SETS[0].mu(a, b, c)
            for cnt in range(1, len(self.SETS)):
                ret = min(ret, self.SETS[cnt].mu(a, b, c))
            return ret

class Maximum(_basev):
    """A class used to agregate rules or calculate the Maximum of two given sets.

    Attributes
    ----------
    A : FuzzySet
        The first FuzzySet.
    B : FuzzySet
        The second FuzzySet.     
    term : str
        The term used to identify the resulting set.
    symmetric : bool
        If True, a, b, c points will be applied to A and B in order to get the mu's.
        
    Methods
    -------
    mu(a, b)
        Used to calculate the Mamdani's Rule mu value for a and b,
        it cant handle FuzzySets and FuzzyRelations,
        see method documentation.
        
    Inherited Methods
    -----------------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (a,b) from matrices or vectors m1 and m2. (FuzzyBase._basev)   
    """
      
    def __init__(self, SETS, term=''):
        """Use to construct the object.
        
        This class can't mix a FuzzySet and a FuzzyRelation.
        
        Parameters
        ----------
        A : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.
        B : FuzzySet
            A FuzzySet or FuzzyRelation instance defining a mu function.             
        term : str
            The term used to identify the resulting set.
        """        
        self.SETS = SETS
        self.term = term
            
    def mu(self, a, b, c='', d=''):
        """Use to calculate the Mamdani's mu for two FuzzySets or two FuzzyRelations.
        
        Parameters
        ----------
        a: float
            A number within the A related universe to calculate mu value.
        b: float
            A number within the B related universe to calculate mu value .           
        
        Returns
        -------
            The resulting set mu(a,b). Where mu(a,b) is the membership degree of the pair.
        """             
        if c=='':
            ret = self.SETS[0].mu(a, b)
            for cnt in range(1, len(self.SETS)):
                ret = max(ret, self.SETS[cnt].mu(a, b))
            return ret
        else:
            ret = self.SETS[0].mu(a, b, c)
            for cnt in range(1, len(self.SETS)):
                ret = max(ret, self.SETS[cnt].mu(a, b, c))
            return ret

