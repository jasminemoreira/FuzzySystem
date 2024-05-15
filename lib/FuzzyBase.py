# -*- coding: utf-8 -*-
"""
Fuzzy base module.

This module is intended to provide the base for sets and operations.

Created on Tue Nov 28 21:10:19 2023

@author: Jasmine Moreira
"""
import numpy as np

class _base:
    """A base class used to create basic sets and operations.

    Methods
    -------
    vmu(x)
        Use to calculte the mu for each element of a list
    """
    
    def vmu(self, X):
        """Use to calculate the mu's for a list of X.
        
        Parameters
        ----------
        X : [float]
            A list of floats to calculate the mu
        Returns
        -------
            A list of the calculated mu for each X element 
        """
        ret = []
        idx = 0
        ant = 0

        for x in X: 
            ret.append(self.mu(x))
            
        return np.array(ret)


    
    def cardinality(self, U):
        """Use to calculate the cardinality of the set for the universe U.
        
        Parameters
        ----------
        U : Universe
            A Universe instance defining the interval or the points to calculate cardinality
        Returns
        -------
            The calculated cardinality for the set in the universe U 
        """ 
        return sum(self.vmu(U.getDenseU()))   


class _basev:
    """A base class used to create vectorized operations.

    Methods
    -------
    vmu(m1,m2)
        Use to calculate the mu's for a pair (x,y) from matrices or vectors m1 and m2.
    """
    
    def vmu(self,m1, m2):
        """Use to calculate the mu's for a pair (x,y) from matrices or vectors m1 and m2.
        
        Parameters
        ----------
        m1 : float[][]
            A matrix of floats to calculate the mu (x)
        m2 : float[][]
            A matrix of floats to calculate the mu (y)
            
        Returns
        -------
            A list of the calculated mu for each (x,y) pair 
        """          
        if np.shape(m1) != np.shape(m2):
            return "matrices must have the same dimensions"       
        result = np.zeros_like(m1, dtype=float)     
        for i in range(np.shape(m1)[0]):
            for j in range(np.shape(m1)[1]):
                result[i][j] = self.mu(m1[i][j],m2[i][j])
        return result  