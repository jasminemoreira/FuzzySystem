# -*- coding: utf-8 -*-
"""
Fuzzy Structures Module.

This module is intended to provide the basic fuzzy structures 

Created on Sun Oct 15 11:55:02 2023

@author: Jasmine Moreira
"""
import time
from copy import copy
import numpy as np
import pandas as pd
from FuzzyOperations import Possibility, Union, Intersection, epsilon
from FuzzyBasicSets import SingletonSet, FuzzySet
from FuzzyTriangularNorms import T_Norms
from FuzzyRelations import MamdaniConjunction, LarsenConjunction, MaxDisjunction, LukasiewiczImplication, GodelImplication, KleeneImplication
            
class FuzzyPartition:
    """A class to create a fuzzy partition.
    
    A fuzzy partition is a colection many fuzzy sets, representing a linguistic variable.

    Attributes
    ----------
    sets : [FuzzySet]
        A list of fuzzy sets in the same universe.
    name : str
        The name of the linguistic variable.
        
    Methods
    -------
    add_sets(fuzzy_sets)
        Used to add sets to the partition.      
    """
    
    def __init__(self, sets=None, term=''):
        """Use to construct the object.
        
        Parameters
        ----------
        sets : [FuzzySet]
            A list of fuzzy sets in the same universe.
        name : str
            The name of the linguistic variable.
        """
        self.term = term 
        self.sets = {}
        if sets != None:
            self.add_sets(sets)
            
    def __getitem__(self, term):
        """Use add get a member set by its term.
        
        Parameters
        ----------
        term: str
            The term to reference the member set.         
        """        
        return self.sets[term]
                  
    def add_sets(self, fuzzy_sets):
        """Use add sets to the partition.
        
        Parameters
        ----------
        fuzzy_sets: [FuzzySet]
            A list containing the fuzzy sets to be added.         
        """          
        for fset in fuzzy_sets:
            self.sets[fset.term] = fset
            setattr(self, fset.term, fset)
             
            

class FuzzySystem:
    
    def __init__(self, UNIV, ACR, name=''):       
        self.UNIV = UNIV
        self.ACR = ACR
        self.reset();
        
    def reset(self):
        self.ICONC = []
        self.FCONC = None
        self.CONC = None
        self.MoM = 0
        self.Centroid = 0
        
      
    def inference(self, INP, aggregator=9, semantics="Mamdani", defuzzyfication="centroid", precision_mode=False):
            
        self.reset();
        ACR_DF = pd.DataFrame(self.ACR)
        ANT = ACR_DF.iloc[:,0:-1].values.tolist()
        CONS = ACR_DF.iloc[:,-1:].values.tolist()

        UNIVA = self.UNIV[:-1]
        UNIC = self.UNIV[-1]

        # mapeamento dos consequentes no universo
        for con in CONS:
            UNIC.mapSet(con)
            
        # matching das entradas com antecedentes da regra
        MATCH = np.zeros([len(ANT),len(ANT[0])], dtype=float)
        lin = 0
        for Ant in ANT:
            col = 0
            for Set, Inp, Uni in zip(Ant, INP, UNIVA):
                MATCH[lin,col] = Possibility(Inp, Set, Uni, precision_mode=precision_mode).value
                col +=1
            lin += 1

        # agregação dos antecedentes (pode ser por qualquer norma T)
        t_norm = T_Norms()
        norm = getattr(t_norm, 't'+str(aggregator))
        AGANT = np.zeros([len(ANT)], dtype=float)
        lin = 0
        for Match in MATCH:
            aggval = Match[0]
            for cntr in range(1,len(Match)):
                aggval = norm(aggval,Match[cntr],pt=1)           
            AGANT[lin] = aggval 
            lin += 1
    
        # conclusão individual
        ICONC = []
        for idx in range(0, len(AGANT)):
            cons = CONS[idx][0]
            agant = AGANT[idx]          
            conc = self._set_factory(cons, agant, semantics)
            conc.term = 'R'+str(idx+1)+": "+cons.term
            ICONC.append(conc)

        self.ICONC = ICONC
        self.FCONC = Union(ICONC, 'Conclusion')
        self.UNIC = UNIC        
        
        # defuzzificação          
        if defuzzyfication=="Centroid":
            self.calcCentroid()
        elif defuzzyfication=="MoM":        
            self.calcMoM()
 
            
    def _set_factory(self, cons, agant, semantics):
        if semantics == "Mamdani":
            conc = FuzzySet(lambda a: min(cons.mu(a), agant))
        elif semantics == "Larsen":
            conc = FuzzySet(lambda a: cons.mu(a) * agant)
        elif semantics == "MaxDisjunction":
            conc = FuzzySet(lambda a: max(agant, cons.mu(a)))
        elif semantics == "Lukasiewicz":
            conc = FuzzySet(lambda a: min(1, (1-agant)+cons.mu(a)))
        elif semantics == "Godel":
            conc = FuzzySet(lambda a: 1 if agant < cons.mu(a) else cons.mu(a))
        elif semantics == "Kleene":              
            conc = FuzzySet(lambda a: max(1-agant, cons.mu(a)))
        else:
            raise ValueError('Semantics not found.')
        return conc        

        
    def calcCentroid(self):
        uval = self.UNIC.getU()
        res  = self.FCONC.vmu(uval)
        prodmucx = np.matmul(res, uval) 
        summu = np.sum(res) 
        self.Centroid = prodmucx/summu if summu > 0 else 0
        return self.Centroid
        
    def calcMoM(self):
        uval = self.UNIC.getU()
        res  = self.FCONC.vmu(uval)
        maxidx_first = res.argmax()
        maxidx_last = len(res)-res[::-1].argmax()-1
        self.MoM = uval[int((maxidx_last+maxidx_first)/2)] 
        return self.MoM         
        
    def getConclusionAsPartition(self):
        self.calcCentroid()
        self.calcMoM()
        MoMSet = SingletonSet(self.MoM, epsilon, self.UNIC.term+' (MoM)')
        CentroidSet = SingletonSet(self.Centroid, epsilon, self.UNIC.term+' (Centroid)')        
        self.UNIC.mapSet(MoMSet)
        self.UNIC.mapSet(CentroidSet)   
        return FuzzyPartition([self.FCONC, MoMSet, CentroidSet])
    