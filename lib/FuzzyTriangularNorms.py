# -*- coding: utf-8 -*-
"""
Fuzzy triangular norms and co-norms

Created on Sun Oct 22 19:09:41 2023

@author: Jasmine Moreira
"""
class T_Norms:   
    def __init__(self):
        self
        
    def t1(self, a, b, pt):
        return 1/(1+(((1-a)/a)**pt+((1-b)/b)**pt)**(1/pt))

    def t2(self, a, b, pt):
        return max(0,(1+pt)*(a+b-1)-pt*a*b)
    
    def t3(self, a, b, pt):
        return 1-min(1,((1-a)**pt+(1-b)**pt)**(1/pt))
    
    def t4(self, a, b, pt=0):
        return a*b
        
    def t5(self, a, b, pt):
        return a*b/(pt+(1-pt)*(a+b-a*b))
        
    def t6(self, a, b, pt):
        return 1/(1/a**pt+1/b**pt-1)**(1/pt)
    
    def t7(self, a, b, pt):
        return (max(0,a**pt+b**pt-1))**(1/pt)
        
    def t8(self, a, b, pt=0):
        return a if b==1 else b if a==1 else 0
        
    def t9(self, a, b, pt):
        return min(a,b)

class S_Norms:
    def s1(self, a, b, pt):
        return 1/(1+((a/(1-a))**pt+(b/(1-b))**pt)**(1/pt))

    def s2(self, a, b, pt):
        return min(1,a+b-pt*a*b)
    
    def s3(self, a, b, pt):
        return min(1,(a**pt+b**pt)**(1/pt))
    
    def s4(self, a, b, pt=0):
        return a+b-a*b
        
    def s5(self, a, b, pt):
        return (a+b-a*b-(1-pt)*a*b)/(1-(1-pt)*a*b)
        
    def s6(self, a, b, pt):
        return 1-1/(1/(1-a)**pt+1/(1-b)**pt-1)**(1/pt)
    
    def s7(self, a, b, pt):
        return 1-max(0,((1-a)**pt+(1-b)**pt-1)**(1/pt))
        
    def s8(self, a, b, pt=0):
        return a if b==0 else b if a==0 else 1
        
    def s9(self, a, b, pt):
        return max(a,b)
    
    
    