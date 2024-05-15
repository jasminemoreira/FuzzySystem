# -*- coding: utf-8 -*-
"""
Fuzzy plot functions.

Created on Mon Oct 23 09:03:04 2023

@author: Jasmine Moreira
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from FuzzySystem import FuzzyPartition
from FuzzyBasicSets import SingletonSet, DiscreteSet
from FuzzyOperations import Possibility

def plotMatrix(X, Y, A):
    """Use to print the A FuzzyRelation mu values matrix for X, Y universes.
    
    Parameters
    ----------
    X : Universe
        A universe to provide x values for mesh generation and mu calculation.
    Y : Universe
        A universe to provide y values for mesh generation and mu calculation.
    A : FuzzyRelation, TNorm, SNorm 
        A two sets fuzzy operation wich can handle an (x,y) pair.
    """    
    pd.set_option('display.width',1000)
    MX, MY = np.meshgrid(X.getU(), Y.getU())   
    pd.set_option('display.precision', 2)
    df = pd.DataFrame(A.vmu(MX, MY), columns=X.getU(), index=Y.getU())
    return df.style.format(precision=2) \
      .background_gradient(cmap ='jet', axis=None, vmin=0, vmax=1) \
      .format_index(precision=2, axis=0) \
      .format_index(precision=2, axis=1) 
      
def plot4D(X, Y, Z, RULE, title='', mask=0.01, alpha=0.2):  
    """Use to plot the a fuzzy rule for X, Y and Z universes.
    
    Parameters
    ----------
    X : Universe
        A universe to provide x values for the rule mu calculation.
    Y : Universe
        A universe to provide y values for the rule and mu calculation.
    Z : Universe
        A universe to provide z values for the rule and mu calculation.    
    RULE : FuzzyRelation, FuzzyRule 
        A rule wich can manipulate x, y and z coordinates in order to produce a mu value.
    title : str
        The plot title.
    mask : float
        The minimum mu value for plotting (points with mu below this value won't be shown).
    alpha: float
        A transparency factor for all the points in order to allow internal regions visualization.
    """
    data = getMatrixFromRule(X, Y, Z, RULE)   
    ticks = 5
    xticks_idx = [int(item) for item in list(range(0,X.npoints,int(X.npoints/ticks)))]
    xticks_labels = np.round(X.getU(),2)[xticks_idx]
    yticks_idx = [int(item) for item in list(range(0,Y.npoints,int(Y.npoints/ticks)))]
    yticks_labels = np.round(Y.getU(),2)[yticks_idx]    
    zticks_idx = [int(item) for item in list(range(0,Z.npoints,int(Z.npoints/ticks)))]
    zticks_labels = np.round(Z.getU(),2)[zticks_idx]
    
    plt.figure()
    ax = plt.axes(projection='3d')
    mask = data > mask
    idx = np.arange(int(np.prod(data.shape)))
    x, y, z = np.unravel_index(idx, data.shape)
    ax.scatter(x, y, z, c=data.flatten(), s=20.0 * mask, edgecolor="face", 
               alpha=alpha, marker="o", cmap='jet', linewidth=0, vmin=0, vmax=1)
    ax.set_xticks(xticks_idx, labels=xticks_labels)
    ax.set_yticks(yticks_idx, labels=yticks_labels)
    ax.set_zticks(zticks_idx, labels=zticks_labels)
    ax.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
    ax.set_title(title)
    ax.set_xlabel(X.term, labelpad=5)
    ax.set_ylabel(Y.term, labelpad=10)
    ax.set_zlabel(Z.term, labelpad=5)
    plt.show()
    
def getMatrixFromRule(X, Y, Z, RULE):
    """Use to get the mu matrix for a rule, using X, Y and Z universes.
    
    Parameters
    ----------
    X : Universe
        A universe to provide x values for the rule mu calculation.
    Y : Universe
        A universe to provide y values for the rule and mu calculation.
    Z : Universe
        A universe to provide z values for the rule and mu calculation.    
    RULE : FuzzyRelation, FuzzyRule 
        A rule wich can manipulate x, y and z coordinates in order to produce a mu value.
    """
    rule = RULE.mu
    data = np.zeros((X.npoints, Y.npoints, Z.npoints), dtype=float)
    xu = X.getU()
    yu = Y.getU()
    zu = Z.getU()
    for x in range(0,len(xu)):
        for y in range(0,len(yu)):
            for z in range(0,len(zu)):
                data[x, y, z] = rule(xu[x], yu[y], zu[z])
    return data

def plotMesh(X, Y, A, title=''):
    """Use to create a 3D plot for X,Y meshes and A FuzzyRelation.
    
    Parameters
    ----------
    X : Universe
        A universe to provide x values for mesh generation and mu calculation.
    Y : Universe
        A universe to provide y values for mesh generation and mu calculation.
    A : FuzzyRelation, TNorm, SNorm 
        A two sets fuzzy operation wich can handle an (x,y) pair.
    title : str
        The plot title.    
    """
    xvals = X.getU()
    yvals = Y.getU()
    xmg,ymg = np.meshgrid(xvals, yvals)
    zmg = A.vmu(xmg,ymg)
    plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(xmg, ymg, zmg, cmap='jet', edgecolor='green')
    ax.set_title(title) 
    ax.set_xlabel(X.term, labelpad=5)
    ax.set_ylabel(Y.term, labelpad=10)
    ax.set_zlabel(A.term, labelpad=5)
    plt.show()

def plot(X,SETS, title='', colors='', show_overlaid=False):
    """Use to create a plot for one or more sets in the X universe.
    
    Parameters
    ----------
    X : Universe
        A universe to provide x values for plot generation and mu calculation.
    SETS : [] FuzzySet, FuzzyOperation...
        A list of fuzzy sets or fuzzy operations.
    title : str
        The plot title.
    colors : [str]
        A list of colors to be applied.
    show_overlaid: bool
        Change the size of each set plot in order to allow overlaid points visualization.
    """
    plt.figure()
    legends = []
    sets = []
    if type(SETS) is FuzzyPartition:
        if title == '':
            title = SETS.term
        for key in list(SETS.sets.keys()):
            sets.append(SETS.sets[key])
        size = len(sets)
    elif type(SETS) is list:
        sets = SETS
        size = len(sets)  
    else:
        sets = [SETS]
        if colors == '':
            colors=['tab:blue']
        else:
            colors = [colors]
        size = 1 

    if colors=='':
        colors = __generate_tableau_palette(size)
        
    times = 1
    zorder = 1000
    for set_,color in zip(sets, colors):
        X.mapSet(set_)
        Y, s = __getY(set_, X)
        if show_overlaid or s==0:
            s = 7
        if np.nanmax(Y) > 0:
            legends.append(set_.term)
            plt.scatter(X.getU(), Y, c=color,s=int(s**times), zorder=zorder)
            if show_overlaid:
                times += 0.3
                zorder -= 1
    plt.title(title)
    if not '' in legends:
        plt.legend(legends)
    plt.show()

def __getY(set_,X):
    s = 0
    if type(set_) is DiscreteSet or X.discrete==True:
        s = 15   
    Y = pd.DataFrame(set_.vmu(X.getU()))   
    if type(set_) is SingletonSet:
        Y.replace(0, np.nan, inplace=True)
        s = 50
    elif type(set_) is Possibility:
        Y.replace(0, np.nan, inplace=True)
        s = 50
    return Y,s

def __generate_tableau_palette(N):
    colors = ['blue', 'orange','red', 'green', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    palette = [colors[i % len(colors)] for i in range(N)]
    return palette