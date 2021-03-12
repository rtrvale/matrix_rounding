# -*- coding: utf-8 -*-
"""
Spyder Editor

DESCRIPTION:
    This module provides functions to force a 2D array to have desired row 
    and column sums. Its basic usage is to do automatic rounding of an array
    in Excel, in order to force it to have a given vector of row sums and 
    column sums. The entries of the input array are expected to be positive.
    
    For example, the 3x3 array [[1.2, 3.4, 2.4], [3.9, 4.0, 2.1], [7.9, 1.6, 0.5]]
    is to be rounded to an array of integers in such a way that the row sums of 
    the rounded array are [7, 10, 10] and the column sums are [13, 9, 5]. One 
    possible solution is [[1, 3, 3], [4, 4, 2], [8, 2, 0]].
    
    This problem, in which each entry is rounded up or down to the nearest integer,
    can be solved with the smooth_and_round_array function. 
    
    For some applications we needed to round further than the nearest integer.
    For this purpose, steps of the RAS algorithm were added in the 
    rescale_smooth_and_round_array function.
    
EXAMPLES:
    
X = np.array([[1.2, 3.4, 2.4], [3.9, 4.0, 2.1], [7.9, 1.6, 0.5]])
row_sums = np.array([7, 10, 10])
col_sums = np.array([13, 9, 5])

A = smooth_and_round_array(X, row_sums, col_sums)
A.sum(axis=0) # equal to col_sums
A.sum(axis=1) # equal to row_sums

X[0, 0] = 5
B = rescale_smooth_and_round_array(X, row_sums, col_sums)
B[1].sum(axis=0) # equal to col_sums
B[1].sum(axis=1) # equal to row_sums
    
"""

import numpy as np
import pandas as pd
import os
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow


def prod_target_matrix(row_sums, col_sums):
    """
    organises the row and column sums into a matrix so that the max flow
    algorithm can be applied.
    
    Parameters
    ----------
    row_sums : list of ints with length (M)
    col_sums : list of ints with length (N)
    
    Returns
    -------
    array : ndarray with M x N
    """
        
    m = len(row_sums)
    n = len(col_sums)
    matrix = []
    
    matrix += [[0] + list(row_sums) + [0] * n + [0]]
    
    for i in range(m): 
        matrix += [[0] * (m + 1) + [1] * n + [0]]
    
    for i in range(n): 
        matrix += [[0] * (m + n + 1) + [col_sums[i]]]
    
    matrix += [[0] * (m + n + 2)]
    
    return matrix


def solve_for_row_and_column_sums(row_sums, col_sums):
    """
    applies the max flow algorithm for matrix rounding
    
    """
    
    m = len(row_sums)    
    n = len(col_sums)
    
    # produce a sparse matrix object and solve
    graph = csr_matrix(prod_target_matrix(row_sums, col_sums))
    
    solution = maximum_flow(graph, 0, m+n+1).residual.toarray()[1:(m+1), (m+1):(m+n+1)]
    
    return solution

def smooth_and_round_array(array, row_sums, col_sums):
    '''
    round off an array to get desired row and column sums
  
    array :    np.array or pandas dataframe; not rounded
    row_sums : np.array or list, desired row sums
    col_sums : np.array or list, desired column sums
    '''

    # check that the row and column sums are consistent
    assert (np.array(row_sums).sum() - np.array(col_sums).sum() < 0.01)

    # get the row and column sums for the fractional parts
    round_row_sums = (row_sums - np.floor(array).sum(axis=1)).astype('int')
    round_col_sums = (col_sums - np.floor(array).sum(axis=0)).astype('int')

    # solve for a binary matrix
    rounded_solution = solve_for_row_and_column_sums(list(round_row_sums), list(round_col_sums))

    # return the integer part plus the solution for the fractional part
    solution = np.floor(array) + rounded_solution

    return solution

  
def rescale_smooth_and_round_array(org_array, row_sums, col_sums, max_try=100):
    '''
    make an array have the correct row and column sums, including rescaling
    steps if required. This function is a wrapper for the
    smooth_and_round_array function and includes iterations of the RAS 
    algorithm to get close enough to a solution for 
    smooth_and_round_array to work

    Parameters
    ----------
    array : numpy array
    row_sums : list
    col_sums : list

    Returns
    -------
    tuple (solved, solution)
    
    solved : Boolean; True if solution was reached
    solution : numpy array of integers

    '''
    
    # the numbers don't matter, as long as tol < 1
    tol = 0.5
    err = 1
    count = 0    
     
    # generate initial solution
    array = np.copy(org_array)
    
    # generate initial solution
    solution = smooth_and_round_array(array, row_sums, col_sums)
    
    # calculate current error
    err = np.abs(solution.sum(axis=0) - col_sums).sum() + np.abs(solution.sum(axis=1) - row_sums).sum()
    
    while (err > tol) and (count < max_try):
                 
        # rescale array to have the correct column sums
        array = np.array(pd.DataFrame(array).divide(array.sum(axis=0))) * np.array(col_sums)
         
        # rescale array to have the correct row sums
        array = array.transpose()
        array = np.array(pd.DataFrame(array).divide(array.sum(axis=0))) * np.array(row_sums)
        array = array.transpose()
         
        # generate new solution
        solution = smooth_and_round_array(array, row_sums, col_sums)
        
        # calculate current error
        err = np.abs(solution.sum(axis=0) - col_sums).sum() + np.abs(solution.sum(axis=1) - row_sums).sum()
        
        count += 1
        
    solved = True if err <= tol else False
    
    return solved, solution
