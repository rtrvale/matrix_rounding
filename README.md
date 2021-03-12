# matrix_rounding
Code for rounding a matrix to nearest integer, preserving row and column sums

A piece of code for rounding entries of a matrix to the nearest integer in such a way that row and column sums are preserved. This turns out to be more difficult than it looks, and the solution calculates the maximum flow through an associated graph. I am saving the code here in case I ever need to use it in the future.

In the original application, more severe modification of the matrix entries was sometimes needed, so there is a second version of the code in which the matrix is rescaled using the RAS algorithm in between attempts to use the matrix rounding algorithm.
