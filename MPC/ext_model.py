import numpy as np

def extmodel_(Am, Bm, Cm):
    """
        This creates an augemented matrix A, B and C based;
    """
    r_a, c_a = Am.shape
    r_c, c_c = (Cm @ Am).shape
    
    r = r_a + r_c
    c = c_a + 1
    
    A = np.zeros((r, c))
    A[-1, -1] = 1
    
    A[:, :c-1] = np.vstack((Am, Cm @ Am))
    B = np.vstack((Bm, Cm @ Bm))
    
    row, col = A.shape
    C = np.zeros((1, col))
    C[0, -1] = 1
    
    return A, B, C