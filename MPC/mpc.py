import numpy as np
from typing import Tuple

def mpcgain(
    Ap: np.ndarray,
    Bp: np.ndarray,
    Cp: np.ndarray,
    Nc: int,
    Np: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute MPC gain matrices for a discrete-time linear system.

    Parameters:
    - Ap (np.ndarray): State matrix A
    - Bp (np.ndarray): Input matrix B
    - Cp (np.ndarray): Output matrix C
    - Nc (int): Control horizon
    - Np (int): Prediction horizon

    Return value:
    - Phi_Phi (np.ndarray)
    - Phi_F (np.ndarray)
    - Phi_R (np.ndarray)
    - A_e (np.ndarray): Augmented system matrix
    - B_e (np.ndarray): Augmented input matrix
    - C_e (np.ndarray): Augmented output matrix
    """
    m1, n1 = Cp.shape
    n1, n_in = Bp.shape
    
    A_e = np.eye(n1 + m1, n1 + m1)
    A_e[:n1, :n1] = Ap
    A_e[n1:n1 + m1, :n1] = Cp @ Ap
    
    B_e = np.zeros((n1 + m1, n_in))
    B_e[:n1, :] = Bp
    B_e[n1:n1 + m1, :] = Cp @ Bp
    
    C_e = np.zeros((m1, n1 + m1))
    C_e[:, n1:n1 + m1] = np.eye(m1)
    
    n = n1 + m1
    h = np.zeros((Np, n))
    F = np.zeros((Np, n))
    
    h[0, :] = C_e
    F[0, :] = C_e @ A_e
    
    for kk in range(1, Np):
        h[kk, :] = h[kk - 1, :] @ A_e
        F[kk, :] = F[kk - 1, :] @ A_e
    
    v = h @ B_e
    
    Phi = np.zeros((Np, Nc))
    Phi[:, 0] = v[:, 0]
    
    for i in range(1, Nc):
        Phi[:, i] = np.hstack((np.zeros(i), v[:Np - i, 0]))
    
    BarRs = np.ones((Np, 1))
    
    Phi_Phi = Phi.T @ Phi
    Phi_F = Phi.T @ F
    Phi_R = Phi.T @ BarRs
    
    return Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e


Ap = np.array([[0.8]])  # 1x1 matrix
Bp = np.array([[0.1]])  # 1x1 matrix
Cp = np.array([[1]])    # 1x1 matrix
Nc = 4
Np = 10


Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e = mpcgain(Ap, Bp, Cp, Nc, Np)


print("Phi_Phi:\n", Phi_Phi)
print("Phi_F:\n", Phi_F)
print("Phi_R:\n", Phi_R)
print("A_e:\n", A_e)
print("B_e:\n", B_e)
print("C_e:\n", C_e)