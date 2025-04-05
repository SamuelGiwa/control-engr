import numpy as np
from mpc import mpcgain  

# Define system matrices
Ap = np.array([[1, 1], 
               [0, 1]])

Bp = np.array([[0.5], 
               [1]])  

Cp = np.array([[1, 0]])  

Dp = 0  

Np = 20
Nc = 4

# Get MPC gain matrices
Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e = mpcgain(Ap, Bp, Cp, Nc, Np)

# Dimensions
n, n_in = B_e.shape  # e.g. 2, 1

# Initialize states
xm = np.array([[0], [0]])  # Ensure xm is a 2x1 column vector
Xf = np.zeros((n, 1))

# Simulation parameters
N_sim = 100
r = np.ones((N_sim, 1))  # reference trajectory
u = 0
y = 0

# Preallocate storage
u1 = np.zeros((N_sim, 1))
y1 = np.zeros((N_sim, 1))

# Simulation loop
for kk in range(N_sim):
    # Compute optimal control increment
    DeltaU = np.linalg.inv(Phi_Phi + 0.1 * np.eye(Nc)) @ (Phi_R @ r[kk] - Phi_F @ Xf)
    deltau = DeltaU[0, 0]

    # Update control input
    u = u + deltau

    # Save inputs
    u1[kk] = u

    # Save current state
    xm_old = xm.copy()

    # Plant state update
    xm_new = Ap @ xm + Bp * u  # Ensure xm is a 2x1 column vector, Ap is 2x2, Bp is 2x1

    # Output computation
    y = Cp @ xm_new  # y will be a scalar, since Cp is 1x2 and xm_new is 2x1

    # Save output
    y1[kk] = y

    # Ensure y is reshaped as a 2x1 matrix to match the shape of xm
    y_reshaped = np.array([[y], [0]])  # Make y a 2x1 vector, append 0 for correct shape

    # Update the state vector (combine xm_new with the reshaped output)
    xm = np.vstack((xm_new, y_reshaped))  # Stack the new state and output

    # Update the extended state vector Xf
    Xf = np.vstack((xm - xm_old, y_reshaped))  # Stack the difference and output

# Optional: plot
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(y1, label='Output y')
plt.plot(u1, label='Control input u')
plt.plot(r, '--', label='Reference r')
plt.xlabel('Time Step')
plt.title('MPC Response')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()