# Model Predictive Control System Design and Implementation

This project implements a **Model Predictive Control (MPC)** system in Python, based on the methodology and examples from the book:

📘 **"Model Predictive Control System Design and Implementation Using MATLAB"**  
✍️ *by Liuping Wang*

## 📌 Author

**Written by:** Samuel Giwa Boluwatife

🔗 [LinkedIn](http://www.linkedin.com/in/samuelgiwa)  

---

## 🛠️ Description

This project demonstrates how to design and simulate a simple MPC controller for a discrete-time linear system using Python and NumPy. The implementation includes:

- State-space modeling of the plant
- Construction of MPC gain matrices
- Simulation of control behavior over a prediction horizon
- Visualization of control input and system output

---

## 🧮 System Dynamics

The controlled system is defined in state-space form:

```python
Ap = np.array([[1, 1],
               [0, 1]])
Bp = np.array([[0.5],
               [1]])
Cp = np.array([[1, 0]])