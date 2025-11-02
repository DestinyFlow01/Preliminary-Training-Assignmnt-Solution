# Preliminary Training Assignment Solution

Solution for the preliminary training before entering the real code.

---

## Tasks

1. **Approximating boundary layer parameters**  
   Compute displacement thickness, momentum thickness, drag coefficient, etc. for various velocity profiles (Pohlhausen, Schlichting, Majdalani–Xuan, etc.) and compare them to Blasius results for flow over a flat plate.  
   *Focus:* syntax, numerical integration/derivatives, and design thinking.

2. **Solving the original Blasius equation numerically and obtain the corresponding boundary layer parameters**  
   This problem focuses on numerical derivatives, ODE solvers on semi‑infinite domains and domains with singularities, inverse problems, optimization algorithm, code design decision and code efficiency (not limited to space and time complexity, but these two are essential).  

   The breakdown of methods:
   - **a. Shooting method**  
   - **b. Pair of IVPs**  
   - **c. Crocco's transformation**
   - **d. Solving energy equation for flow over flat plate with $Ec \neq 0$**

   For all of these methods, you will need to determine the Blasius constant. You may take it directly from *Frank M. White*, but it is better to start with a random initial value of σ and optimize it until the result matches the known value with at least **10⁻¹⁰ % accuracy**.

---

## Note

- Task 1 may be implemented inefficiently (no parallelization or optimization required).  
- Starting from Task 2, **you are required to make the code as efficient as possible** in almost any possible aspect, since the experience on making these solvers will be used intensively in the real code.
