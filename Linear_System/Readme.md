# Variable description for Ax = b
```
mat_M: A
vec_x: x
vec_b: b
T0: translation vectors for subrange algorithm
Matrix_size: dimension of the square matrix & size of vectors
```

## Qubit variables for xi
```
vector x = (x1, x2, ---, xn)
x1 = q1 + 2q2+ 4q3 + --- + 2^(n-1)qn - q(n+1) - 2q(n+2) - 4q(n+3) + --- + 2^(n-1)q2n
Qubits: n in xi
```

## QUBO matrix
```
QM: QUBO matrix
2*Qubits*Matrix_size: Dimension of QUBO matrix
T: qubit variable vector
```

## Results
```
Abs_Min: Absolute minimum value -(vec_b*)*vec_b (negative constant term comming from |Ax-b|^2) 
(T*)*QM*T: One of annealing results
vec_Anl: Absolute minimum vector of T
Annealing: Absolute minimum value
sol_x: real value calculated by the qubit vector vec_Anl
Ax_b: a vector A(sol_x) - b
```
