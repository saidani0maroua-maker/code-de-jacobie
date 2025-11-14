import numpy as np

def fill_2D_table(n):
    tab2D = np.zeros((n, n), dtype=int)
    print(f"\nEnter the values of the {n}×{n} matrix A ({n*n} total elements):\n")
    for i in range(n):
        for j in range(n):
            tab2D[i, j] = int(input(f"A[{i}][{j}]: "))
    return tab2D

def fill_1D_table(n):
    tab1D = np.zeros(n, dtype=int)
    print(f"\nNow enter the values of vector b ({n} elements):\n")
    for i in range(n):
        tab1D[i] = int(input(f"b[{i}]: "))
    return tab1D


def checking(tab2D, n):
    sums = np.zeros(n)
    node = 1
    for i in range(n):
        sum_rest = 0
        for j in range(n):
            if i != j:
                sum_rest += abs(tab2D[i, j])
        sums[i] = sum_rest
        pivot = abs(tab2D[i, i])
        print(f"\nRow {i+1}: pivot = {pivot}, sum_rest = {sum_rest}")
        if pivot <= sum_rest:
            node = 0
    return node

def jacobi_all_x(A, b, n, iterations):
    x = np.zeros(n)
    all_x = np.zeros((iterations, n),dtype=float)
    for k in range(iterations):
        x_new = np.zeros(n)
        for i in range(n):
            sum_rest = 0
            for j in range(n):
                if j != i:
                    sum_rest += A[i][j] * x[j]
            x_new[i] = (b[i] - sum_rest) / A[i][i]
        x = x_new.copy()
        all_x[k] = x.copy()
    return all_x

n = int(input("Enter n: "))
A = fill_2D_table(n)
B = fill_1D_table(n)
OK = checking(A, n)

if OK == 0:
    print("We can't use the Jacobi method because the matrix is not diagonally dominant.")
else:
    iterations = int(input("Enter number of iterations: "))
    X_tab = jacobi_all_x(A, B, n, iterations)
    print("\nJacobi results for each iteration:\n", X_tab)
