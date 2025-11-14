#include <iostream>
#include <cmath>
using namespace std;

const int N= 10;  // Maximum size of the matrix/vector
const double EPS = 0.1; // Convergence tolerance (epsilon)

// Function to fill matrix A and vector b
bool fillMatrix(double A[N][N], double b[N], int &n) {
    // Ask user for number of equations
    do {
        cout << "Enter the number of equations (n): ";
        cin >> n;
        if (n <= 0 || n > N)
            cout << "Invalid n. Please enter a number between 1 and " << N << ".\n";
    } while (n <= 0 || n > N);

    // Fill matrix A
    cout << "Enter elements of matrix A row by row:\n";
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cin >> A[i][j];

    // Fill vector b
    cout << "Enter elements of vector b:\n";
    for (int i = 0; i < n; ++i)
        cin >> b[i];

    return true;
}

//  Make matrix diagonally dominant
bool makeDiagonallyDominant(double A[N][N], double b[N], int n) {
    for (int i = 0; i < n; ++i) {
        // Sum of non-diagonal elements in row i
        double sum = 0;
        for (int j = 0; j < n; ++j)
            if (j != i) sum += fabs(A[i][j]);

        // If diagonal element is smaller than sum, try to swap with lower rows
        if (fabs(A[i][i]) < sum) {
            bool swapped = false;

            for (int k = i + 1; k < n; ++k) {
                double sum_k = 0;
                for (int j = 0; j < n; ++j)
                    if (j != i) sum_k += fabs(A[k][j]);

                // Check if swapping row k with row i makes it diagonally dominant
                if (fabs(A[k][i]) >= sum_k) {
                    for (int j = 0; j < n; ++j)
                        swap(A[i][j], A[k][j]); // Swap rows in A
                    swap(b[i], b[k]); // Swap corresponding elements in b
                    swapped = true;
                    break;
                }
            }

            if (!swapped) {
                cout << "Cannot make row " << i << " diagonally dominant.\n";
                return false;
            }
        }
    }
    return true;
}

// Jacobi solver
void jacobiSolver(double A[N][N], double b[N], double x[N], int n, int maxIter) {
    double x_old[N];

    cout << "\nStarting Jacobi Iterations \n";

    for (int iter = 1; iter <= maxIter; iter++) {

        // Store current x in x_old
        for (int i = 0; i < n; i++)
            x_old[i] = x[i];

        // Compute new values using Jacobi formula
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++)
                if (j != i)
                    sum += A[i][j] * x_old[j]; // sum(a_ij * x_j^(k))
            x[i] = (b[i] - sum) / A[i][i]; // x_i^(k+1) = (b_i - sum)/a_ii
        }

        // Print current iteration
        cout << "Iteration " << iter << ": ";
        for (int i = 0; i < n; i++)
            cout << x[i] << "  ";
        cout << endl;

        // Check convergence using epsilon
        double maxDiff = 0;
        for (int i = 0; i < n; i++)
            maxDiff = max(maxDiff, fabs(x[i] - x_old[i]));

        if (maxDiff < EPS) {
            cout << "\nConverged at iteration " << iter << " (max difference < " << EPS << ")\n";
            break;
        }
    }
}

// Display vector
void displayVector(double vec[N], int n) {
    cout << "\nFinal solution vector x:\n";
    for (int i = 0; i < n; i++)
        cout << "x[" << i << "] = " << vec[i] << endl;
}

// Main function
int main() {
    double A[N][N], b[N], x[N];
    int n, iterations;

    // Step 1: Fill matrix and vector
    fillMatrix(A, b, n);
    cout << "\nMatrix A and vector b have been filled successfully.\n";

    // Step 2: Make diagonally dominant if possible
    if (!makeDiagonallyDominant(A, b, n))
        cout << " Matrix is not diagonally dominant — Jacobi method may not converge.\n";
else
        cout << "Matrix is diagonally dominant or successfully rearranged.\n";

    // Step 3: Ask user for initial guess vector
    cout << "\nEnter initial guess vector (x⁰):\n";
    for (int i = 0; i < n; i++)
        cin >> x[i];

    // Step 4: Ask user for maximum number of iterations
    cout << "\nEnter maximum number of iterations: ";
    cin >> iterations;

    // Step 5: Solve using Jacobi
    jacobiSolver(A, b, x, n, iterations);

    // Step 6: Display final solution
    displayVector(x, n);

    return 0;
}


