#include <iostream> 
#include <cmath> 
using namespace std; 
 
#define N 10  // Maximum size of the matrix/vector 
 
// Function to fill the matrix A and vector b 
bool fillMatrix(double A[N][N], double b[N], int &n) { 
    do { 
        cout << "Enter the number of equations (n): "; 
        cin >> n; 
        if (n <= 0 || n > N) { 
            cout << "Invalid n. Please enter a number between 1 and " << N << ".\n"; 
        } 
    } while (n <= 0 || n > N); 
 
    cout << "Enter elements of matrix A row by row:\n"; 
    for (int i = 0; i < n; ++i) { 
        cout << "Row " << i << ": "; 
        for (int j = 0; j < n; ++j) { 
            cin >> A[i][j]; 
        } 
    } 
 
    cout << "Enter elements of vector b:\n"; 
    for (int i = 0; i < n; ++i) { 
        cout << "b[" << i << "]: "; 
        cin >> b[i]; 
    } 
 
    return true; 
} 
 
// Function to check and make matrix diagonally dominant if possible 
bool makeDiagonallyDominant(double A[N][N], double b[N], int n) { 
    for (int i = 0; i < n; ++i) { 
        double sum = 0; 
        for (int j = 0; j < n; ++j) { 
            if (j != i) 
                sum += fabs(A[i][j]); 
        } 
 
        // If the diagonal element is smaller than the sum of others in its row 
        if (fabs(A[i][i]) < sum) { 
            bool swapped = false; 
            for (int k = i + 1; k < n; ++k) { 
                double sum_k = 0; 
                for (int j = 0; j < n; ++j) 
                    if (j != i) 
                        sum_k += fabs(A[k][j]); 
 
                if (fabs(A[k][i]) >= sum_k) { 
                    // Swap rows i and k in A 
                    for (int j = 0; j < n; ++j) { 
                        double temp = A[i][j]; 
                        A[i][j] = A[k][j]; 
                        A[k][j] = temp; 
                    } 
                    // Swap elements in b 
                    double tempb = b[i]; 
                    b[i] = b[k]; 
                    b[k] = tempb; 
 
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
void jacobiSolver(double A[N][N], double b[N], double x[N], int n, double epsilon = 10^-1, int maxIter = 1000) { 
    double x_old[N]; 
    int iter = 0; 
 
    while (iter < maxIter) { 
        iter++; 
 
        // Copy current x into x_old 
        for (int i = 0; i < n; i++) 
            x_old[i] = x[i]; 
 
        // Compute new values using Jacobi formula 
        for (int i = 0; i < n; i++) { 
            double sum = 0; 
            for (int j = 0; j < n; j++) 
                if (j != i) 
                    sum += A[i][j] * x_old[j]; 
            x[i] = (b[i] - sum) / A[i][i]; 
        } 
 
        // Check convergence 
        double maxDiff = 0; 
        for (int i = 0; i < n; i++) { 
            double diff = fabs(x[i] - x_old[i]); 
            if (diff > maxDiff) 
                maxDiff = diff; 
        } 
 
        if (maxDiff <= epsilon) { 
            cout << " Converged after " << iter << " iterations.\n"; 
            return; 
        } 
    } 
 
    cout << " Did not converge after " << maxIter << " iterations.\n"; 
} 
 
// Function to display a vector 
void displayVector(double vec[N], int n) { 
    cout << "\nSolution vector x:\n"; 
    for (int i = 0; i < n; i++) 
        cout << "x[" << i << "] = " << vec[i] << endl; 
} 
 
// ---- Main ---- 
int main() { 
    double A[N][N], b[N], x[N]; 
    int n; 
 
    // Step 1: Fill A and b 
    fillMatrix(A, b, n); 
    cout << "\nMatrix A and vector b have been filled successfully.\n"; 
 
    // Step 2: Make diagonally dominant
if (!makeDiagonallyDominant(A, b, n)) 
        cout << " Matrix is not diagonally dominant — Jacobi method may not converge.\n"; 
    else 
        cout << " Matrix is diagonally dominant or successfully rearranged.\n"; 
// Step 3: Ask user for initial guess 
    cout << "\nEnter initial guess vector (x⁰):\n"; 
    for (int i = 0; i < n; i++) { 
        cout << "x[" << i << "]: "; 
        cin >> x[i]; 
    } 
 
    // Step 4: Solve with Jacobi 
    jacobiSolver(A, b, x, n); 
 
    // Step 5: Display result 
    displayVector(x, n); 
 
    return 0; 
}
