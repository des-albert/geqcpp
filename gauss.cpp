#include "geq.h"
#include <algorithm>

using namespace std;

double ** array2d(int, int);
void array2del(double **, int);

void gauss(double *b, double **a, int n) {

    double **aug = array2d(n, n + 1);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }

    // Pivoting

    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            if (abs(aug[i][i]) < abs(aug[k][i])) {
                for (int j = 0; j <= n; j++) {
                    double temp = aug[i][j];
                    aug[i][j] = aug[k][j];
                    aug[k][j] = temp;
                }
            }
        }
    }

    //  Gauss elimination

    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {
            double temp = aug[k][i] / aug[i][i];
            for (int j = 0; j <=n; j++) {
                aug[k][j] = aug[k][j] - temp * aug[i][j];
            }
        }
    }

    // Back substitution

    for (int i = n - 1; i >= 0; i--) {
        b[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            if (j != i) {
                b[i] = b[i] - aug[i][j] * b[j];
            }
            b[i] = b[i]/aug[i][i];
        }
    }

    array2del(aug, n);
}

