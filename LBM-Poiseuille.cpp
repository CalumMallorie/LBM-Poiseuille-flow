#include <iostream>

using namespace std;

void macroVariables() {
    
}

int main() {
    // INITIALISE SIMULATION PARAMETERS:
    const int nx = 5; int ny = 5; // domain size
    const int n_t = 64000; // timesteps
    const double tau = 1.0; // relaxation time
    const double omega = 1/tau; // collision operator
    const double forceDensity = 1.0e-5; // force density
    const int density = 1; // density

    // LATTICE PARAMETERS:
    int nQ = 9; // number of velocities
    const double w[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights
    const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // CONVERGENCE PARAMETERS
    double tol = 1.0e-12;

    // INITIALISE ARRAYS
    double xForce[nx][ny]; double yForce[nx][ny]; // x- and y- forcing
    double velU[nx][ny]; double velV[nx][ny]; // x- and y- velocities
    int rho[nx][ny]; // density
    double fEq[nx][ny][nQ]; double fProp[nx][ny][nQ]; double f[nx][ny][nQ]; // distributions
    
    for (int x = 0; x < nx; x++) { // loop over arrays to set initial values
        for (int y = 0; y < ny; y++) {
            xForce[x][y] = forceDensity; yForce[x][y] = 0; // force only acts in x-direction
            velU[x][y] = 0; velV[x][y] = 0; // initially at rest
            rho[x][y] = 1; // constant density
            for (int i = 0; i < nQ; i++) {
                // assuming density = 1 and zero velocity
                fEq[x][y][i] = w[i] - 0.5*3*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
                fProp[x][y][i] = w[i] - 0.5*3*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
                f[x][y][i] = w[i] - 0.5*3*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
            }
        }
    }

    // MAIN LOOP
    for (int t = 1; t < n_t; t++) {

        // compute macroscopic variables:
        for (int x = 0; x < nx; x++) { // loop over arrays to set initial values, using eq. 6.2
            for (int y = 0; y < ny; y++) {
                rho[x][y] = fEq[x][y][0] + fEq[x][y][1] + fEq[x][y][2] + fEq[x][y][3] + fEq[x][y][4] + fEq[x][y][5] + fEq[x][y][6] + fEq[x][y][7] + fEq[x][y][8];
                velU[x][y] = fProp[x][y][1] + fProp[x][y][5] + fProp[x][y][8] - fProp[x][y][3] - fProp[x][y][6] - fProp[x][y][7] +  xForce[x][y]/2;
                velV[x][y] = fProp[x][y][2] + fProp[x][y][5] + fProp[x][y][6] - fProp[x][y][4] - fProp[x][y][7] - fProp[x][y][8] +  yForce[x][y]/2;
            }
        }

        // check convergence:
        // compute average x velocity
        double sum = 0.0;
        double size = 0.0;
        double avgU;
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                sum += velU[x][y];
                size++;
            }
        }
        avgU = sum / size;



    }
    cout << "simulation complete!";
    return 0;
}