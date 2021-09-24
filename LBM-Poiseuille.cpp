#include <iostream>
#include <fstream>
using namespace std;



int main() {
    // INITIALISE SIMULATION PARAMETERS:
    const int nx = 5; int ny = 25; // domain size
    const int n_t = 64000; // timesteps
    const double tau = 1.; // relaxation time
    const double omega = 1/tau; // collision operator
    const double forceDensity = 1.e-5; // force density
    const int density = 1; // density

    // LATTICE PARAMETERS:
    int nQ = 9; // number of velocities
    const double w[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights
    const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // CONVERGENCE PARAMETERS
    double tolerance = 1.0e-12; // tolerance
    double avgU_old = 100; // initialise stored average velocity for comparison

    // INITIALISE ARRAYS
    double xForce[nx][ny]; double yForce[nx][ny]; // x- and y- forcing
    double S[nx][ny][nQ]; // forcing source term
    double velU[nx][ny]; double velV[nx][ny]; // x- and y- velocities
    int rho[nx][ny]; // density
    double fEq[nx][ny][nQ]; double fProp[nx][ny][nQ]; double f[nx][ny][nQ]; // distributions
    
    for (int x = 0; x < nx; x++) { // loop over arrays to set initial values
        for (int y = 0; y < ny; y++) {
            xForce[x][y] = forceDensity; yForce[x][y] = 0.; // force only acts in x-direction
            velU[x][y] = 0.; velV[x][y] = 0.; // initially at rest
            rho[x][y] = 1.; // constant density
            for (int i = 0; i < nQ; i++) {
                // assuming density = 1 and zero velocity
                fEq[x][y][i] = w[i] - 0.5*3.*w[i]*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
                fProp[x][y][i] = w[i] - 0.5*3.*w[i]*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
                f[x][y][i] = w[i] - 0.5*3.*w[i]*(cx[i]*xForce[x][y] + cy[i]*yForce[x][y]);
                S[nx][ny][nQ] = 0.; // initialsise with zero source term
            }
        }
    }

    // MAIN LOOP
    for (int t = 1; t < n_t; t++) {

        // compute macroscopic variables:
        for (int x = 0; x < nx; x++) { // loop over arrays to update rho, u and v using eq_6.2
            for (int y = 0; y < ny; y++) {
                rho[x][y] = fEq[x][y][0] + fEq[x][y][1] + fEq[x][y][2] + fEq[x][y][3] + fEq[x][y][4] + fEq[x][y][5] + fEq[x][y][6] + fEq[x][y][7] + fEq[x][y][8];
                velU[x][y] = fProp[x][y][1] + fProp[x][y][5] + fProp[x][y][8] - fProp[x][y][3] - fProp[x][y][6] - fProp[x][y][7] +  xForce[x][y]/2;
                velV[x][y] = fProp[x][y][2] + fProp[x][y][5] + fProp[x][y][6] - fProp[x][y][4] - fProp[x][y][7] - fProp[x][y][8] +  yForce[x][y]/2;
            }
        }

        // check convergence:
        double sum = 0.0;
        double size = 0.0;
        double avgU; // average x velocity
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                sum += velU[x][y];
                size++;
            }
        }
        avgU = sum / size;
        if (abs(avgU-avgU_old) < tolerance) {
            cout << "Tolerance criteria reached." << endl;
            break; // stop main loop if tolerance is met
            
        } else {
            avgU_old = avgU; // update average u at t = t-dt
        }
        
        // compute equilibrium distributions using eq_4.38
        for (int x = 0; x < nx; x++) {
            for (int y = 0; y < ny; y++) {
                fEq[x][y][0] = w[0]*(rho[x][y] + 3*(velU[x][y]*cx[0] + velV[x][y]*cy[0]));
                fEq[x][y][1] = w[1]*(rho[x][y] + 3*(velU[x][y]*cx[1] + velV[x][y]*cy[1]));
                fEq[x][y][2] = w[2]*(rho[x][y] + 3*(velU[x][y]*cx[2] + velV[x][y]*cy[2]));
                fEq[x][y][3] = w[3]*(rho[x][y] + 3*(velU[x][y]*cx[3] + velV[x][y]*cy[3]));
                fEq[x][y][4] = w[4]*(rho[x][y] + 3*(velU[x][y]*cx[4] + velV[x][y]*cy[4]));
                fEq[x][y][5] = w[5]*(rho[x][y] + 3*(velU[x][y]*cx[5] + velV[x][y]*cy[5]));
                fEq[x][y][6] = w[6]*(rho[x][y] + 3*(velU[x][y]*cx[6] + velV[x][y]*cy[6]));
                fEq[x][y][7] = w[7]*(rho[x][y] + 3*(velU[x][y]*cx[7] + velV[x][y]*cy[7]));
                fEq[x][y][8] = w[8]*(rho[x][y] + 3*(velU[x][y]*cx[8] + velV[x][y]*cy[8]));
            }
        }

        // compute forcing source term
        for (int x = 0; x < nx; x++) { // loop over arrays to update rho, u and v using eq_6.2
            for (int y = 0; y < ny; y++) {
                S[x][y][0] = (1. - omega/2)*(w[0]*3*(cx[0]*xForce[x][y] + cy[0]*yForce[x][y]));
                S[x][y][1] = (1. - omega/2)*(w[1]*3*(cx[1]*xForce[x][y] + cy[1]*yForce[x][y]));
                S[x][y][2] = (1. - omega/2)*(w[2]*3*(cx[2]*xForce[x][y] + cy[2]*yForce[x][y]));
                S[x][y][3] = (1. - omega/2)*(w[3]*3*(cx[3]*xForce[x][y] + cy[3]*yForce[x][y]));
                S[x][y][4] = (1. - omega/2)*(w[4]*3*(cx[4]*xForce[x][y] + cy[4]*yForce[x][y]));
                S[x][y][5] = (1. - omega/2)*(w[5]*3*(cx[5]*xForce[x][y] + cy[5]*yForce[x][y]));
                S[x][y][6] = (1. - omega/2)*(w[6]*3*(cx[6]*xForce[x][y] + cy[6]*yForce[x][y]));
                S[x][y][7] = (1. - omega/2)*(w[7]*3*(cx[7]*xForce[x][y] + cy[7]*yForce[x][y]));
                S[x][y][8] = (1. - omega/2)*(w[8]*3*(cx[8]*xForce[x][y] + cy[8]*yForce[x][y]));
            }
        }

        // compute collisions 
        for (int x = 0; x < nx; x++) { // loop over arrays to compute post-collision distributions, using eq_3.9 + S[x][y][n]
            for (int y = 0; y < ny; y++) {
                f[x][y][0] = (1. - omega)*fProp[x][y][0] + omega*fEq[x][y][0] + S[x][y][0];
                f[x][y][1] = (1. - omega)*fProp[x][y][1] + omega*fEq[x][y][1] + S[x][y][1];
                f[x][y][2] = (1. - omega)*fProp[x][y][2] + omega*fEq[x][y][2] + S[x][y][2];
                f[x][y][3] = (1. - omega)*fProp[x][y][3] + omega*fEq[x][y][3] + S[x][y][3];
                f[x][y][4] = (1. - omega)*fProp[x][y][4] + omega*fEq[x][y][4] + S[x][y][4];
                f[x][y][5] = (1. - omega)*fProp[x][y][5] + omega*fEq[x][y][5] + S[x][y][5];
                f[x][y][6] = (1. - omega)*fProp[x][y][6] + omega*fEq[x][y][6] + S[x][y][6];
                f[x][y][7] = (1. - omega)*fProp[x][y][7] + omega*fEq[x][y][7] + S[x][y][7];
                f[x][y][8] = (1. - omega)*fProp[x][y][8] + omega*fEq[x][y][8] + S[x][y][8];
            }
        }

        // apply periodic boundary conditions using virtual nodes at x = 0 and x = nx, using eq_5.15
        for (int y = 0; y < ny; y++) {
            // set f(near side) =  f(far side)
            f[0][y][1] = f[nx-2][y][1];
            f[0][y][5] = f[nx-2][y][5];
            f[0][y][8] = f[nx-2][y][8];
            // set f(far side) = f(near side)
            f[nx-2][y][3] = f[0][y][3];
            f[nx-2][y][6] = f[0][y][6];
            f[nx-2][y][7] = f[0][y][7];
        }

        // compute streaming step
        for (int x = 0; x < nx; x++) { 
            for (int y = 1; y < ny-1; y++) {
                fProp[x][y][0]                   = f[x][y][0];
                fProp[(x + 1) % nx][y][1]        = f[x][y][1];
                fProp[x][y+1][2]                 = f[x][y][2];
                fProp[(x - 1 + nx) % nx][y][3]   = f[x][y][3];             
                fProp[x][y-1][4]                 = f[x][y][4];
                fProp[(x + 1) % nx][y+1][5]      = f[x][y][5];
                fProp[(x + 1) % nx][y+1][6]      = f[x][y][6];
                fProp[(x - 1 + nx) % nx][y-1][7] = f[x][y][7];
                fProp[(x + 1) % nx][y-1][8]      = f[x][y][8]; 
                
            }
        }

        // apply bounce-back boundary conditions at top and bottom walls
        for (int x = 0; x < nx; x++) {
            // top wall
            fProp[x][ny-1][4] = f[x][ny-1][2];
            fProp[x][ny-1][7] = f[x][ny-1][5];
            fProp[x][ny-1][8] = f[x][ny-1][6];
            // bottom wall
            fProp[x][0][2] = f[x][0][4];
            fProp[x][0][5] = f[x][0][7];
            fProp[x][0][6] = f[x][0][8];
        }
        cout << avgU << endl;
        

    }


 std::ofstream out("test.csv");

for (auto& row : velU) {
  for (auto col : row)
    out << col <<',';
  out << '\n';

}
   
    
    


    return 0;
}