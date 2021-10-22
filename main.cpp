#include <iostream>
#include <fstream>
#include "LBM.hh"

using namespace std;
using namespace LBM;

int main() {
    // CONVERGENCE PARAMETERS
    const double tolerance = 1.0e-12; // tolerance
    double avgU_old = 100; // initialise stored average velocity for comparison

    // INITIALISE ARRAYS
    double xForce[nx][ny]; double yForce[nx][ny]; // x- and y- forcing
    double S[nx][ny][nQ]; // forcing source term
    double velU[nx][ny]; double velV[nx][ny]; // x- and y- velocities
    double rho[nx][ny]; // density
    double fEq[nx][ny][nQ]; double fProp[nx][ny][nQ]; double f[nx][ny][nQ]; // distributions
    initialiseArrays(rho, fEq, velU, velV, fProp, xForce, yForce, f, S); // set all arrays to initial values

    // MAIN LOOP
    for (int t = 1; t < n_t; t++) {
        computeMacros(rho, fEq, velU, velV, fProp, xForce, yForce);
        computeEquilibrium(fEq, rho, velU, velV);
        computeCollisions(f, fEq, fProp, S);
        computePeriodicBC(f);
        computeStreaming(f, fProp);
        computeBounceBackBC(f, fProp);
        
        // check convergence:
        double avgU = computeAvgVelocity(velU);
        if (abs(avgU_old-avgU) < tolerance) {
            cout << "Tolerance criteria met; stopping simulation" << endl;
            break;
        } else {
            avgU_old = avgU;
        }

        // output average velocity and timestep to console:
        if (t%1000 == 0) {
        cout << "t = " << t << ", " << "avgU = " << avgU << endl;
        }
    }

    // save final velocity 
    std::ofstream out("velocityU_field.csv");
    for (auto& row : velU) {
    for (auto col : row)
        out << col <<',';
    out << '\n';
    }

    return 0;
}



