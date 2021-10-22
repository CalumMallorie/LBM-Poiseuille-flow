namespace LBM
{
    // LBM INPUT PARAMETERS:
    const int n_t = 64000; // max timesteps
    const double tau = 1.; // relaxation time
    const double omega = 1/tau; // collision operator
    const double forceDensity = 1.e-5; // force density
    const int density = 1; // density

    // GLOBAL LATTICE PARAMETERS:
    const int nx = 50; const int ny = 50; // domain size
    const int nQ = 9; // number of velocities
    const double w[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.}; // lattice weights
    const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // FUNCTION DECLARATIONS:
    void computeMacros(double rho[nx][ny], double fEq[nx][ny][nQ], double velU[nx][ny], double velV[nx][ny], double fProp[nx][ny][nQ], double xForce[nx][ny], double yForce[nx][ny]);
    void initialiseArrays(double rho[nx][ny], double fEq[nx][ny][nQ], double velU[nx][ny], double velV[nx][ny],double fProp[nx][ny][nQ], double xForce[nx][ny], double yForce[nx][ny], double f[nx][ny][nQ], double S[nx][ny][nQ]);
    void computeEquilibrium(double fEq[nx][ny][nQ], double rho[nx][ny], double velU[nx][ny], double velV[nx][ny]);
    void computeForcing(double S[nx][ny][nQ], double xForce[nx][ny], double yForce[nx][ny]);
    void computeCollisions(double f[nx][ny][nQ], double fEq[nx][ny][nQ], double fProp[nx][ny][nQ], double S[nx][ny][nQ]);
    void computePeriodicBC(double f[nx][ny][nQ]);
    void computeStreaming(double f[nx][ny][nQ], double fProp[nx][ny][nQ]);
    void computeBounceBackBC(double f[nx][ny][nQ], double fProp[nx][ny][nQ]);
    double computeAvgVelocity(double velU[nx][ny]);

}