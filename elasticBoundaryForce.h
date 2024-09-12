#ifndef ELASTICBOUNDARYFORCE_H
#define ELASTICBOUNDARYFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticBoundaryForce
{
public:
	elasticBoundaryForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticBoundaryForce();
	void computeFboundarty();
	void computeJboundarty();

    VectorXd totalForce;
    
private:
	elasticPlate *plate;
    timeStepper *stepper;
	
	double EA;
    double mu;
    int ndof;
    
    double xk, yk, xkp1, ykp1, l_k;
    double phiBar;

    VectorXd flocal;
    MatrixXd Jss;

    int ind1, ind2;

    Vector2d p, p1;

    VectorXi localDOF;

    VectorXd computeBoundaryForce(double xa, double ya, double xb, double yb, double phiBar);
    MatrixXd computeBoundaryJacobian(double xa, double ya, double xb, double yb, double phiBar);
    VectorXd ListVec(double a1, double a2, double a3, double a4);
    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4);
};

#endif
