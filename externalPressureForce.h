#ifndef EXTERNALPRESSUREFORCE_M
#define EXTERNALPRESSUREFORCE_M

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalPressureForce
{
public:
	externalPressureForce(elasticPlate &m_plate, timeStepper &m_stepper);
	~externalPressureForce();
	void computeFp();
	void computeJp();

	double pressure;

	VectorXd totalForce;
	
private:
	elasticPlate *plate;
	timeStepper *stepper;

    int ind1, ind2;
    Vector2d p1, p2;

    double l_k, rBar1, rBar2;

    Vector2d tangentVec;
    Vector2d normalVec;
    Vector2d force;

    VectorXd flocal;

    VectorXi localDOF;
};

#endif
