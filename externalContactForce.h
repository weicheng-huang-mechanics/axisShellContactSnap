#ifndef EXTERNALCONATCTFORCE_H
#define EXTERNALCONATCTFORCE_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class externalContactForce
{
public:
	externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar, double m_ballRadius);
	~externalContactForce();

	void computeFc();
	void computeJc();

    VectorXi ifContact;
    int ifRebound;

    int contactNode;

    VectorXd totalForce;

private:
	elasticPlate *plate;
    timeStepper *stepper;

    double ballRadius;

    double dEnergydD;
    double d2EnergydD2;

    double stiffness;
    double dBar;

    VectorXi indexArray;
    MatrixXd jacobian;

    Vector2d dDdEdge;
    Matrix2d Id3;

    Matrix2d d2DdEdge2;
    Matrix2d d2EdEdge2;

    Vector2d f;

    int ind;
};

#endif
