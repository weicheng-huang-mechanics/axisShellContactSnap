#include "externalPressureForce.h"

externalPressureForce::externalPressureForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
	stepper = &m_stepper;

	pressure = 0.0;

	localDOF = VectorXi::Zero(4);
	flocal = VectorXd::Zero(4);
}

externalPressureForce::~externalPressureForce()
{
	;
}

void externalPressureForce::computeFp()
{
	//pressure = 100.0;

	totalForce = VectorXd::Zero(plate->ndof);

	for (int k = 0; k < plate->edgeNum; k++)
	{
		ind1 = plate->v_edgeElement[k].nv_1;
		ind2 = plate->v_edgeElement[k].nv_2;

		//p1 = plate->getVertexOld(ind1);
		//p2 = plate->getVertexOld(ind2);

		p1 = plate->getVertex(ind1);
		p2 = plate->getVertex(ind2);

		l_k = plate->v_edgeElement[k].refLength;

		rBar1 = plate->v_edgeElement[k].x_1_start(0);
		rBar2 = plate->v_edgeElement[k].x_2_start(0);

		tangentVec = (p2 - p1) / (p2 - p1).norm();

		normalVec(0) =   tangentVec(1);
		normalVec(1) = - tangentVec(0);

		force = - normalVec * pressure * l_k * M_PI * (rBar1 + rBar2) / 2;

		flocal(0) = force(0);
		flocal(1) = force(1);
		flocal(2) = force(0);
		flocal(3) = force(1);

		//cout << flocal.norm() << endl;

		localDOF(0) = 2 * ind1 + 0;
		localDOF(1) = 2 * ind1 + 1;
		localDOF(2) = 2 * ind2 + 0;
		localDOF(3) = 2 * ind2 + 1;

		for (int i = 0; i < 4; i++)
		{
			stepper->addForce(localDOF(i), - flocal(i));

			totalForce(localDOF(i)) = totalForce(localDOF(i)) + flocal(i);
		}
	}	
}

void externalPressureForce::computeJp()
{
	;
}
