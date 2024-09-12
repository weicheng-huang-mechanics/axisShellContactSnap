#include "externalContactForce.h"

externalContactForce::externalContactForce(elasticPlate &m_plate, timeStepper &m_stepper, 
        double m_stiffness, double m_dBar, double m_ballRadius)
{
	plate = &m_plate;
	stepper = &m_stepper;

    stiffness = m_stiffness;
    dBar = m_dBar;

    ballRadius = m_ballRadius;

    ifContact = VectorXi::Zero(plate->nv - 1);
    ifRebound = 0;

    indexArray = VectorXi::Zero(4);
    jacobian = MatrixXd::Zero(4, 4);

    totalForce = VectorXd::Zero(plate->ndof);

    Id3<<1,0,
         0,1;
}

externalContactForce::~externalContactForce()
{
	;
}

void externalContactForce::computeFc()
{
	Vector2d xEnd = plate->getVertex(plate->nv - 1);

	totalForce = VectorXd::Zero(plate->ndof);

	ifContact = VectorXi::Zero(plate->nv - 1);

	for(int i = 0; i < plate->nv - 1; i++)
	{
		Vector2d xCurrent = plate->getVertex(i);

		double d = (xEnd - xCurrent).norm() - ballRadius;

		//cout << i << " " << d << endl;

		if (d <= dBar)
		{
			ifContact(i) = 1;

			dDdEdge = (xEnd - xCurrent) / (xEnd - xCurrent).norm();

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			f = stiffness * dEnergydD * dDdEdge;

			for (int k = 0; k < 2; k++)
			{
				ind = 2 * i + k;
				stepper->addForce(ind, - f[k]);

				totalForce(ind) = totalForce(ind) + f[k];
			}

			for (int k = 0; k < 2; k++)
			{
				ind = 2 * (plate->nv-1) + k;
				stepper->addForce(ind, f[k]);

				totalForce(ind) = totalForce(ind) - f[k];
			}
		}
	}

	contactNode = plate->nv - 3;

	for(int i = 0; i < plate->nv - 1; i++)
	{
		if (ifContact(i) == 1)
		{
			contactNode = i;
			break;
		} 
	}
}

void externalContactForce::computeJc()
{
	Vector2d xEnd = plate->getVertex(plate->nv - 1);

	for(int i = 0; i < plate->nv - 1; i++)
	{
		indexArray(0) = 2 * (plate->nv - 1) + 0;
		indexArray(1) = 2 * (plate->nv - 1) + 1;
		
		indexArray(2) = 2 * i + 0;
		indexArray(3) = 2 * i + 1;
		
		Vector2d xCurrent = plate->getVertex(i);

		double d = (xEnd - xCurrent).norm() - ballRadius;

		if (d <= dBar)
		{
			dDdEdge = (xEnd - xCurrent) / (xEnd - xCurrent).norm();

			d2DdEdge2 = (Id3 - dDdEdge * dDdEdge.transpose()) / (xEnd - xCurrent).norm();

			dEnergydD = - 2 * (d - dBar) * log(d / dBar) - (d - dBar) * (d - dBar) / d;

			d2EnergydD2 = - 2 * log(d / dBar) - 2 * (d - dBar) / d - 2 * (d - dBar) / d + (d - dBar) * (d - dBar) / (d * d);

			d2EdEdge2 = stiffness * ( d2EnergydD2 * dDdEdge * dDdEdge.transpose() + dEnergydD * d2DdEdge2 );

			jacobian.block(0,0,2,2) =   d2EdEdge2;
			jacobian.block(2,2,2,2) =   d2EdEdge2;
			jacobian.block(2,0,2,2) = - d2EdEdge2;
			jacobian.block(0,2,2,2) = - d2EdEdge2;

			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 4; k++)
				{
					stepper->addJacobian(indexArray(j), indexArray(k), jacobian(k,j));
				}
			}

		}
	}
}
