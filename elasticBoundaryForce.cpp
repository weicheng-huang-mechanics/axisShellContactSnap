#include "elasticBoundaryForce.h"
#include <iostream>

elasticBoundaryForce::elasticBoundaryForce(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;
	
	Jss.setZero(4, 4);
	flocal = VectorXd::Zero(4);
	EA = plate->EA;
	mu = plate->Possion;

	localDOF = VectorXi::Zero(4);
}

elasticBoundaryForce::~elasticBoundaryForce()
{
	;
}

void elasticBoundaryForce::computeFboundarty()
{
	totalForce = VectorXd::Zero(plate->ndof);

  ind1 = plate->v_edgeElement[plate->edgeNum - 1].nv_1;
  ind2 = plate->v_edgeElement[plate->edgeNum - 1].nv_2;

  phiBar = plate->v_edgeElement[plate->edgeNum - 1].phiBar;

  p = plate->getVertex(ind1);
  p1 = plate->getVertex(ind2);

  xk = p[0];
  yk = p[1];
  xkp1 = p1[0];
  ykp1 = p1[1];

  l_k = plate->v_edgeElement[plate->edgeNum - 1].refLength;

  flocal = computeBoundaryForce(xk, yk, xkp1, ykp1, phiBar);

  //cout << (xkp1 - xk) / (p1 - p).norm() << " " << tBarx << " ";
  //cout << (ykp1 - yk) / (p1 - p).norm() << " " << tBary << " ";

  flocal = - flocal * plate->EA * 1.0;

  localDOF = VectorXi::Zero(4);

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

void elasticBoundaryForce::computeJboundarty()
{
  ind1 = plate->v_edgeElement[plate->edgeNum - 1].nv_1;
  ind2 = plate->v_edgeElement[plate->edgeNum - 1].nv_2;

  phiBar = plate->v_edgeElement[plate->edgeNum - 1].phiBar;

  p = plate->getVertex(ind1);
  p1 = plate->getVertex(ind2);

  xk = p[0];
  yk = p[1];
  xkp1 = p1[0];
  ykp1 = p1[1];

  l_k = plate->v_edgeElement[plate->edgeNum - 1].refLength;

  localDOF = VectorXi::Zero(4);

  localDOF(0) = 2 * ind1 + 0;
  localDOF(1) = 2 * ind1 + 1;
  localDOF(2) = 2 * ind2 + 0;
  localDOF(3) = 2 * ind2 + 1;

  Jss = computeBoundaryJacobian(xk, yk, xkp1, ykp1, phiBar);
  Jss = Jss * plate->EA * 1.0;

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      stepper->addJacobian(localDOF(i), localDOF(j), Jss(i,j));
    }
  }
}

VectorXd elasticBoundaryForce::computeBoundaryForce(double xa, double ya, double xb, double yb, double phiBar)
{
  VectorXd vecResult;

  vecResult = ListVec((2*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
    (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
   (-2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
    ((-xa + xb)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
   (-2*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
    (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
   (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
    ((-xa + xb)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))));

  return vecResult;
}


MatrixXd elasticBoundaryForce::computeBoundaryJacobian(double xa, double ya, double xb, double yb, double phiBar)
{
  MatrixXd matResult;

  matResult = ListMat(ListVec((2*pow(-ya + yb,2))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,3)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,5)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (-2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (-2*pow(-ya + yb,2))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,3)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,5)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2)))),
   ListVec((-2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    2/(pow(-xa + xb,2)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)),
    (2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    -2/(pow(-xa + xb,2)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2))),
   ListVec((-2*pow(-ya + yb,2))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,3)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,5)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (2*pow(-ya + yb,2))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,3)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,5)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    (-2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2)))),
   ListVec((2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    -2/(pow(-xa + xb,2)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)),
    (-2*(-ya + yb))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) + 
     (4*pow(-ya + yb,2)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,4)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (2*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,2)*(1 + pow(-ya + yb,2)/pow(-xa + xb,2))),
    2/(pow(-xa + xb,2)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2)) - 
     (4*(-ya + yb)*(-phiBar + atan((-ya + yb)/(-xa + xb))))/
      (pow(-xa + xb,3)*pow(1 + pow(-ya + yb,2)/pow(-xa + xb,2),2))));

  return matResult;
}

VectorXd elasticBoundaryForce::ListVec(double a1, double a2, double a3, double a4)
{
  VectorXd vecResult;

  vecResult.setZero(4, 1);

  vecResult(0) = a1;
  vecResult(1) = a2;
  vecResult(2) = a3;
  vecResult(3) = a4;

  return vecResult;
}

MatrixXd elasticBoundaryForce::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, VectorXd a4)
{
  MatrixXd matResult;

  matResult.setZero(4, 4);

  matResult.col(0) = a1;
  matResult.col(1) = a2;
  matResult.col(2) = a3;
  matResult.col(3) = a4;

  return matResult;
}