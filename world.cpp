#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");		

	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");

	thickness = m_inputData.GetScalarOpt("thickness");
	Possion = m_inputData.GetScalarOpt("Possion");

	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
	nv = m_inputData.GetIntOpt("nv");

	totalCompress = m_inputData.GetScalarOpt("totalCompress");
	inputPressure = m_inputData.GetScalarOpt("inputPressure");
	delta = m_inputData.GetScalarOpt("delta");

	height = m_inputData.GetScalarOpt("height");
	stiffness = m_inputData.GetScalarOpt("stiffness");
	dBar = m_inputData.GetScalarOpt("dBar");

	ballRadius = m_inputData.GetScalarOpt("ballRadius");
	massSpeed = m_inputData.GetScalarOpt("massSpeed");

	currentCompress = 0.0;
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER";
    //name << "_nv_" << nv;
    name << "_totalCompress_" << totalCompress;
    name << "_ballRadius_" << ballRadius;
    name << "_delta_" << delta;
    name << "_thickness_" << thickness;
    //name << "_dBar_" << dBar;
    //name << "_stiffness_" << stiffness / 1e6;
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	//if (timeStep == Nstep)
	{
		//cout << nv << endl;

		computeReactionForce();

		Vector2d xBall = plate->getVertex(plate->nv - 1);
		Vector2d xEnd = plate->getVertex(plate->nv - 2);

		Vector2d x1 = plate->getVertex(nv - 1);
		Vector2d x2 = plate->getVertex(nv - 2);
		Vector2d x3 = plate->getVertex(nv - 3);

		Vector2d t1 = (x1 - x2) / (x1 - x2).norm();
		Vector2d t2 = (x2 - x3) / (x2 - x3).norm();

		double deltaL = ( (x2 - x3).norm() + (x1 - x2).norm() ) / 2;

		double curvatureLocal1 = 2 * ( t1(0) * t2(1) - t1(1) * t2(0) ) / ( (1 + t1.dot(t2) ) * deltaL );

		//double deltaL = plate->v_bendingElement[2].voroniLength;

		//cout << plate->v_bendingElement[2].voroniLength << " " << ( (x2 - x3).norm() + (x1 - x2).norm() ) / 2 << endl;

		//outfile << currentCompress << " " << curvatureLocal << " " << ( xBall(1) - xEnd(1) ) - ballRadius << " " << reactionForce(plate->ndof - 1) << endl;

		Vector2d x11 = plate->getVertex(nv);
		Vector2d x22 = plate->getVertex(nv - 1);
		Vector2d x33 = plate->getVertex(nv - 2);

		x11(0) = 0.0;

		Vector2d t11 = (x11 - x22) / (x11 - x22).norm();
		Vector2d t22 = (x22 - x33) / (x22 - x33).norm();

		double deltaLL = ( (x22 - x33).norm() + (x11 - x22).norm() ) / 2;

		double curvatureLocal2 = 2 * ( t11(0) * t22(1) - t11(1) * t22(0) ) / ( (1 + t11.dot(t22) ) * deltaLL );


		//outfile << currentCompress << " " << (curvatureLocal1 + curvatureLocal2) / 2  << " " << ( xBall(1) - xEnd(1) ) - ballRadius << " " << reactionForce(plate->ndof - 1) << endl;


		/*

		double curvatureLocal = 0.0;
		int temp = 0;

		//for (int i = m_externalContactForce->contactNode; i < plate->nv - 2; i++)
		{
			Vector2d x1 = plate->getVertex(m_externalContactForce->contactNode);
			Vector2d x2 = plate->getVertex(plate->nv - 2);
			Vector2d x3 = x1;
			x3(0) = - x1(0);
			x2(0) = 0.0;

			Vector2d t1 = (x1 - x2) / (x1 - x2).norm();
			Vector2d t2 = (x2 - x3) / (x2 - x3).norm();

			double deltaL = ( (x2 - x3).norm() + (x1 - x2).norm() ) / 2;

			double curvatureLocal1 = 2 * ( t1(0) * t2(1) - t1(1) * t2(0) ) / ( (1 + t1.dot(t2) ) * deltaL );

			curvatureLocal = curvatureLocal + curvatureLocal1;

			temp = temp + 1;
		}

		//cout << temp << endl;

		curvatureLocal = curvatureLocal / temp;

		outfile << currentCompress << " " << curvatureLocal  << " " << ( xBall(1) - xEnd(1) ) - ballRadius << " " << reactionForce(plate->ndof - 1) << endl;

		*/

		//outfile << currentCompress << " " << ( xBall(1) - xEnd(1) ) - ballRadius << " " << reactionForce(plate->ndof - 1) << endl;

		//outfile << currentCompress << " " << reactionForce(plate->ndof - 1) << endl;

		//outfile << currentTime << " " << xEnd(1) << endl;

		//cout << plate->nv << endl;

		/*

		if ( abs(reactionForce(plate->ndof - 1)) < 1.0 )
		{
			for (int i = 0; i < nv; i++)
			{
				Vector2d xCurrent = plate->getVertex(i);
				outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << endl;
			}
		}

		*/

		for (int i = 0; i < plate->nv; i++)
		{
			Vector2d xCurrent = plate->getVertex(i);
			outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << endl;
		}

	
		for (int i = 0; i < nv; i++)
		{
			/*
			Vector2d x1 = plate->getVertex(i);
			Vector2d x2 = plate->getVertex(i + 1);

			Vector2d tangent = (x2 - x1) / (x2 - x1).norm();
			Vector2d surfaceNormal;
			surfaceNormal(0) =   tangent(1);
			surfaceNormal(1) = - tangent(0);

			Vector2d xMid = (x2 + x1) / 2;

			double lineA;
			double lineB;

			if ( abs(surfaceNormal(0)) > 1e-5)
			{
				lineA = surfaceNormal(1) / surfaceNormal(0);
				lineB = xMid(1) - lineA * xMid(0);
			}
			else
			{
				lineA = surfaceNormal(1) / (surfaceNormal(0)+1e-5);
				lineB = xMid(1) - lineA * xMid(0);
			}

			Vector2d yPos;
			yPos(0) = 0.0;
			yPos(1) = lineB;

			Vector2d xEnd = plate->getVertex(nv - 1);

			double distance1 = (xMid - yPos).norm();
			double distance2 = (xEnd - yPos).norm();

			outfile << currentTime << " " << i << " " << distance1 << " " << distance2 << endl;

			*/

			/*
			Vector2d x1 = plate->getVertex(i);
			Vector2d x2 = plate->getVertex(i + 1);
			Vector2d x3 = plate->getVertex(i + 2);

			x2 = plate->getVertex(nv);
			x2(0) = 0.0;

			x3 = x1;
			x3(0) = - x1(0);

			Vector2d t1 = (x1 - x2) / (x1 - x2).norm();
			Vector2d t2 = (x2 - x3) / (x2 - x3).norm();

			double deltaL = ( (x2 - x3).norm() + (x1 - x2).norm() ) / 2;

			double curvatureLocal = 2 * ( t1(0) * t2(1) - t1(1) * t2(0) ) / ( (1 + t1.dot(t2) ) * deltaL );

			outfile << currentTime << " " << i << " " << curvatureLocal << endl;

			*/

			//Vector2d xCurrent = plate->getVertex(i);
			//outfile << currentTime << " " << xCurrent(0) << " " << xCurrent(1) << endl;
		}

	}
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = new elasticPlate(YoungM, density, thickness, Possion, 
		deltaTime, nv, delta, height, massSpeed);

	plateBoundaryCondition();

	plate->setup();

	stepper = new timeStepper(*plate);

	// set up force
	m_inertialForce = new inertialForce(*plate, *stepper);
	m_gravityForce = new externalGravityForce(*plate, *stepper, gVector);
	m_dampingForce = new dampingForce(*plate, *stepper, viscosity);
	m_elasticStretchingForce = new elasticStretchingForce(*plate, *stepper);
	m_elasticBendingForce = new elasticBendingForce(*plate, *stepper);
	//m_externalPressureForce = new externalPressureForce(*plate, *stepper);
	m_elasticBoundaryForce = new elasticBoundaryForce(*plate, *stepper);
	m_externalContactForce = new externalContactForce(*plate, *stepper, stiffness, dBar, ballRadius);

	plate->updateTimeStep();

	// set up first jacobian
	m_inertialForce->setFirstJacobian();
	m_dampingForce->setFirstJacobian();
	m_elasticStretchingForce->setFirstJacobian();
	m_elasticBendingForce->setFirstJacobian();

	stepper->first_time_PARDISO_setup();

	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	Vector2d xStart = plate->getVertex(0);
	plate->setOneBoundaryCondition(xStart(1), 0, 1);

	Vector2d xEnd1 = plate->getVertex(plate->nv - 2);
	plate->setOneBoundaryCondition(xEnd1(0), plate->nv - 2, 0);

	Vector2d xEnd2 = plate->getVertex(plate->nv - 3);
	plate->setOneBoundaryCondition(xEnd2(0), plate->nv - 3, 0);


	Vector2d xBall = plate->getVertex(plate->nv - 1);
	plate->setOneBoundaryCondition(xBall(0), plate->nv - 1, 0);
	plate->setOneBoundaryCondition(xBall(1), plate->nv - 1, 1);
}

void world::updateTimeStep()
{

	if (currentCompress < totalCompress)
	{
		Vector2d xBall = plate->getVertex(plate->nv - 1);
		xBall(1) = xBall(1) - deltaTime * massSpeed;
		plate->setOneBoundaryCondition(xBall(0), plate->nv - 1, 0);
		plate->setOneBoundaryCondition(xBall(1), plate->nv - 1, 1);

		currentCompress = currentCompress + deltaTime * massSpeed;
	}

	bool goodSolved = false;

	while (goodSolved == false)
	{
		// Start with a trial solution for our solution x
		plate->updateGuess(); // x = x0 + u * dt

		updateEachStep();

		goodSolved = true;
	}

	plate->updateTimeStep();

	if (render) 
	{
		cout << "time: " << currentTime << endl;
	}

	currentTime += deltaTime;
		
	timeStep++;
}

void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;

	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		m_inertialForce->computeFi();
		m_dampingForce->computeFd();

		//m_gravityForce->computeFg();
		m_externalContactForce->computeFc();

		m_elasticStretchingForce->computeFs();
		m_elasticBendingForce->computeFb();
		//m_externalPressureForce->computeFp();
		m_elasticBoundaryForce->computeFboundarty();
	
		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		//cout << normf << " ";

		normf = 0.0;
		
		if (solved == false)
		{
			m_inertialForce->computeJi();
			//m_gravityForce->computeJg();
			m_dampingForce->computeJd();

			m_elasticStretchingForce->computeJs();
			m_elasticBendingForce->computeJb();
			//m_externalPressureForce->computeJp();
			m_elasticBoundaryForce->computeJboundarty();
			m_externalContactForce->computeJc();
			
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			timeStep = Nstep;
			break;
		}
	}

	//cout << endl;

	if (render)
	{
		cout << "iter " << iter << endl;
	}

	for (int i = 0; i < plate->edgeNum; i++)
	{
		double e111 = plate->v_edgeElement[i].edgeLength;
		double e222 = plate->v_edgeElement[i].refLength;

		if ( abs(e111 - e222) / abs(e222) > 0.5 )
		{
			timeStep = Nstep;
			break;
		}
	}
}

int world::simulationRunning()
{

	if (timeStep < Nstep)
	{
		return 1;
	}
	else
	{
		return - 1;
	}
}

Vector2d world::getScaledCoordinate(int i, int j)
{
	Vector2d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

Vector2d world::getScaledPos(int i)
{
	Vector2d xCurrent;

	xCurrent = plate->getVertex(i) * scaleRendering;

	return xCurrent;
}

int world::getNv()
{
	return plate->nv;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}

void world::computeReactionForce()
{
	reactionForce = VectorXd::Zero(plate->ndof);

	m_elasticStretchingForce->computeFs();
	m_elasticBendingForce->computeFb();
	//m_externalPressureForce->computeFp();
	m_elasticBoundaryForce->computeFboundarty();
	m_externalContactForce->computeFc();

	reactionForce = - m_externalContactForce->totalForce - m_elasticBoundaryForce->totalForce - m_elasticStretchingForce->totalForce - m_elasticBendingForce->totalForce;
}

double world::getBallRadius()
{
	return ballRadius * scaleRendering;
}

int world::getIfContact(int i)
{
	return m_externalContactForce->ifContact(i);
}

int world::getVectex()
{
	return plate->nv - 1;
}