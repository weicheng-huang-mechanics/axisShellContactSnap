#include "elasticPlate.h"

elasticPlate::elasticPlate(double m_YoungM, double m_density, double m_thickness, 
		double m_Possion, double m_dt, int m_nv, double m_delta, double m_height, double m_massSpeed)
{
	YoungM = m_YoungM;
	density = m_density;
	thickness = m_thickness;
	Possion = m_Possion;
	dt = m_dt;
	nv = m_nv;
	delta = m_delta;
	height = m_height;

	beta0 = 8.83 / 180 * M_PI;

	massSpeed = m_massSpeed;

	EA = YoungM * thickness / (1 - Possion * Possion);
	EI = YoungM * thickness * thickness * thickness / ( 12 * (1 - Possion * Possion) );

	setupGeometry();

	ndof = 2 * nv;
	x = VectorXd::Zero(ndof);
	x0 = VectorXd::Zero(ndof);
	u = VectorXd::Zero(ndof);

	u(ndof - 1) = - massSpeed;

	for (int i = 0; i < nv; i++)
	{
		x(2 * i + 0) = v_nodes[i](0);
		x(2 * i + 1) = v_nodes[i](1);
	}
	x0 = x;
	x00 = x;

	x = x0 + u * dt;

	//cout << x << endl;

	refHeight = x(ndof - 1) - x(1);

	computeEdge();
	computeBending();

	setupMass();

	//set up constraint map
	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

elasticPlate::~elasticPlate()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;
}

void elasticPlate::setup()
{
	ncons = 0;
    for (int i=0; i < ndof; i++)
    {
		if (isConstrained[i] > 0)
		{
			ncons++;
		}
	}
	uncons = ndof - ncons;

	unconstrainedMap = new int[uncons]; // maps xUncons to x
	fullToUnconsMap = new int[ndof];
	setupMap();
}

void elasticPlate::setupMap()
{
	int c = 0;
	for (int i=0; i < ndof; i++)
	{
		if (isConstrained[i] == 0)
		{
			unconstrainedMap[c] = i;
			fullToUnconsMap[i] = c;
			c++;
		}
	}
}

void elasticPlate::setupMass()
{
	massArray = VectorXd::Zero(ndof);

	double deltaMass;

	int index1;
	int index2;

	double localRadius;

	for (int i = 0; i < edgeNum; i++)
	{
		localRadius = (v_edgeElement[i].x_1_start(0) + v_edgeElement[i].x_2_start(0)) / 2;

		deltaMass = 2 * M_PI * localRadius * v_edgeElement[i].refLength * thickness * density / 2;

		index1 = v_edgeElement[i].nv_1;
		index2 = v_edgeElement[i].nv_2;

		massArray(2 * index1 + 0) = massArray(2 * index1 + 0) + deltaMass;
		massArray(2 * index1 + 1) = massArray(2 * index1 + 1) + deltaMass;
	
		massArray(2 * index2 + 0) = massArray(2 * index2 + 0) + deltaMass;
		massArray(2 * index2 + 1) = massArray(2 * index2 + 1) + deltaMass;

		if ( i == edgeNum - 1 )
		{
			massArray(2 * index1 + 0) = 10 * massArray(2 * index1 + 0);
			massArray(2 * index1 + 1) = 10 * massArray(2 * index1 + 1);
			massArray(2 * index2 + 0) = 10 * massArray(2 * index2 + 0);
			massArray(2 * index2 + 1) = 10 * massArray(2 * index2 + 1);

		}
	}

	massArray(2 * (nv-1) + 0) = 1.0;
	massArray(2 * (nv-1) + 1) = 1.0;
}

int elasticPlate::getIfConstrained(int k)
{
	return isConstrained[k];
}

void elasticPlate::setVertexBoundaryCondition(Vector2d position, int k)
{
	isConstrained[2 * k + 0] = 1;
	isConstrained[2 * k + 1] = 1;
	
	// Store in the constrained dof vector
	x(2 * k + 0) = position(0);
	x(2 * k + 1) = position(1);
}

void elasticPlate::setOneBoundaryCondition(double position, int i, int k)
{
	isConstrained[2 * i + k] = 1;

	x(2 * i + k) = position;
}

Vector2d elasticPlate::getVertex(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x(2 * i + 0);
	xCurrent(1) = x(2 * i + 1);
	
	return xCurrent;
}

Vector2d elasticPlate::getVertexOld(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x0(2 * i + 0);
	xCurrent(1) = x0(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVertexInitial(int i)
{
	Vector2d xCurrent;

	xCurrent(0) = x00(2 * i + 0);
	xCurrent(1) = x00(2 * i + 1);

	return xCurrent;
}

Vector2d elasticPlate::getVelocity(int i)
{
	Vector2d uCurrent;

	uCurrent(0) = u(2 * i + 0);
	uCurrent(1) = u(2 * i + 1);
	
	return uCurrent;
}

void elasticPlate::updateEdgePair()
{
	for (int i = 0; i < edgeNum; i++)
	{
		v_edgeElement[i].x_1 = getVertex(v_edgeElement[i].nv_1);
		v_edgeElement[i].x_2 = getVertex(v_edgeElement[i].nv_2);
		v_edgeElement[i].edgeLength = (v_edgeElement[i].x_1 - v_edgeElement[i].x_2).norm();
	}
}

void elasticPlate::updateBendingPair()
{
	for (int i = 0; i < bendingNum; i++)
	{
		v_bendingElement[i].x_1 = getVertex(v_bendingElement[i].nv_1);
		v_bendingElement[i].x_2 = getVertex(v_bendingElement[i].nv_2);
		v_bendingElement[i].x_3 = getVertex(v_bendingElement[i].nv_3);

		v_bendingElement[i].e_1 = v_bendingElement[i].x_2 - v_bendingElement[i].x_1;
		v_bendingElement[i].e_2 = v_bendingElement[i].x_3 - v_bendingElement[i].x_2;

		v_bendingElement[i].norm_1 =  v_bendingElement[i].e_1.norm();
		v_bendingElement[i].norm_2 =  v_bendingElement[i].e_2.norm();

		v_bendingElement[i].t_1 = v_bendingElement[i].e_1 / v_bendingElement[i].norm_1;
		v_bendingElement[i].t_2 = v_bendingElement[i].e_2 / v_bendingElement[i].norm_2;
	}
}

void elasticPlate::prepareForIteration()
{
	updateEdgePair();
	updateBendingPair();
}

void elasticPlate::updateTimeStep()
{
	prepareForIteration();

	// compute velocity
	u = (x - x0) / dt;

	// update x
	x0 = x;
}

void elasticPlate::updateGuess()
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] = x[unconstrainedMap[c]] + u[unconstrainedMap[c]] * dt;
	}
}

void elasticPlate::updateNewtonMethod(VectorXd m_motion)
{
	for (int c=0; c < uncons; c++)
	{
		x[unconstrainedMap[c]] -= m_motion[c];
	}
}

void elasticPlate::setupGeometry()
{
	v_nodes.clear();
    edge.clear();
    bending.clear();
    constraint.clear();

    double thetaSpace = (M_PI / 2.0000001 - delta);

    double deltaTheta = thetaSpace / ( nv );

    for (int i = 0; i < nv; i++)
    {
    	Vector2d xCurrent;

    	xCurrent(0) = cos(deltaTheta * (i-1) + delta);
    	xCurrent(1) = sin(deltaTheta * (i-1) + delta);

    	//xCurrent(1) = - xCurrent(1);

    	cout << xCurrent << endl;

    	v_nodes.push_back(xCurrent);
    }

    Vector2d xEnd = v_nodes[nv - 1];
    Vector2d xCurrent;
    xCurrent(0) = - xEnd(0);
    xCurrent(1) = xEnd(1);
    v_nodes.push_back(xCurrent);
    nv = nv + 1;


    for (int i = 0; i < nv - 1; i++)
    {
    	Vector2i edgeCurrent;

    	edgeCurrent(0) = i;
    	edgeCurrent(1) = i + 1;

    	edge.push_back(edgeCurrent);
    }

    for (int i = 1; i < nv - 1; i++)
    {
    	Vector3i bendingCurrent;

    	bendingCurrent(0) = i - 1;
    	bendingCurrent(1) = i; 
    	bendingCurrent(2) = i + 1; 

    	bending.push_back(bendingCurrent);
    }


    Vector2d xBall;
    xBall(0) = 0.0;
    xBall(1) = height;
    v_nodes.push_back(xBall);
    nv = nv + 1;

    /*

    ifstream inFile1;
	inFile1.open("inputdata/nodesInput.txt");
	double a, b, c;
	nv = 0;
	while(inFile1 >> a >> b)
	{
		Vector2d xCurrent;

		xCurrent(0) = a;
		xCurrent(1) = b;

		nv = nv + 1;

		v_nodes.push_back(xCurrent);
	}
	inFile1.close();

	ifstream inFile2;
	inFile2.open("inputdata/edgeInput.txt");
	int d, e;
	while(inFile2 >> d >> e)
	{
		Vector2i edgeCurrent;

    	edgeCurrent(0) = d;
    	edgeCurrent(1) = e;

    	edge.push_back(edgeCurrent);
	}
	inFile2.close();

	ifstream inFile3;
	inFile3.open("inputdata/bendingInput.txt");
	int f, g, h;
	while(inFile3 >> f >> g >> h)
	{
		Vector3i bendingCurrent;

    	bendingCurrent(0) = f;
    	bendingCurrent(1) = g; 
    	bendingCurrent(2) = h; 

    	bending.push_back(bendingCurrent);
	}
	inFile3.close();

	*/
}

void elasticPlate::computeEdge()
{
	edgeNum = 0;
	v_edgeElement.clear();

	for (int i = 0; i < edge.size(); i++)
	{
		Vector2i edgeCurrent = edge[i];

		edgeElement m_edgeElement;

		m_edgeElement.nv_1 = edgeCurrent(0);
		m_edgeElement.nv_2 = edgeCurrent(1);

		m_edgeElement.x_1_start = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2_start = getVertex(m_edgeElement.nv_2);

		m_edgeElement.x_1 = getVertex(m_edgeElement.nv_1);
		m_edgeElement.x_2 = getVertex(m_edgeElement.nv_2);

		m_edgeElement.refLength = (m_edgeElement.x_2- m_edgeElement.x_1).norm();
		m_edgeElement.edgeLength = m_edgeElement.refLength;

		Vector2d tangent = (m_edgeElement.x_2 - m_edgeElement.x_1) / (m_edgeElement.x_2 - m_edgeElement.x_1).norm();

		m_edgeElement.phiBar = atan( tangent(1) / tangent(0) );

		m_edgeElement.arrayNum = VectorXi::Zero(4);

		m_edgeElement.arrayNum(0) = 2 * m_edgeElement.nv_1 + 0;
		m_edgeElement.arrayNum(1) = 2 * m_edgeElement.nv_1 + 1;
		
		m_edgeElement.arrayNum(2) = 2 * m_edgeElement.nv_2 + 0;
		m_edgeElement.arrayNum(3) = 2 * m_edgeElement.nv_2 + 1;
		
		v_edgeElement.push_back(m_edgeElement);

		edgeNum = edgeNum + 1;
	}
}

void elasticPlate::computeBending()
{
	bendingNum = 0;
	v_bendingElement.clear();

	for (int i = 0; i < bending.size(); i++)
	{
		Vector3i bendingCurrent = bending[i];

		bendingElement m_bendingElement;

		m_bendingElement.nv_1 = bendingCurrent(0);
		m_bendingElement.nv_2 = bendingCurrent(1);
		m_bendingElement.nv_3 = bendingCurrent(2);

		m_bendingElement.x_1_start = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2_start = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3_start = getVertex(m_bendingElement.nv_3);

		m_bendingElement.x_1 = getVertex(m_bendingElement.nv_1);
		m_bendingElement.x_2 = getVertex(m_bendingElement.nv_2);
		m_bendingElement.x_3 = getVertex(m_bendingElement.nv_3);

		m_bendingElement.e_1 = m_bendingElement.x_2 - m_bendingElement.x_1;
		m_bendingElement.e_2 = m_bendingElement.x_3 - m_bendingElement.x_2;

		m_bendingElement.norm_1 =  m_bendingElement.e_1.norm();
		m_bendingElement.norm_2 =  m_bendingElement.e_2.norm();

		m_bendingElement.voroniLength = (m_bendingElement.norm_1 + m_bendingElement.norm_2) / 2;

		m_bendingElement.t_1 = m_bendingElement.e_1 / m_bendingElement.norm_1;
		m_bendingElement.t_2 = m_bendingElement.e_2 / m_bendingElement.norm_2;

		m_bendingElement.t = (m_bendingElement.t_1 + m_bendingElement.t_2) / 2;

		m_bendingElement.nBar1 = 2 * (m_bendingElement.t_1(0) * m_bendingElement.t_2(1) - m_bendingElement.t_1(1) * m_bendingElement.t_2(0) ) / ( (1 + m_bendingElement.t_1.dot(m_bendingElement.t_2)) * m_bendingElement.voroniLength );
		m_bendingElement.nBar2 = m_bendingElement.t(1) / ( m_bendingElement.t.norm() * m_bendingElement.x_2_start(0) );

		//cout << m_bendingElement.nBar1 << " " << m_bendingElement.nBar2 << endl;

		m_bendingElement.arrayNum = VectorXi::Zero(6);

		m_bendingElement.arrayNum(0) = 2 * m_bendingElement.nv_1 + 0;
		m_bendingElement.arrayNum(1) = 2 * m_bendingElement.nv_1 + 1;
		
		m_bendingElement.arrayNum(2) = 2 * m_bendingElement.nv_2 + 0;
		m_bendingElement.arrayNum(3) = 2 * m_bendingElement.nv_2 + 1;
		
		m_bendingElement.arrayNum(4) = 2 * m_bendingElement.nv_3 + 0;
		m_bendingElement.arrayNum(5) = 2 * m_bendingElement.nv_3 + 1;

		v_bendingElement.push_back(m_bendingElement);

		bendingNum = bendingNum + 1;
	}
}

void elasticPlate::resetMap()
{
	delete isConstrained;
	delete unconstrainedMap;
	delete fullToUnconsMap;

	isConstrained = new int[ndof];
    for (int i=0; i < ndof; i++)
    {
		isConstrained[i] = 0;
    }
}

void elasticPlate::setVelocity()
{
	for (int i = 0; i < nv; i++)
	{
		//u(2 * i + 1) = - height;
	}

	//x = x0 + u * dt;
}