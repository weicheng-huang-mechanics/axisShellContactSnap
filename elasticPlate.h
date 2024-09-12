#ifndef ELASTICPLATE_H
#define ELASTICPLATE_H

#include "eigenIncludes.h"
#include <fstream>

struct edgeElement
{
	int nv_1;
	int nv_2;

	Vector2d x_1;
	Vector2d x_2;

	double phiBar;

	Vector2d x_1_start;
	Vector2d x_2_start;

	double refLength;
	double edgeLength;

	VectorXi arrayNum;
};

struct bendingElement
{
	int nv_1;
	int nv_2;
	int nv_3;

	VectorXi arrayNum;

	Vector2d x_1;
	Vector2d x_2;
	Vector2d x_3;

	Vector2d x_1_start;
	Vector2d x_2_start;
	Vector2d x_3_start;

	Vector2d e_1;
	Vector2d e_2;

	double norm_1;
	double norm_2;

	Vector2d t_1;
	Vector2d t_2;

	Vector2d t;

	double nBar1;
	double nBar2;

	double voroniLength;

	double EI_local;
};

class elasticPlate
{
	public:
	elasticPlate(double m_YoungM, double m_density, double m_thickness, 
		double m_Possion, double m_dt, int m_nv, double m_delta, double m_height, double m_massSpeed);
	~elasticPlate();

	double YoungM;
	double thickness;
	double Possion;
	double dt;
	double density;
	double delta;

	double massSpeed;

	double height;

	double beta0;

	Vector2d getVertex(int i);
	Vector2d getVertexOld(int i);
	Vector2d getVertexInitial(int i);
	Vector2d getVelocity(int i);

	VectorXd x;
	VectorXd x0;
	VectorXd x00;
	VectorXd u;

	std::vector<Vector2d> v_nodes;
    std::vector<Vector2i> edge;
    std::vector<Vector3i> bending;
    std::vector<int> constraint;

	std::vector<edgeElement> v_edgeElement;
	std::vector<bendingElement> v_bendingElement;

	int temp;

	int nv;
	int edgeNum;
	int bendingNum;

	int ndof;
	int uncons;
	int ncons;

	void setupGeometry();

	void setVertexBoundaryCondition(Vector2d position, int k);
	void setOneBoundaryCondition(double position, int i, int k);

	void computeEdge();
	void computeBending();

	void updateEdgePair();
	void updateBendingPair();

	// boundary conditions
	int* isConstrained;
	int getIfConstrained(int k);
	int* unconstrainedMap;
	int* fullToUnconsMap;
	void setup();
	void setupMap();
	void setVelocity();

	void updateTimeStep();
	void updateGuess();
	void updateNewtonMethod(VectorXd m_motion);
	void prepareForIteration();

	void resetMap();

	VectorXd massArray;
	void setupMass();
	VectorXi boundaryIndex;

	double EA;
	double EI;

	double refHeight;

	private:
};

#endif
