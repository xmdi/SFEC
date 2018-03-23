#include <vector>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Dense>

struct hexCell
{
	float a,b,o,d,L;
};

struct mesh
{
	Eigen::MatrixXf nodes,elems,BC,NF;
};


void generateHexCellMesh(hexCell cell, float elmSize, mesh geo)
{
	int d=roundf(cell.d/elmSize); // elements along depth
	int e=roundf(cell.b*cell.L/elmSize);  // elements across cell wall thickness
	float stumpLength=cell.L/2*(cell.a-cell.b/sqrt(3)); // length of top and bot stumps
	int f=roundf(stumpLength)/elmSize); // elements along vertical wall semi-length
	float inclinedLength=cell.L*(1-cell.b/sqrt(3)); // length of inclined walls
	int g=roundf(inclinedLength/elmSize); // elements along the inclined wall length
	int node=0;

	// bottom stump
	float w=cell.b*cell.L/e;
	float h=stumpLength/f;
	float t=cell.d*cell.L/d;
	for (int k=0; k<=d; k++)
		for (int j=0; j<=f; j++)
			for (int i=0; i<=e; i++)
				geo.nodes.row(node++)<<i*w,j*h,-k*t;
	
	
	
	// bottom left
	h=inclinedLength/g;
	for (int k=0; k<=d; k++)
		for (int j=0; j<=g; j++)
			for (int i=0; i<=e; i++)
				geo.nodes.row(node++)<<i*w*sin(cell.o)-j*h*cos(cell.o),stumpLength+i*w*cos(cell.o)+j*h*sin(cell.o),-k*t;
	// bottom right
	for (int k=0; k<=d; k++)
		for (int j=0; j<=g; j++)
			for (int i=0; i<=e; i++) // TO-DO!!!!! pre-evaluate some of these expressions to save runtime!
				geo.nodes.row(node++)<<cell.b*cell.L/2+i*w*sin(cell.o)+j*h*cos(cell.o),cell.b*cell.L*sqrt(3)/2+stumpLength-i*w*cos(cell.o)+j*h*sin(cell.o),-k*t;
	// top left
	for (int k=0; k<=d; k++)
		for (int j=0; j<=g; j++)
			for (int i=0; i<=e; i++)
				geo.nodes.row(node++)<<i*w*sin(cell.o)-j*h*cos(cell.o),stumpLength+i*w*cos(cell.o)+j*h*sin(cell.o),-k*t;
	// top right
	for (int k=0; k<=d; k++)
		for (int j=0; j<=g; j++)
			for (int i=0; i<=e; i++)
				geo.nodes.row(node++)<<i*w*sin(cell.o)-j*h*cos(cell.o),stumpLength+i*w*cos(cell.o)+j*h*sin(cell.o),-k*t;
	
	



}
