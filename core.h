#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <algorithm>

struct prop
{
	Eigen::Matrix<float,6,6,Eigen::DontAlign> D;
	prop(float E, float v)
	{
		D<<1-v,v,v,0,0,0,
		   v,1-v,v,0,0,0,
		   v,v,1-v,0,0,0,
		   0,0,0,.5-v,0,0,
		   0,0,0,0,.5-v,0,
		   0,0,0,0,0,.5-v;
		D*=E/((1+v)*(1-2*v)); 
	}
};

struct mesh
{
	Eigen::MatrixXf nodes,elems,BC,NF;
};

struct math
{	
	Eigen::MatrixXf K; 
	Eigen::VectorXf d,f;
	math(int DOF)	
	{
		K=Eigen::MatrixXf::Zero(DOF,DOF);
		f=Eigen::VectorXf::Zero(DOF);
		d=Eigen::VectorXf::Zero(DOF);
	}
};

float dN1dr(Eigen::Vector3f x){return -(1-x(1))*(1-x(2))/8;}
float dN1ds(Eigen::Vector3f x){return -(1-x(0))*(1-x(2))/8;}
float dN1dt(Eigen::Vector3f x){return -(1-x(0))*(1-x(1))/8;}
float dN2dr(Eigen::Vector3f x){return -(1+x(1))*(1-x(2))/8;}
float dN2ds(Eigen::Vector3f x){return (1-x(0))*(1-x(2))/8;}
float dN2dt(Eigen::Vector3f x){return -(1-x(0))*(1+x(1))/8;}
float dN3dr(Eigen::Vector3f x){return (1+x(1))*(1-x(2))/8;}
float dN3ds(Eigen::Vector3f x){return (1+x(0))*(1-x(2))/8;}
float dN3dt(Eigen::Vector3f x){return -(1+x(0))*(1+x(1))/8;}
float dN4dr(Eigen::Vector3f x){return (1-x(1))*(1-x(2))/8;}
float dN4ds(Eigen::Vector3f x){return -(1+x(0))*(1-x(2))/8;}
float dN4dt(Eigen::Vector3f x){return -(1+x(0))*(1-x(1))/8;}
float dN5dr(Eigen::Vector3f x){return -(1-x(1))*(1+x(2))/8;}
float dN5ds(Eigen::Vector3f x){return -(1-x(0))*(1+x(2))/8;}
float dN5dt(Eigen::Vector3f x){return (1-x(0))*(1-x(1))/8;}
float dN6dr(Eigen::Vector3f x){return -(1+x(1))*(1+x(2))/8;}
float dN6ds(Eigen::Vector3f x){return (1-x(0))*(1+x(2))/8;}
float dN6dt(Eigen::Vector3f x){return (1-x(0))*(1+x(1))/8;}
float dN7dr(Eigen::Vector3f x){return (1+x(1))*(1+x(2))/8;}
float dN7ds(Eigen::Vector3f x){return (1+x(0))*(1+x(2))/8;}
float dN7dt(Eigen::Vector3f x){return (1+x(0))*(1+x(1))/8;}
float dN8dr(Eigen::Vector3f x){return (1-x(1))*(1+x(2))/8;}
float dN8ds(Eigen::Vector3f x){return -(1+x(0))*(1+x(2))/8;}
float dN8dt(Eigen::Vector3f x){return (1+x(0))*(1-x(1))/8;}

auto elementalHex8(prop &mat, auto &X)
{
//	N1=(1-r)(1-s)(1-t)/8
//		dN1/dr=-(1-s)(1-t)/8
//		dN1/ds=-(1-r)(1-t)/8
//		dN1/dt=-(1-r)(1-s)/8
//	N2=(1-r)(1+s)(1-t)/8
//		dN2/dr=-(1+s)(1-t)/8
//		dN2/ds=(1-r)(1-t)/8
//		dN2/dt=-(1-r)(1+s)/8
//	N3=(1+r)(1+s)(1-t)/8
//		dN3/dr=(1+s)(1-t)/8
//		dN3/ds=(1+r)(1-t)/8
//		dN3/dt=-(1+r)(1+s)/8
//	N4=(1+r)(1-s)(1-t)/8
//		dN4/dr=(1-s)(1-t)/8
//		dN4/ds=-(1+r)(1-t)/8
//		dN4/dt=-(1+r)(1-s)/8
//	N5=(1-r)(1-s)(1+t)/8
//		dN5/dr=-(1-s)(1+t)/8
//		dN5/ds=-(1-r)(1+t)/8
//		dN5/dt=(1-r)(1-s)/8
//	N6=(1-r)(1+s)(1+t)/8
//		dN6/dr=-(1+s)(1+t)/8
//		dN6/ds=(1-r)(1+t)/8
//		dN6/dt=(1-r)(1+s)/8
//	N7=(1+r)(1+s)(1+t)/8
//		dN7/dr=(1+s)(1+t)/8
//		dN7/ds=(1+r)(1+t)/8
//		dN7/dt=(1+r)(1+s)/8
//	N8=(1+r)(1-s)(1+t)/8
//		dN8/dr=(1-s)(1+t)/8
//		dN8/ds=-(1+r)(1+t)/8
//		dN8/dt=(1+r)(1-s)/8

	Eigen::MatrixXf intpts(8,3);
	intpts<<1/sqrt(3),1/sqrt(3),1/sqrt(3),
		1/sqrt(3),1/sqrt(3),-1/sqrt(3),
		1/sqrt(3),-1/sqrt(3),1/sqrt(3),
		1/sqrt(3),-1/sqrt(3),-1/sqrt(3),
		-1/sqrt(3),1/sqrt(3),1/sqrt(3),
		-1/sqrt(3),1/sqrt(3),-1/sqrt(3),
		-1/sqrt(3),-1/sqrt(3),1/sqrt(3),
		-1/sqrt(3),-1/sqrt(3),-1/sqrt(3);
	
	Eigen::MatrixXf Kel=Eigen::MatrixXf::Zero(24,24);

	for (int i=0;i<intpts.rows();i++)
	{
		float J11=dN1dr(intpts.row(i))*X(0,0)+
			dN2dr(intpts.row(i))*X(1,0)+
			dN3dr(intpts.row(i))*X(2,0)+
			dN4dr(intpts.row(i))*X(3,0)+
			dN5dr(intpts.row(i))*X(4,0)+
			dN6dr(intpts.row(i))*X(5,0)+
			dN7dr(intpts.row(i))*X(6,0)+
			dN8dr(intpts.row(i))*X(7,0);
		float J21=dN1ds(intpts.row(i))*X(0,0)+
			dN2ds(intpts.row(i))*X(1,0)+
			dN3ds(intpts.row(i))*X(2,0)+
			dN4ds(intpts.row(i))*X(3,0)+
			dN5ds(intpts.row(i))*X(4,0)+
			dN6ds(intpts.row(i))*X(5,0)+
			dN7ds(intpts.row(i))*X(6,0)+
			dN8ds(intpts.row(i))*X(7,0);
		float J31=dN1dt(intpts.row(i))*X(0,0)+
			dN2dt(intpts.row(i))*X(1,0)+
			dN3dt(intpts.row(i))*X(2,0)+
			dN4dt(intpts.row(i))*X(3,0)+
			dN5dt(intpts.row(i))*X(4,0)+
			dN6dt(intpts.row(i))*X(5,0)+
			dN7dt(intpts.row(i))*X(6,0)+
			dN8dt(intpts.row(i))*X(7,0);
		float J12=dN1dr(intpts.row(i))*X(0,1)+
			dN2dr(intpts.row(i))*X(1,1)+
			dN3dr(intpts.row(i))*X(2,1)+
			dN4dr(intpts.row(i))*X(3,1)+
			dN5dr(intpts.row(i))*X(4,1)+
			dN6dr(intpts.row(i))*X(5,1)+
			dN7dr(intpts.row(i))*X(6,1)+
			dN8dr(intpts.row(i))*X(7,1);
		float J22=dN1ds(intpts.row(i))*X(0,1)+
			dN2ds(intpts.row(i))*X(1,1)+
			dN3ds(intpts.row(i))*X(2,1)+
			dN4ds(intpts.row(i))*X(3,1)+
			dN5ds(intpts.row(i))*X(4,1)+
			dN6ds(intpts.row(i))*X(5,1)+
			dN7ds(intpts.row(i))*X(6,1)+
			dN8ds(intpts.row(i))*X(7,1);
		float J32=dN1dt(intpts.row(i))*X(0,1)+
			dN2dt(intpts.row(i))*X(1,1)+
			dN3dt(intpts.row(i))*X(2,1)+
			dN4dt(intpts.row(i))*X(3,1)+
			dN5dt(intpts.row(i))*X(4,1)+
			dN6dt(intpts.row(i))*X(5,1)+
			dN7dt(intpts.row(i))*X(6,1)+
			dN8dt(intpts.row(i))*X(7,1);
		float J13=dN1dr(intpts.row(i))*X(0,2)+
			dN2dr(intpts.row(i))*X(1,2)+
			dN3dr(intpts.row(i))*X(2,2)+
			dN4dr(intpts.row(i))*X(3,2)+
			dN5dr(intpts.row(i))*X(4,2)+
			dN6dr(intpts.row(i))*X(5,2)+
			dN7dr(intpts.row(i))*X(6,2)+
			dN8dr(intpts.row(i))*X(7,2);
		float J23=dN1ds(intpts.row(i))*X(0,2)+
			dN2ds(intpts.row(i))*X(1,2)+
			dN3ds(intpts.row(i))*X(2,2)+
			dN4ds(intpts.row(i))*X(3,2)+
			dN5ds(intpts.row(i))*X(4,2)+
			dN6ds(intpts.row(i))*X(5,2)+
			dN7ds(intpts.row(i))*X(6,2)+
			dN8ds(intpts.row(i))*X(7,2);
		float J33=dN1dt(intpts.row(i))*X(0,2)+
			dN2dt(intpts.row(i))*X(1,2)+
			dN3dt(intpts.row(i))*X(2,2)+
			dN4dt(intpts.row(i))*X(3,2)+
			dN5dt(intpts.row(i))*X(4,2)+
			dN6dt(intpts.row(i))*X(5,2)+
			dN7dt(intpts.row(i))*X(6,2)+
			dN8dt(intpts.row(i))*X(7,2);
		
		float Ji11=J22*J33-J23*J32;
		float Ji33=J11*J22-J12*J21;
		float Ji23=J31*J12-J32*J11;
		float Ji21=J32*J13-J12*J33;
		float Ji13=J21*J22-J31*J22;
		float Ji22=J33*J11-J31*J13;
		float Ji12=J23*J31-J21*J33;
		float Ji31=J13*J23-J13*J22;
		float Ji32=J13*J21-J23*J11;
		float detJ=J11*Ji11+J12*Ji21+J13*Ji31;
			
		Eigen::Matrix3f J;
		J<<J11,J12,J13,
			J21,J22,J23,
			J31,J32,J33;

		Eigen::Matrix3f Ji;
		Ji=J.inverse();

		Eigen::MatrixXf B=Eigen::MatrixXf::Zero(6,24);	
		Eigen::Vector3f temp;
		Eigen::Vector3f temp2;

		// 1
		temp<<dN1dr(intpts.row(i)),dN1ds(intpts.row(i)),dN1dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,0)=temp2(0);
		B(3,1)=temp2(0);
		B(5,2)=temp2(0);
		B(1,1)=temp2(1);
		B(3,0)=temp2(1);
		B(4,2)=temp2(1);
		B(2,2)=temp2(2);
		B(4,1)=temp2(2);
		B(5,0)=temp2(2);
		// 2
		temp<<dN2dr(intpts.row(i)),dN2ds(intpts.row(i)),dN2dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,3)=temp2(0);
		B(3,4)=temp2(0);
		B(5,5)=temp2(0);
		B(1,4)=temp2(1);
		B(3,3)=temp2(1);
		B(4,5)=temp2(1);
		B(2,5)=temp2(2);
		B(4,4)=temp2(2);
		B(5,3)=temp2(2);
		// 3
		temp<<dN3dr(intpts.row(i)),dN3ds(intpts.row(i)),dN3dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,6)=temp2(0);
		B(3,7)=temp2(0);
		B(5,8)=temp2(0);
		B(1,7)=temp2(1);
		B(3,6)=temp2(1);
		B(4,8)=temp2(1);
		B(2,8)=temp2(2);
		B(4,7)=temp2(2);
		B(5,6)=temp2(2);
		// 4
		temp<<dN4dr(intpts.row(i)),dN4ds(intpts.row(i)),dN4dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,9)=temp2(0);
		B(3,10)=temp2(0);
		B(5,11)=temp2(0);
		B(1,10)=temp2(1);
		B(3,9)=temp2(1);
		B(4,11)=temp2(1);
		B(2,11)=temp2(2);
		B(4,10)=temp2(2);
		B(5,9)=temp2(2);
		// 5
		temp<<dN5dr(intpts.row(i)),dN5ds(intpts.row(i)),dN5dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,12)=temp2(0);
		B(3,13)=temp2(0);
		B(5,14)=temp2(0);
		B(1,13)=temp2(1);
		B(3,12)=temp2(1);
		B(4,14)=temp2(1);
		B(2,14)=temp2(2);
		B(4,13)=temp2(2);
		B(5,12)=temp2(2);
		// 6
		temp<<dN6dr(intpts.row(i)),dN6ds(intpts.row(i)),dN6dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,15)=temp2(0);
		B(3,16)=temp2(0);
		B(5,17)=temp2(0);
		B(1,16)=temp2(1);
		B(3,15)=temp2(1);
		B(4,17)=temp2(1);
		B(2,17)=temp2(2);
		B(4,16)=temp2(2);
		B(5,15)=temp2(2);
		// 7
		temp<<dN7dr(intpts.row(i)),dN7ds(intpts.row(i)),dN7dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,18)=temp2(0);
		B(3,19)=temp2(0);
		B(5,20)=temp2(0);
		B(1,19)=temp2(1);
		B(3,18)=temp2(1);
		B(4,20)=temp2(1);
		B(2,20)=temp2(2);
		B(4,19)=temp2(2);
		B(5,18)=temp2(2);
		// 8
		temp<<dN8dr(intpts.row(i)),dN8ds(intpts.row(i)),dN8dt(intpts.row(i));
		temp2=Ji*temp;
		B(0,21)=temp2(0);
		B(3,22)=temp2(0);
		B(5,23)=temp2(0);
		B(1,22)=temp2(1);
		B(3,21)=temp2(1);
		B(4,23)=temp2(1);
		B(2,23)=temp2(2);
		B(4,22)=temp2(2);
		B(5,21)=temp2(2);
		Kel-=B.transpose()*mat.D*B*detJ;	
	}
	return Kel;
}

void applyBoundaryConditions(math &eqn, mesh &geo)
{
	for (int i=0; i<geo.NF.rows(); i++)
		eqn.f(3*geo.NF(i,0)+geo.NF(i,1))=geo.NF(i,2); // apply nodal force
	for (int i=0; i<geo.BC.rows(); i++)
	{
		int m=3*geo.BC(i,0)+geo.BC(i,1); // constrained DOF
		eqn.d(m)=geo.BC(i,2); // apply boundary condition, line not required, recomputed later
		for (int j=0; j<eqn.d.rows(); j++) // loop over DOFS
			if (j-m) // zero out rows and columns of K, move terms to f
			{
				eqn.f(j)-=eqn.K(j,m)*eqn.d(m);
				eqn.K(m,j)=0;
				eqn.K(j,m)=0;
			}
			else // put 1 along the "diagonal" and manually set f to match d
			{
				eqn.f(m)=eqn.d(m);
				eqn.K(m,m)=1;
			}
	}
}

void globalHex8(prop &mat, math &eqn, mesh &geo)
{
	Eigen::MatrixXf X=Eigen::MatrixXf::Zero(8,3);
	Eigen::MatrixXf Kel=Eigen::MatrixXf::Zero(24,24);
	for (int i=0; i<geo.elems.rows(); i++)
	{
		X<<geo.nodes.row((geo.elems(i,0))),
		   geo.nodes.row((geo.elems(i,1))),
		   geo.nodes.row((geo.elems(i,2))),
		   geo.nodes.row((geo.elems(i,3))),
		   geo.nodes.row((geo.elems(i,4))),
		   geo.nodes.row((geo.elems(i,5))),
		   geo.nodes.row((geo.elems(i,6))),
		   geo.nodes.row((geo.elems(i,7)));
		Kel=elementalHex8(mat,X);
		for (int j=0; j<8; j++) 
			for (int k=0; k<8; k++)
			{		
				int m=geo.elems(i,j);
				int n=geo.elems(i,k);
				eqn.K(3*m,3*n)+=Kel(3*j,3*k);
				eqn.K(3*m+1,3*n)+=Kel(3*j+1,3*k);
				eqn.K(3*m+2,3*n)+=Kel(3*j+2,3*k);
				eqn.K(3*m,3*n+1)+=Kel(3*j,3*k+1);
				eqn.K(3*m+1,3*n+1)+=Kel(3*j+1,3*k+1);
				eqn.K(3*m+2,3*n+1)+=Kel(3*j+2,3*k+1);
				eqn.K(3*m,3*n+2)+=Kel(3*j,3*k+2);
				eqn.K(3*m+1,3*n+2)+=Kel(3*j+1,3*k+2);
				eqn.K(3*m+2,3*n+2)+=Kel(3*j+2,3*k+2);
			}
	}
}

void writeData(math &eqn, mesh &geo, std::string filename)
{
	int w[]={0, 1, 5, 1, 2, 6, 2, 3, 7, 3, 0, 4, 5, 6, 7, 4};
	std::ofstream out0(("dat/"+filename+"_0.dat").c_str());
	for (int i=0; i<geo.elems.rows(); i++)
	{	
		for (int k=0; k<16; k++)
			out0<<geo.nodes(geo.elems(i,w[k]),0)<<"  "
				<<geo.nodes(geo.elems(i,w[k]),1)<<"  "
				<<geo.nodes(geo.elems(i,w[k]),2)<<"\n";
		out0<<"\n\n";
	}
	std::ofstream out1(("dat/"+filename+"_1.dat").c_str());
	for (int i=0; i<geo.elems.rows(); i++)
	{
		for (int k=0; k<16; k++)
			out1<<geo.nodes(geo.elems(i,w[k]),0)+eqn.d(3*geo.elems(i,w[k]))<<"  "
				<<geo.nodes(geo.elems(i,w[k]),1)+eqn.d(3*geo.elems(i,w[k])+1)<<"  "
				<<geo.nodes(geo.elems(i,w[k]),2)+eqn.d(3*geo.elems(i,w[k])+2)<<"\n";
		out1<<"\n\n";
	}
}
