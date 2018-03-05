#include "core.h"

void testModel_1()
{
	std::chrono::high_resolution_clock::time_point t0=std::chrono::high_resolution_clock::now(); // start timer

	prop mat(32,.3333333); // material property structure, Young's Modulus & Poisson's Ratio
	mesh geo; // model mesh structure
	math eqn(24); // matrix equation structure, number of DOF

	geo.nodes=Eigen::MatrixXf::Zero(8,3); // node location list
	geo.nodes<<-1,-1,-1,
	   1,-1,-1,
	   1,1,-1,
	   -1,1,-1,
	   -1,-1,1,
	   1,-1,1,
	   1,1,1,
	   -1,1,1;
	geo.elems=Eigen::MatrixXf::Zero(1,8); // nodes per element list
	geo.elems<<0,1,2,3,4,5,6,7;
		geo.BC=Eigen::MatrixXf::Zero(7,3); // boundary conditions
	geo.BC<<0,0,0,
		0,1,0,	
		0,2,0,	
		3,0,0,
		3,2,0,
	    4,0,0,
	    7,0,0;
	geo.NF=Eigen::MatrixXf::Zero(4,3); // nodal forces
	geo.NF<<1,0,10,
		2,0,10,
	    5,0,10,
	    6,0,10;

	globalHex8(mat,eqn,geo); // assemble stiffness matrix from 8-node brick elements

	applyBoundaryConditions(eqn,geo); // assign BCs and NFs

	eqn.d=eqn.K.llt().solve(eqn.f); // solve matrix problem for DOFs

	writeData(eqn,geo,"testModel_1"); // write output to dat files

	std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now(); // stop timer
	double dif=std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count(); // elapsed time

	// print relevant results
	std::cout<<"TEST1 ::: DOF: "<<eqn.K.rows()<<",   Relative Error: "<<(eqn.K*eqn.d-eqn.f).norm()/eqn.f.norm()<<",   Elapsed Time: "<<dif/1e9<<" seconds\n";
}

void testModel_2()
{
	std::chrono::high_resolution_clock::time_point t0=std::chrono::high_resolution_clock::now(); // start timer

	prop mat(32,.3333333); // material property structure, Young's Modulus & Poisson's Ratio
	mesh geo; // model mesh structure
	math eqn(36); // matrix equation structure, number of DOF

	geo.nodes=Eigen::MatrixXf::Zero(12,3); // node location list
	geo.nodes<<-1,-1,-1,
	   1,-1,-1,
	   1,1,-1,
	   -1,1,-1,
	   -1,-1,1,
	   1,-1,1,
	   1,1,1,
	   -1,1,1,
	   -1,-1,3,
	   1,-1,3,
	   1,1,3,
	   -1,1,3;
	geo.elems=Eigen::MatrixXf::Zero(2,8); // nodes per element list
	geo.elems<<0,1,2,3,4,5,6,7,
		4,5,6,7,8,9,10,11;
	geo.BC=Eigen::MatrixXf::Zero(9,3); // boundary conditions
	geo.BC<<0,0,0,
		0,1,0,	
		0,2,0,	
		3,0,0,
		3,2,0,
	    4,0,0,
	    7,0,0,
	    8,0,0,
	    11,0,0;
	geo.NF=Eigen::MatrixXf::Zero(6,3); // nodal forces
	geo.NF<<1,0,10,
		2,0,10,
	    5,0,20,
	    6,0,20,
	    9,0,10,
	    10,0,10;

	globalHex8(mat,eqn,geo); // assemble stiffness matrix from 8-node brick elements

	applyBoundaryConditions(eqn,geo); // assign BCs and NFs

	eqn.d=eqn.K.llt().solve(eqn.f); // solve matrix problem for DOFs

	writeData(eqn,geo,"testModel_2"); // write output to dat files

	std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now(); // stop timer
	double dif=std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count(); // elapsed time

	// print relevant results
	std::cout<<"TEST2 ::: DOF: "<<eqn.K.rows()<<",   Relative Error: "<<(eqn.K*eqn.d-eqn.f).norm()/eqn.f.norm()<<",   Elapsed Time: "<<dif/1e9<<" seconds\n";
}

void testModel_3()
{
	std::chrono::high_resolution_clock::time_point t0=std::chrono::high_resolution_clock::now(); // start timer

	prop mat(32,.3333333); // material property structure, Young's Modulus & Poisson's Ratio
	mesh geo; // model mesh structure
	math eqn(96); // matrix equation structure, number of DOF

	geo.nodes=Eigen::MatrixXf::Zero(32,3); // node location list
	geo.nodes<<0,0,0, 
		1,0,0,
		1,1,0,
		0,1,0,
		0,0,1, 
		1,0,1,
		1,1,1,
		0,1,1,
		0,0,2, 
		1,0,2,
		1,1,2,
		0,1,2,
		0,0,3, 
		1,0,3,
		1,1,3,
		0,1,3,
	    4,0,0, 
		5,0,0,
		5,1,0,
		4,1,0,
		4,0,1, 
		5,0,1,
		5,1,1,
		4,1,1,
		4,0,2, 
		5,0,2,
		5,1,2,
		4,1,2,
		4,0,3, 
		5,0,3,
		5,1,3,
		4,1,3;
	geo.elems=Eigen::MatrixXf::Zero(6,8); // nodes per element list
	geo.elems<<0,1,2,3,4,5,6,7,
			4,5,6,7,8,9,10,11,
			8,9,10,11,12,13,14,15,
			16,17,18,19,20,21,22,23,
			20,21,22,23,24,25,26,27,
			24,25,26,27,28,29,30,31;
	
	geo.BC=Eigen::MatrixXf::Zero(24,3); // boundary conditions
	geo.BC<<0,0,0,
		0,1,0,	
		0,2,0,	
		1,0,0,
		1,1,0,
	    1,2,0,
	    2,0,0,
	    2,1,0,
	    2,2,0,
	    3,0,0,
	    3,1,0,
	    3,2,0,
		16,0,0,
		16,1,0,	
		16,2,0,	
		17,0,0,
		17,1,0,
	    17,2,0,
	    18,0,0,
	    18,1,0,
	    18,2,0,
	    19,0,0,
	    19,1,0,
	    19,2,0;
	geo.NF=Eigen::MatrixXf::Zero(4,3); // nodal forces
	geo.NF<<13,0,.25,
			14,0,.25,
			28,0,-.25,
			31,0,-.25;

	globalHex8(mat,eqn,geo); // assemble stiffness matrix from 8-node brick elements

	applyBoundaryConditions(eqn,geo); // assign BCs and NFs

	eqn.d=eqn.K.llt().solve(eqn.f); // solve matrix problem for DOFs

	writeData(eqn,geo,"testModel_3"); // write output to dat files

	std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now(); // stop timer
	double dif=std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count(); // elapsed time

	// print relevant results
	std::cout<<"TEST3 ::: DOF: "<<eqn.K.rows()<<",   Relative Error: "<<(eqn.K*eqn.d-eqn.f).norm()/eqn.f.norm()<<",   Elapsed Time: "<<dif/1e9<<" seconds\n";
}

int main()
{
	testModel_1(); // run test 1, single element
	testModel_2(); // run test 2, two elements
	testModel_3(); // run test 3, disjointed bodies
	return 0;
}
