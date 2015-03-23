/* * cut_cell_integration.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: afonsoal*/



#include "cut_cell_3D.h"
#include <fstream>
#include <iostream>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <iostream>
#include <math.h> // log(), pow
#include "SetPolyhedron.h"
using namespace dealii;

// This function will pass all the relevant objects (fe, fe_values, quadrature)
// to private, in order to be accessible to all the methods
// For some strange reason, I could not pass FEValues the same way I passed
// quadrature_formula & fe (to private) . WHY????
// However, I can pass individually for each member function.
//template <int dim>
cut_cell_3D::cut_cell_3D(FEValues<3> const & fe_values,
		FE_Q<3> const& fe,
		QGauss<3> const &quadrature_formula,
		QGauss<3> const &face_quadrature_formula
		/*,int dim*/)
//:   fe1(fe), quadrature_formula1(quadrature_formula),
//		face_quadrature_formula1(face_quadrature_formula)
/*,fe_values1(fe_values)*/

{
//	jacobian_inverse.reinit(2,2);
//	coefficients.reinit(4,4);
//	dofs_per_cell = fe1.dofs_per_cell;
//	n_q_points    = quadrature_formula1.size();
//	for (unsigned int i=0; i<2; ++i)
//		for (unsigned int j=0; j<2; ++j)
//		{
//			jacobian_inverse(i,j) = fe_values.inverse_jacobian(0)[i][j];
////			std::cout << "jacobian_inverse(i,j): "
////			<< fe_values.inverse_jacobian(0)[i][j] << "\n";
//		}
//
//
//	for (unsigned int i = 0;i<4;++i)
//		for (unsigned int j = 0;j<4;++j)
//		{
//			if (i==j || (i==0 && j==3))
//				coefficients(i,j) = 1;
//			else if((i==0 && (j==1 || j==2)) ||
//					(j==3 && (i==1 || i==2)) )
//				coefficients(i,j) = -1;
//			else coefficients(i,j) = 0;
//		}
//	cell_quadrature_weights = quadrature_formula1.get_weights();
//
//	face_quadrature_weights = face_quadrature_formula1.get_weights();
//	face_quadrature_points = face_quadrature_formula1.get_points();
//	n_face_q_points = face_quadrature_formula1.size();
//
//	jacobian_determinant = fe_values.JxW(0)/
//			(cell_quadrature_weights[0] /**cell_quadrature_weights[1]*/ );
////	 cell_quadrature_weights is already wi*wi (= 0.25*0.25 = 0.625)
////	std::cout << "cell_quadrature_weights[0]: " << cell_quadrature_weights[0] << "\n";
////				std::cout << "face_quadrature_weights[0]: " << face_quadrature_weights[0] << "\n";
//	//			std::cout << "fe_values.JxW (0)1: " << fe_values.JxW(0) << "\n";
//	//			std::cout << "fe_values.shape_value (0,0) " << fe_values.shape_value (0,0) << "\n";
//


}

double cut_cell_3D::return_rhs_face_integration (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i,
		const double face_length)
{
//	Vector<double> exponents_matrix_x(4);
//	Vector<double> exponents_matrix_y(4);
//	exponents_matrix_x(1) = 1;
//	exponents_matrix_x(3) = 1;
//	exponents_matrix_y(2) = 1;
//	exponents_matrix_y(3) = 1;
//
//	double xt;
//	double yt;
//	double final_face_integration = 0;
//	for (unsigned int q_point = 0;q_point<face_quadrature_points.size();++q_point)
//	{
//		xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Point<2> row_sum;
//		double face_integration_scalar = 0;
//		for (unsigned int i = 0;i<4;++i) {
//			Point<2> new_Xt;
//			// First row, representing d phi_i/dx * d phi_j/dx
//			new_Xt[0] = pow(xt,(exponents_matrix_x(i)+1))/(2*(1+exponents_matrix_x(i)))
//																																																													* pow(yt,exponents_matrix_y(i));
//			new_Xt[1] = pow(yt,(exponents_matrix_y(i)+1))/(2*(1+exponents_matrix_y(i)))
//																																																													* pow(xt,exponents_matrix_x(i));
//			new_Xt *= coefficients(dof_i,i);
//			row_sum += new_Xt;
//		}
//
//		face_integration_scalar = row_sum*normal; // point*point = scalar (Is this the correct form?)
//		face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];
//		face_integration_scalar *= face_length;
//		final_face_integration += face_integration_scalar;
//	} // end quadrature sum
//	return final_face_integration;
	return 0;
}

int cut_cell_3D::GetDelta (const int var1, const int var2)
{
	// Return the value of the "delta" multiplier, which acts to select the appropriate multiplication
	// of phi_i, phi_j.
	int delta;

	if (var1 == var2)	delta = 1;
	else 				delta = 0;

	return delta;
}




int cut_cell_3D::GetACoeff (const int index, const int dof, const int var)
{
	FullMatrix c_coeff(8,8);
	Vector</*double*/int> a_coeff(7);
	a_coeff(0) = delta(0,var)*c_coeff(1,dof)+delta(1,var)*c_coeff(2,dof)+delta(2,var)*c_coeff(3,dof);
	a_coeff(1) = delta(1,var)*c_cpeff(4,dof)+delta(2,var)*c_coeff(6,dof);
	a_coeff(2) = delta(0,var)*c_cpeff(4,dof)+delta(2,var)*c_coeff(5,dof);
	a_coeff(3) = delta(1,var)*c_cpeff(5,dof)+delta(1,var)*c_coeff(6,dof);
	a_coeff(4) = delta(0,var)*c_cpeff(7,dof);
	a_coeff(5) = delta(1,var)*c_cpeff(7,dof);
	a_coeff(6) = delta(2,var)*c_cpeff(7,dof);
	return a_coeff(index)
}
// Returns face integration for the stiffness (A) matrix.

double cut_cell_3D::return_face_integration (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,
		const int dof_i, const int dof_j, const double face_length)
{
}
Point<2> cut_cell_3D::CompHVector(const double xt, const double yt, const int m, const int n)
{

	Vector<double> k_p(23);

	k_p(0) = GetACoeff(0,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(1) = GetACoeff(0,dof_i,m)*GetACoeff(1,dof_j,n)+GetACoeff(1,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(2) = GetACoeff(0,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(1,dof_i,m)*GetACoeff(2,dof_j,n)
	+GetACoeff(2,dof_i,m)*GetACoeff(1,dof_j,n)+GetACoeff(6,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(3) = GetACoeff(1,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(3,dof_i,m)*GetACoeff(6,dof_j,n)+
			GetACoeff(4,dof_i,m)*GetACoeff(1,dof_j,n)+GetACoeff(5,dof_i,m)*GetACoeff(2,dof_j,n)+
			GetACoeff(6,dof_i,m)*GetACoeff(3,dof_j,n);
	k_p(4) = GetACoeff(0,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(1,dof_i,m)*GetACoeff(3,dof_j,n)+
			GetACoeff(3,dof_i,m)*GetACoeff(1,dof_j,n)+GetACoeff(5,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(5) = GetACoeff(0,dof_i,m)*GetACoeff(2,dof_j,n)+GetACoeff(2,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(6) = GetACoeff(6,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(2,dof_i,m)*GetACoeff(3,dof_j,n)+
				GetACoeff(3,dof_i,m)*GetACoeff(2,dof_j,n)+GetACoeff(4,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(7) = GetACoeff(0,dof_i,m)*GetACoeff(3,dof_j,n)+GetACoeff(3,dof_i,m)*GetACoeff(0,dof_j,n);
	k_p(8) = GetACoeff(2,dof_i,m)*GetACoeff(2,dof_j,n);
	k_p(9) = GetACoeff(6,dof_i,m)*GetACoeff(1,dof_j,n)+GetACoeff(1,dof_i,m)*GetACoeff(7,dof_j,n);
	k_p(10) = GetACoeff(5,dof_i,m)*GetACoeff(6,dof_j,n)+GetACoeff(6,dof_i,m)*GetACoeff(5,dof_j,n);
	k_p(11) = GetACoeff(6,dof_i,m)*GetACoeff(6,dof_j,n);
	k_p(12) = GetACoeff(1,dof_i,m)*GetACoeff(5,dof_j,n)+GetACoeff(5,dof_i,m)*GetACoeff(1,dof_j,n);
	k_p(13) = GetACoeff(5,dof_i,m)*GetACoeff(5,dof_j,n);
	k_p(14) = GetACoeff(2,dof_i,m)*GetACoeff(2,dof_j,n);
	k_p(15) = GetACoeff(2,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(4,dof_i,m)*GetACoeff(2,dof_j,n);
	k_p(16) = GetACoeff(2,dof_i,m)*GetACoeff(6,dof_j,n)+GetACoeff(6,dof_i,m)*GetACoeff(2,dof_j,n);
	k_p(17) = GetACoeff(3,dof_i,m)*GetACoeff(3,dof_j,n);
	k_p(18) = GetACoeff(3,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(4,dof_i,m)*GetACoeff(3,dof_j,n);
	k_p(19) = GetACoeff(3,dof_i,m)*GetACoeff(5,dof_j,n)+GetACoeff(5,dof_i,m)*GetACoeff(3,dof_j,n);
	k_p(20) = GetACoeff(4,dof_i,m)*GetACoeff(4,dof_j,n);
	k_p(21) = GetACoeff(5,dof_i,m)*GetACoeff(4,dof_j,n)+GetACoeff(4,dof_i,m)*GetACoeff(5,dof_j,n);
	k_p(22) = GetACoeff(4,dof_i,m)*GetACoeff(6,dof_j,n)+GetACoeff(6,dof_i,m)*GetACoeff(4,dof_j,n);

	std::vector<std::vector<int> > alfa(23,3);
//	FullMatrix<double> alfa(23,3);
	int alfa[][3] = {{0,0,0},{1,0,0},{1,1,0},
					  {1,1,1},{1,0,1},{0,1,0},
					  {0,1,1},{0,0,1},{2,0,0},
					  {2,1,0},{2,1,1},{2,2,0},
					  {2,0,1},{2,0,2},{0,2,0},
					  {0,2,1},{1,2,0},{0,0,2},
					  {0,1,2},{1,0,2},{0,2,2},
					  {1,1,2},{1,2,1}};

	Point<2> H;
	Point<2> sum_H;
	// normal vector of the plane F (before projection) = (n_A,n_B,n_C)
	for (unsigned int p = 0; i<22; ++i) {
		beta_1_x= alfa[p][0]+1;
		beta_1_y= alfa[p][0];
		beta_2_x= alfa[p][1];
		beta_2_y= alfa[p][1]+1;
		beta_1_z= alfa[p][0];
		beta_2_z= alfa[p][1];
		if (alfa[p][2] == 0) {
			H(0) = pow(A,beta_1_x+1)*pow(B,beta_2_x)*n_A/(beta_1_x+1)
					+pow(A,beta_1_y+1)*pow(B,beta_2_y)*n_B/(beta_1_y+1)
					+pow(A,beta_1_z+2)*pow(B,beta_2_z)*n_A*n_C/(beta_1_z+2)
					+pow(A,beta_1_z+1)*pow(B,beta_2_z+1)*n_B*n_C/(beta_1_z+1)
					+pow(A,beta_1_z+1)*pow(B,beta_2_z)*n_C/(beta_1_z+1);
			H(1) = H(2) = 0;
		}
		else if (alfa[p][2] == 1) {

			H(0) = pow(A,beta_1_x+2)*pow(B,beta_2_x)*n_A*n_A/(beta_1_x+2)
				+pow(A,beta_1_x+1)*pow(B,beta_2_x+1)*n_A*n_B/(beta_1_x+1)
				+pow(A,beta_1_x+1)*pow(B,beta_2_x)*n_A*w/(beta_1_x+1)
				+pow(A,beta_1_y+2)*pow(B,beta_2_y)*n_A*n_B/(beta_1_y+2)
				+pow(A,beta_1_y+1)*pow(B,beta_2_y+1)*n_B*n_B/(beta_1_y+1)
				+pow(A,beta_1_y+1)*pow(B,beta_2_y)*n_B*w/(beta_1_y+1)
				+pow(A,beta_1_z+3)*pow(B,beta_2_z)*n_A*n_A*n_C/(beta_1_z+3)
				+pow(A,beta_1_z+2)*( pow(B,beta_2_z+1)*2*n_A*n_B*n_C
						+ pow(B,beta_2_z)*2*n_A*n_C*w ) /(beta_1_z+2)
				+pow(A,beta_1_z+1)*( pow(B,beta_2_z+2)*n_B*n_B*n_C
						+ pow(B,beta_2_z)*w*w*n_C+ 2*pow(B,beta_2_z+1)*2*n_B*n_C*w ) /(beta_1_z+1);
			H(1) = H(2) = 0;
		}
		else if (alfa[p][2] == 2) {
			H(0) = 	n_A*pow(A,beta_1_x+3)*pow(B,beta_2_x)*n_A*n_A/(beta_1_x+3)+
					n_A*pow(A,beta_1_x+2)*( pow(B,beta_2_x+1)*2*n_A*n_B
							+ pow(B,beta_2_x)*2*n_A*w )/(beta_1_x+2)
					+n_A*pow(A,beta_1_x+1)*( pow(B,beta_2_x+2)*n_B*n_B +
									pow(B,beta_2_x)*w*w + pow(B,beta_2_x+1)*2*n_B*w )/(beta_1_x+1)
					+// y(beta)-part
					n_B*pow(A,beta_1_y+3)*pow(B,beta_2_y)*n_A*n_A/(beta_1_y+3)+
					n_B*pow(A,beta_1_y+2)*( pow(B,beta_2_y+1)*2*n_A*n_B
							+ pow(B,beta_2_y)*2*n_A*w )/(beta_1_y+2)
					+n_B*pow(A,beta_1_y+1)*( pow(B,beta_2_y+2)*n_B*n_B +
									pow(B,beta_2_y)*w*w + pow(B,beta_2_y+1)*2*n_B*w )/(beta_1_y+1)
					+ // z(gamma)-part
					n_C*pow(A,beta_1_z+4)*pow(B,beta_2_z)*n_A*n_A*n_A/(beta_1_y+4)
					+n_C*pow(A,beta_1_z+3)*( 3*n_A*n_B*pow(B,beta_2_z+1)
							+ pow(B,beta_2_z)*n_A*n_A*w )/(beta_1_y+3)
					+n_C*pow(A,beta_1_z+2)*(
							3*w*w*n_A*pow(B,beta_2_z) +
							+3*n_A*pow(B,beta_2_z+1) +
							+6*n_A*n_B*w*pow(B,beta_2_z+1) ) /(beta_1_z+2)
					+n_C*pow(A,beta_1_z+1)*(
							n_A*n_A*n_A*pow(B,beta_2_z+3) +
							+w*w*w*pow(B,beta_2_z) +
							+3*n_B*n_B*pow(B,beta_2_z+2)
							+3*n_B*w*w*pow(B,beta_2_z+1) ) /(beta_1_z+1) ;
		} // end if alfa 3 == ...
//		H*=k_p[p];
		sum_H+=H*k_p[p];
	} // end loop over p=0,22;
	return sum_H;
}

double cut_cell_3D::CompVolumeIntegral ()
{
	double nx, ny, nz;

	xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
	yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
	std::vector<Point<2> > row_sum(4);
	Point<2> new_Xt;
	int jacobian_determinant = 1;
	for (unsigned int l = 0;l<3;++l) {
		for (unsigned int m = 0;m<3;++m) {
			for (unsigned int n = 0;n<3;++n) {

//				for (unsigned int p = 0;p<22;++p) { // Already into integral_H... important because
				// the integral formula changes depending on the p! (alfa[][2], to be exact)
				// Get inside a face of the polyhedron!
					for (unsigned int it_face = 0;it_face < p->numFaces; ++it_face) {
						// Get information about this face belonging to polyhedron p.
						f = &p->faces[it_face];

						// Projection: choose where the face will be projected; maximize nz (n_gamma)
					    nx = fabs(f->norm[X]);
					    ny = fabs(f->norm[Y]);
					    nz = fabs(f->norm[Z]);			// A => alfa, B=> beta, C=>gamma
					    if (nx > ny && nx > nz) C = X;  // X = 0, Y = 1, Z = 3 (defined in the beginning)
					    else C = (ny > nz) ? Y : Z;
					    A = (C + 1) % 3;
					    B = (A + 1) % 3;

					    // Get inside a LINE of the FACE!
						for (unsigned int it_line = 0; it_line<f->_lines;++it_line) {
//							&line = plane(it_line);
							l = &f->lines[it_line];
							l->X0 =
							// Relevant info: line_normal, X0 and X1 (unit coordinates of the line
							// on plane)
							for (unsigned int q_point = 0;q_point<n_q_points /*line quad.,2*/ ;++q_point) {
								xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
								yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
								Point<2> H_t;
								H_t = GetHVector(xt,yt,m,n);
//								new_Xt[1] = for the time being, equal to 0, I may change later
								// change these line_... to something like line.normal_vector(),
								// line.length(),  etc.
								sum_integral += new_Xt*l->norm;
								sum_integral *= l->length*face_quadrature_weights[q_point];
							}
						}
					}
//				 	sum_integral_p += k_p[p]*
//				}
				integral=*J(l,n)*J(l,m);
				total_sum +=integral;
			}
		}
	}
	total_sum *= jacobian_determinant;



			// First row, representing d phi_i/dx * d phi_j/dx
			new_Xt[0] = pow(xt,(exponents_matrix_x(i,j)+1))/(2*(1+exponents_matrix_x(i,j)))
																						* pow(yt,exponents_matrix_y(i,j));
			new_Xt[1] = pow(yt,(exponents_matrix_y(i,j)+1))/(2*(1+exponents_matrix_y(i,j)))
																						* pow(xt,exponents_matrix_x(i,j));
			new_Xt *= new_matrix_coefficients(i,j);
			row_sum[i] += new_Xt;
			// Now just need to multiply each row by its specific
			// Jacobian (J11*J11,etc) and add all together. This will result in a point.
			//Then, multiply this point by the normal vector; by the det(J); by the
			// length of the face; The result in the complete face integration. Repeat
			// for each face and add all together to result in the Aij integration.
		}

	Point<2> face_integration;
	double face_integration_scalar = 0;
	for (unsigned int i = 0;i<4;++i)
	{
		face_integration += (jacobian_vector(i)+jacobian_vector(i+4))
													*row_sum[i]; // scalar*point = point
	}

	face_integration_scalar = face_integration*normal; // point*point = scalar
	face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];

	face_integration_scalar *= face_length;
	//		std::cout << "face_length2: " << face_length << "\n";
	//		face_length is ok;
	// X0 and X1 are ok;
	// normal is ok;
	// face_quadrature_weights are ok;
	// jacobian_determinant is ok;
	// This is the final form for the integration of this specific face.
	// Now all face_integration are added for each face to yield Aij.
	final_face_integration += face_integration_scalar;





//	//Create the Jacobian vector, that will multiply each instance of d phi_i/dx, etc.
//	Vector<double> jacobian_vector(8);
//	jacobian_vector[0] = jacobian_inverse(0,0)*jacobian_inverse(0,0); // J11*J11
//	jacobian_vector[1] = jacobian_inverse(0,0)*jacobian_inverse(0,1); // J11*J12
//	jacobian_vector[2] = jacobian_inverse(0,1)*jacobian_inverse(0,0); // J12*J11
//	jacobian_vector[3] = jacobian_inverse(0,1)*jacobian_inverse(0,1); // J12*J12
//	jacobian_vector[4] = jacobian_inverse(1,0)*jacobian_inverse(1,0); // J21*J21
//	jacobian_vector[5] = jacobian_inverse(1,0)*jacobian_inverse(1,1); // J21*J22
//	jacobian_vector[6] = jacobian_inverse(1,1)*jacobian_inverse(1,0); // J22*J21
//	jacobian_vector[7] = jacobian_inverse(1,1)*jacobian_inverse(1,1); // J22*J22
//
//	// Matrices representing the exponents that raise the (x,y) term in the polynomials.
//	// Here, I created two 4x4 matrices, with the exponents of x and y each.
//	// Maybe it is better to implement this as an object of the type Table <2 <Point<2>>
//
//	FullMatrix<double> exponents_matrix_x(4,4);
//	FullMatrix<double> exponents_matrix_y(4,4);
//
//	exponents_matrix_x(1,1) = 1;
//	exponents_matrix_x(1,3) = 1;
//	exponents_matrix_x(2,1) = 1;
//	exponents_matrix_x(2,3) = 1;
//	exponents_matrix_x(3,1) = 1;
//	exponents_matrix_x(3,3) = 2;
//
//	exponents_matrix_y(0,2) = 1;
//	exponents_matrix_y(0,3) = 2;
//	exponents_matrix_y(1,2) = 1;
//	exponents_matrix_y(1,3) = 1;
//	exponents_matrix_y(2,2) = 1;
//	exponents_matrix_y(2,3) = 1;
//
//	// This is the Matrix representing the coefficients that multiply the terms x,y
//	// in the new face integration formula. This will change for each dof, therefore
//	// it needs to be evaluated every change in dof (i,j). They depend on dof_i,dof_j
//	FullMatrix<double> new_matrix_coefficients(4,4);
//
//	// Need to find some way to define this better.
//	// Column 0
//	new_matrix_coefficients(0,0) = coefficients(dof_i,1)*coefficients(dof_j,1);
//	new_matrix_coefficients(1,0) = coefficients(dof_i,1)*coefficients(dof_j,2);
//	new_matrix_coefficients(2,0) = coefficients(dof_i,2)*coefficients(dof_j,1);
//	new_matrix_coefficients(3,0) = coefficients(dof_i,2)*coefficients(dof_j,2);
//	// Column 1 etc.
//	new_matrix_coefficients(0,1) = 0;
//	new_matrix_coefficients(1,1) = coefficients(dof_i,1)*coefficients(dof_j,3);
//	new_matrix_coefficients(2,1) = coefficients(dof_i,3)*coefficients(dof_j,1);
//	new_matrix_coefficients(3,1) = coefficients(dof_i,2)*coefficients(dof_j,3)
//																		+ coefficients(dof_i,3)*coefficients(dof_j,2);
//	new_matrix_coefficients(0,2) = coefficients(dof_i,1)*coefficients(dof_j,3)
//																		+ coefficients(dof_i,3)*coefficients(dof_j,1);
//	new_matrix_coefficients(1,2) = coefficients(dof_i,3)*coefficients(dof_j,2);
//	new_matrix_coefficients(2,2) = coefficients(dof_i,2)*coefficients(dof_j,3);
//	new_matrix_coefficients(3,2) = 0;
//
//	new_matrix_coefficients(0,3) = coefficients(dof_i,3)*coefficients(dof_j,3);
//	new_matrix_coefficients(1,3) = coefficients(dof_i,3)*coefficients(dof_j,3);
//	new_matrix_coefficients(2,3) = coefficients(dof_i,3)*coefficients(dof_j,3);
//	new_matrix_coefficients(3,3) = coefficients(dof_i,3)*coefficients(dof_j,3);
//
//
//	double xt;
//	double yt;
//	double final_face_integration = 0;
//	for (unsigned int q_point = 0;q_point<face_quadrature_points.size();++q_point)
//	{
//		// I assume that the points are in the right direction, ie, X0->X1,
//		// due to the change in return_normal_vector.
//		xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//		std::vector<Point<2> > row_sum(4);
//		Point<2> new_Xt;
//		for (unsigned int i = 0;i<4;++i)
//			for (unsigned int j = 0;j<4;++j)
//			{
//				// First row, representing d phi_i/dx * d phi_j/dx
//				new_Xt[0] = pow(xt,(exponents_matrix_x(i,j)+1))/(2*(1+exponents_matrix_x(i,j)))
//																			* pow(yt,exponents_matrix_y(i,j));
//				new_Xt[1] = pow(yt,(exponents_matrix_y(i,j)+1))/(2*(1+exponents_matrix_y(i,j)))
//																			* pow(xt,exponents_matrix_x(i,j));
//				new_Xt *= new_matrix_coefficients(i,j);
//				row_sum[i] += new_Xt;
//				// Now just need to multiply each row by its specific
//				// Jacobian (J11*J11,etc) and add all together. This will result in a point.
//				//Then, multiply this point by the normal vector; by the det(J); by the
//				// length of the face; The result in the complete face integration. Repeat
//				// for each face and add all together to result in the Aij integration.
//			}
//		Point<2> face_integration;
//		double face_integration_scalar = 0;
//		for (unsigned int i = 0;i<4;++i)
//		{
//			face_integration += (jacobian_vector(i)+jacobian_vector(i+4))
//										*row_sum[i]; // scalar*point = point
//		}
//
//		face_integration_scalar = face_integration*normal; // point*point = scalar (Is this the correct form?)
//		face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];
//
//		face_integration_scalar *= face_length;
//		//		std::cout << "face_length2: " << face_length << "\n";
//		//		face_length is ok;
//		// X0 and X1 are ok;
//		// normal is ok;
//		// face_quadrature_weights are ok;
//		// jacobian_determinant is ok;
//		// This is the final form for the integration of this specific face.
//		// Now all face_integration are added for each face to yield Aij.
//		final_face_integration += face_integration_scalar;
//
//	} // end quadrature sum
//
//	return final_face_integration;
		return 0.0;
}

double cut_cell_3D::getTermC (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const Point<2> &face_normal_vector,
		const int dof_i, const int dof_j)
// Checked if input are really from the new face integration; OK
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//		//		std::cout << "new_Xt: " << new_Xt << "\n";
//		// This will be J⁻¹*fe.shape_grad(dof_i j ,new_Xt)
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		Vector<double> J_fe_shape_grad_j(2);
//		tmp = 0;
//		fe1.shape_grad(dof_j,new_Xt).unroll(tmp);
//		jacobian_inverse.vmult(J_fe_shape_grad_j,tmp);
//		// Checked, multiplication is ok
//		face_integration +=
//				-((J_fe_shape_grad_i[0]*face_normal_vector[0]
//				                                           + J_fe_shape_grad_i[1]*face_normal_vector[1])
//						*fe1.shape_value(dof_j,new_Xt)*face_quadrature_weights[q_point]
//						                                                       *face_length
//			          +
//			          (J_fe_shape_grad_j[0]*face_normal_vector[0]
//			                         + J_fe_shape_grad_j[1]*face_normal_vector[1])
//
//			     *fe1.shape_value(dof_i,new_Xt)*face_quadrature_weights[q_point]*face_length
//				);
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::getTermD (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const int dof_i, const int dof_j, const double alfa)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//				(  fe1.shape_value(dof_i,new_Xt)
//						*fe1.shape_value(dof_j,new_Xt)
//						*face_quadrature_weights[q_point]  )
//						*alfa*face_length ;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

// Before, I used a getTermJ where every local_dof variable was multiplied by two, and it could be used
// only for the ubulk assembly. This one can assemble both usurface and ubulk, with the difference that
// the local_dof_K and local_dof_K_neighbor must be set as if they were a 8 dof cell, ie,
// EXTENDED_local_dof_K   = {0,2,4,6,-1,-1} for usurface
// EXTENDED_local_dof_K_p = {1,3,5,7,-1,-1} for ubulk
// same for neighbor.

double cut_cell_3D::getTermJ (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor,const FEValuesExtractors::Scalar uvariable)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell[uvariable].gradient /*.shape_grad*/(/*2**/local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell[uvariable].gradient/* .shape_grad*/ (/*2**/local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF i");
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell[uvariable].gradient(/*2**/local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell[uvariable].gradient(/*2**/local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF j");
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}

//	return j;
	return 0;

}
// Calculate J term when the neighbor cell is inside
double cut_cell_3D::getTermJ_mixed (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	const FEValuesExtractors::Scalar pressure /*usurface */(1);
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF i");
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF j");
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}
//
//	return j;
	return 0;

}

double cut_cell_3D::getTermJ_OneVar (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values.shape_grad(local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values.shape_grad(local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else (assert(0));
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values.shape_grad(local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values.shape_grad(local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else (assert(0));
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}
//	return j;
	return 0;

}

double cut_cell_3D::getTermD2 (const Point<2> &X0, const Point<2> &X1,
		const Point<2> &face_normal_vector,	const int dof_i, const double alfa,
		const double g_D, const double face_length)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//		//		std::cout << "new_Xt: " << new_Xt << "\n";
//		// This will be J⁻¹*fe.shape_grad(dof_i j ,new_Xt)
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		tmp = 0;
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad_i = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		// Checked, multiplication is ok
//
//		face_integration +=
//				(	alfa*fe1.shape_value(dof_i,new_Xt)
//						-
//						(J_fe_shape_grad_i[0]*face_normal_vector[0]
//					   + J_fe_shape_grad_i[1]*face_normal_vector[1])	)
//						*g_D*face_quadrature_weights[q_point]*face_length;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::mass_matrix (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,
		const int dof_i, const int dof_j, const double face_length)
{
//	std::vector<int> exponents_matrix_x =  {0,1,2,0,0,1,2,1,2 };
//	std::vector<int> exponents_matrix_y =  {0,0,0,1,2,1,1,2,2 };
//
//	std::vector<double> new_coefficients =
//	{
//			coefficients(dof_i,0)*coefficients(dof_j,0),
//			coefficients(dof_i,0)*coefficients(dof_j,1)+coefficients(dof_i,1)*coefficients(dof_j,0),
//			coefficients(dof_i,1)*coefficients(dof_j,1),
//			coefficients(dof_i,0)*coefficients(dof_j,2)+coefficients(dof_i,2)*coefficients(dof_j,0),
//			coefficients(dof_i,2)*coefficients(dof_j,2),
//
//			coefficients(dof_i,0)*coefficients(dof_j,3)+coefficients(dof_i,1)*coefficients(dof_j,2) +
//			coefficients(dof_i,2)*coefficients(dof_j,1)+coefficients(dof_i,3)*coefficients(dof_j,0),
//
//			coefficients(dof_i,1)*coefficients(dof_j,3)+coefficients(dof_i,3)*coefficients(dof_j,1),
//			coefficients(dof_i,2)*coefficients(dof_j,3)+coefficients(dof_i,3)*coefficients(dof_j,2),
//			coefficients(dof_i,3)*coefficients(dof_j,3) };
//
//	double xt;
//	double yt;
//	double final_face_integration = 0;
//	for (unsigned int q_point = 0;q_point<face_quadrature_points.size();++q_point)
//	{
//		xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Point<2> row_sum;
//		double face_integration_scalar = 0;
//		for (unsigned int i = 0;i<9;++i) {
//			Point<2> new_Xt;
//			// First row, representing d phi_i/dx * d phi_j/dx
//			new_Xt[0] = pow(xt,(exponents_matrix_x[i]+1))/(2*(1+exponents_matrix_x[i]))
//									* pow(yt,exponents_matrix_y[i]);
//			new_Xt[1] = pow(yt,(exponents_matrix_y[i]+1))/(2*(1+exponents_matrix_y[i]))
//									* pow(xt,exponents_matrix_x[i]);
//
//			new_Xt *= new_coefficients[i];
//			row_sum += new_Xt;
//		}
//
//		face_integration_scalar = row_sum*normal; // point*point = scalar (Is this the correct form?)
//		face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];
//		face_integration_scalar *= face_length;
//		final_face_integration += face_integration_scalar;
//	} // end quadrature sum
//	return final_face_integration;
	return 0;
}


double cut_cell_3D::getTermConstraintBoundary (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const int dof_i, const int dof_j)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//				(  fe1.shape_value(dof_i,new_Xt)
//						*face_quadrature_weights[q_point]  )
//						*face_length ;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}


double cut_cell_3D::getTermBeltramiBoundary (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
		const double face_length)
{
//	FullMatrix<double> identity_matrix (IdentityMatrix(2));
//	FullMatrix<double> beltrami_operator(2,2);
//	Point<2> beltrami_0;
//	Point<2> beltrami_1;
//	beltrami_operator(0,0) = -normal[0] * normal[0];
//	beltrami_operator(0,1) = -normal[1] * normal[0];
//	beltrami_operator(1,0) = -normal[0] * normal[1];
//	beltrami_operator(1,1) = -normal[1] * normal[1];
//	beltrami_operator.add(1,identity_matrix);
//
//	beltrami_0(0) = beltrami_operator(0,0);
//	beltrami_0(1) = beltrami_operator(0,1);
//	beltrami_1(0) = beltrami_operator(1,0);
//	beltrami_1(1) = beltrami_operator(1,1);
//
//
//
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		Vector<double> J_fe_shape_grad_j(2);
//		tmp = 0;
//		fe1.shape_grad(dof_j,new_Xt).unroll(tmp);
//		jacobian_inverse.vmult(J_fe_shape_grad_j,tmp);
//
//
////		 Ignoring Beltrami operator ("=0 for 2D")
//
////		face_integration +=
////						(       (  J_fe_shape_grad_i[0]
////								 * J_fe_shape_grad_j[0] )
////								 +
////								 ( J_fe_shape_grad_i[1]
////		                         * J_fe_shape_grad_j[1] ) )
////				*face_quadrature_weights[q_point]*face_length;
//
//
////		Considering v2 of Beltrami Operator. Here n(x)n is not equal to zero. The results obviously differ
////		from the above.
//		tmp = 0;
//		beltrami_operator.vmult(tmp,J_fe_shape_grad_i);
//		J_fe_shape_grad_i = tmp;
//		tmp = 0;
//		beltrami_operator.vmult(tmp,J_fe_shape_grad_j);
//		J_fe_shape_grad_j = tmp;
//
//		face_integration +=
//				(       (  J_fe_shape_grad_i[0]
//				                             * J_fe_shape_grad_j[0] )
//						+
//						( J_fe_shape_grad_i[1]
//						                    * J_fe_shape_grad_j[1] ) )
//						                    *face_quadrature_weights[q_point]*face_length;
//
//
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::getTermBeltramiBoundaryRHS (const Point<2> &X0,
		const Point<2> &X1,	const int dof_i,
		const double face_length, const double fs)
{
//	double face_integration = 0;
//		for (int q_point = 0;q_point<n_face_q_points;++q_point)
//		{
//			Point<2> new_Xt;
//			new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//			new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//			face_integration +=
//					( fs*fe1.shape_value(dof_i,new_Xt) )
//					*face_quadrature_weights[q_point]*face_length;
//
//		} // end quadrature sum
//		return face_integration;
	return 0;
}

double cut_cell_3D::constraintVector (const Point<2> &X0,
		const Point<2> &X1,	const int dof_i,
		const double face_length)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=  ( fe1.shape_value(dof_i,new_Xt) )
//								*face_quadrature_weights[q_point]*face_length;
//
//	} // end quadrature sum
//	return face_integration;
	return 0;

}

double cut_cell_3D::getTermCoupling (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
		const double face_length,
		const int multiplier_u_surface_i,
		const int multiplier_u_surface_j,
		const int multiplier_u_bulk_i,
		const int multiplier_u_bulk_j, const double b_B, const double b_S)
{

//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//		 (  	b_B*b_B*
//				multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
//				 	 * multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
//				-
//				b_B*b_S*
//				multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
//					 * multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt)
//				-
//				b_B*b_S
//					 * multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt)*multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
//				+
//				b_S*b_S*
//				multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt)
//					 * multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt) )
//					*face_quadrature_weights[q_point]*face_length;
//
//		// This is equivalent to the above
////				( (multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
////						-  multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt) )
////						*
////						(multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
////								-  multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt)) )
////						*face_quadrature_weights[q_point]*face_length;
//
//	} // end quadrature sum
//	return face_integration;
	return 0;

}

