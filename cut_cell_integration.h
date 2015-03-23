/*
 * cut_cell_integration.h
 *
 *  Created on: Nov 5, 2014
 *      Author: afonsoal*/


#ifndef CUT_CELL_INTEGRATION_H_
#define CUT_CELL_INTEGRATION_H_
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
/*If you want to declare Obj1 and Obj2 in your .h file, add extern in the .h file like so:

extern SA Obj1, Obj2;
but you should declare the objects in a .cpp file in your project:

main.cpp

SA Obj1, Obj2;

1. Add a forward declaration of library_class_1 in header.h before using it.

class library_class_1;
2. Make sure to #include the .h file that defines library_class_1 in header.cc.*/


using namespace dealii;

class cut_cell_integration
{
private:

//	FEValues<2> fe_values;
//	FEValues<2> fe_values2;
	FE_Q<2> fe1;
	QGauss<2> quadrature_formula1;
	QGauss<1> face_quadrature_formula1;
	int   dofs_per_cell;
	int   n_q_points;
	int 	n_face_q_points;
	FullMatrix<double> jacobian_inverse;
	FullMatrix<double> coefficients;
	std::vector<double > cell_quadrature_weights;
	std::vector<double > face_quadrature_weights;
	std::vector<Point<1> > face_quadrature_points;
	double jacobian_determinant;
	int dim;


//	FE_Q<2>& fe1;
//	const unsigned int   n_q_points    = quadrature_formula.size();
//	const unsigned int n_face_q_points = face_quadrature_formula.size();

public:

	cut_cell_integration(FEValues<2> const & fe_values,
			FE_Q<2> const& fe,
			QGauss<2> const &quadrature_formula,
			QGauss<1> const &face_quadrature_formula);

	double return_face_integration (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const int dof_j, const double face_length);

	double return_rhs_face_integration (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i,
			const double face_length);

	double getTermC (const Point<2> &X0, const Point<2> &X1,
			const double face_length, const Point<2> &face_normal_vector,
			const int dof_i, const int dof_j);

	double getTermD (const Point<2> &X0, const Point<2> &X1,
				const double face_length, const int dof_i,
				const int dof_j, const double alfa);

	double getTermJ (FEFaceValues<2> const & fe_face_values,
			FEFaceValues<2> const & fe_face_values_neighborCell,
			const int dof_i, const int dof_j,
			const std::vector<int> & local_dof_K,
					const std::vector<int> & local_dof_K_neighbor,
					const FEValuesExtractors::Scalar uvariable );

	double getTermD2(const Point<2> &X0, const Point<2> &X1,
			const Point<2> &face_normal_vector, const int dof_i, const double alfa,
			const double g_D, const double face_length);
	double mass_matrix (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const int dof_j, const double face_length);

	double CompConstraintUbulk (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const double face_length);

	double CompMassMatrixSurface (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i, const int dof_j,
			const double face_length);
	double getTermConstraintBoundary (const Point<2> &X0, const Point<2> &X1,
			const double face_length, const int dof_i, const int dof_j);
	double getTermBeltramiBoundary (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
			const double face_length);

	double getTermBeltramiBoundaryRHS (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i,
			const double face_length, const double fs);
	double constraintVector (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i,
			const double face_length);

	double getTermCoupling (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
			const double face_length, const int corrector_u_i, const int corrector_u_j,
			const int corrector_p_i, const int corrector_p_j, const double b_B, const double b_S);

	double getTermJ_mixed (FEFaceValues<2> const & /*NULL_*/fe_face_values,
				FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
				/*const*/ int dof_i, /*const */int dof_j,
				const std::vector<int> & local_dof_K,
				const std::vector<int> & local_dof_K_neighbor);

	double getTermJ_OneVar (FEFaceValues<2> const & /*NULL_*/fe_face_values,
			FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
			/*const*/ int dof_i, /*const */int dof_j,
			const std::vector<int> & local_dof_K,
			const std::vector<int> & local_dof_K_neighbor);

	double getTermN2 (const Point<2> &X0, const Point<2> &X1,
			const Point<2> &face_normal_vector,	const int dof_i,
			const double g_N, const double face_length, const double gamma_N, const double cell_diameter);
	double getTermNlhs (const Point<2> &X0, const Point<2> &X1,
			const Point<2> &face_normal_vector,	const int dof_i, const int dof_j,
			const double face_length, const double gamma_N, const double cell_diameter);
	double CompMatrix_kC (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i,
			std::vector<double> &u_bulk_at_nodes,const double face_length);

	double CompMassConservation_face (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const double face_length, std::vector<double> &u_bulk_at_nodes);

	double CompMassConservation_boundary (const Point<2> &X0,
			const Point<2> &X1, std::vector<double> &u_at_nodes,const double face_length);
};

#endif


