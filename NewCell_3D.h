/*
 * NewCell_3D.h
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#ifndef NEWCELL_3D_H_
#define NEWCELL_3D_H_

#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include "NewFace_3D.h"
#include <math.h> // fabs
using namespace dealii;

//template<int dim>
class NewCell_3D
{
private:
//	int face_index;
	int vertices_index;


//	struct NewFace {
//		Point<2> X0;
//		Point<2> X1;
//		Point<2> normal_vector;
//		double real_face_length;
//		double unit_face_length;
//		int face_index;
//		bool is_boundary_face;
//	};

	std::vector<Point<3> > vertices;

	Vector<double> _levelset;
	struct real_faces_info{
		int face_it;
		bool is_stabilization_face;
	};
//	std::vector<std::pair< int, bool>>real_face_info;
	void CompFaceArea();
public:
	std::vector<real_faces_info> real_faces_vector;
	Point<3> cell_centroid;
	std::vector<NewFace_3D> Obj_VectorNewFace;
//	NewFace_3D Obj_NewFace;
	bool is_surface_cell;
	/*hp::*/DoFHandler<3>::active_cell_iterator cell_index;
	int number_of_faces;
	void InputNewRealFaceInfo (int _face_it);
	void SetNewRealFaceStabilization ( int _face_it, bool _is_stabilization_face);
//	void InputNewFace();
	void InputNewFace (NewFace_3D &CutFace);
	NewCell_3D(bool _is_surface_cell);
	void SetIndex(/*hp::*/DoFHandler<3>::active_cell_iterator _index);
	void CompBoundaryFace(NewFace_3D &CutFaceBoundary);
	///// NEED TO REDO.


	void OutputVertices();
	void CompCellCentroid();
	void ReorderAllVertices();
	void OrganizeVertices();

	void Compute(); // Change name later

};


#endif /* NEWCELL_H_ */
