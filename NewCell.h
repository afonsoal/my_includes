/*
 * NewCell.h
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#ifndef NEWCELL_H_
#define NEWCELL_H_

#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
using namespace dealii;

//template<int dim>
class NewCell
{
private:
	int face_index;
	int vertices_index;

	struct NewFace {
		Point<2> X0;
		Point<2> X1;
		Point<2> normal_vector;
		double face_length;
//		int face_index; // MIGRATED TO NEW NewCell CLASS, OUTSIDE my_includes
	};

	std::vector<Point<2> > vertices;
	Point<2> centroid;
	Vector<double> _levelset;
	hp::DoFHandler<2>::active_cell_iterator cell_index;

public:
//	std::vector<NewFace> Obj_VectorNewFace; // MIGRATED TO NEW NewCell CLASS, OUTSIDE my_includes
	NewFace ObjNewFace[5];
	NewCell(const Vector<double> & _levelset);
	NewCell();
	void setCoordinates (const int i, const Point<2> &X0,const Point<2> &X1);
	Point<2> getCentroid();
	void setVertices();
	void outputVertices();
	Point<2> getCutFaceNormal();
	double distance (const Point<2> &X0, const Point<2>& X1);
	Point<2> getIntersection(const Point<2> &X0, const Point<2>& X1,
//			const Vector<double> & levelset_,
			const int k0,const int k1);
	void SetIndex(hp::DoFHandler<2>::active_cell_iterator _index);
	bool is_surface_cell;

	Point<2> X0_CutFace;
	Point<2> X1_CutFace;
	double real_face_length_CutFace;
	void SetCutFacePoints(Point<2> X0, Point<2> X1);
	void SetCutFaceLength(double FaceLength);
};


#endif /* NEWCELL_H_ */
