/*
 * NewCell.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#include "NewCell.h"
#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>


using namespace dealii;

//template <int 2>
//NewCell::NewCell(/*const Vector<double> & levelset*/) /*: levelset_ (levelset)*/
NewCell::NewCell(const Vector<double> & _levelset): _levelset (_levelset)
{}
NewCell::NewCell() {}

//template <int 2>
void NewCell::setCoordinates (const int i, const Point<2> &X0,const Point<2> &X1)
{
	face_index = i;
	ObjNewFace[face_index].X0 = X0;
	ObjNewFace[face_index].X1 = X1;
	ObjNewFace[face_index].face_length = distance(X0,X1);
}

void NewCell::SetIndex(hp::DoFHandler<2>::active_cell_iterator _index)
{
	cell_index = _index;
}

void NewCell::SetCutFacePoints(Point<2> X0, Point<2> X1)
{
	X0_CutFace = X0;
	X1_CutFace = X1;
}

void NewCell::SetCutFaceLength(double FaceLength)
{
	real_face_length_CutFace = FaceLength;
}

//template <int dim>
void NewCell::setVertices()
{
//	NewFace Obj_NewFace; // MIGRATED TO NEW NewCell CLASS, OUTSIDE my_includes
//	Obj_VectorNewFace.push_back(Obj_NewFace);
	if (std::find(vertices.begin(), vertices.end(), ObjNewFace[face_index].X0)
						== vertices.end())
	{
		vertices.push_back(ObjNewFace[face_index].X0);
		vertices_index++;
	}

	if (std::find(vertices.begin(), vertices.end(), ObjNewFace[face_index].X1)
						== vertices.end())
	{

		vertices.push_back(ObjNewFace[face_index].X1);
		vertices_index++;
	}

}
//template <int dim>
void NewCell ::outputVertices()
{
	std::cout << "vertices.end()[0]" << vertices.end()[0] << "\n";
	std::cout << "vertices.end()[1]" << vertices.end()[1] << "\n";
	std::cout << "vertices.size()" << vertices.size() << "\n";
	for (unsigned int i = 0; i<vertices.size(); ++i)
		std::cout << "vertices[" << i << "]: "<< vertices[i] << "\n";

}

//template <int dim>
Point <2> NewCell::getCentroid()
{
	assert(vertices.size()>2);
//	if (vertices.size()>2)
//	{
		double X_centroid = 0;
		double Y_centroid = 0;
		for (unsigned int i = 0; i<vertices.size(); ++i)
		{
			X_centroid += vertices[i][0];
			Y_centroid += vertices[i][1];
		}
		Point<2> Centroid (X_centroid/vertices.size(),Y_centroid/vertices.size());
		return Centroid;
//	}

}
Point<2> NewCell::getCutFaceNormal()
{
	// Calculates the normal of the new cut face. For the other faces, use the usual
	// fe_face_values.normal_vector(q_point)

	Point<2> cut_face_X0 = ObjNewFace[face_index].X0;
	Point<2> cut_face_X1 = ObjNewFace[face_index].X1;

	double dy = (cut_face_X1[1] - cut_face_X0[1]);
	double dx = (cut_face_X1[0] - cut_face_X0[0]);
	double sqrt_dxdy = sqrt(dx*dx+dy*dy);

	dy =  dy / sqrt_dxdy;
	dx =  dx / sqrt_dxdy;


	// If this is the right normal vector, it means that the line is going from X0 to X1!
	Point<2> n1(dy,-dx);
	// If this is the right normal vector, it means that the line is going from X1 to X0!
	// We can rearrange this Vector so that it becomes oriented from X0 to X1. This
	// will fix the parametric equation in return_face_integration.
	Point<2> n2(-dy,dx);


	Point<2> aux_n1 = (n1+(ObjNewFace[face_index].X0+ObjNewFace[face_index].X1)/2);
	Point<2> aux_n2 = (n2+(ObjNewFace[face_index].X0+ObjNewFace[face_index].X1)/2);
	if (vertices.size()>2)
		Point<2> centroid = getCentroid();
	else assert(false && "vertices.size()<= 2");

	Point<2> normal;
	Point<2> Xtemp;

	if ( distance(aux_n1, centroid) >
			distance(aux_n2, centroid) )
		normal = n1;

	else
	{
		normal= n2;
		// Reorder X0, X1
		Xtemp = ObjNewFace[face_index].X0;
		ObjNewFace[face_index].X0 = ObjNewFace[face_index].X1;
		ObjNewFace[face_index].X1 = Xtemp;
	}
	return normal;
}


double NewCell::distance (const Point<2> &X0, const Point<2> &X1){

	double dy = (X1[1] - X0[1]);
	double dx = (X1[0] - X0[0]);
	double sqrt_dxdy = sqrt(dx*dx+dy*dy);
	return sqrt_dxdy;

}

// the intersection is based on the levelset values on adjacent nodes!
// Therefore the approximation levelset = 1 outside, = 0 on the boundary and = -1
// will never work for this case.
Point<2> NewCell::getIntersection(const Point<2> &X0, const Point<2>& X1,
		// const Vector<double>  & levelset_,
		const int k0,const int k1)
//(Vector<double> const & levelset)
//		std::vector<types::global_dof_index> global_dof_index_X0,
//		std::vector<types::global_dof_index> global_dof_index_X1){
{
	double px = 0;
//	std::cout << "levelset_.size()  " << _levelset.size() << "\n";
//	std::cout << "IN k0: " << k0 << "\n";
//	std::cout << "IN k1: " << k1 << "\n";
//	std::cout << "IN levelset_[k0]: " << _levelset[k0] << "\n";
//	std::cout << "IN levelset_[k1]: " << _levelset[k1] << "\n";
	Point<2> intersection = (px-_levelset[k0])*(X1-X0)/
			(_levelset[k1]-_levelset[k0])+X0;
	return intersection;
}
