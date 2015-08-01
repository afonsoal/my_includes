/*
 * NewCell_3D.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#include "NewCell_3D.h"
#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include "NewFace_3D.h"
#include <math.h> // fabs

using namespace dealii;

//template <int 2>
//NewCell_3D::NewCell_3D(/*const Vector<double> & levelset*/) /*: levelset_ (levelset)*/
NewCell_3D::NewCell_3D(	bool _is_surface_cell)
{
	is_surface_cell = _is_surface_cell;
	if (is_surface_cell)
		number_of_faces = 0;

	real_faces_vector.reserve(GeometryInfo<3>::faces_per_cell);
}
void NewCell_3D::InputNewRealFaceInfo(int _face_it)
{
	real_faces_info obj_real_face;
	obj_real_face.face_it = _face_it;
	obj_real_face.is_stabilization_face = false /*_is_stabilization_face*/;
	real_faces_vector.push_back(obj_real_face);
}

void NewCell_3D::SetNewRealFaceStabilization ( int _face_it, bool _is_stabilization_face)
{
	real_faces_vector[_face_it].is_stabilization_face = _is_stabilization_face;
}
void NewCell_3D::InputNewFace (NewFace_3D &CutFace)
{
//	std::cout<< "Call to InputNewFace \n";
//	CutFace.face_index = number_of_faces;

	CutFace.SetFaceIndex(number_of_faces);
//	CutFace.CompFaceArea(); // moved to ReorderAllVertices
	++number_of_faces;
	Obj_VectorNewFace.push_back(CutFace);
}


void NewCell_3D::SetIndex(/*hp::*/DoFHandler<3>::active_cell_iterator _index)
{
	cell_index = _index;
}

void NewCell_3D::CompBoundaryFace(NewFace_3D &CutFaceBoundary)
{
//	std::cout << "Call to CompBoundaryFace \n";
//	std::cout << "number_of_faces: " << number_of_faces << "\n";
//	std::cout << "number_of_lines: " << CutFaceBoundary.number_of_lines << "\n";
//	NewFace_3D CutFace;
	int line_no = 0;
	bool add_1_more_line = false;
	for (int face_it = 0; face_it < number_of_faces; ++face_it) // Here the CutBoundaryFace is not counted yet.
		for (int line_it = 0; line_it <  Obj_VectorNewFace[face_it].number_of_lines/*CutFaceBoundary.number_of_lines*/; ++line_it)
	{
//			std::cout << "face_it: " << face_it << " line_it: " << line_it << " is_boundary_face: " << Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].is_boundary_face << " end" << std::endl;
			if (Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].is_boundary_line)
			{
//				 OLD: USE REAL unit X0, X1. Should use the real!
				Point<3> X0 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0;
				Point<3> X1 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1;
				// NEW: USE mapped unit X0, X1
//				Point<3> X0 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0_unit;
//				Point<3> X1 = Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1_unit;
				for (unsigned int coord = 0; coord<3; ++coord)
				{
//					std::cout << "fabs(pow(10,-17): " << fabs(pow(10,-17)) << std::endl;
//					std::cout << "(fabs(pow(10,-17) < pow(10,-10)): " << (fabs(pow(10,-17)) < pow(10,-10)) << std::endl;
					if (fabs(X0[coord]) < pow(10,-10))
					{
//						std::cout << "(fabs(X0[coord]) < pow(10,-10)) \n";
							X0[coord] = 0.0;
					}
					if (fabs(X1[coord]) < pow(10,-10))
					{
//						std::cout << "(fabs(X1[coord]) < pow(10,-10)) \n";
							X1[coord] = 0.0;
					}
				}
//				 std::cout << "CutFaceBoundary.Obj_VectorNewLine[line_it].X0: " << X0 << std::endl;
//				 std::cout << "CutFaceBoundary.Obj_VectorNewLine[line_it].X1: " << X1 << std::endl;
				// SetCoordinates inputs newLine
				// Up to this point, I already know how many lines CutFaceBoundary has, because I sum 1 line every time I find a cut line.
				// (in the main.cpp) (why??). That's why here I do not need to add_1_more_line.
				 CutFaceBoundary.SetCoordinates(line_no,X0,X1,true /*is boundary (line)*/,add_1_more_line);
				line_no++;
			}
	}
	CutFaceBoundary.is_boundary_face = true;
//	CutFaceBoundary.CompAllLineNormals();
}

void NewCell_3D ::OutputVertices()
{
	std::cout << "vertices.end()[0]" << vertices.end()[0] << "\n";
	std::cout << "vertices.end()[1]" << vertices.end()[1] << "\n";
	std::cout << "vertices.size()" << vertices.size() << "\n";
	for (unsigned int i = 0; i<vertices.size(); ++i)
		std::cout << "vertices[" << i << "]: "<< vertices[i] << "\n";
}

// Create a vector of vertices with non-repeated vertices in order to calculate the centroid of the cell.
void NewCell_3D ::OrganizeVertices()
{
	for (int face_it = 0; face_it < number_of_faces; ++face_it)
	{
		for (int line_it = 0; line_it < Obj_VectorNewFace[face_it].number_of_lines; ++line_it)
		{
			if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0)
			== vertices.end())
			{
				vertices.push_back(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X0);
				vertices_index++;
			}

			if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1)
			== vertices.end())
			{
				vertices.push_back(Obj_VectorNewFace[face_it].Obj_VectorNewLine[line_it].X1);
				vertices_index++;
			}
		}
	}
}

// Compute the geometric center of the cell.
void NewCell_3D::CompCellCentroid()
{
	OrganizeVertices();
//	assert(vertices.size()>2);
	double X_centroid = 0;
	double Y_centroid = 0;
	double Z_centroid = 0;
	for (unsigned int i = 0; i<vertices.size(); ++i)
	{
		X_centroid += vertices[i][0];
		Y_centroid += vertices[i][1];
		Z_centroid += vertices[i][2];
	}
	cell_centroid[0] = X_centroid/vertices.size();
	cell_centroid[1] = Y_centroid/vertices.size();
	cell_centroid[2] = Z_centroid/vertices.size();
}
// Need to be called AFTER setting up normal vectors!
void NewCell_3D::ReorderAllVertices()
{
	for (int face_it = 0; face_it < number_of_faces; ++face_it)
	{
		Obj_VectorNewFace[face_it].ReorderAllLineVertices();
		Obj_VectorNewFace[face_it].SortVertices();
		Obj_VectorNewFace[face_it].CompFaceArea();
	}
}

