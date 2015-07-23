/*
 * NewFace_3D.cpp
 *
 *  Created on: Nov 16, 2014
 *      Author: afonsoal
 */

#include "NewFace_3D.h"
#include <fstream>
#include <iostream>
//#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <fenv.h> // Catch NaN


using namespace dealii;

//template <int 2>
//Constructor for new faces formed by cutting old faces, not the Boundary face itself (which is completely new)
NewFace_3D::NewFace_3D (const Vector<double> & _levelset, const MappingQ1<3> & _mapping, DoFHandler<3>::active_cell_iterator _cell_index)
:
		levelset (_levelset),
		mapping (_mapping),
		cell_index (_cell_index)

{
	is_boundary_face = false;
	number_of_lines = 0;
}
//Constructor of the Cut Boundary face.
NewFace_3D::NewFace_3D(const MappingQ1<3> & _mapping, DoFHandler<3>::active_cell_iterator _cell_index)
	:
		mapping (_mapping),
		cell_index (_cell_index)
{
	is_boundary_face = true;
	number_of_lines = 0;
}

// Set coordinates of the line. Input information of the line index, if is a boundary line and the real face length of the face. The unit face length is input separately.
void NewFace_3D::SetCoordinates (const int _line_index, const Point<3> &X0,const Point<3> &X1, bool _is_boundary_line, bool _sum_1_more_line)
{
	InputNewLine(_sum_1_more_line);
	//	face_no = i;
	Obj_VectorNewLine[_line_index].line_index =  _line_index;

//	_is_boundary_line identifies that this is the new cut LINE.
	Obj_VectorNewLine[_line_index].is_boundary_line = _is_boundary_line;

//	line_index = i;
	Obj_VectorNewLine[_line_index].X0 = X0;
	Obj_VectorNewLine[_line_index].X1 = X1;
	Obj_VectorNewLine[_line_index].real_face_length = distance(X0,X1);

	SetVertices(_line_index);
	SetUnitLineLength(_line_index,X0,X1);
	// Can only compute the normal vector after having input all the faces... and only the cut face!
//	if (_is_boundary_line) {
//		std::cout << "face_no: " << face_no << std::endl;
//		CompCutFaceNormal();s
//	}
}
void NewFace_3D::InputNewLine (bool _sum_1_more_line)
{
	if (_sum_1_more_line)
		++number_of_lines;

	Obj_VectorNewLine.push_back(Obj_NewLine);
}
void NewFace_3D::Add1Line()
{
//	std::cout << "ADD 1 LINE TO CutFaceBoundary \n";
	++number_of_lines;
}
void NewFace_3D::SetFaceIndex(int _face_index)
{
	face_index = _face_index;
}

void NewFace_3D::SetVertices(const int _line_index)
{
	if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewLine[_line_index].X0)
	== vertices.end())
	{
		vertices.push_back(Obj_VectorNewLine[_line_index].X0);
		vertices_index++;
	}

	if (std::find(vertices.begin(), vertices.end(), Obj_VectorNewLine[_line_index].X1)
	== vertices.end())
	{
		vertices.push_back(Obj_VectorNewLine[_line_index].X1);
		vertices_index++;
	}

}
void NewFace_3D ::OutputVertices()
{
	std::cout << "vertices.end()[0]" << vertices.end()[0] << "\n";
	std::cout << "vertices.end()[1]" << vertices.end()[1] << "\n";
	std::cout << "vertices.size()" << vertices.size() << "\n";
	for (unsigned int i = 0; i<vertices.size(); ++i)
		std::cout << "vertices[" << i << "]: "<< vertices[i] << "\n";
}

/*Point <3>*/void NewFace_3D::CompFaceCentroid()
{
	assert(vertices.size()>2);
		double X_centroid = 0;
		double Y_centroid = 0;
		double Z_centroid = 0;
		for (unsigned int i = 0; i<vertices.size(); ++i)
		{
			X_centroid += vertices[i][0];
			Y_centroid += vertices[i][1];
			Z_centroid += vertices[i][2];
		}
		Point<3> Centroid (X_centroid/vertices.size(),Y_centroid/vertices.size(),Z_centroid/vertices.size());
//		return Centroid;
		// I am not sure if I should do this, do I lose any accuracy with this?
//		for (unsigned int coord = 0; coord<3; ++coord)
//			if (fabs(Centroid[coord]) < pow(10,-10))
//				Centroid[coord] = 0.0;
		face_centroid = Centroid;

}

double NewFace_3D::distance (const Point<3> &X0, const Point<3> &X1){

	double dx = (X1[0] - X0[0]);
	double dy = (X1[1] - X0[1]);
	double dz = (X1[2] - X0[2]);
	double sqrt_dxdydz = sqrt(dx*dx+dy*dy+dz*dz);

	return sqrt_dxdydz;

}

// the intersection is based on the levelset values on adjacent nodes!
Point<3> NewFace_3D::GetIntersection(const Point<3> &X0, const Point<3>& X1,const int k0,const int k1)
{
	double px = 0;
	Point<3> intersection = (px-levelset[k0])*(X1-X0)/
			(levelset[k1]-levelset[k0])+X0;
	return intersection;
}

void NewFace_3D:: SetFaceNormal(Point<3> normal)
{
	face_normal_vector = normal;
}

// Calculate the normal vector of the CUT FACE.
//	 Calculates the normal of the new cut FACe. For the other faces, use the usual
//	// fe_face_values.normal_vector(q_point)
void NewFace_3D:: CompCutFaceNormal(Point<3> cell_centroid)
{
//	The face_index of the CutFace is always (number of faces - 1), because it is the last face to be input.

//	The normal vector from the cut face can come from the cross product of any two vectors of the polygon describing the cut face.
// These two vectors will come from three different points of the polygon.
//  The result is two possible vectors, one pointing inwards and the other outwards the cell (this is the one I want)
//	First, select any two lines (here [0] and [1]) and create two vectors

	// NOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPENOPE
//	This is a bad idea that cost me a couple of hours! Choosing any two lines can result in taking parallel lines, so that
//	the cross product is zero!!
	Point<3> vector_0, vector_1, normal_vector_0;
	for (unsigned  int line_it = 0; line_it < number_of_lines-1; ++line_it)
	{
		/*Point<3> */vector_0 = Obj_VectorNewLine[line_it].X0 - Obj_VectorNewLine[line_it].X1;
		/*Point<3> */vector_1 = Obj_VectorNewLine[line_it+1].X0 - Obj_VectorNewLine[line_it+1].X1;
		/*Point<3> */normal_vector_0 = CrossProduct(vector_0,vector_1);
		double length = sqrt(normal_vector_0[0]*normal_vector_0[0]+normal_vector_0[1]*normal_vector_0[1]+normal_vector_0[2]*normal_vector_0[2]);
		if (length-pow(10,-10) > 0.0)
			break;
	}

	// If this is the right normal vector, it means that the line is going from P0 to P1 to P2
//	/*Point<3> */normal_vector_0 = CrossProduct(vector_0,vector_1);
	double length = sqrt(normal_vector_0[0]*normal_vector_0[0]+normal_vector_0[1]*normal_vector_0[1]+normal_vector_0[2]*normal_vector_0[2]);
	normal_vector_0 = normal_vector_0/length;

	// If this is the right normal vector, it means that the line is going from P0 to P2 to P1!
	Point<3> normal_vector_1 = -normal_vector_0;

	// aux_nX represents a Point given by the vector normal_vector_X starting on the FACE centroid.
	// For the normal_vector_X to be the real normal_vector, this Point aux_nX needs to be further from the CELL centroid than its counter part, which points inwards the cell.

	Point<3> aux_n0 =  normal_vector_0+face_centroid;
	Point<3> aux_n1 =  normal_vector_1+face_centroid;

	Point<3> normal;

	if ( distance(aux_n0, cell_centroid) > distance(aux_n1, cell_centroid) )
		normal = normal_vector_0;

	else
		normal= normal_vector_1;
	if (!(normal == normal))
	{
		std::cout << "vector_0: " << vector_0 << std::endl;
		std::cout << "vector_1: " << vector_1 << std::endl;
		std::cout << "length: " << length << std::endl;
		std::cout << "normal_vector_0: " << normal_vector_0 << std::endl;
		std::cout << "normal_vector_1 " << normal_vector_1 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X0: " << Obj_VectorNewLine[0].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[0].X1: " << Obj_VectorNewLine[0].X1 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X0: " << Obj_VectorNewLine[1].X0 << std::endl;
		std::cout << "Obj_VectorNewLine[1].X1: " << Obj_VectorNewLine[1].X1 << std::endl;
		std::cout << "NORMAL: " << normal << std::endl;
		std::cout << "face_index: " << face_index << std::endl;
		std::cout << "is_boundary_face: " << is_boundary_face << std::endl;
		std::cout << "face_centroid: " << face_centroid << std::endl;
		std::cout << "cell_centroid: " << cell_centroid << std::endl;
		assert(normal == normal); // Catch NaN
	}

	SetFaceNormal(normal);
//	face_normal_vector = normal;
}

// Take the cross product of two vectors a and b.
Point<3> NewFace_3D::CrossProduct(Point<3> a, Point<3> b)
{
	Point<3> Product;

        //Cross product formula
	Product[0] = (a[1] * b[2]) - (a[2] * b[1]);
	Product[1] = (a[2] * b[0]) - (a[0] * b[2]);
	Product[2] = (a[0] * b[1]) - (a[1] * b[0]);
	return Product;
}

/////////////////
// Feed it a line (A<->B) and a point, and it will give you the nearest point on the line to the target point.
// A and B are two points which define the line, Centroid is "the point."
// This can be called only after inputting all the lines of the cut face, because I need the Face Centroid.
// Adapted from http://arstechnica.com/civis/viewtopic.php?t=149128
Point<3> NewFace_3D::CompLineNormal (Point<3> A, Point<3> B	)
{

//double u = ((face_centroid[0]-A[0])*(B[0]-B[1])) + ((face_centroid[1] - A[1]) * (B[1] - A[1])) + ((face_centroid[2] - A[2]) * (B[2] - A[2]));
//var u = ((pX.x - p1.x) * (p2.x - p1.x)) + ((pX.y - p1.y) * (p2.y - p1.y)) + ((pX.z - p1.z) * (p2.z - p1.z))
double u = ((face_centroid[0] - A[0]) * (B[0] - A[0])) + ((face_centroid[1] - A[1]) * (B[1] - A[1])) + ((face_centroid[2] - A[2]) * (B[2] - A[2]));
double dist = distance(A,B);
u = u/(dist*dist);

Point<3> intersection;
// This is the point of intersection of the line AB with the normal line passing through the face_centroid.
intersection[0] = A[0] + u * (B[0] - A[0]);
intersection[1] = A[1] + u * (B[1] - A[1]);
intersection[2] = A[2] + u * (B[2] - A[2]);

Point<3> normal_vector;
//The returning point and the original point will define your line. It will be perpendicular by merit of being the shortest line from the point to the target line.
normal_vector = (intersection-face_centroid)/distance(intersection,face_centroid);

for (unsigned int coord = 0; coord<3; ++coord)
	if (fabs(normal_vector[coord]) < pow(10,-10))
		normal_vector[coord] = 0.0;

return normal_vector;
}

void NewFace_3D::CompAllLineNormals ()
{
	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
		Obj_VectorNewLine[line_it].normal_vector = CompLineNormal(Obj_VectorNewLine[line_it].X0,Obj_VectorNewLine[line_it].X1);
}

void NewFace_3D::ReorderAllLineVertices()
{
	for (int line_it = 0; line_it < number_of_lines ; ++line_it)
	{
		Point<3> X0 = Obj_VectorNewLine[line_it].X0;
		Point<3> X1 = Obj_VectorNewLine[line_it].X1;
		Point<3> X0_X1 = (X1 - X0)/distance(X0,X1); // Vector going from X0 to X1 (X0->X1)
		Point<3> CrossProduct_X0X1 = CrossProduct(X0_X1,Obj_VectorNewLine[line_it].normal_vector);
		CrossProduct_X0X1 = CrossProduct_X0X1 / distance(CrossProduct_X0X1, Point<3>(0,0,0));
//		if (CrossProduct(X1_X0,Obj_VectorNewLine[line_it].normal_vector) ==  face_normal_vector)
//
		bool change_order = true;
//		This is equivalent to
//		!(CrossProduct_X1X0 == face_normal_vector), but considering a tolerance.
//		If at least one coordinate of CrossProduct_X1X0 is different than its counterpart on normal_vector, they are different.
				for (unsigned int coord = 0; coord<3; ++coord)
					if ( !( (fabs(CrossProduct_X0X1[coord] - face_normal_vector[coord] ) < pow(10,-10) ) ) )
						change_order = false;

		if (change_order)
//			if (0)
		{
//			std::cout << "Change order of points! \n";
//			std::cout << "cell_index: \n" << cell_index << std::endl;
//			std::cout << "Obj_VectorNewLine[line_it].global_line_index: " << Obj_VectorNewLine[line_it].global_line_index << std::endl;
			Obj_VectorNewLine[line_it].X0 = X1;
			Obj_VectorNewLine[line_it].X1 = X0;
		}
	}
}

//Set the unit length of the line. I could do this automatically inside NewFace_3D class if I input the mapping given by deal.ii, but I'd rather not for now (why? just to save "space"?)
// As an alternative, I could call this function from SetCoordinates, and ask the user to input unit_X0_1 in SetCoordinates function call.
void NewFace_3D::SetUnitLineLength(const int _line_index,const Point<3> & real_X0,const Point<3> &real_X1)
{
	Obj_VectorNewLine[_line_index].unit_face_length = distance (mapping.transform_real_to_unit_cell(cell_index,real_X0),mapping.transform_real_to_unit_cell(cell_index,real_X1));
}

void NewFace_3D::SetGlobalLineIndex(int _line_index, int _global_line_index)
{
	Obj_VectorNewLine[_line_index].global_line_index = _global_line_index;
}

/////////////
// Compute projection variables. Based on SetPolyhedron.cpp .h and on Mirtich's work (1999)
void NewFace_3D::CompProjectionVars() // Change name later
{
//	feenableexcept(FE_INVALID | FE_OVERFLOW); // Catch NaN
	CompProjectionVars_was_called = true;
	double nx = fabs(face_normal_vector[X]);
	double ny = fabs(face_normal_vector[Y]);
	double nz = fabs(face_normal_vector[Z]);			// A => alfa, B=> beta, C=>gamma
	if (nx > ny && nx > nz) gamma_ = X;  // X = 0, Y = 1, Z = 3 (defined in the beginning)
	else gamma_ = (ny > nz) ? Y : Z;

	alfa = (gamma_ + 1) % 3;
	beta = (alfa + 1) % 3;

	  w= - face_normal_vector[X] * vertices [0][X]
	         - face_normal_vector[Y] * vertices[0][Y]
	         - face_normal_vector[Z] * vertices[0][Z];

		Point<2> centroid_projected;
		centroid_projected[0] =face_centroid[alfa];
		centroid_projected[1] =face_centroid[beta];

	  for (int line_it = 0; line_it < number_of_lines ; ++line_it)
	  {
		  Obj_VectorNewLine[line_it].X_0_projection[alfa] = Obj_VectorNewLine[line_it].X0[alfa];
		  Obj_VectorNewLine[line_it].X_0_projection[beta] = Obj_VectorNewLine[line_it].X0[beta];
		  Obj_VectorNewLine[line_it].X_0_projection[gamma_] = (Obj_VectorNewLine[line_it].X0[alfa]*face_normal_vector[alfa]+
				  Obj_VectorNewLine[line_it].X0[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);

		  Obj_VectorNewLine[line_it].X_1_projection[alfa] = Obj_VectorNewLine[line_it].X1[alfa];
		  Obj_VectorNewLine[line_it].X_1_projection[beta] = Obj_VectorNewLine[line_it].X1[beta];
		  Obj_VectorNewLine[line_it].X_1_projection[gamma_] = (Obj_VectorNewLine[line_it].X1[alfa]*face_normal_vector[alfa]+
				  Obj_VectorNewLine[line_it].X1[beta]*face_normal_vector[beta]+w)/(-face_normal_vector[gamma_]);

		  //	X_projection[initial point,final point][X,Y coordinate]
		  double dy = Obj_VectorNewLine[line_it].X_1_projection[beta]-Obj_VectorNewLine[line_it].X_0_projection[beta];
		  double dx = Obj_VectorNewLine[line_it].X_1_projection[alfa]-Obj_VectorNewLine[line_it].X_0_projection[alfa];
		  double length = sqrt(dx*dx+dy*dy);

		  dy =  dy / length;
		  dx =  dx / length;

		  Obj_VectorNewLine[line_it].length_projection = length; // Is this the same as unit_length?

		  // Procedure to find the Normal vector of the projected line (normal_projection).
//		  I believe I could do this inside the method which implements the normal vector of the face.
			// If this is the right normal vector, it means that the line is going from X0 to X1!
			Point<2> n1(dy,-dx);
			// If this is the right normal vector, it means that the line is going from X1 to X0!
			// We can rearrange this Vector so that it becomes oriented from X0 to X1. This
			// will fix the parametric equation in return_face_integration.
			Point<2> n2(-dy,dx);

			// Exctract only the relevant (alfa, beta) components of the X_projection points.
			Point<2> aux_n00;
			Point<2> aux_n01;
			aux_n00[0] =Obj_VectorNewLine[line_it].X_0_projection[alfa];
			aux_n00[1] = Obj_VectorNewLine[line_it].X_0_projection[beta];

			aux_n01[0] =Obj_VectorNewLine[line_it]. X_1_projection[alfa];
			aux_n01[1] = Obj_VectorNewLine[line_it].X_1_projection[beta];

			// aux_ni represents a point beginning on the middle of the line and being pointed by the normal
			// vector. If the resulting point points outwards the plane, this is the right normal vector.
			Point<2> aux_n1;
			Point<2> aux_n2;
			aux_n1 = (n1+(aux_n00+aux_n01)/2);
		//	aux_n1 = (n1+(l->X_projection[0]+l->X_projection[1])/2);
		//	Point<2> aux_n2 = (n2+(l->X_projection[0]+l->X_projection[1])/2);
			aux_n2 = (n2+(aux_n00+aux_n01)/2);

			Point<3> Xtemp;

			if ( aux_n1.distance(centroid_projected) > aux_n2.distance(centroid_projected) )
				Obj_VectorNewLine[line_it].normal_projection = n1;

			else
				Obj_VectorNewLine[line_it].normal_projection = n2;


	  } // end loop Lines
}
