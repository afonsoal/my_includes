/** hyper_ball_hack.cpp
 *
 *  Created on: Nov 20, 2014
 *      Author: afonsoal*/

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <iostream>

#include "hyper_ball_hack.h"

using namespace dealii;

hyper_ball_hack::hyper_ball_hack() {}
void hyper_ball_hack::reinit()
{
	vector_new_cell_struct.clear();
	new_vertices_vector.clear();
	new_cell_counter = 0;
	new_vertices_index = 0;
}
void hyper_ball_hack::set_variables(const std::vector<Point<2> > & _support_points,
		const int _dofs_per_cell)
{
	support_points = _support_points;
	dofs_per_cell = _dofs_per_cell;
}

void hyper_ball_hack::set_new_mesh
		(std::vector<types::global_dof_index> cell_global_dof_indices,
				const bool _cell_is_boundary)
{
	vector_new_cell_struct.push_back(new_cell_struct());
	vector_new_cell_struct[new_cell_counter].new_vertices_cell.resize(4);
	// Set the coordinates for new nodes, if they haven't been assigned yet.
	// Set also the global index of the vertex (node).

	for (int i=0;i<dofs_per_cell;++i)	{
		int j = cell_global_dof_indices[i];

		// Set the coordinates for all the vertices of this cell.
		vector_new_cell_struct[new_cell_counter].new_vertices_cell[i] =
				support_points[j];
		// Identify if cell is boundary(true) or inside (false).
		vector_new_cell_struct[new_cell_counter].cell_is_boundary = _cell_is_boundary;

		// Verify if this vertex has already been assigned to the vector of global
		// vertices indices. In this vector, the vertices cannot repeat themselves.
		std::vector<Point<2> >::iterator it;
		it = std::find(new_vertices_vector.begin(), new_vertices_vector.end(),
				support_points[j]);

		// Vertex corresponding to support_points[j] was not found; must assign it
		// to the global vector of nodal coordinates and global vector of nodal indices
		if (it == new_vertices_vector.end())
		{
			// Global vectors
			new_vertices_vector.push_back(support_points[j]);
			new_vertices_vector_index.push_back(new_vertices_index);

			// Local (cell) vector of vertices; must hold the same vertex numbering
			// as the global vector
			vector_new_cell_struct[new_cell_counter].new_cells.vertices[i] =
					new_vertices_index;

			new_vertices_index++;
		}
		// support_points[j] was already assigned
		else {
			// Location in new_vertice_vector where support_points[j] was assigned;
			// This is the global index of the vertex
			auto pos = it - new_vertices_vector.begin();
			// Assign local cell vector of vertices with the global node number
			// corresponding to the position where it was last found in the
			// global vector.

			vector_new_cell_struct[new_cell_counter].new_cells.vertices[i] = pos;
		}
	} // end for i<dofs_per_cell
	new_cell_counter++;
}

void hyper_ball_hack::create_new_triangulation (Triangulation<2> &tria)
{
	std::vector<CellData<2> > cells_data (
			vector_new_cell_struct.size(), CellData<2>());
	for (unsigned int i=0; i< vector_new_cell_struct.size(); ++i) {
		for (unsigned int j=0; j< 4; ++j) {

			cells_data[i].vertices[j] = vector_new_cell_struct[i].new_cells.vertices[j];
		}
	}

	tria.create_triangulation
	( new_vertices_vector, cells_data, SubCellData());

	//	GridReordering<2>::reorder_cells(cells); // necessary?
}
