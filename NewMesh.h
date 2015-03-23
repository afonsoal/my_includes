/*
 * NewMesh.h
 *
 *  Created on: Nov 20, 2014
 *      Author: afonsoal
 */

#ifndef NEWMESH_H_
#define NEWMESH_H_

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <map>

using namespace dealii;
class NewMesh
{
private:

	int dofs_per_cell;
	int new_cell_counter;
	int new_vertices_index;

	// Structure containing the information relevant to create a new triangulation.
	struct new_cell_struct
	{
		std::vector <Point<2>> new_vertices_cell; // size = 4; will repeat between cells
		CellData<2> new_cells; // 4 vertices
		bool cell_is_boundary;
		int custom_cell_index;
	};
	std::vector <Point<2>> new_vertices_vector;
	std::vector <int> new_vertices_vector_index;
//	enum
//	{
//	not_surface_domain_id,
//	surface_domain_id
//	};

public:
	std::vector<new_cell_struct> vector_new_cell_struct;
	std::vector<Point<2> > support_points;
	NewMesh();
	void reinit();
	void set_variables(const std::vector<Point<2> > & _support_points,
		const int _dofs_per_cell);

	void set_new_mesh (std::vector<types::global_dof_index> cell_global_dof_indices,
			const bool cell_is_boundary );
	void create_new_triangulation (Triangulation<2> &tria);
	void set_triangulation_indices (Triangulation<2> &tria);

};

#endif

