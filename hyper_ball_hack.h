/*
 * hyper_ball_hack.h
 *
 *  Created on: Nov 20, 2014
 *      Author: afonsoal
 */

#ifndef HYPER_BALL_HACK_H_
#define HYPER_BALL_HACK_H_

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <map>

using namespace dealii;
class hyper_ball_hack
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
	};
	std::vector <Point<2>> new_vertices_vector;
	std::vector <int> new_vertices_vector_index;

public:
	std::vector<new_cell_struct> vector_new_cell_struct;
	std::vector<Point<2> > support_points;
	hyper_ball_hack();
	void reinit();
	void set_variables(const std::vector<Point<2> > & _support_points,
		const int _dofs_per_cell);

	void set_new_mesh (std::vector<types::global_dof_index> cell_global_dof_indices,
			const bool cell_is_boundary );
	void create_new_triangulation (Triangulation<2> &tria);

};

#endif

