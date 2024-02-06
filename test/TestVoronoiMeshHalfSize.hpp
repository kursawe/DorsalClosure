#ifndef TESTVORONOIMESHHALFSIZE_HPP_
#define TESTVORONOIMESHHALFSIZE_HPP_

#include <cxxtest/TestSuite.h>
#include "MyVoronoiVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "ModifiedVoronoiVertexMeshGenerator.hpp"
#include "PlaneStickyBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "FarhadifarForce.hpp"
#include "AbstractMesh.hpp"
#include "VoronoiTargetAreaModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "MyApoptoticCellKiller.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Debug.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OutputFileHandler.hpp"
#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "MyCellKiller.hpp"
#include "ChasteSerialization.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include <boost/serialization/base_object.hpp>
#include "FakePetscSetup.hpp"
#include "RandomNumberGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"

class TestVoronoiMeshHalfSize : public AbstractCellBasedTestSuite
{
public:
    void TestVoronoiMeshHalfSizeV1()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(8,22,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
		MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(CellLabel, p_label);
    
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);

			p_model->SetStemCellG1Duration(48);
            p_model->SetTransitCellG1Duration(60);
            p_model->SetMaxTransitGenerations(3);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                               (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            //no proliferation
			p_cell->SetBirthTime(birth_time);
            //p_cell->SetCellProliferativeType(p_diff_type);

            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

			if (this_location(1) > 5.9 && this_location(1) < 34.5)
			{
				cell_iter->AddCellProperty(p_label);
                cell_iter->SetApoptosisTime(DBL_MAX);	
                cell_iter->SetCellProliferativeType(p_diff_type);		
			}
            else
            {
                cell_iter->SetCellProliferativeType(p_stem_type);
                //cell_iter->SetCellProliferativeType(p_diff_type);
            }
			if (this_location(0) < -0.2 && (this_location(1) < 5.5 || this_location(1) > 35.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 27.4 && (this_location(1) < 5.5 || this_location(1) > 35.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 27.0 && this_location(1) > 34.0 && this_location(1) < 37.0)
			{
				cell_iter->Kill();
			}
		}

        OffLatticeSimulation<2> simulator(cell_population);
		
        simulator.SetOutputDirectory("TestHalfMesh");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1);//0.01
        simulator.SetEndTime(10.0);  //300
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(VoronoiTargetAreaModifier<2>, p_growth_modifier);
		//double lec_area = sqrt(3.0)*3.625*3.625/2.0;
		//p_growth_modifier.SetLecArea(lec_area);
        simulator.AddSimulationModifier(p_growth_modifier);


		double top_height = 39.8 ;
		//double bottom_height = 3.0/25.0;
		double bottom_height = 0.2;
		double right_bound = 27.204975;
		double left_bound = 0.0;
		unsigned num_nodes = cell_population.GetNumNodes();

		//iterate over vertices in the cell population 
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
			std::set<unsigned> containing_element_indices = p_this_node->rGetContainingElementIndices();
			bool hb_neighbour = false;
			bool lec_neighbour = false;
			int elem_with_node = p_this_node->GetNumContainingElements();
			if (elem_with_node == 1)
			{
				if (node_position[1] < 6.5 && node_position[1] > 4.5 )
				{
					node_position[1] -= 0.7;
				}
				else if (node_position[1] > 33.5 && node_position[1] < 35.5)
				{
					node_position[1] += 0.7;
				}
			}

			for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
             iter != containing_element_indices.end();
             iter++)
        	{
				//CellPtr p_current_cell = this->GetCellUsingLocationIndex(*iter);
				if (cell_population.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>())
				{
					hb_neighbour = true;
				}
				else
				{
					lec_neighbour = true;
				}
			}
			if (hb_neighbour == true && lec_neighbour == true)
			{
				if (node_position[1] < 6.5 && node_position[1] > 4.5 )
				{
					node_position[1] -= 0.7;
				}
				else if (node_position[1] > 33.5 && node_position[1] < 35.5)
				{
					node_position[1] += 0.7;
				}
			}
		}
        

		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			if ( node_position[1] > 39.8 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < 0.2 )
			{
				assert( p_this_node->IsBoundaryNode() );
			    node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < 0.01) 
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[0] > 27.0 && node_position[1] < 35.5)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 27.0 && node_position[1] > 36.8)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 26.8 && node_position[1] > 34.5 && node_position[1] < 36.8)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
		}

		// BOTTOM
		c_vector<double,2> bottom_point = zero_vector<double>(2);
		bottom_point(0) = 0.0;
		bottom_point(1) = bottom_height;

		c_vector<double,2> bottom_normal = zero_vector<double>(2);
		bottom_normal(1) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_lower_boundary_condition, (&cell_population, bottom_point, bottom_normal) );
		simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_condition);

		// TOP
		c_vector<double,2> top_point = zero_vector<double>(2);
		top_point(0) = 0.0;
		top_point(1) = top_height;

		c_vector<double,2> top_normal = zero_vector<double>(2);
		top_normal(1) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_top_boundary_condition, (&cell_population, top_point, top_normal));
		simulator.AddCellPopulationBoundaryCondition(p_top_boundary_condition);

		// LEFT
		c_vector<double,2> left_point = zero_vector<double>(2);
		left_point(0) = left_bound;
		left_point(1) = 0.0;

		c_vector<double,2> left_normal = zero_vector<double>(2);
		left_normal(0) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_left_boundary_condition, (&cell_population, left_point, left_normal));
		simulator.AddCellPopulationBoundaryCondition(p_left_boundary_condition);

		// RIGHT
		c_vector<double,2> right_point = zero_vector<double>(2);
		right_point(0) = right_bound;
		right_point(1) = 0.0;

		c_vector<double,2> right_normal = zero_vector<double>(2);
		right_normal(0) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_right_boundary_condition, (&cell_population, right_point, right_normal));
		simulator.AddCellPopulationBoundaryCondition(p_right_boundary_condition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.01));
        simulator.AddCellKiller(p_cell_killer);
	
	    cell_population.RemoveDeadCells();
		cell_population.Update();

        simulator.Solve();
    }
	void xTestVoronoiMeshHalfSizeV1_2()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(8,22,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
		MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(CellLabel, p_label);
    
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);

			p_model->SetStemCellG1Duration(48);
            p_model->SetTransitCellG1Duration(60);
            p_model->SetMaxTransitGenerations(6);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                               (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            //no proliferation
			p_cell->SetBirthTime(birth_time);
            //p_cell->SetCellProliferativeType(p_diff_type);

            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

			if (this_location(1) > 5.9 && this_location(1) < 34.5)
			{
				cell_iter->AddCellProperty(p_label);
                cell_iter->SetApoptosisTime(DBL_MAX);	
                cell_iter->SetCellProliferativeType(p_diff_type);		
			}
            else
            {
                cell_iter->SetCellProliferativeType(p_stem_type);
                //cell_iter->SetCellProliferativeType(p_diff_type);
            }
			if (this_location(0) < -0.2 && (this_location(1) < 5.5 || this_location(1) > 35.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 27.4 && (this_location(1) < 5.5 || this_location(1) > 35.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 27.0 && this_location(1) > 34.0 && this_location(1) < 37.0)
			{
				cell_iter->Kill();
			}
		}

        OffLatticeSimulation<2> simulator(cell_population);
		
        simulator.SetOutputDirectory("TestHalfMesh_2");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1);//0.01
        simulator.SetEndTime(10.0);  //300//160
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(VoronoiTargetAreaModifier<2>, p_growth_modifier);
		//double lec_area = sqrt(3.0)*3.625*3.625/2.0;
		//p_growth_modifier.SetLecArea(lec_area);
        simulator.AddSimulationModifier(p_growth_modifier);


		double top_height = 39.8 ;
		//double bottom_height = 3.0/25.0;
		double bottom_height = 0.2;
		double right_bound = 27.62;
		double left_bound = 0.0;
		unsigned num_nodes = cell_population.GetNumNodes();

		//iterate over vertices in the cell population 
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
			std::set<unsigned> containing_element_indices = p_this_node->rGetContainingElementIndices();
			bool hb_neighbour = false;
			bool lec_neighbour = false;
			int elem_with_node = p_this_node->GetNumContainingElements();
			if (elem_with_node == 1)
			{
				if (node_position[1] < 6.5 && node_position[1] > 4.5 )
				{
					node_position[1] -= 0.7;
				}
				else if (node_position[1] > 33.5 && node_position[1] < 35.5)
				{
					node_position[1] += 0.7;
				}
			}

			for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
             iter != containing_element_indices.end();
             iter++)
        	{
				//CellPtr p_current_cell = this->GetCellUsingLocationIndex(*iter);
				if (cell_population.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>())
				{
					hb_neighbour = true;
				}
				else
				{
					lec_neighbour = true;
				}
			}
			if (hb_neighbour == true && lec_neighbour == true)
			{
				if (node_position[1] < 6.5 && node_position[1] > 4.5 )
				{
					node_position[1] -= 0.7;
				}
				else if (node_position[1] > 33.5 && node_position[1] < 35.5)
				{
					node_position[1] += 0.7;
				}
			}

			//if (node_position[1] > 6.5 && node_position[1] < 7.5 && node_position[0] > -0.1 && node_position[0] < 3.7)
			//{
			//	PRINT_VARIABLE(node_position[0]);
			//}
		}
        

		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			if ( node_position[1] > 39.8 )
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < 0.2 )
			{
				//assert( p_this_node->IsBoundaryNode() );
			    node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < 0.01) 
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[0] > 27.0 && node_position[1] < 35.5)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 27.0 && node_position[1] > 36.8)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 26.8 && node_position[1] > 34.5 && node_position[1] < 36.8)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
		}

		// BOTTOM
		c_vector<double,2> bottom_point = zero_vector<double>(2);
		bottom_point(0) = 0.0;
		bottom_point(1) = bottom_height;

		c_vector<double,2> bottom_normal = zero_vector<double>(2);
		bottom_normal(1) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_lower_boundary_condition, (&cell_population, bottom_point, bottom_normal) );
		simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_condition);

		// TOP
		c_vector<double,2> top_point = zero_vector<double>(2);
		top_point(0) = 0.0;
		top_point(1) = top_height;

		c_vector<double,2> top_normal = zero_vector<double>(2);
		top_normal(1) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_top_boundary_condition, (&cell_population, top_point, top_normal));
		simulator.AddCellPopulationBoundaryCondition(p_top_boundary_condition);

		// LEFT
		c_vector<double,2> left_point = zero_vector<double>(2);
		left_point(0) = left_bound;
		left_point(1) = 0.0;

		c_vector<double,2> left_normal = zero_vector<double>(2);
		left_normal(0) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_left_boundary_condition, (&cell_population, left_point, left_normal));
		simulator.AddCellPopulationBoundaryCondition(p_left_boundary_condition);

		// RIGHT
		c_vector<double,2> right_point = zero_vector<double>(2);
		right_point(0) = right_bound;
		right_point(1) = 0.0;

		c_vector<double,2> right_normal = zero_vector<double>(2);
		right_normal(0) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_right_boundary_condition, (&cell_population, right_point, right_normal));
		simulator.AddCellPopulationBoundaryCondition(p_right_boundary_condition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.008));
        simulator.AddCellKiller(p_cell_killer);
	
	    cell_population.RemoveDeadCells();
		cell_population.Update();

        simulator.Solve();
    }
	void xTestVoronoiSmallerMesh()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(6,20,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
		MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(CellLabel, p_label);
    
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);

			p_model->SetStemCellG1Duration(48);
            p_model->SetTransitCellG1Duration(60);
            p_model->SetMaxTransitGenerations(6);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                               (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            //no proliferation
			p_cell->SetBirthTime(birth_time);
            //p_cell->SetCellProliferativeType(p_diff_type);

            cells.push_back(p_cell);
        }


	
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

			if (this_location(1) > 8.0 && this_location(1) < 42.0)
			{
				cell_iter->AddCellProperty(p_label);
                cell_iter->SetApoptosisTime(DBL_MAX);	
                cell_iter->SetCellProliferativeType(p_diff_type);		
			}
            else
            {
                cell_iter->SetCellProliferativeType(p_stem_type);
                //cell_iter->SetCellProliferativeType(p_diff_type);
            }
			if (this_location(0) < -0.2 && (this_location(1) < 8.0 || this_location(1) > 42.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 29.8 && (this_location(1) < 8.0 || this_location(1) > 42.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 29.0 && this_location(1) > 41.0 && this_location(1) < 45.5)
			{
				cell_iter->Kill();
			}
		}

        OffLatticeSimulation<2> simulator(cell_population);
		
        simulator.SetOutputDirectory("TestSmallerMesh");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1);//0.01
        simulator.SetEndTime(10.0);  //300//30
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(VoronoiTargetAreaModifier<2>, p_growth_modifier);
		//double lec_area = sqrt(3.0)*3.625*3.625/2.0;
		//p_growth_modifier.SetLecArea(lec_area);
        simulator.AddSimulationModifier(p_growth_modifier);

		//MARK;
		double top_height = 49.5 ;
		//double bottom_height = 3.0/25.0;
		double bottom_height = 0.4;
		double right_bound = 29.5; //29.5
		double left_bound = 0.0;
		unsigned num_nodes = cell_population.GetNumNodes();


		//double current_a = 0.3726703612;
		//double desired_a = sqrt(2.0/(3.0*sqrt(3.0)));
		//double node_multiplier = desired_a/current_a;

		//iterate over vertices in the cell population 
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			//node_position[0] *= node_multiplier;
			//node_position[1] *= node_multiplier;
			std::set<unsigned> containing_element_indices = p_this_node->rGetContainingElementIndices();
			bool hb_neighbour = false;
			bool lec_neighbour = false;
			int elem_with_node = p_this_node->GetNumContainingElements();
			if (elem_with_node == 1)
			{
				if (node_position[1] < 10.0 && node_position[1] > 6.5 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 40.5 && node_position[1] < 43.5)
				{
					node_position[1] += 1.0;
				}
			}

			for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
             iter != containing_element_indices.end();
             iter++)
        	{
				//CellPtr p_current_cell = this->GetCellUsingLocationIndex(*iter);
				if (cell_population.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>())
				{
					hb_neighbour = true;
				}
				else
				{
					lec_neighbour = true;
				}
			}
			if (hb_neighbour == true && lec_neighbour == true)
			{
				if (node_position[1] < 10.0 && node_position[1] > 6.5 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 40.5 && node_position[1] < 43.5)
				{
					node_position[1] += 1.0;
				}
			}

			//if (node_position[1] > 6.5 && node_position[1] < 7.5 && node_position[0] > -0.1 && node_position[0] < 3.7)
			//{
			//	PRINT_VARIABLE(node_position[0]);
			//}
		}
        
		
		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
			if (node_position[1] > 1.0 && node_position[1] < 1.5 && node_position[0] > -1.0 && node_position[0] < 1.0)
			{
				//PRINT_VARIABLE(node_position[0]);
				//PRINT_VARIABLE(node_position[1]);
			}
			if ( node_position[1] > 49.8 )
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < 0.2 )
			{
				//assert( p_this_node->IsBoundaryNode() );
			    node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < 0.0) 
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[0] > 29.5 && node_position[1] > 7.0 && node_position[1] < 41.0)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 29.5 && node_position[1] < 5.0)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 29.54 && node_position[1] < 6.0 && node_position[1] > 5.8)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
				//PRINT_VARIABLE(node_position[0]);
			}
			else if ( node_position[0] > 29.5 && node_position[1] > 45.0)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 28.9 && node_position[1] > 42.3 && node_position[1] < 45.5)
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
		}

		// BOTTOM
		c_vector<double,2> bottom_point = zero_vector<double>(2);
		bottom_point(0) = 0.0;
		bottom_point(1) = bottom_height;

		c_vector<double,2> bottom_normal = zero_vector<double>(2);
		bottom_normal(1) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_lower_boundary_condition, (&cell_population, bottom_point, bottom_normal) );
		simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_condition);

		// TOP
		c_vector<double,2> top_point = zero_vector<double>(2);
		top_point(0) = 0.0;
		top_point(1) = top_height;

		c_vector<double,2> top_normal = zero_vector<double>(2);
		top_normal(1) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_top_boundary_condition, (&cell_population, top_point, top_normal));
		simulator.AddCellPopulationBoundaryCondition(p_top_boundary_condition);

		// LEFT
		c_vector<double,2> left_point = zero_vector<double>(2);
		left_point(0) = left_bound;
		left_point(1) = 0.0;

		c_vector<double,2> left_normal = zero_vector<double>(2);
		left_normal(0) = -1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_left_boundary_condition, (&cell_population, left_point, left_normal));
		simulator.AddCellPopulationBoundaryCondition(p_left_boundary_condition);

		// RIGHT
		c_vector<double,2> right_point = zero_vector<double>(2);
		right_point(0) = right_bound;
		right_point(1) = 0.0;

		c_vector<double,2> right_normal = zero_vector<double>(2);
		right_normal(0) = 1.0;

		MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_right_boundary_condition, (&cell_population, right_point, right_normal));
		simulator.AddCellPopulationBoundaryCondition(p_right_boundary_condition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.008));
        simulator.AddCellKiller(p_cell_killer);
	
	    cell_population.RemoveDeadCells();
		cell_population.Update();
        simulator.Solve();
    }
};

#endif /*TESTVORONOIMESHHALFSIZEHPP_*/