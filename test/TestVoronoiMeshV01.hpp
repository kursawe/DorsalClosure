#ifndef TESTVORONOIMESHV01_HPP_
#define TESTVORONOIMESHV01_HPP_

#include <cxxtest/TestSuite.h>
#include "MyVoronoiVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "ModifiedVoronoiVertexMeshGenerator.hpp"
#include "PlaneStickyBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "VoronoiDataWriter.hpp"
#include "FarhadifarForce.hpp"
#include "AbstractMesh.hpp"
#include "MyTargetAreaModifier.hpp"
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

class TestVoronoiMeshV01 : public AbstractCellBasedTestSuite
{
public:
    //trying to use voronoi mesh 
    void xTestVoronoiMesh()
    {
		//multiply y by 0.7 in myvoronoivertexmeshgenerator
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        MyVoronoiVertexMeshGenerator generator(6,10,3,0,3.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(CellLabel, p_label);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		//for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		//	cell_iter != cell_population.End();
		//	++cell_iter)
	    //{
		//   c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		//    if (this_location(1) < 3.0 && this_location(1) > 9.0)

		//    {
        //    	cell_iter->AddCellProperty(p_label);
                
        //    }
	    //}


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVoronoiMesh");
        simulator.SetSamplingTimestepMultiple(160);
		simulator.SetDt(0.01);
        simulator.SetEndTime(10.0);
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

		double top_height = 17.3*0.7;
		double bottom_height = 0.0;
		double right_bound = 9.1;
		double left_bound = -0.5;

		unsigned num_nodes = cell_population.GetNumNodes();

		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			if ( node_position[1] > 17.2*0.7 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			else if ( node_position[1] < 0.2*0.7 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < -0.4 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[0] > 9.1 )
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
	
	    
        simulator.Solve();
    }
	 void xTestVoronoiMeshBoundaries()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        MyVoronoiVertexMeshGenerator generator(6,10,3,0,3.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
            if (cell_iter<36 || cell_iter>59)
            //histoblasts
            {
                //proliferation 
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(48);
                p_model->SetTransitCellG1Duration(5);
                p_model->SetMaxTransitGenerations(3);

                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);
            }
            else
            //LECs
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                //no proliferation
                p_cell->SetCellProliferativeType(p_diff_type);

                p_cell->AddCellProperty(p_label);
                
                //p_cell->SetApoptosisTime(DOUBLE_UNSET); 
                //p_cell->SetApoptosisTime(20);  

                cells.push_back(p_cell);

            }
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);


		//for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
		//	cell_iter != cell_population.End();
		//	++cell_iter)
	    //{
		//    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		//    if (this_location(1) < 1.8 && this_location(1) > 7)

		//    {
                //cell_iter->ReadyToDivide();
                
        //    }
	    //}

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVoronoiMeshBoundaries");
        simulator.SetSamplingTimestepMultiple(10); //georgia 160
		simulator.SetDt(0.1);//georgia 0.01
        simulator.SetEndTime(10.0); 	//97 //georgia 90
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
		//p_force->SetBoundaryLineTensionParameter(0.3);      //0.12 default 
        //p_force->SetPerimeterContractilityParameter(0.01);  //0.04 default
        //p_force->SetAreaElasticityParameter(0.8);            //1.0 default 
        //p_force->SetLineTensionParameter(0.06);  
        simulator.AddForce(p_force);
       
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        //MARK;
	    
        // BOTTOM
	    c_vector<double,2> bottom_point = zero_vector<double>(2);
	    bottom_point(1) = 0.0;

	    c_vector<double,2> bottom_normal = zero_vector<double>(2);
	    bottom_normal(1) = -1.0;

	    MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_lower_boundary_condition, (&cell_population, bottom_point, bottom_normal) );
	    simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_condition);

	    // TOP
	    c_vector<double,2> top_point = zero_vector<double>(2);
	    top_point(1) = 17.0;

	    c_vector<double,2> top_normal = zero_vector<double>(2);
	    top_normal(1) = 1.0;

	    MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_top_boundary_condition, (&cell_population, top_point, top_normal));
	    simulator.AddCellPopulationBoundaryCondition(p_top_boundary_condition);

	    // LEFT
		c_vector<double,2> left_point = zero_vector<double>(2);
	    left_point(0) = 0.0;

	    c_vector<double,2> left_normal = zero_vector<double>(2);
	    left_normal(0) = -1.0;

	    MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_left_boundary_condition, (&cell_population, left_point, left_normal));
	    simulator.AddCellPopulationBoundaryCondition(p_left_boundary_condition);

	    // RIGHT
	    c_vector<double,2> right_point = zero_vector<double>(2);
	    right_point(0) = 10.5;

	    c_vector<double,2> right_normal = zero_vector<double>(2);
	    right_normal(0) = 1.0;

	    MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_right_boundary_condition, (&cell_population, right_point, right_normal));
	    simulator.AddCellPopulationBoundaryCondition(p_right_boundary_condition);
        
		//MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.046));
        //simulator.AddCellKiller(p_cell_killer);
		
		simulator.Solve();
    }
	void TestLabels()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(16,22,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
		MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		//MAKE_PTR(StemCellProliferativeType, p_stem_type);
		MAKE_PTR(CellLabel, p_label);
    
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);

			//p_model->SetStemCellG1Duration(48);
            //p_model->SetTransitCellG1Duration(5);
            //p_model->SetMaxTransitGenerations(3);

            //double birth_time = - RandomNumberGenerator::Instance()->ranf() *
             //                   (  p_model->GetStemCellG1Duration()
              //                   + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            //no proliferation
			//p_cell->SetBirthTime(birth_time);
            p_cell->SetCellProliferativeType(p_diff_type);

            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

			if (this_location(1) < 7.0 || this_location(1) > 48.0)
			{
				cell_iter->AddCellProperty(p_label);			
			}
			if (this_location(1) > 47.0 && (this_location(0) < 0.0 || this_location(0) > 78.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(1) < 7.0 && (this_location(0) < 0.0 || this_location(0) > 78.0))
			{
				cell_iter->Kill();
			}
			else if (this_location(1) > 4.5 && this_location(1) < 5.5 && this_location(0) > 77.5)
			{
				cell_iter->Kill();
			}
			else if (this_location(1) > 48.0 && this_location(1) < 49.5 && this_location(0) > 77.0)
			{
				cell_iter->Kill();
			}
			else if (this_location(1) > 49.5 && this_location(1) < 50.5 && this_location(0) > 77.0)
			{
				cell_iter->Kill();
			}
		}

        OffLatticeSimulation<2> simulator(cell_population);
		
        simulator.SetOutputDirectory("TestLabelledCells");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1); //0.01
        simulator.SetEndTime(5.0); //100
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

		double top_height = 55.0 ;
		//double bottom_height = 3.0/25.0;
		double bottom_height = 0.0;
		double right_bound = 77.5;
		double left_bound = 0.0;
		unsigned num_nodes = cell_population.GetNumNodes();

		//iterate over vertices in the cell population 
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
			std::set<unsigned> containing_element_indices = p_this_node->rGetContainingElementIndices();
			//std::set<unsigned> containing_element_indices = this->GetNode(node_index)->rGetContainingElementIndices();
			//iterate over the elements which contain this node 
			//for (typename Node<SPACE_DIM>::ContainingElementIterator elem_iter = p_this_node->ContainingElementsBegin();
            // elem_iter != p_this_node->ContainingElementsEnd();
            // ++elem_iter)
			bool hb_neighbour = false;
			bool lec_neighbour = false;
			int elem_with_node = p_this_node->GetNumContainingElements();
			if (elem_with_node == 1)
			{
				if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 46.5 && node_position[1] < 49.0)
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
			//unsigned nearest_index = p_mesh->GetNearestNodeIndex(node_position);
			//Node<2>* p_nearest_node = cell_population.GetNode(nearest_index);
			//c_vector<double, 2>& nearest_node_position = p_nearest_node->rGetModifiableLocation();
			//double dist_in_y = nearest_node_position[1] - node_position[1];  
			if (hb_neighbour == true && lec_neighbour == true)
			{
				//HOW_MANY_TIMES_HERE("lec hb boundary");
				//if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				//{
				//	if (fabs(dist_in_y) < 0.1)
				//	{
				//		if (dist_in_y < 0)
				//		{
				//			MARK;
				//			if ((nearest_node_position[0] - node_position[0]) > 0)
				//			{
				//				node_position[0] -= 1.0;
				//				MARK;
				//			}
				//			else
				//			{
				//				node_position[0] += 1.0;
				//			}
				//			node_position[1] -= 0.5;
				//		}
				//		else
				//		{
				//			node_position[1] -= 1.0;
				//		}
				//	}
				//	else
				//	{
				//		node_position[1] -= 1.0;
				//	}
				//}
				if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 46.5 && node_position[1] < 49.0)
				{
					node_position[1] += 1.0;
				}
			}
		}
        

		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
//
			if ( node_position[1] > 55.0 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < 0.2 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < 0.01) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[0] > 77.49 && node_position[1] < 48.01)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 77.49 && node_position[1] > 49.5)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[0] > 76.9 && node_position[1] > 47.5 && node_position[1] < 51.0)
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			//else if ( node_position[0] > 76.9 && node_position[1] > 48.0 && node_position[1] < 50.0)
			//{
			//	assert( p_this_node->IsBoundaryNode() );
			//	node_position[0] = right_bound + 0.001;
			//}
//			else if ( node_position[1] > 49.79 && node_position[1] < 55.0 && node_position[0] < -0.49) //10, 46
//			{
//				assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = left_bound - 0.001;
//			}
//			else if ( node_position[1] > 49.0 && node_position[1] < 49.79 && node_position[0] < -0.03) //10, 46
//			{
				//assert( p_this_node->IsBoundaryNode() );
				//PRINT_VARIABLE(node_position[0]);
//				node_position[0] = left_bound - 0.001;
//			}
//			else if ( node_position[1] > -0.2 && node_position[1] < 7.0 && node_position[0] < -0.49) //10, 46
//			{
//				assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = left_bound - 0.001;
//			}
//			else if ( node_position[1] > 7.0 && node_position[1] < 7.6 && node_position[0] < 0.01) //10, 46
//			{
				//assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = left_bound - 0.001;
//			}

//			else if ( node_position[1] > 6.5 && node_position[1] < 47.0 && node_position[0] > 77.49 ) //10, 46
//			{
//				assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = right_bound + 0.001;
//			}
//			if (fabs(node_position[0] - 77.5) < 0.01 && fabs(node_position[1] - 47.8556) < 0.01)
//			{
//				node_position[0] = 78.0;
//				node_position[1] = 47.3968;
//			}
//			if ( node_position[1] > 48.0 && node_position[1] < 54.9 && node_position[0] > 78.9 ) //10, 46
//			{
//				assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = right_bound + 0.001;
//			}
//			else if ( node_position[1] > 47.0 && node_position[1] < 48.0 && node_position[0] > 77.1 ) //10, 46
//			{
//				MARK;
//				node_position[0] = right_bound + 0.001;
//			}
//			if (fabs(node_position[0] -  78.5) < 0.01 && fabs(node_position[1] - 5.34049) < 0.01)
//			{
//				node_position[0] = 77.53;
//				node_position[1] = 5.46055;
//			}
//			else if (node_position[1] < 5.5 && node_position[0] > 78.49) //10, 46
//			{
//				assert( p_this_node->IsBoundaryNode() );
//				node_position[0] = right_bound + 0.001;
//			}
			//else if (node_position[1] < 5.5 && node_position[1] > 5.0 && node_position[0] > 78.2)
			//{
				//PRINT_VARIABLE(node_position[0]);
				//PRINT_VARIABLE(node_position[1]);
			//	assert( p_this_node->IsBoundaryNode() );
			//	node_position[0] = right_bound + 0.001;
			//}
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
	
	    cell_population.RemoveDeadCells();
		cell_population.Update();

        simulator.Solve();
    }
	void xTestVoronoiMeshShape()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(16,22,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(CellLabel, p_label);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVoronoiMeshShape");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1); //0.01
        simulator.SetEndTime(10.0); //10
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

		//double top_height = 12.0 - (3.0/25.0) ;
		//double bottom_height = 3.0/25.0;
		//double right_bound = ;
		//double left_bound = ;

		//unsigned num_nodes = cell_population.GetNumNodes();

		// Iterate over vertices in the cell population
		//for (unsigned node_index=0; node_index<num_nodes; node_index++)
		//{
		//	Node<2>* p_this_node = cell_population.GetNode(node_index);
		//	c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

		//	if ( node_position[1] > 11.9 )
		//	{
		//		assert( p_this_node->IsBoundaryNode() );
		//		node_position[1] = top_height + 0.001;
		//	}
		//	else if ( node_position[1] < 0.1 )
		//	{
		//		assert( p_this_node->IsBoundaryNode() );
		//		node_position[1] = bottom_height - 0.001;
		//	}

		//	if ( node_position[0] <  )
		//	{
		//		assert( p_this_node->IsBoundaryNode() );
		//		node_position[0] = left_bound - 0.001;
		//	}
		//	else if ( node_position[0] >  )
		//	{
		//		assert( p_this_node->IsBoundaryNode() );
		//		node_position[0] = right_bound + 0.001;
		//	}
		//}

		// BOTTOM
		//c_vector<double,2> bottom_point = zero_vector<double>(2);
		//bottom_point(0) = 0.0;
		//bottom_point(1) = bottom_height;

		//c_vector<double,2> bottom_normal = zero_vector<double>(2);
		//bottom_normal(1) = -1.0;

		//MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_lower_boundary_condition, (&cell_population, bottom_point, bottom_normal) );
		//simulator.AddCellPopulationBoundaryCondition(p_lower_boundary_condition);

		// TOP
		//c_vector<double,2> top_point = zero_vector<double>(2);
		//top_point(0) = 0.0;
		//top_point(1) = top_height;

		//c_vector<double,2> top_normal = zero_vector<double>(2);
		//top_normal(1) = 1.0;

		//MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_top_boundary_condition, (&cell_population, top_point, top_normal));
		//simulator.AddCellPopulationBoundaryCondition(p_top_boundary_condition);

		// LEFT
		//c_vector<double,2> left_point = zero_vector<double>(2);
		//left_point(0) = left_bound;
		//left_point(1) = 0.0;

		//c_vector<double,2> left_normal = zero_vector<double>(2);
		//left_normal(0) = -1.0;

		//MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_left_boundary_condition, (&cell_population, left_point, left_normal));
		//simulator.AddCellPopulationBoundaryCondition(p_left_boundary_condition);

		// RIGHT
		//c_vector<double,2> right_point = zero_vector<double>(2);
		//right_point(0) = right_bound;
		//right_point(1) = 0.0;

		//c_vector<double,2> right_normal = zero_vector<double>(2);
		//right_normal(0) = 1.0;

		//MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_right_boundary_condition, (&cell_population, right_point, right_normal));
		//simulator.AddCellPopulationBoundaryCondition(p_right_boundary_condition);
	
	    
        simulator.Solve();
    }
	void xTestLECDeath()
    {
        //cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
        ModifiedVoronoiVertexMeshGenerator generator(16,22,6,0,1.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        //CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
		MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
		MAKE_PTR(CellLabel, p_label);
    
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

		
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			//unsigned cell_row;
			//cell_row = cell_iter%22;
			//PRINT_VARIABLE(cell_row);			
            //if (cell_iter < 32)
            //{
				//lecs
                //no proliferation 
            //    UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            //    p_model->SetDimension(2);
            //    CellPtr p_cell(new Cell(p_state, p_model));
            //    p_cell->SetCellProliferativeType(p_diff_type);
			//	p_cell->AddCellProperty(p_label);
        
            //    cells.push_back(p_cell);
            //}
            //else
            //hbs
            //{
        	//	UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            //    p_model->SetDimension(2);
            //    CellPtr p_cell(new Cell(p_state, p_model));
                //proliferation
            //    p_cell->SetCellProliferativeType(p_diff_type);
				//p_cell->SetCellProliferativeType(p_stem_type);
				//p_model->SetStemCellG1Duration(48);
                //p_model->SetTransitCellG1Duration(5);
                //p_model->SetMaxTransitGenerations(3);

                //double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                //                    (  p_model->GetStemCellG1Duration()
                //                    + p_model->GetSG2MDuration() );

                //p_cell->SetBirthTime(birth_time);

            //    cells.push_back(p_cell);

            //}
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            //no proliferation
            p_cell->SetCellProliferativeType(p_diff_type);

            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

			if (this_location(1) > 7.0 && this_location(1) < 48.0)
			{
				cell_iter->AddCellProperty(p_label);
				//cell_iter->SetApoptosisTime(DOUBLE_UNSET);
						
			}
		}

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestLecDeath");
        simulator.SetSamplingTimestepMultiple(10);//160
		simulator.SetDt(0.1);//0.01
        simulator.SetEndTime(5.0);//50
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

		double top_height = 54.8 ;
		//double bottom_height = 3.0/25.0;
		double bottom_height = 0.2;
		double right_bound = 78.75 ;
		double left_bound = -1.25;
		unsigned num_nodes = cell_population.GetNumNodes();

		//iterate over vertices in the cell population 
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();
			std::set<unsigned> containing_element_indices = p_this_node->rGetContainingElementIndices();
			//std::set<unsigned> containing_element_indices = this->GetNode(node_index)->rGetContainingElementIndices();
			//iterate over the elements which contain this node 
			//for (typename Node<SPACE_DIM>::ContainingElementIterator elem_iter = p_this_node->ContainingElementsBegin();
            // elem_iter != p_this_node->ContainingElementsEnd();
            // ++elem_iter)
			bool hb_neighbour = false;
			bool lec_neighbour = false;
			int elem_with_node = p_this_node->GetNumContainingElements();
			if (elem_with_node == 1)
			{
				if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 46.5 && node_position[1] < 49.0)
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
			//unsigned nearest_index = p_mesh->GetNearestNodeIndex(node_position);
			//Node<2>* p_nearest_node = cell_population.GetNode(nearest_index);
			//c_vector<double, 2>& nearest_node_position = p_nearest_node->rGetModifiableLocation();
			//double dist_in_y = nearest_node_position[1] - node_position[1];  
			if (hb_neighbour == true && lec_neighbour == true)
			{
				//HOW_MANY_TIMES_HERE("lec hb boundary");
				//if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				//{
				//	if (fabs(dist_in_y) < 0.1)
				//	{
				//		if (dist_in_y < 0)
				//		{
				//			MARK;
				//			if ((nearest_node_position[0] - node_position[0]) > 0)
				//			{
				//				node_position[0] -= 1.0;
				//				MARK;
				//			}
				//			else
				//			{
				//				node_position[0] += 1.0;
				//			}
				//			node_position[1] -= 0.5;
				//		}
				//		else
				//		{
				//			node_position[1] -= 1.0;
				//		}
				//	}
				//	else
				//	{
				//		node_position[1] -= 1.0;
				//	}
				//}
				if (node_position[1] < 9.0 && node_position[1] > 6.4 )
				{
					node_position[1] -= 1.0;
				}
				else if (node_position[1] > 46.5 && node_position[1] < 49.0)
				{
					node_position[1] += 1.0;
				}
			}
		}
        

		// Iterate over vertices in the cell population
		for (unsigned node_index=0; node_index<num_nodes; node_index++)
		{
			Node<2>* p_this_node = cell_population.GetNode(node_index);
			c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			if ( node_position[1] > 55.0 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < 0.2 )
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[1] > 7.5 && node_position[1] < 49.0 && node_position[0] < 0.001) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[1] > 49.79 && node_position[1] < 55.0 && node_position[0] < -0.49) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[1] > 49.0 && node_position[1] < 49.79 && node_position[0] < -0.03) //10, 46
			{
				//assert( p_this_node->IsBoundaryNode() );
				//PRINT_VARIABLE(node_position[0]);
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[1] > -0.2 && node_position[1] < 7.0 && node_position[0] < -0.49) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}
			else if ( node_position[1] > 7.0 && node_position[1] < 7.6 && node_position[0] < 0.01) //10, 46
			{
				//assert( p_this_node->IsBoundaryNode() );
				node_position[0] = left_bound - 0.001;
			}

			else if ( node_position[1] > 6.5 && node_position[1] < 47.0 && node_position[0] > 77.49 ) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			if (fabs(node_position[0] - 77.5) < 0.01 && fabs(node_position[1] - 47.8556) < 0.01)
			{
				node_position[0] = 78.0;
				node_position[1] = 47.3968;
			}
			if ( node_position[1] > 48.0 && node_position[1] < 54.9 && node_position[0] > 78.9 ) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			else if ( node_position[1] > 47.0 && node_position[1] < 48.0 && node_position[0] > 77.1 ) //10, 46
			{
				MARK;
				node_position[0] = right_bound + 0.001;
			}
			if (fabs(node_position[0] -  78.5) < 0.01 && fabs(node_position[1] - 5.34049) < 0.01)
			{
				node_position[0] = 77.53;
				node_position[1] = 5.46055;
			}
			else if (node_position[1] < 5.5 && node_position[0] > 78.49) //10, 46
			{
				assert( p_this_node->IsBoundaryNode() );
				node_position[0] = right_bound + 0.001;
			}
			//else if (node_position[1] < 5.5 && node_position[1] > 5.0 && node_position[0] > 78.2)
			//{
				//PRINT_VARIABLE(node_position[0]);
				//PRINT_VARIABLE(node_position[1]);
			//	assert( p_this_node->IsBoundaryNode() );
			//	node_position[0] = right_bound + 0.001;
			//}
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
	
	    //MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.02));
        //simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }
};

#endif /*TESTVORONOIMESHV01_HPP_*/

    
	    

	



    
        


       

