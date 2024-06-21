#ifndef TESTONECELL_HPP_
#define TESTONECELLHPP_

//INCLUDE NECCESSARY FILES

// Default Chaste files (needed for everything)
#include "AbstractMesh.hpp" //Chaste file, Abstract base class for all meshes.
#include <boost/serialization/base_object.hpp> //Chaste file for boost
#include <boost/shared_ptr.hpp> //Chaste file, allows you to point to things
#include "ChasteSerialization.hpp" //Chaste file, This header is a wrapper including some of the Boost serialization library headers, along with a couple of standard C++ headers required to fix bugs in Boost.
#include "CheckpointArchiveTypes.hpp" //Chaste file, Intended for use by SerializationExportWrapper.hpp, by tests of archiving and code that needs to directly create an archive itself
#include <cxxtest/TestSuite.h> //Default needed for all tests, includes header files
#include "Debug.hpp" //Chaste file that allows you to debug and add in checks
#include "FakePetscSetup.hpp" //Chaste file that enforces running the file on only one process
#include "OutputFileHandler.hpp" //Chaste file, This file abstracts stuff that needs to be done when creating output files for tests.
#include "PetscSetupAndFinalize.hpp" // Chaste file This file is designed to be included by any test suites that use PETSc. It controls the PETSc initialisation and finalisation.
#include "SimulationTime.hpp" //Chaste file, Simulation time object stores the simulation time.
#include "SmartPointers.hpp" //Chaste file,  Includes the Boost shared_ptr smart pointer, and defines some useful macros to save typing when using it.

// Cell based model Chaste default files (needed for vertex models)
#include "AbstractCellBasedTestSuite.hpp" //Chaste file for running cell based models
#include "AbstractCellPopulation.hpp" //Chaste file, An abstract facade class encapsulating a cell population. Stores info about cells
#include "AbstractCellProperty.hpp" //Chaste file, Each cell has a collection of cell properties, which may express such concepts as mutation states, inherited labels, or per-cell data.
#include "ApoptoticCellProperty.hpp"//Chaste file, The ApoptoticCellProperty object keeps track of the number of cells that are apoptotic, as well as what colour should be used by the visualizer to display cells that are apoptotic.
#include "CellBasedSimulationArchiver.hpp" //Chaste file, CellBasedSimulationArchiver handles the checkpointing (saving and loading) of all the various AbstractCellBasedSimulation objects. 
#include "CellsGenerator.hpp" //Chaste file, A helper class for generating a vector of cells for a given mesh.
#include "CellLabel.hpp" //Chaste file, Each Cell owns a CellPropertyCollection, which may include a shared pointer to an object of this type. When a Cell that is labelled divides, the daughter cells are both labelled.
#include "MutableVertexMesh.hpp" // Chaste file to create a mutable vertex-based mesh class, which inherits from VertexMesh and allows for local remeshing.
#include "OffLatticeSimulation.hpp" //Chaste file, Run an off-lattice 2D or 3D cell-based simulation using an off-lattice cell population.
#include "RandomNumberGenerator.hpp" //Chaste file, A special singleton class allowing one to generate different types of random number
#include "VertexBasedCellPopulation.hpp" //Chaste file that defines vertex based cell popualtions
#include "CellIdWriter.hpp"
//Cell based model Chaste choices
#include "AbstractSimpleGenerationalCellCycleModel.hpp" //Chaste file for cell cycle model
#include "DifferentiatedCellProliferativeType.hpp" //Chaste file for defining differentiated type
#include "FarhadifarForce.hpp" //Chaste file, A force class for use in Vertex-based simulations.
#include "TransitCellProliferativeType.hpp" //Chaste file for defining proliferative types
#include "UniformG1GenerationalCellCycleModel.hpp" //Chaste file that defines cell cycle model
#include "WildTypeCellMutationState.hpp" //Chaste file, Subclass of AbstractCellMutationState defining a 'wild type' mutation state.

// Georgia specified files
#include "FionaVoronoiVertexMeshGenerator.hpp" // Georgia File that generates cell populations and initial condition
//#include "VoronoiVertexMeshGenerator.hpp"
//#include "ModifiedVoronoiVertexMeshGenerator.hpp" // Georgia File that generates cell populations and initial condition
//#include "MyApoptoticCellKiller.hpp" // Georgia File for defining cell death
#include "MyCellKiller.hpp" // Georgia File for defining cell death
#include "PlaneStickyBoundaryCondition.hpp" // Georgia File that imposes boundary condition
#include "VoronoiTargetAreaModifier.hpp" // Georgia file that updates target area modifier


//Fiona added files
#include "VolumeTrackingModifier.hpp"
//#include "FionaTrackGen.hpp"
//#include "CellProliferativePhasesCountWriter.hpp"
#include "FionaGenWriter.hpp"
//#include "CellAncestorWriter.hpp"

// CODE ACTUALLY STARTS BELOW
class TestOneCell : public AbstractCellBasedTestSuite
{
public:
// Test Starts Here
	void TestOneCellTest1()
    {
        //Generate Mesh cells across, cells up, rows of histoblasts on the bottom, number of relaxation steps, target area
		//double hb_edge_length = sqrt(2.0)/(sqrt(3.0*sqrt(3.0)));

		FionaVoronoiVertexMeshGenerator generator(13,30,10,0,1.0); // GEORGIA FILE
		//VoronoiVertexMeshGenerator generator(13,30,0,1.0); // GEORGIA FILE
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh(); // Define as a mutable 2D vertex mesh
        
		//Create cells to populate the mesh
        std::vector<CellPtr> cells;
        
		// Make pointers/short hand to point to each of these propoerties of cells
		MAKE_PTR(WildTypeCellMutationState, p_state); // Allows us to track the mutation state of cells
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type); // Differentiated cells do not divide
		MAKE_PTR(TransitCellProliferativeType, p_transit_type); // Transit cells can divide
		MAKE_PTR(CellLabel, p_label); // LECs are labelled cells
		

		// For each cell
		for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
    	{
			// Assign the cell cycle model to be a uniform G1 cell cycle
			UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;

			// Set the spatial dimension of each cell to be 2D
            p_model->SetDimension(2);

			// Set length of time that histoblast transit cells remains in G1
            p_model->SetTransitCellG1Duration(80);// This should be approx 8 times more than other 3 phases default=80


			// Set number of generations/divisions? How many generations a transit cell lives for before becoming fully differentiated.
            p_model->SetMaxTransitGenerations(1); // (1 results in target area of histoblasts being roughly right at end)

			// Birth time is the time that the cell was born (hours), can be negative if cell already alive befoore simulation
			// Set birth time to be random number (0,1) mutiplied by stem cell life cycle hours. OTHER 3 phases approx 10 time units.
            double birth_time = - RandomNumberGenerator::Instance()-> ranf() *
                               (  p_model->GetTransitCellG1Duration()
                                 + p_model->GetSG2MDuration() ); 
								 
			// Update cell pointers
            CellPtr p_cell(new Cell(p_state, p_model));
            //Set birth time to be value above
			p_cell->SetBirthTime(birth_time);
			
            // Add data to the vector (add cells to vector of cell data)
            cells.push_back(p_cell);

			

        }

		// Define a vertex based cell popualtion on the mesh using the cells above
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

		for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
				cell_iter != cell_population.End();
				++cell_iter)
		{
			c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);
			

			// Assign cells to LECs if in this area
			 //if (this_location(1) > 6.0 && this_location(1) < 59.0) 
			 if (this_location(1) > 15.0 && this_location(1) < 21.0 && this_location(0)>5.0 && this_location(0)<6)
			 {
			 	cell_iter->AddCellProperty(p_label); // Label the LECS 
                 cell_iter->SetApoptosisTime(DBL_MAX); // Set apoptosis time to be infinte to start with or cell dies too quickly leaving gaps	
                 cell_iter->SetCellProliferativeType(p_diff_type); // Set to be differentiated/non-proliferative		
			 }
			// // Assign other cells (histoblasys) 
             else
             {
                 //cell_iter->SetCellProliferativeType(p_transit_type);
				
                 cell_iter->SetCellProliferativeType(p_diff_type);
             }
			 // If whole cell outside the boundaries then kill the cells
			if (this_location(0) < -0.2 && (this_location(1) < 4.5 || this_location(1) > 60.5)) // all hbs on the LHS
			{
				cell_iter->Kill();
			}
			else if (this_location(0) > 26.5 && (this_location(1) < 4.5 || this_location(1) > 60.5)) // all hbs on the RHS
			{
			   	cell_iter->Kill();
			}
			else if (this_location(0) < -0.2 && this_location(1) > 4.0 && this_location(1) < 61.0) // specific hb/lec boundary cell
			{
			   	cell_iter->Kill();
			}
			else if (this_location(0) > 26.7 && this_location(1) > 4.0 && this_location(1) < 61.0) // specific hb/lec boundary cell
			{
			   	cell_iter->Kill();
			}
			
			else if (this_location(0) < 0 && this_location(1) > 60.5 && this_location(1) < 61.5)  // specific hb/lec boundary cell
			{
			   	cell_iter->Kill();
			}

			//if (this_location(0) < -1.5 || this_location(1) <-1.5) // all hbs on the LHS
			//{
			//	cell_iter->Kill();
			//}

			
			
			
			
		}
///////////////////// BELOW NOT PART OF INITIAL CONDITION /////////////////////
        OffLatticeSimulation<2> simulator(cell_population);
		
		// Set timestep details for simulatinos
        simulator.SetOutputDirectory("TestOneCell");
        simulator.SetSamplingTimestepMultiple(100); // 100 means each data value plotted is order 1 time unit, 1000 means each data value plotted in order 10 time units
		simulator.SetDt(0.01);
        simulator.SetEndTime(1000.0);
        
		// Set the force to be used by cells to be a Farhadifar Force (CHASTE DEFINED). Forces defined here but not implemented until area modifier?
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
	   
	   	// Update the target area of cells (GEORGIA DEFINED)
        MAKE_PTR(VoronoiTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

		// Include Volume Tracker - allows us to visualise/plot cell volumes/areas in paraview
		MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
	

		 // Set c and y boundary points
		 double top_height = 64.9 ; // 65 -0.2 so that top histoblasts are halfed
		 double bottom_height = 0.1;  //0 +0.2 so that bottom histoblasts are halfed
		 double right_bound = 26.5; 
		 double left_bound = 0.0;
		 unsigned num_nodes = cell_population.GetNumNodes();

		
		
		//iterate over vertices in the cell population, this section modifies the histoblast cells on the boundary of the LECs
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
		 			if (node_position[1] > 4.0 && node_position[1] < 6.0 )
		 			{
		 				node_position[1] -= 0.8;//0.85;//1.3;//1.35;//1.82; //-0.5  0.8 t=164, 0.9 t=57, 0.7 t=162
		 			}
		 			else if (node_position[1] > 59.0 && node_position[1] < 61.0)
		 			{
		 				node_position[1] += 0.8;//0.85;//1.3;//1.35;//1.82; 
		 			}
		 		}
				
				
		 		for (std::set<unsigned>::iterator iter = containing_element_indices.begin();
		 		iter != containing_element_indices.end();
		 		iter++)
		 		{
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
		  		if (node_position[1] > 4.0 && node_position[1] < 6.0 ) 
		  		{
		  			node_position[1] -= 0.8;//0.85;//1.3;//1.35;//1.82; 
		  		}
		  		else if (node_position[1] > 59.0 && node_position[1] < 61.0) 
		  		{
		  			node_position[1] += 0.8;//0.85;//1.3;//1.35;//1.82; 
		  		}
		  	}

			if (node_position[1] > 59.0 && node_position[1] < 61.0)
		 			{
		 				node_position[1] += 0.2;//1.25;//1.2;//0.45;//0.44;   0 = t=104, 0.3 t=56, 0.2 t=164, 0.1 t=151
						
		 			}
					else if (node_position[1] > 4.0 && node_position[1] < 5.0)
		 			{
		 				node_position[1] -= 0.2;//1.2;//0.45;//0.44; 
						
		 			}


						if (node_position[1] > 60.0 && node_position[1] < 60.8 && node_position[0]<0.01)
		 			{
		 				node_position[1] -= 1.0;//1.25;//1.2;//0.45;//0.44;   0 = t=104, 0.3 t=56, 0.2 t=164, 0.1 t=151
						
		 			}
					
			
		
		// //  	//PRINT_VARIABLE(node_position)
			}
        
		

		
		  // Iterate over vertices in the cell population and ensure cells near the boundary are set as boundary nodes
		 for (unsigned node_index=0; node_index<num_nodes; node_index++)
		  {
		  	Node<2>* p_this_node = cell_population.GetNode(node_index);
		  	c_vector<double, 2>& node_position = p_this_node->rGetModifiableLocation();

			if ( node_position[1] > top_height)
			{
				p_this_node->SetAsBoundaryNode(true); 
				node_position[1] = top_height + 0.001;
			}
			if ( node_position[1] < bottom_height )// no difference if set as 0.4
			{ 
				p_this_node->SetAsBoundaryNode(true);
			    node_position[1] = bottom_height - 0.001;
			}

			if ( node_position[0] < left_bound+0.01) //needs to be 0.01
			{
				p_this_node->SetAsBoundaryNode(true); 
				node_position[0] = left_bound - 0.001;
			}

			if ( node_position[0] > right_bound-0.01) //needs to be 0.01
			{
				p_this_node->SetAsBoundaryNode(true); 
				node_position[0] = right_bound + 0.001;
			}
			
		
		
		  }

		  // Boundary conditions
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

		 // // LEFT
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

		// // //Only form of cell death is from MyCellKiller. The input probability is the probability that in an hour's worth of trying, the cell killer will have successfully killed a given cell.
         MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.008));//default = 0.008
         simulator.AddCellKiller(p_cell_killer);
	
		//cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
		cell_population.AddCellWriter<FionaGenWriter>();
		cell_population.AddCellWriter<CellIdWriter>();
		
		
		// Update cell populations
	    cell_population.RemoveDeadCells();
		cell_population.Update();
		
	
        simulator.Solve();

		
	
    }
		
    	



};
// END OF FILE
#endif /*TESTONECELL_HPP_*/