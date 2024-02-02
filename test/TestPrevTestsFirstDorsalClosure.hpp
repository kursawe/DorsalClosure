#ifndef TESTFIRSTDORSALCLOSUREHELLO_HPP_
#define TESTFIRSTDORSALCLOSUREHELLO_HPP_

//#include "FirstDorsalClosure.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"

#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
//#include "TargetAreaLinearGrowthModifier.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbstractMesh.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
//#include "StemCellProliferativeType.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include <cxxtest/TestSuite.h>
#include "Debug.hpp"
#include "WildTypeCellMutationState.hpp"
//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
#include "CheckpointArchiveTypes.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "AbstractCellPopulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"
#include "OutputFileHandler.hpp"
#include "MutableVertexMesh.hpp"
#include "FixedVertexBasedDivisionRule.hpp"
//#include "GammaG1CellCycleModel.hpp"
#include "RandomCellKiller.hpp"
#include "MyCellKiller.hpp"

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "FakePetscSetup.hpp"

//#include "RandomNumberGenerator.hpp"

class TestFirstDorsalClosure : public AbstractCellBasedTestSuite
{
public:
    //test whether LECs are labelled - show up as different colour in paraview 
    void xTestCellLabel()
    {
        
        HoneycombVertexMeshGenerator generator(6, 6);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        MAKE_PTR(CellLabel, p_label);
	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 3.6)
		    {
			    cell_iter->AddCellProperty(p_label);
                
	        }

	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_label");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);
        
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 4;
        uppery_normal(1) = 1.0;

        simulator.Solve();
    }
    //test whether the cells are bounded by a 'box'
    void xTestBoundary()
    {
        
        HoneycombVertexMeshGenerator generator(6, 6);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);    
         
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_boundary");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 6;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.5;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.5;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);
        
        simulator.Solve();
        MARK;
    }
    //test labelled LECs and boundary together 
    void xTestLabelandBoundary()
    {
        
        HoneycombVertexMeshGenerator generator(6, 6);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        MAKE_PTR(CellLabel, p_label);
	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 3.6)
		    {
			    cell_iter->AddCellProperty(p_label);
                
	        }

	    }
        
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.SetOutputDirectory("TestFirstDorsalClosure_labelandboundary"); 
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 5.5;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);
         
        simulator.Solve();
        MARK;
    }
    //test hbs dividing and box boundary 
    void xTestProliferationandBoundary()
    {
        
        HoneycombVertexMeshGenerator generator(6, 6);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        MAKE_PTR(CellLabel, p_label);
	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 3.6)
		    {
			    cell_iter->AddCellProperty(p_label);
                
	        }
            else
            {
                //cell_iter->AddCellProperty(p_transit_type);
                //only one round of division 
                cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
				cell_iter->SetCellProliferativeType(p_transit_type);
            }

	    }
        
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.SetOutputDirectory("TestFirstDorsalClosure_prolifandboundary");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);
         
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 5.5;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);
        
        simulator.Solve();
    }
    //test LECs dying and box boundary 
    void xTestApoptosisandBoundary()
    {
        
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        MAKE_PTR(CellLabel, p_label);
	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7.0)
		    {
			    cell_iter->AddCellProperty(p_label);
                cell_iter->StartApoptosis(false);

	        }
            else
            {
                cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
				cell_iter->SetCellProliferativeType(p_transit_type);
            }

	    }



        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_apoptosisandboundary");
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(15.0);
        
 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        simulator.Solve();
    }
    //test more than one round of cell division 
    void xTestIncreasedProliferation()
    {
        
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements();
			++cell_iter)
	    {
		    UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);

	    }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_increasedproliferation");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(50.0);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();


    }
    //test hbs proliferating and LEC death 
    void xTestProliferationAndDeath()
    {
        
        HoneycombVertexMeshGenerator generator(6, 7);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);

        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        //VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        //c_vector<double,2> cell_division_axis;
        //cell_division_axis[0] = 1.0;
        //cell_division_axis[1] = 0.0;
        //MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));

        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>29)
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(0.005);
                p_model->SetTransitCellG1Duration(0.005);

                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);
            }
            else
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->SetApoptosisTime(20);
                cells.push_back(p_cell);

            }
        }
	    

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 4.4)
		    {
			    //cell_iter->AddCellProperty(p_label);
                cell_iter->StartApoptosis(false);
                
	        }

	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_proliferationanddeath");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(25.0);
        
 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 6.4;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        simulator.Solve();
    }

    //test changing parameters to deal with hourglass shape 
    void TestBoundaryTension()
    {
        //setting up the mesh that the elements (cells) lie on
        //6 cells across, 10 up 
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        //creating cells
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);

        //trying to label cells as apoptotic
        //boost::shared_ptr<AbstractCellProperty> p_apoptotic(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
        
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        //trying to have directed proliferation
        //c_vector<double,2> cell_division_axis;
        //cell_division_axis[0] = 1.0;
        //cell_division_axis[1] = 0.0;
        //MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));

        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>47)
            //histoblasts
            {
                //proliferation 
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(0.05);
                p_model->SetTransitCellG1Duration(0.05);

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

                //trying to label lecs as apoptotic
                //p_cell->AddCellProperty(p_apoptotic);
                //p_cell->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());

                p_cell->SetApoptosisTime(80);

                //trying to set target area to 0 instead of using SetApoptosisTime
                p_cell->GetCellData()->SetItem("target area", 0.0);

                cells.push_back(p_cell);

            }
        }
	    
        //creating cell population 
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //trying to have cells die randomly 
        //RandomCellKiller<2> random_cell_killer(&cell_population, 0.5);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7)
            //LECs apoptosis
		    {
			    cell_iter->AddCellProperty(p_label);
                cell_iter->StartApoptosis(true);
                
	        }

	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_boundarytension");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(0.5);
        simulator.SetDt(0.1);

        //changed parameters are to try and deal with 'hourglass' shape 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(0.3);      //0.12 default 
        p_force->SetPerimeterContractilityParameter(0.01);  //0.04 default
        p_force->SetAreaElasticityParameter(0.8);            //1.0 default 
        p_force->SetLineTensionParameter(0.06);               //0.12 default
        simulator.AddForce(p_force);
       
       //required for Farhadifar force 
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        simulator.Solve();
    }
    //test for proper shape of tissue, apoptosis of LECs, hb proliferation and box boundary
    void xTestClosure()
    {
        
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);

        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        //VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        //c_vector<double,2> cell_division_axis;
        //cell_division_axis[0] = 1.0;
        //cell_division_axis[1] = 0.0;
        //MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));

        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>47)
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(0.05);
                p_model->SetTransitCellG1Duration(0.05);

                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);
            }
            else
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->SetApoptosisTime(80);
                //p_cell->GetCellData()->SetItem("target area", 0.0);
                cells.push_back(p_cell);

            }
        }
	    

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        //RandomCellKiller<2> random_cell_killer(&cell_population, 0.5);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7)
		    {
			    //cell_iter->AddCellProperty(p_label);
                cell_iter->StartApoptosis(true);
                
	        }

	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_closure");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(50.0);
        
 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        simulator.Solve();
    }
    //trying to understand the proliferation code
    void xTestUnderstandingProliferation()
    {
        
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(CellLabel, p_label);

        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>47)
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(5);
                //p_model->SetTransitCellG1Duration(5);

                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);
            }
            else
            {
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                //The spatial dimension (1, 2 or 3) needs to be set on the cell-cycle model before it is passed to the cell.
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                /*
                * We also alter the default cell-cycle times.
                */
                p_model->SetStemCellG1Duration(0.05);
                p_model->SetTransitCellG1Duration(0.05);

                /*
                * We now define a random birth time, chosen from [-T,0], where
                * T = t,,1,, + t,,2,,, where t,,1,, is a parameter representing the G,,1,, duration
                * of a 'stem' cell, and t,,2,, is the basic S+G,,2,,+M phases duration...
                */
                double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() );

                /*
                * ...then we set the birth time and push the cell back into the vector
                * of cells.
                */
                p_cell->SetBirthTime(birth_time);
                cells.push_back(p_cell);

            }
        }
         VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7)
		    {
			    cell_iter->AddCellProperty(p_label);

                
	        }

	    }
        
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestFirstDorsalClosure_understandingproliferation");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(50.0);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();


    }
    void xTestRandomLecDeaths()
    {
        //setting up the mesh that the elements (cells) lie on
        //6 cells across, 10 up 
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        //creating cells
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);

        //trying to label cells as apoptotic
        //boost::shared_ptr<AbstractCellProperty> p_apoptotic(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
   
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        //trying to have directed proliferation
        //c_vector<double,2> cell_division_axis;
        //cell_division_axis[0] = 1.0;
        //cell_division_axis[1] = 0.0;
        //MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));


        //deathtype used to either labels the lec so it dies actively or else it dies passively
        unsigned deathtype;
        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>47)
            //histoblasts
            {
                //proliferation 
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(0.05);
                p_model->SetTransitCellG1Duration(0.05);

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

                //p_label used to allow for random deaths within the lecs
                //instead of just passive constriction from the proliferating hbs
                deathtype = RandomNumberGenerator::Instance()->randMod(101);
                //PRINT_VARIABLE(deathtype);
                if (deathtype<41)
                {
                    p_cell->AddCellProperty(p_label);
                }
                else
                {
                    //MARK;
                    p_cell->SetApoptosisTime(80);
                }

                //trying to label lecs as apoptotic
                //p_cell->AddCellProperty(p_apoptotic);
                //p_cell->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());


                //trying to set target area to 0 instead of using SetApoptosisTime
                //p_cell->GetCellData()->SetItem("target area", 0.0);

                cells.push_back(p_cell);

            }
        }
	    
        //creating cell population 
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //trying to have cells die randomly 
        //RandomCellKiller<2> random_cell_killer(&cell_population, 0.5);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7)
            //LECs apoptosis
		    {
                //p_label used to allow for random deaths within the lecs
                //instead of just passive constriction from the proliferating hbs
                //deathtype = RandomNumberGenerator::Instance()->randMod(4);
                //PRINT_VARIABLE(deathtype);
                if (!cell_iter->template HasCellProperty<CellLabel>())
                {
                    //MARK;
                    cell_iter->StartApoptosis(true);
                }
                
	        }

	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("closurev0_randomlecdeaths");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(57.0);
        
        //changed parameters are to try and deal with 'hourglass' shape 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(0.3);      //0.12 default 
        p_force->SetPerimeterContractilityParameter(0.01);  //0.04 default
        p_force->SetAreaElasticityParameter(0.8);            //1.0 default 
        p_force->SetLineTensionParameter(0.06);               //0.12 default
        simulator.AddForce(p_force);
       
       //required for Farhadifar force 
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.045));
        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }
    void xTestRandomDeath()
    {
        //setting up the mesh that the elements (cells) lie on
        //6 cells across, 10 up 
        HoneycombVertexMeshGenerator generator(6, 10);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        //creating cells
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);

        //trying to label cells as apoptotic
        //boost::shared_ptr<AbstractCellProperty> p_apoptotic(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
   
        //cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        //trying to have directed proliferation
        //c_vector<double,2> cell_division_axis;
        //cell_division_axis[0] = 1.0;
        //cell_division_axis[1] = 0.0;
        //MAKE_PTR_ARGS(FixedVertexBasedDivisionRule<2>, p_div_rule, (cell_division_axis));


        for (unsigned cell_iter=0; cell_iter<p_mesh->GetNumElements(); ++cell_iter)
        {
            if (cell_iter<12 || cell_iter>47)
            //histoblasts
            {
                //proliferation 
                UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel;
                p_model->SetDimension(2);
                CellPtr p_cell(new Cell(p_state, p_model));
                p_cell->SetCellProliferativeType(p_stem_type);

                p_model->SetStemCellG1Duration(0.05);
                p_model->SetTransitCellG1Duration(0.05);

                double birth_time = - 2* (RandomNumberGenerator::Instance()->ranf() *
                                    (  p_model->GetStemCellG1Duration()
                                    + p_model->GetSG2MDuration() ));

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
                p_cell->SetApoptosisTime(10);

                //trying to label lecs as apoptotic
                //p_cell->AddCellProperty(p_apoptotic);
                //p_cell->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());


                //trying to set target area to 0 instead of using SetApoptosisTime
                //p_cell->GetCellData()->SetItem("target area", 0.0);

                //p_cell->GetCellData()->SetItem("target area", 0.25);
                cells.push_back(p_cell);

            }
        }
	    
        //creating cell population 
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        //trying to have cells die randomly 
        //RandomCellKiller<2> random_cell_killer(&cell_population, 0.5);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) > 1.8 && this_location(1) < 7)
            //LECs apoptosis
		    {
                //p_label used to allow for random deaths within the lecs
                //instead of just passive constriction from the proliferating hbs
                //deathtype = RandomNumberGenerator::Instance()->randMod(4);
                //PRINT_VARIABLE(deathtype);
                if (!cell_iter->template HasCellProperty<CellLabel>())
                {
                    //MARK;
                    cell_iter->StartApoptosis(true);
                }
                
	        }
            else
            {
                cell_iter->ReadyToDivide();
            }
	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("closurev0_randomdeaths");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(80.0);
        
        //changed parameters are to try and deal with 'hourglass' shape 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(0.3);      //0.12 default 
        p_force->SetPerimeterContractilityParameter(0.01);  //0.04 default
        p_force->SetAreaElasticityParameter(0.8);            //1.0 default 
        p_force->SetLineTensionParameter(0.06);               //0.12 default
        simulator.AddForce(p_force);
       
       //required for Farhadifar force 
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.045));
        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }

};
#endif /*TESTFIRSTDORSALCLOSUREHELLO_HPP_*/