#ifndef TESTNOVOIDS_HPP_
#define TESTNOVOIDS_HPP_

//#include "FirstDorsalClosure.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "MyTargetAreaModifier.hpp"
#include "MyApoptoticCellKiller.hpp"
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
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
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

#include <TimeStepper.hpp>

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "FakePetscSetup.hpp"
#include "PlaneStickyBoundaryCondition.hpp"
#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"


class TestNoVoids : public AbstractCellBasedTestSuite
{
public:
//all cell removals are via a T2 swap
    void TestVoids()
    {
        //setting up the mesh that the elements (cells) lie on
        //6 cells across, 10 up 
        HoneycombVertexMeshGenerator generator(6, 10, false, 0.01, 0.001, 0.5*sqrt(3.0));

        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        //creating cells
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);


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
                
                p_cell->SetApoptosisTime(DOUBLE_UNSET); 
                //p_cell->SetApoptosisTime(20);  

                cells.push_back(p_cell);

            }
        }
	    
        //creating cell population 
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) < 1.8 && this_location(1) > 7)
            //LECs apoptosis
		    {
                //cell_iter->ReadyToDivide();
                
            }
	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("test_void");
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
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        //uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        //lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        //leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        //rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.046));
        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }
//trying to improve size of hbs 
    void xTestSize()
    {

        
        //setting up the mesh that the elements (cells) lie on
        //6 cells across, 10 up 
        HoneycombVertexMeshGenerator generator(6, 10, false, 0.01, 0.001, 0.5*sqrt(3.0));

        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();
        
        //creating cells
        std::vector<CellPtr> cells;
        
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(CellLabel, p_label);


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

                p_model->SetStemCellG1Duration(48);
                p_model->SetTransitCellG1Duration(5);
                p_model->SetMaxTransitGenerations(3); 
                //SG2MDuration = 10
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
                
                p_cell->SetApoptosisTime(DOUBLE_UNSET); 
                //p_cell->SetApoptosisTime(20);  

                cells.push_back(p_cell);

            }
        }
	    
        //creating cell population 
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			cell_iter != cell_population.End();
			++cell_iter)
	    {
		    c_vector<double, 2> this_location = cell_population.GetLocationOfCellCentre(*cell_iter);

		    if (this_location(1) < 1.8 && this_location(1) > 7)
            //LECs apoptosis
		    {
                //cell_iter->ReadyToDivide();
                
	        }
	    }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("test_size");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(97.0);
        
        //changed parameters are to try and deal with 'hourglass' shape 
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(0.3);      //0.12 default 
        p_force->SetPerimeterContractilityParameter(0.01);  //0.04 default
        p_force->SetAreaElasticityParameter(0.8);            //1.0 default 
        p_force->SetLineTensionParameter(0.06);               //0.12 default
        simulator.AddForce(p_force);
       
       //required for Farhadifar force 
        MAKE_PTR(MyTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        //upper boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> uppery_point = zero_vector<double>(2);
        c_vector<double,2> uppery_normal = zero_vector<double>(2);
        uppery_point(1) = 9.0;
        uppery_normal(1) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, uppery_p_boundarycondition, (&cell_population, uppery_point, uppery_normal));
        //uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

        //lower boundary
       
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> lowery_point = zero_vector<double>(2);
        c_vector<double,2> lowery_normal = zero_vector<double>(2);
        lowery_point(1) = 0.0;
        lowery_normal(1) = -1.0;
        
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, lowery_p_boundarycondition, (&cell_population, lowery_point, lowery_normal));
        //lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

        //left boundary

        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> leftx_point = zero_vector<double>(2);
        c_vector<double,2> leftx_normal = zero_vector<double>(2);
        leftx_point(0) = 0.0;
        leftx_normal(0) = -1.0;
    
        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, leftx_p_boundarycondition, (&cell_population, leftx_point, leftx_normal));
        //leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

        //right boundary
    
        //define point on the plane boundary and a normal to the plane 
        c_vector<double,2> rightx_point = zero_vector<double>(2);
        c_vector<double,2> rightx_normal = zero_vector<double>(2);
        rightx_point(0) = 6.5;
        rightx_normal(0) = 1.0;

        //make a pointer to a PlaneBoundaryCondition, pass in point and normal
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, rightx_p_boundarycondition, (&cell_population, rightx_point, rightx_normal));
        //rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);

        MAKE_PTR_ARGS(MyCellKiller<2>, p_cell_killer, (&cell_population, 0.046));
        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();
    }
};

#endif /*TESTNOVOIDS_HPP_*/