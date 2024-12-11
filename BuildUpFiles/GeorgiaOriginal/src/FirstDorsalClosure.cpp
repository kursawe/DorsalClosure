#include "FirstDorsalClosure.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include <CellId.hpp>
//#include "VertexModelDataWriter.hpp"
#include <CellLabel.hpp>
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbstractMesh.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellBasedSimulationArchiver.hpp"

FirstDorsalClosure::FirstDorsalClosure()
	: mInitialTargetArea(1.0)
{
}

FirstDorsalClosure::~FirstDorsalClosure()
{
}

void FirstDorsalClosure::SetupSimulation()
{
	SetupSingletons();

    CreateCellPopulation();

    // make the simulation

	// Note: we need to initialise the boost pointers in this way - if we use the macro it doesn't work because it won't overwrite the old pointer.
	// Remember that the pointer is being initialised to NULL in the constructor automatically for us, some neat boost magic.
    mpSimulator.reset(new OffLatticeSimulation<2>(*mpCellPopulation));

	mpSimulator->SetSamplingTimestepMultiple(200);

	// Make the Farhadifar force
	MAKE_PTR(FarhadifarForce<2>, p_force);

	// before passing the force to the simulation
	mpSimulator->AddForce(p_force);

    MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
    mpSimulator->AddSimulationModifier(p_growth_modifier);

	// set the time step
	mpSimulator->SetDt(0.01);

	mpSimulator->SetEndTime(40.0);
}

void FirstDorsalClosure::SetupSingletons()
{
	// Set up what the test suite would do
	SimulationTime::Instance()->SetStartTime(0.0);
	CellPropertyRegistry::Instance()->Clear();
	CellId::ResetMaxCellId();
}

void FirstDorsalClosure::DestroySingletons()
{
	// this is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
}

MutableVertexMesh<2,2>* FirstDorsalClosure::CreateMesh()
{
    mpMeshGenerator.reset( new HoneycombVertexMeshGenerator(6, 6) );
    MutableVertexMesh<2,2>* p_mesh = mpMeshGenerator->GetMesh();
    return p_mesh;
}

void FirstDorsalClosure::CreateCellPopulation()
{
	std::vector<CellPtr> cells;

    MutableVertexMesh<2,2>* p_mesh;

	p_mesh = this->CreateMesh();
	
	MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
	CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;

	cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

    mpCellPopulation.reset(new VertexBasedCellPopulation<2>(*p_mesh, cells));

	// Here we define the lecs
	MAKE_PTR(CellLabel, p_label);
	for (typename AbstractCellPopulation<2>::Iterator cell_iter = mpCellPopulation->Begin();
			cell_iter != mpCellPopulation->End();
			++cell_iter)
	{
		c_vector<double, 2> this_location = mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

		if (this_location(1) > 1.8 && this_location(1) < 3.6)
		{
			cell_iter->AddCellProperty(p_label);
	    }

	}

	//lets make a writer that outputs polygon classes and other vertex based data and add it to the population
	//mpCellPopulation->AddCellWriter<VertexModelDataWriter>();
}

void FirstDorsalClosure::ApplyBoundaries()
{
	this->MakeAndApplyBoundaryConditions();
}

void FirstDorsalClosure::MakeAndApplyBoundaryConditions()
{

	// BOTTOM

    c_vector<double,2> lowery_point = zero_vector<double>(2);
    c_vector<double,2> lowery_normal = zero_vector<double>(2);
    lowery_point(1) = 0.5;
    lowery_normal(1) = -1.0;
    
    //make a pointer to a PlaneBoundaryCondition, pass in point and normal
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, lowery_p_boundarycondition, (mpCellPopulation.get(), lowery_point, lowery_normal));
    lowery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
    mpSimulator->AddCellPopulationBoundaryCondition(lowery_p_boundarycondition);

	// TOP
	
    c_vector<double,2> uppery_point = zero_vector<double>(2);
    c_vector<double,2> uppery_normal = zero_vector<double>(2);
    uppery_point(1) = 4;
    uppery_normal(1) = 1.0;

    //make a pointer to a PlaneBoundaryCondition, pass in point and normal
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, uppery_p_boundarycondition, (mpCellPopulation.get(), uppery_point, uppery_normal));
    uppery_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
    mpSimulator->AddCellPopulationBoundaryCondition(uppery_p_boundarycondition);

	// LEFT
	
    c_vector<double,2> leftx_point = zero_vector<double>(2);
    c_vector<double,2> leftx_normal = zero_vector<double>(2);
    leftx_point(0) = 0.5;
    leftx_normal(0) = -1.0;
    
    //make a pointer to a PlaneBoundaryCondition, pass in point and normal
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, leftx_p_boundarycondition, (mpCellPopulation.get(), leftx_point, leftx_normal));
    leftx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
    mpSimulator->AddCellPopulationBoundaryCondition(leftx_p_boundarycondition);

	// RIGHT
	
    c_vector<double,2> rightx_point = zero_vector<double>(2);
    c_vector<double,2> rightx_normal = zero_vector<double>(2);
    rightx_point(0) = 6;
    rightx_normal(0) = 1.0;

    //make a pointer to a PlaneBoundaryCondition, pass in point and normal
    MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, rightx_p_boundarycondition, (mpCellPopulation.get(), rightx_point, rightx_normal));
    rightx_p_boundarycondition->SetUseJiggledNodesOnPlane(true);
    mpSimulator->AddCellPopulationBoundaryCondition(rightx_p_boundarycondition);


}
/*
void FirstDorsalClosure::AddSimulationModifier(boost::shared_ptr<AbstractCellBasedSimulationModifier<2,2> > pSimulationModifier)
{
	mpSimulator->AddSimulationModifier(pSimulationModifier);
}
*/
boost::shared_ptr< OffLatticeSimulation<2> > FirstDorsalClosure::GetSimulation()
{
	return mpSimulator;
}

void FirstDorsalClosure::SetOutputDirectory(std::string outputDirectory)
{
	mpSimulator->SetOutputDirectory(outputDirectory);
}

void FirstDorsalClosure::RunAndDestroySimulation()
{
	mpSimulator->Solve();
	CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(mpSimulator.get());
	DestroySimulation();
}

void FirstDorsalClosure::DestroySimulation()
{
	// Delete all the simulation data from memory before running next mutant.
	mpSimulator.reset();
	mpCellPopulation.reset();
	mpMeshGenerator.reset();
	DestroySingletons();
}
/*
void FirstDorsalClosure::SetRandomInitialConditionLogical(bool randomInitialConditionLogical)
{
	mApplyRandomInitialCondition = randomInitialConditionLogical;
}

bool FirstDorsalClosure::GetRandomInitialConditionLogical()
{
	return mApplyRandomInitialCondition;
}
*/
void FirstDorsalClosure::SetInitialTargetArea(double initialTargetArea)
{
	mInitialTargetArea = initialTargetArea;
}

double FirstDorsalClosure::GetInitialTargetArea()
{
	return mInitialTargetArea;
}
