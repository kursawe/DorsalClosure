#ifndef FIRSTDORSALCLOSURE_HPP_
#define FIRSTDORSALCLOSURE_HPP_

#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"

class FirstDorsalClosure
{

protected:

	/**
	 * The name that describes this mutant. Will be used by the scenario to set the simulation name.
	 */
	//std::string mMutantName;


	/**
	 * The initial target area of cells in the tissue.
	 */
	double mInitialTargetArea;

	/**
	 * The mesh generator. It should only be deleted once the simulation is over, otherwise we loose the mesh.
	 * This is a workaround, instead we should make it optional for the generator not to delete the mesh and delete
	 * it ourselves.
	 */
	boost::shared_ptr< HoneycombVertexMeshGenerator > mpMeshGenerator;


	/**
	 * The population.
	 */
	boost::shared_ptr< VertexBasedCellPopulation<2> > mpCellPopulation;

	/**
	 * The simulation that we would like to run for this mutant
	 */
	boost::shared_ptr< OffLatticeSimulation<2> > mpSimulator;

    /**
     * Helper method for RunAndDestroySimulation(), frees memory for this simulation.
     */
    void DestroySimulation();

    /**
     * Helper method for CreateCellPopulation. Creates mesh.
     *
     */
    virtual MutableVertexMesh<2,2>* CreateMesh();

    /**
     * Helper method for SetupSimulation. Designs the mesh and sets up the cell population
     *
     * @returns the population with cells correctly set up and the mesh correctly defined.
     */
    void CreateCellPopulation();

    /**
     * Helper method for CreateCellPopulation.
     *
     * Makes all the boundaries straight lines and apply sticky boundary conditions.
     */
    void ApplyBoundaries();

    /**
     * Helper method for ApplyStraightBoundaries.
     *
     * Adds the boundary condition that keeps the boundaries straight.
     */
    void MakeAndApplyBoundaryConditions();

    /**
     * Helper method for SetupSimulation. Executes any setup code that is mutant specific.
     *
     * @returns the population with cells correctly set up and the mesh correctly defined.
     */
    //virtual void MutantSpecificSimulationSetup()=0;

public:

    /**
     * Constructor. Does nothing.
     */
    FirstDorsalClosure();

    /**
     * Destructor.
     */
    virtual ~FirstDorsalClosure();

    /**
     * Method to be called by parent program to set up the simulation.
     */
    void SetupSimulation();

    /**
     * Add a SimulationModifier to be used by this mutant. This is intended for the scenario to further specify the simulation.
     *
     * @param pSimulationModifier pointer to a SimulationModifier
     */
    //void AddSimulationModifier(boost::shared_ptr<AbstractCellBasedSimulationModifier<2,2> > pSimulationModifier);

    /**
     * Set the output directory of the mutant, to be done by the scenario.
     *
     * @param outputDirectory the output directory to use
     */
    void SetOutputDirectory(std::string outputDirectory);

    /**
     * Helper method for SetupSimulation(), sets up the singletons for the mutant simulation.
     * This would usually be done by the setup method of the abstract cell based test suite.
     */
    static void SetupSingletons();

    /**
     * Helper method for RunAndDestroySimulation(), destroys the singletons from the mutant simulation.
     * This would usually be done by the tearDown method of the abstract cell based test suite.
     */
    static void DestroySingletons();

    /**
     * Get the simulation for this mutant to be able to modify it.
     *
     * @returns pointer to the simulation modifier for this mutant.
     */
    boost::shared_ptr< OffLatticeSimulation<2> > GetSimulation();

    /**
     * Run the simulation for this mutant.
     */
    void RunAndDestroySimulation();

    /**
     * The name of the mutant is used for the output directory structure.
     *
     * @returns the name of the mutant
     */
    std::string GetName();

    /**
     * Set the way we do boundary conditions
     *
     * @params straightBoundaryLogical
     */
    //void SetStraightBoundaryLogical(bool straightBoundaryLogical);

    /**
     * Get the way we do boundary conditions
     *
     * @returns mApplyStraightBoundaries
     */
    //bool GetStraightBoundaryLogical();

    /**
     * Set the way we generate initial conditions
     *
     * @params randomInitialConditionLogical
     */
    //void SetRandomInitialConditionLogical(bool randomInitialConditionLogical);

    /**
     * Get the way we generate initial conditions
     *
     * @returns mApplyRandomInitialCondition
     */
    //bool GetRandomInitialConditionLogical();

    /**
     * Set the initial target area of cells in the tissue.
     *
     * @params initialTargetArea
     */
    void SetInitialTargetArea(double initialTargetArea);

    /**
     * Get the initial target area of cells in the tissue.
     *
     * @returns mInitialTargetArea
     */
    double GetInitialTargetArea();


};

#endif /*FIRSTDORSALCLOSURE*/