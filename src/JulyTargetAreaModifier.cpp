#include "JulyTargetAreaModifier.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Debug.hpp"
#include "CellLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "Debug.hpp"
#include "FarhadifarForce.hpp"
template<unsigned DIM>
JulyTargetAreaModifier<DIM>::JulyTargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>(),
      mGrowthDuration(DOUBLE_UNSET),
      mApoptosisDuration(20.0), //20
      mLecArea(75.0)      //  
{
}

template<unsigned DIM>
JulyTargetAreaModifier<DIM>::~JulyTargetAreaModifier()
{
}


template<unsigned DIM>
void JulyTargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)

{
    double cell_target_area;
    double lec_area = mLecArea;

    if (pCell->IsDead())
    {
         PRINT_VARIABLE(pCell->GetCellId());
         if (pCell->HasCellProperty<CellLabel>())
            {
                TRACE("LEC dead")
            }
        else
            {
            TRACE("Histoblast dead")
            }
    }
   
    
    // Get target area A of a healthy cell in S, G2 or M phase
    if (pCell->HasCellProperty<CellLabel>())
    {
        //lecs are labelled, and are given larger target area
        cell_target_area = lec_area;
    }
    else
    {
        //hbs given target area = 1
        cell_target_area = 1*this->mReferenceTargetArea; 
    }
    double growth_duration = mGrowthDuration;
    if (growth_duration == DOUBLE_UNSET)
    {
        if (dynamic_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel()) == nullptr)
        {
            EXCEPTION("If SetGrowthDuration() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");
        }
        AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel());

        growth_duration = p_model->GetG1Duration();

        // If the cell is differentiated then its G1 duration is infinite
        if (growth_duration == DBL_MAX)
        {
          //  MARK;
            // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
            growth_duration = p_model->GetTransitCellG1Duration();      // = 2 
        }
    }
    double apoptosis_duration = mApoptosisDuration;
  if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        //PRINT_2_VARIABLES(pCell->GetCellId(), pCell->GetStartOfApoptosisTime());
        // Age of cell when apoptosis begins is given by pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime()
        //maybe - cell target area needs to be adjusted here to allow for the cells which have not reached the target area before they are labelled as apoptotic 
        //but in TestClosureV0 it is longer than 2 time units as there are other forces opposing the growth to the target area size//
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < growth_duration)
        {
            //MARK;
            cell_target_area *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/growth_duration);

      

        }

        

        //current time minus time when the cell became apoptotic 
        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();
        cell_target_area *= 1.0 - (0.5*(time_spent_apoptotic/apoptosis_duration));
        //area is a positive quantity 

        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
        //unsigned num_lecs = CellPropertyRegistry::Instance()->Get<CellLabel>()->GetCellCount();
		//PRINT_VARIABLE(num_lecs);
    }
    else
    {

        if (pCell->HasCellProperty<CellLabel>())
        {
            cell_target_area=lec_area;
        }
        else
        {


        double cell_age = pCell->GetAge();

        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
        if (cell_age < growth_duration)
        {
            cell_target_area *= 0.5*(1 + cell_age/growth_duration);
        }
        else
        {
            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_area = 0.5*this->mReferenceTargetArea;
            }
        }
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
   
}

template<unsigned DIM>
double JulyTargetAreaModifier<DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned DIM>
void JulyTargetAreaModifier<DIM>::SetGrowthDuration(double growthDuration)
{
    assert(growthDuration >= 0.0);
    mGrowthDuration = growthDuration;
}

template<unsigned DIM>
double JulyTargetAreaModifier<DIM>::GetApoptosisDuration()
{
    return mApoptosisDuration;
}

template<unsigned DIM>
void JulyTargetAreaModifier<DIM>::SetApoptosisDuration(double apoptosisDuration)
{
    assert(apoptosisDuration >= 0.0);
    mApoptosisDuration = apoptosisDuration;
}

template<unsigned DIM>
double JulyTargetAreaModifier<DIM>::GetLecArea()
{
    return mLecArea;
}

template<unsigned DIM>
void JulyTargetAreaModifier<DIM>::SetLecArea(double lecArea)
{
    assert(lecArea >= 0.0);
    mLecArea = lecArea;
}
template<unsigned DIM>
void JulyTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    *rParamsFile << "\t\t\t<ApoptosisDuration>" << mApoptosisDuration << "</ApoptosisDuration>\n";
    *rParamsFile << "\t\t\t<LecArea>" << mLecArea << "</LecArea>\n";

    // Next, call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class JulyTargetAreaModifier<1>;
template class JulyTargetAreaModifier<2>;
template class JulyTargetAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(JulyTargetAreaModifier)