#include "FionaTargetAreaModifier.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Debug.hpp"
#include "LecLabel.hpp"
#include "EarlyDeathLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FionaTargetAreaModifier<DIM>::FionaTargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>(),
      mGrowthDuration(DOUBLE_UNSET),
      mApoptosisDuration(20.0), //20
      mLecArea(75.0),      //  This is correct for actual area of LECs at start
      mLecShrinkingScale(1.0), //default 1
      mScaleEarlyCellDeath(1.0)
{
}

template<unsigned DIM>
FionaTargetAreaModifier<DIM>::~FionaTargetAreaModifier()
{
}


template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)

{
    double cell_target_area;
    double lec_area = mLecArea;

    // Get target area A of a healthy cell in S, G2 or M phase
    if (pCell->HasCellProperty<LecLabel>())
    {
        //lecs are labelled, and are given larger target area
        cell_target_area = lec_area;
    }
    else
    {
        cell_target_area = 1;  //ResultsNew=1
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
            // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
            growth_duration = p_model->GetTransitCellG1Duration(); 
        }
    }

        


 
  if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
      
        double apoptosis_duration = mApoptosisDuration;

        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();

        if (pCell->HasCellProperty<EarlyDeathLabel>())
        {
            cell_target_area *= 1.0 - (mScaleEarlyCellDeath*(mLecShrinkingScale)*(time_spent_apoptotic/apoptosis_duration));
        }
        else
        {
            cell_target_area *= 1.0 - ((mLecShrinkingScale)*(time_spent_apoptotic/apoptosis_duration));
        
        }


        
       
        
        
        //area is a positive quantity 
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    
    }
    else
    {

        if (pCell->HasCellProperty<LecLabel>())
        {
           
            
            cell_target_area = mLecArea;
            
    
        }
        else
        {
            double cell_age = pCell->GetAge();

        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
        if (cell_age < growth_duration)
        {
            cell_target_area = 0.5*(1 + cell_age/growth_duration);
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
                cell_target_area = 0.5;
            }
            else
            {
                cell_target_area =1;
            }
        }
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
double FionaTargetAreaModifier<DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::SetGrowthDuration(double growthDuration)
{
    assert(growthDuration >= 0.0);
    mGrowthDuration = growthDuration;
}

template<unsigned DIM>
double FionaTargetAreaModifier<DIM>::GetLecShrinkingScale()
{
    return mLecShrinkingScale;
}

template<unsigned DIM>
double FionaTargetAreaModifier<DIM>::GetScaleEarlyCellDeath()
{
    return mScaleEarlyCellDeath;
}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::SetLecShrinkingScale(double LecShrinkingScale)
{
    mLecShrinkingScale = LecShrinkingScale;
}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::SetScaleEarlyCellDeath(double ScaleEarlyCellDeath)
{
    mScaleEarlyCellDeath = ScaleEarlyCellDeath;
}

template<unsigned DIM>
double FionaTargetAreaModifier<DIM>::GetApoptosisDuration()
{
    return mApoptosisDuration;
}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::SetApoptosisDuration(double apoptosisDuration)
{
    assert(apoptosisDuration >= 0.0);
    mApoptosisDuration = apoptosisDuration;
}

template<unsigned DIM>
double FionaTargetAreaModifier<DIM>::GetLecArea()
{
    return mLecArea;
}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::SetLecArea(double lecArea)
{
    assert(lecArea >= 0.0);
    mLecArea = lecArea;
}
template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    *rParamsFile << "\t\t\t<ApoptosisDuration>" << mApoptosisDuration << "</ApoptosisDuration>\n";
    *rParamsFile << "\t\t\t<LecArea>" << mLecArea << "</LecArea>\n";
    *rParamsFile << "\t\t\t<LecShrinkingScale>" << mLecShrinkingScale << "</LecShrinkingScale>\n";
    *rParamsFile << "\t\t\t<ScaleEarlyCellDeath>" << mScaleEarlyCellDeath << "</ScaleEarlyCellDeath>\n";
    


    // Next, call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FionaTargetAreaModifier<1>;
template class FionaTargetAreaModifier<2>;
template class FionaTargetAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FionaTargetAreaModifier)