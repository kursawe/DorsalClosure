/*

Copyright (c) 2005-2022, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "FionaTargetAreaModifier.hpp"
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
FionaTargetAreaModifier<DIM>::FionaTargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>(),
      mGrowthDuration(DOUBLE_UNSET),
      mApoptosisDuration(20.0), //20
      mLecArea(38)      //  12.5 if TestFionaElongated, 32 for TestFionaDorsalClosure, 37.5 for TestInitialConditionAlt
{
}

template<unsigned DIM>
FionaTargetAreaModifier<DIM>::~FionaTargetAreaModifier()
{
}

//template<unsigned DIM>
//double MyTargetAreaModifier<DIM>::GetCurrentCellArea(AbstractCellPopulation<DIM>& rCellPopulation)
//{
//    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
//    unsigned num_elements = p_cell_population->GetNumElements();
//    std::vector<double> element_areas(num_elements);

//    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
//         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
//         ++elem_iter)
//    {
//        unsigned elem_index = elem_iter->GetIndex();
//        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
//    }
//    return element_areas
//}

template<unsigned DIM>
void FionaTargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    double cell_target_area;
    double lec_area = mLecArea;

    
    
    // if (pCell->GetCellId()==273)
    // {
    //     pCell->
    // }



	 //VertexBasedCellPopulation<SPACE_DIM>* p_CellPopulation = static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);
	 //VertexElement<SPACE_DIM, SPACE_DIM>* this_element = p_CellPopulation->GetElementCorrespondingToCell(*cell_iter);
	 //unsigned number_nodes = this_element->GetNumNodes();
	 //for (unsigned node_index=0; node_index < number_nodes; node_index++)
	 //{
	 //	PRINT_VECTOR(this_element->GetNode(node_index)->rGetLocation());

//	 }

    // Get target area A of a healthy cell in S, G2 or M phase
    if (pCell->HasCellProperty<CellLabel>())
    {
        //lecs are labelled, and are given larger target area
        
        cell_target_area = lec_area;
       
    }
    else
    {
        //hbs given target area = 1
        cell_target_area = 1*this->mReferenceTargetArea;  //ResultsNew=1
        //cell_target_area *= 1.0;
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


    if (pCell->HasCellProperty<CellLabel>())
    {

        //cell_target_area *= 0.9;

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
        cell_target_area *= 1.0 - (1*(time_spent_apoptotic/apoptosis_duration));
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