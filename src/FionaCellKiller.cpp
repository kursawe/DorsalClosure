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
//the only change made (from randomcellkiller) is that only cells given a celllabel are labelled for apoptosis
#include "FionaCellKiller.hpp"
#include "LecLabel.hpp"
#include "EarlyDeathLabel.hpp"
#include "CellPropertyRegistry.hpp"
#include "Debug.hpp"
#include "ApoptoticCellProperty.hpp"
#include "VertexBasedCellPopulation.hpp"


template<unsigned DIM>
FionaCellKiller<DIM>::FionaCellKiller(AbstractCellPopulation<DIM>* pCellPopulation, double ProbabilityShrinking, double ProbabilityEarlyCellDeath, double ProbabilityBiasLecDeath)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mProbabilityShrinking(ProbabilityShrinking),
          mProbabilityEarlyCellDeath(ProbabilityEarlyCellDeath),
          mProbabilityBiasLecDeath(ProbabilityBiasLecDeath)
{
    if ((mProbabilityShrinking<0) || (mProbabilityShrinking>1))
    {
        EXCEPTION("Probability of shrinking must be between zero and one");
    }

    if ((mProbabilityEarlyCellDeath<0) || (mProbabilityEarlyCellDeath>1))
    {
        EXCEPTION("Probability of early cell death must be between zero and one");
    }
}

template<unsigned DIM>
double FionaCellKiller<DIM>::GetShrinkingProbability() const
{
    return mProbabilityShrinking;
}

template<unsigned DIM>
double FionaCellKiller<DIM>::GetEarlyCellDeathProbability() const
{
    return mProbabilityEarlyCellDeath;
}

template<unsigned DIM>
double FionaCellKiller<DIM>::GetBiasLecDeathProbability() const
{
    return mProbabilityBiasLecDeath;
}

template<unsigned DIM>
void FionaCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
    /*
     * We assume a constant time step and that there are an integer number (n = 1/dt)
     * of time steps per hour. We also assume that this method is called every time step
     * and that the probabilities of dying at different times are independent.
     *
     * Let q=mProbabilityOfDeathInAnHour and p="probability of death in a given time step".
     *
     * Probability of not dying in an hour:
     * (1-q) = (1-p)^n = (1-p)^(1/dt).
     *
     * Rearranging for p:
     * p = 1 - (1-q)^dt.
     */
    double shrinking_this_timestep = 1.0 - pow((1.0 - mProbabilityShrinking), SimulationTime::Instance()->GetTimeStep());
    double earlycelldeath_this_timestep = 1.0 - pow((1.0 - mProbabilityEarlyCellDeath), SimulationTime::Instance()->GetTimeStep());
    
    if (!pCell->HasApoptosisBegun() &&
        RandomNumberGenerator::Instance()->ranf() < shrinking_this_timestep)
    {
        // Mark the cell as apoptotic and store removal information if required.
        this->mpCellPopulation->StartApoptosisOnCell(pCell, "FionaCellKiller");
    }

    MutableVertexMesh<DIM, DIM>& vertex_mesh = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation)->rGetMesh();

    // Get the element index corresponding to this cell
    unsigned elem_index = this->mpCellPopulation->GetLocationIndexUsingCell(pCell);

    // Get the set of neighbouring element indices
    std::set<unsigned> neighbouring_elem_indices = vertex_mesh.GetNeighbouringElementIndices(elem_index);

    // Check if any of the corresponding cells have the CellLabel property
    unsigned num_hb_neighbours = 0;
    for (std::set<unsigned>::iterator elem_iter = neighbouring_elem_indices.begin();
    elem_iter != neighbouring_elem_indices.end();
    ++elem_iter)
    {
        if (this->mpCellPopulation->GetCellUsingLocationIndex(*elem_iter)->template HasCellProperty<LecLabel>())
        {}
        else
        {
            num_hb_neighbours++;
        }                
    }

    // ...and if none do, then kill this cell
    if (num_hb_neighbours > 0)
    {
        // Mark the cell as killed and store removal information if required.
        if (!pCell->HasApoptosisBegun() && RandomNumberGenerator::Instance()->ranf() < (mProbabilityBiasLecDeath)*earlycelldeath_this_timestep)
        {
            // Mark the cell as apoptotic and store removal information if required.
            this->mpCellPopulation->StartApoptosisOnCell(pCell, "FionaCellKiller");
            boost::shared_ptr<AbstractCellProperty> ed_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EarlyDeathLabel>();
            pCell->AddCellProperty(ed_label);
        }
        else if (pCell->HasApoptosisBegun() && RandomNumberGenerator::Instance()->ranf() < (mProbabilityBiasLecDeath)*earlycelldeath_this_timestep)
        {
            boost::shared_ptr<AbstractCellProperty> ed_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EarlyDeathLabel>();
            pCell->AddCellProperty(ed_label);
        }
    }
    else
    {
        // Mark the cell as killed and store removal information if required.
        if (!pCell->HasApoptosisBegun() && RandomNumberGenerator::Instance()->ranf() < ((1.0-mProbabilityBiasLecDeath)*earlycelldeath_this_timestep))
        {
            // Mark the cell as apoptotic and store removal information if required.
            this->mpCellPopulation->StartApoptosisOnCell(pCell, "FionaCellKiller");
            boost::shared_ptr<AbstractCellProperty> ed_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EarlyDeathLabel>();
            pCell->AddCellProperty(ed_label);
        }
        else if (pCell->HasApoptosisBegun() && RandomNumberGenerator::Instance()->ranf() < ((1.0-mProbabilityBiasLecDeath)*earlycelldeath_this_timestep))
        {
            boost::shared_ptr<AbstractCellProperty> ed_label = pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EarlyDeathLabel>();
            pCell->AddCellProperty(ed_label);
        }
    }
    

}

template<unsigned DIM>
void FionaCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        //lecs are labelled, hbs are not. so here only the labelled lecs are able to die 
        if (cell_iter->template HasCellProperty<LecLabel>())
        {
            //MARK;
            CheckAndLabelSingleCellForApoptosis(*cell_iter);
        }
    }
}

template<unsigned DIM>
void FionaCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProbabilityShrinking>" << mProbabilityShrinking << "</ProbabilityShrinking>\n";
    *rParamsFile << "\t\t\t<ProbabilityEarlyCellDeath>" << mProbabilityShrinking << "</ProbabilityEarlyCellDeath>\n";
    *rParamsFile << "\t\t\t<ProbabilityBiasLecDeath>" << mProbabilityBiasLecDeath << "</ProbabilityBiasLecDeath>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class FionaCellKiller<1>;
template class FionaCellKiller<2>;
template class FionaCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FionaCellKiller)