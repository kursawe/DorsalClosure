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

#ifndef FIONATARGETAREAMODIFIER_HPP_
#define FIONATARGETAREAMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractTargetAreaModifier.hpp"

/**
 * A target area modifier class in which the target area of a cell grows linearly, up to
 * mReferenceTargetArea, over a prescribed duration.
 *
 * If used with a phase-based cell-cycle model (such as FixedG1GenerationalCellCycleModel),
 * the target area of a cell increases linearly from the value 0.5*mReferenceTargetArea
 * up to mReferenceTargetArea over the course of the cell's G1 phase.
 *
 * If used with a non-phase-based cell-cycle model, the target area of a cell increases
 * linearly from the value 0.5*mReferenceTargetArea up to mReferenceTargetArea while the
 * cell's age is less than mGrowthDuration.
 *
 * Here mReferenceTargetArea, mGrowthDuration and mApoptosisTime are settable member variables. The default
 * value of mReferenceTargetArea is 1.0, the default value of mGrowthDuration is DOUBLE_UNSET and the default 
 * value of mApoptosisDuration is 20.0.
 *
 * Note that if mGrowthDuration is set by the user, then this value is are used to prescribe
 * target area growth as described earlier, regardless of whether a phase-based cell-cycle
 * model is present.
 * 
 * mApoptosisDuration defines how quickly the cell constricts i.e. how quickly it's target area approaches 0
 */
template<unsigned DIM>
class FionaTargetAreaModifier : public AbstractTargetAreaModifier<DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTargetAreaModifier<DIM> >(*this);
        archive & mGrowthDuration & mApoptosisDuration & mLecArea & mLecShrinkingScale & mScaleEarlyCellDeath;
    }

    /**
     * The duration that a cell's target area to increase from 0.5*mReferenceTargetArea
     * to mReferenceTargetArea at the start of its cell cycle. Defaults to DOUBLE_UNSET.
     * If this variable is set using SetGrowthDuration(), then it is used regardless of
     * whether a phase-based cell-cycle model is used.
     */
    double mGrowthDuration;

    double mApoptosisDuration;

    double mLecArea;

    double mLecShrinkingScale;

    double mScaleEarlyCellDeath;

public:

    /**
     * Default constructor.
     */
    FionaTargetAreaModifier();

    /**
     * Destructor.
     */
    virtual ~FionaTargetAreaModifier();

    /**
     * 
     * @param rCellPopulation reference to the cell population
    */
  // virtual void GetCurrentCellArea(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden UpdateTargetAreaOfCell() method.
     *
     * @param pCell pointer to the cell
     */
    virtual void UpdateTargetAreaOfCell(const CellPtr pCell);

    /**
     * @return #mGrowthDuration
     */
    double GetGrowthDuration();

    /**
     * Set #mGrowthDuration.
     *
     * @param growthDuration the new value of #mGrowthDuration
     */
    void SetGrowthDuration(double growthDuration);


    /**
     * @return #mLecShrinkingScale
     */
    double GetLecShrinkingScale();

     /**
     * @return #mScaleEarlyCellDeath
     */
    double GetScaleEarlyCellDeath();



    /**
     * Set #mLecShrinkingScale.
     *
     * @param LecShrinkingScale the new value of #mGrowthDuration
     */
    void SetLecShrinkingScale(double LecShrinkingScale);

        /**
     * Set #mScaleEarlyCellDeath.
     *
     * @param ScaleEarlyCellDeath the new value of #mGrowthDuration
     */
    void SetScaleEarlyCellDeath(double ScaleEarlyCellDeath);

    /**
     * @return #mApoptosisDuration
     */
    double GetApoptosisDuration();

    /**
     * Set #mApoptosisDuration.
     *
     * @param apoptosisDuration the new value of #mApoptosisDuration
     */
    void SetApoptosisDuration(double apoptosisDuration);
    
    /**
     * @return #mLecArea
     */
    double GetLecArea();

    /**
     * Set #mLecArea.
     *
     * @param lecArea the new value of #mLecArea
     */
    void SetLecArea(double lecArea);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FionaTargetAreaModifier)

#endif /*FIONATARGETAREAMODIFIER_HPP_*/