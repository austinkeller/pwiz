//
// PrecursorMaskCodecTest.cpp
//
//
// Original author: Austin Keller <atkeller .@. uw.edu>
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// you may not use this file except in compliance with the License. 
// You may obtain a copy of the License at 
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software 
// distributed under the License is distributed on an "AS IS" BASIS, 
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and 
// limitations under the License.
//

#include "pwiz/analysis/demux/PrecursorMaskCodec.hpp"
#include "pwiz/data/msdata/examples.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/unit.hpp"

using namespace pwiz::util;
using namespace pwiz::analysis;
using namespace pwiz::msdata;

class PrecursorMaskCodecTest {

public:

    void Run()
    {
        SetUp();
        ExceptionTest();
        TearDown();
    }

protected:

    virtual void SetUp()
    {
    }

    void TearDown()
    {
    }

    void PrecursorMaskCodecDummyInitialize(SpectrumList_const_ptr slPtr, bool variableFill = false)
    {
        PrecursorMaskCodec pmc(slPtr, variableFill);
    }

    void ExceptionTest()
    {
        // Remember which spectra correspond to what states
        const int MS1_INDEX_0 = 0;
        const int MS2_INDEX_0 = 1;
        const int MS1_INDEX_1 = 2;
        const int MS2_INDEX_1 = 3;
        const int MS1_INDEX_2 = 4;

        // Generate test data
        MSDataPtr msd = boost::make_shared<MSData>();
        examples::initializeTiny(*msd);
        auto spectrumListPtr = msd->run.spectrumListPtr;

        // This should fail to read through the example dataset because there are too few spectra to interpret an acquisition scheme
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, false), runtime_error,
            "Could not determine demultiplexing scheme. Too few spectra to determine the number of precursor windows.");
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, true), runtime_error,
            "Could not determine demultiplexing scheme. Too few spectra to determine the number of precursor windows.");

        /* The number of precursors for each spectrum cannot change between spectra and cannot be empty (Mathematically, this might be
        * allowable for demultiplexing but the code makes some simplifying assumptions that rely on these two being true)
        */

        // Make the number of precursors for the first MS2 spectrum equal to zero
        msd = boost::make_shared<MSData>();
        examples::initializeTiny(*msd);
        spectrumListPtr = msd->run.spectrumListPtr;
        spectrumListPtr->spectrum(MS2_INDEX_0)->precursors.clear();

        // This should fail to read through the example dataset because the first MS2 spectrum is empty
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, false), runtime_error,
            "MS2 spectrum is missing precursor information.");
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, true), runtime_error,
            "MS2 spectrum is missing precursor information.");

        // Remove all MS2 spectra
        msd = boost::make_shared<MSData>();
        examples::initializeTiny(*msd);
        spectrumListPtr = msd->run.spectrumListPtr;
        shared_ptr<SpectrumListSimple> spectrumListSimple(new SpectrumListSimple);
        spectrumListSimple->dp = boost::make_shared<DataProcessing>(*spectrumListPtr->dataProcessingPtr());
        spectrumListSimple->spectra.push_back(spectrumListPtr->spectrum(MS1_INDEX_0));
        spectrumListSimple->spectra.push_back(spectrumListPtr->spectrum(MS1_INDEX_1));
        spectrumListSimple->spectra.push_back(spectrumListPtr->spectrum(MS1_INDEX_2));
        msd->run.spectrumListPtr = spectrumListSimple;
        spectrumListPtr = msd->run.spectrumListPtr;

        // This should fail because there are no MS2 spectra in the list
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, false), runtime_error,
            "No MS2 scans found for this experiment.");
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, true), runtime_error,
            "No MS2 scans found for this experiment.");

        // Make the number of precursors > 1 but allow the number of precursors to vary
        msd = boost::make_shared<MSData>();
        examples::initializeTiny(*msd);
        spectrumListPtr = msd->run.spectrumListPtr;
        Spectrum& s20 = *spectrumListPtr->spectrum(MS2_INDEX_1);
        s20.precursors.resize(2);
        Precursor& precursor = s20.precursors.back();
        precursor.spectrumID = s20.id;
        precursor.isolationWindow.set(MS_isolation_window_target_m_z, 455.3, MS_m_z);
        precursor.isolationWindow.set(MS_isolation_window_lower_offset, .5, MS_m_z);
        precursor.isolationWindow.set(MS_isolation_window_upper_offset, .5, MS_m_z);
        precursor.selectedIons.resize(1);
        precursor.selectedIons[0].set(MS_selected_ion_m_z, 455.34, MS_m_z);
        precursor.selectedIons[0].set(MS_peak_intensity, 120505, MS_number_of_detector_counts);
        precursor.selectedIons[0].set(MS_charge_state, 2);
        precursor.activation.set(MS_collision_induced_dissociation);
        precursor.activation.set(MS_collision_energy, 35.00, UO_electronvolt);

        // This should fail to read through the example dataset because the number of precursors changes between MS2 spectra
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, false), runtime_error,
            "Precursor sizes are varying between individual MS2 scans. Cannot infer demultiplexing scheme.");
        unit_assert_throws_what(PrecursorMaskCodecDummyInitialize(spectrumListPtr, true), runtime_error,
            "Precursor sizes are varying between individual MS2 scans. Cannot infer demultiplexing scheme.");
    }

    // Test non-multiplexed spectra. This should simply return the input spectra.
    void BasicDIATest()
    {
        //TODO
    }

    void OverlapTest()
    {
        //TODO
    }

    void MSXTest()
    {
        //TODO
    }

    void GapsBetweenPrecursorsTest()
    {
        //TODO
    }

    void VariableFillTest()
    {
        //TODO
    }

    void VaryingPrecursorWidthTest()
    {
        //TODO
    }

    void OverlapMSXTest()
    {
        //TODO
    }

    void WatersSONARTest()
    {
        //TODO
    }

};

int main(int argc, char* argv[])
{
    TEST_PROLOG(argc, argv)

        try
    {
        PrecursorMaskCodecTest tester;
        tester.Run();
    }
    catch (exception& e)
    {
        TEST_FAILED(e.what())
    }
    catch (...)
    {
        TEST_FAILED("Caught unknown exception.")
    }

    TEST_EPILOG
}