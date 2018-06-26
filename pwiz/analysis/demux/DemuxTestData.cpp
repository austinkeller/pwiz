//
// DemuxTestData.cpp
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

#include "pwiz/analysis/demux/DemuxTestData.hpp"
#include <pwiz/data/msdata/MSData.hpp>
#include <pwiz/utility/misc/Export.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <boost/smart_ptr/make_shared.hpp>

namespace pwiz
{
namespace analysis
{
namespace test
{
    using namespace pwiz::msdata;

    void generateAcquisitionScheme(AcquisitionScheme& acquisitionScheme, SimulatedDemuxParams params = SimulatedDemuxParams())
    {

    }


    void initializeMSDDemuxMeta(msdata::MSData& msd)
    {
        msd.id = "urn:lsid:psidev.info:mzML.instanceDocuments.simdemux.pwiz";

        // cvList

        msd.cvs = defaultCVList();

        // fileDescription

        FileContent& fc = msd.fileDescription.fileContent;
        fc.set(MS_MS1_spectrum);
        fc.set(MS_MSn_spectrum);

        SourceFilePtr sfp(new SourceFile);
        sfp->id = "DEMUX_RAW";
        sfp->name = "simulated_demux.raw";
        sfp->location = "file://D:/data/Exp01";
        sfp->set(MS_Thermo_nativeID_format);
        sfp->set(MS_Thermo_RAW_format);
        sfp->set(MS_SHA_1, "1234567890123456789012345678901234567890");
        msd.fileDescription.sourceFilePtrs.push_back(sfp);

        msd.fileDescription.contacts.resize(1);

        // paramGroupList

        ParamGroupPtr pg1(new ParamGroup);
        pg1->id = "CommonMS1SpectrumParams";
        pg1->set(MS_Exactive);
        pg1->set(MS_instrument_serial_number, "Exactive Series slot #1");
        msd.paramGroupPtrs.push_back(pg1);

        // instrumentConfigurationList

        InstrumentConfigurationPtr instrumentConfigurationPtr(new InstrumentConfiguration("IC1"));
        instrumentConfigurationPtr->set(MS_Exactive);
        instrumentConfigurationPtr->componentList.push_back(Component(MS_nanoelectrospray, 1));
        instrumentConfigurationPtr->componentList.push_back(Component(MS_nanospray_inlet, 1));
        instrumentConfigurationPtr->componentList.push_back(Component(MS_orbitrap, 2));
        instrumentConfigurationPtr->componentList.push_back(Component(MS_inductive_detector, 3));

        SoftwarePtr softwareXcalibur(new Software);
        softwareXcalibur->id = "Xcalibur";
        softwareXcalibur->set(MS_Xcalibur);
        softwareXcalibur->version = "2.0-148000/2.0.2.1461";
        instrumentConfigurationPtr->softwarePtr = softwareXcalibur;

        msd.instrumentConfigurationPtrs.push_back(instrumentConfigurationPtr);

        // softwareList

        SoftwarePtr softwarepwiz(new Software);
        softwarepwiz->id = "pwiz";
        softwarepwiz->set(MS_pwiz);
        softwarepwiz->version = "3.0.3505";

        msd.softwarePtrs.push_back(softwarepwiz);
        msd.softwarePtrs.push_back(softwareXcalibur);

        // dataProcessingList

        // (none)

        // run

        msd.run.id = "Experiment 1";
        msd.run.defaultInstrumentConfigurationPtr = instrumentConfigurationPtr;
        msd.run.startTimeStamp = "2007-06-27T15:23:45.00035";
        msd.run.defaultSourceFilePtr = sfp;
    }


    void initializeMSDDemux(msdata::MSData& msd, SimulatedDemuxParams params = SimulatedDemuxParams())
    {
        initializeMSDDemuxMeta(msd);

        InstrumentConfigurationPtr instrumentConfigurationPtr(msd.run.defaultInstrumentConfigurationPtr);
        ParamGroupPtr pg1(msd.paramGroupPtrs.front());

        boost::shared_ptr<SpectrumListSimple> spectrumList(new SpectrumListSimple);
        msd.run.spectrumListPtr = spectrumList;

        // Create acquisition scheme
        SimpleAcquisitionScheme acquisitionScheme;

        // Create elution scheme
        RegularSineElutionScheme elutionScheme;

        SimulatedMassSpec massSpec;
        massSpec.initialize(acquisitionScheme, elutionScheme);

        // Create spectra
        size_t scanNum = 0;
        while (true)
        {
            auto spectrum = massSpec.nextScan();
            if (!spectrum)
                break;
            if (spectrum->msLevel() == 1)
            {
                // Write ms1 spectrum
                spectrumList->spectra.push_back(boost::make_shared<Spectrum>());
                Spectrum& ms1 = *spectrumList->spectra.back();
                boost::format scanfmt("scan=%1%");
                scanfmt % scanNum;
                ms1.id = scanfmt.str();
                ms1.index = scanNum;

                ms1.set(MS_ms_level, 1);

                ms1.set(MS_centroid_spectrum);
                ms1.set(MS_lowest_observed_m_z, 400.39, MS_m_z);
                ms1.set(MS_highest_observed_m_z, 1795.56, MS_m_z);
                ms1.set(MS_base_peak_m_z, 445.347, MS_m_z);
                ms1.set(MS_base_peak_intensity, 120053, MS_number_of_detector_counts);
                ms1.set(MS_total_ion_current, 1.66755e+007);

                ms1.paramGroupPtrs.push_back(pg1);
                ms1.scanList.scans.push_back(Scan());
                ms1.scanList.set(MS_no_combination);
                Scan& ms1scan = ms1.scanList.scans.back();
                ms1scan.instrumentConfigurationPtr = instrumentConfigurationPtr;
                ms1scan.set(MS_scan_start_time, 5.890500, UO_minute);
                ms1scan.set(MS_filter_string, "+ c NSI Full ms [ 400.00-1800.00]");
                ms1scan.set(MS_preset_scan_configuration, 3);
                ms1scan.scanWindows.resize(1);
                ScanWindow& window = ms1.scanList.scans.back().scanWindows.front();
                window.set(MS_scan_window_lower_limit, 400.000000, MS_m_z);
                window.set(MS_scan_window_upper_limit, 1800.000000, MS_m_z);

                BinaryDataArrayPtr ms1_mz(new BinaryDataArray);
                ms1_mz->set(MS_m_z_array, "", MS_m_z);
                ms1_mz->data.resize(15);
                for (int i = 0; i<15; i++)
                    ms1_mz->data[i] = i;

                BinaryDataArrayPtr ms1_intensity(new BinaryDataArray);
                ms1_intensity->set(MS_intensity_array, "", MS_number_of_detector_counts);
                ms1_intensity->data.resize(15);
                for (int i = 0; i<15; i++)
                    ms1_intensity->data[i] = 15 - i;

                ms1.binaryDataArrayPtrs.push_back(ms1_mz);
                ms1.binaryDataArrayPtrs.push_back(ms1_intensity);
                ms1.defaultArrayLength = ms1_mz->data.size();

                // Increment scan index
                ++scanNum;
            }
            else if (spectrum->msLevel() == 2)
            {
                auto ms2scan = boost::dynamic_pointer_cast<MS2Scan>(spectrum);
                if (!ms2scan)
                    throw logic_error("Failed downcast, unknown scan type");
                ms2scan->

                // Increment scan index
                ++scanNum;
            }
                

        }

   

        // chromatograms

        shared_ptr<ChromatogramListSimple> chromatogramList(new ChromatogramListSimple);
        msd.run.chromatogramListPtr = chromatogramList;

        chromatogramList->chromatograms.push_back(ChromatogramPtr(new Chromatogram));
        chromatogramList->chromatograms.push_back(ChromatogramPtr(new Chromatogram));

        Chromatogram& tic = *chromatogramList->chromatograms[0];
        tic.id = "tic";
        tic.index = 0;
        tic.defaultArrayLength = 15;
        tic.set(MS_total_ion_current_chromatogram);

        BinaryDataArrayPtr tic_time(new BinaryDataArray);
        tic_time->set(MS_time_array, "", UO_second);
        tic_time->data.resize(15);
        for (int i = 0; i<15; i++)
            tic_time->data[i] = i;

        BinaryDataArrayPtr tic_intensity(new BinaryDataArray);
        tic_intensity->set(MS_intensity_array, "", MS_number_of_detector_counts);
        tic_intensity->data.resize(15);
        for (int i = 0; i<15; i++)
            tic_intensity->data[i] = 15 - i;

        tic.binaryDataArrayPtrs.push_back(tic_time);
        tic.binaryDataArrayPtrs.push_back(tic_intensity);

        Chromatogram& sic = *chromatogramList->chromatograms[1];
        sic.id = "sic";
        sic.index = 1;
        sic.defaultArrayLength = 10;
        sic.set(MS_selected_ion_current_chromatogram);

        sic.precursor.isolationWindow.set(MS_isolation_window_target_m_z, 456.7, MS_m_z);
        sic.precursor.activation.set(MS_CID);
        sic.product.isolationWindow.set(MS_isolation_window_target_m_z, 678.9, MS_m_z);

        BinaryDataArrayPtr sic_time(new BinaryDataArray);
        sic_time->set(MS_time_array, "", UO_second);
        sic_time->data.resize(10);
        for (int i = 0; i<10; i++)
            sic_time->data[i] = i;

        BinaryDataArrayPtr sic_intensity(new BinaryDataArray);
        sic_intensity->set(MS_intensity_array, "", MS_number_of_detector_counts);
        sic_intensity->data.resize(10);
        for (int i = 0; i<10; i++)
            sic_intensity->data[i] = 10 - i;

        sic.binaryDataArrayPtrs.push_back(sic_time);
        sic.binaryDataArrayPtrs.push_back(sic_intensity);
    }

    double normalPDF(double x, double mu = 0.0, double sigma = 1.0)
    {
        return 1.0 / (sqrt(2.0 * M_PI) * sigma) * exp(-pow(x - mu, 2) / (2.0 * pow(sigma, 2.0)));
    }

    SimpleAcquisitionScheme::SimpleAcquisitionScheme(size_t ms2ScansPerCycle, double startPrecursorMz, double endPrecursorMz, double startProductMz, double endProductMz, int randomSeed)
    {
        _scans.push_back(MS1Scan());
        for (size_t scanNum = 0; scanNum < ms2ScansPerCycle; ++scanNum)
        {
            MS2Scan ms2scan = MS2Scan();
            vector<size_t> demuxIndices;
            Precursor precursor;

            ms2scan.setPrecursors(demuxIndices);
            _scans.push_back(ms2scan);
        }
    }

    size_t SimpleAcquisitionScheme::numScans() const
    {
        return _scans.size();
    }

    size_t SimpleAcquisitionScheme::numPrecursors() const
    {
    }

    boost::shared_ptr<IScanEvent> SimpleAcquisitionScheme::scan(size_t scanNum) const
    {
    }

    SimpleAnalyte::SimpleAnalyte(int randomSeed, double startPrecursorMz, double endPrecursorMz, double startFragmentMz, double endFragmentMz) :
        _randomSeed(randomSeed)
    {
        const size_t numFragments = 5;
        std::mt19937 gen;
        gen.seed(randomSeed);
        std::uniform_real_distribution<double> precursorDist(startPrecursorMz, endPrecursorMz);
        _precursorMz = precursorDist(gen);

        std::uniform_real_distribution<double> fragmentMzDist(startFragmentMz, endFragmentMz);
        std::uniform_real_distribution<double> fragmentRelIntensityDist(0.0, 1.0);
        for (size_t i = 0; i < numFragments; ++i)
        {
            _fragmentMzs->push_back(fragmentMzDist(gen));
            _fragmentRelintensities->push_back(fragmentRelIntensityDist(gen));
        }

        std::sort(_fragmentRelintensities->begin(), _fragmentRelintensities->end(),
            [this](size_t i, size_t j)
            { return (this->_fragmentMzs->at(i) < this->_fragmentMzs->at(j)); }
        );
        std::sort(_fragmentMzs->begin(), _fragmentMzs->end());
    }

    double SimpleAnalyte::precursorMz() const
    {
        return _precursorMz;
    }

    boost::shared_ptr<vector<double>> SimpleAnalyte::fragmentMzs() const
    {
        return _fragmentMzs;
    }

    boost::shared_ptr<vector<double>> SimpleAnalyte::fragmentRelIntensities() const
    {
        return _fragmentRelintensities;
    }

    RegularSineElutionScheme::RegularSineElutionScheme() :
        _sigmaInTime(1.0),
        _timeBetweenPeaks(5.0),
        _timeBetweenLongOscillations(50.0)
    {
    }

    vector<pair<size_t, double>> RegularSineElutionScheme::indexedAnalyteIntensity(double time) const
    {
        vector<pair<size_t, double>> indexedAnalyteIntensities;
        indexedAnalyteIntensities.push_back(pair<size_t, double>(peakIndex(time), intensity(time)));
        return std::move(indexedAnalyteIntensities);
    }

    boost::shared_ptr<IAnalyte> RegularSineElutionScheme::analyte(size_t index) const
    {
        auto it = _analyteCache.find(index);
        if (it == _analyteCache.end())
        {
            _analyteCache[index];
            _analyteCache[index] = boost::make_shared<SimpleAnalyte>(index);
            it = _analyteCache.find(index);
        }
        return it->second;
    }

    void RegularSineElutionScheme::setSigma(double sigmaInTime)
    {
        _sigmaInTime = sigmaInTime;
    }

    void RegularSineElutionScheme::setPeriod(double timeBetweenPeaks)
    {
        _timeBetweenPeaks = timeBetweenPeaks;
    }

    void RegularSineElutionScheme::setSinePeriod(double timeBetweenLongOscillations)
    {
        _timeBetweenLongOscillations = timeBetweenLongOscillations;
    }

    double RegularSineElutionScheme::intensity(double time) const
    {
        double sineMinimumOffset = 1.5;
        double sineScalingFactor = sineMinimumOffset + sin(time * 2.0 * M_PI / _timeBetweenLongOscillations);
        double mu = _timeBetweenPeaks / 2.0;
        double timeWithinGaussianPeak = fmod(time, _timeBetweenPeaks);
        double normalDistValue = normalPDF(timeWithinGaussianPeak, mu, _sigmaInTime);
        return sineScalingFactor * normalDistValue;
    }

    size_t RegularSineElutionScheme::peakIndex(double time) const
    {
        return (size_t)floorl(time / _timeBetweenPeaks);
    }

    SimulatedSpectrum::SimulatedSpectrum()
    {
    }

    int SimulatedSpectrum::msLevel() const
    {
        return _scan->mslevel();
    }

    double SimulatedSpectrum::precursorMz() const
    {
        
    }

    boost::shared_ptr<vector<double>> SimulatedSpectrum::mzs() const
    {
    }

    boost::shared_ptr<vector<double>> SimulatedSpectrum::intensities() const
    {
    }

    boost::shared_ptr<IScanEvent> SimulatedSpectrum::scan() const
    {
    }

    SimulatedMassSpec::SimulatedMassSpec() :
        _minRunDuration(100.0),
        _scanRate(20.0),
        _currentScanNum(0)
    {
    }

    boost::shared_ptr<SimulatedSpectrum> SimulatedMassSpec::nextScan()
    {
        auto spectrum = boost::shared_ptr<SimulatedSpectrum>();
        assert(!spectrum);
        if (_minRunDuration < _currentScanNum * _scanRate)
        {
            spectrum.reset(new SimulatedSpectrum());
            ++_currentScanNum;
        }
        return spectrum;
    }

    void SimulatedMassSpec::initialize(const IAcquisitionScheme& acquisitionScheme, const IElutionScheme& elutionScheme)
    {
        _acquisitionScheme = boost::make_shared<IAcquisitionScheme>(acquisitionScheme);
        _elutionScheme = boost::make_shared<IElutionScheme>(elutionScheme);
    }
} // namespace test
} // namespace analysis
} // namespace pwiz