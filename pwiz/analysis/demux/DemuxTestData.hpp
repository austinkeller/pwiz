//
// DemuxTestData.hpp
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

#ifndef _DEMUXTESTDATA_HPP
#define _DEMUXTESTDATA_HPP

#include "pwiz/utility/misc/Std.hpp"
#include <boost/smart_ptr/shared_ptr.hpp>


namespace pwiz
{
namespace analysis
{
namespace test
{

    struct SimulatedDemuxParams
    {
        size_t numPrecursors = 3;
        size_t numOverlaps = 1;
        size_t numCycles = 10;
        size_t numDemux = 9;
        bool overlapOnly = false;
        int randomDataSeed = 0;
    };

    class IScanEvent
    {
    public:
        IScanEvent();
        virtual ~IScanEvent();
        virtual int mslevel() const = 0;
    };

    class MS1Scan : public IScanEvent
    {
    public:
        MS1Scan();
        virtual ~MS1Scan();

        int mslevel() const final
        {
            return 1;
        }
    };

    class MS2Scan : public IScanEvent
    {
    public:
        MS2Scan();
        virtual ~MS2Scan();

        int mslevel() const final
        {
            return 2;
        }

        size_t numPrecursors() const;
        msdata::Precursor getPrecursor(size_t index) const;
        void setPrecursors(const vector<pair<double, double>>& mzCentersAndWidths);

    private:
        vector<msdata::Precursor> _precursors;
    };

    typedef vector<IScanEvent> AcquisitionScheme;

    class IAcquisitionScheme
    {
    public:
        virtual ~IAcquisitionScheme() = 0 {}

        virtual size_t numScans() const = 0;
        virtual size_t numPrecursors() const = 0;
        virtual boost::shared_ptr<IScanEvent> scan(size_t scanNum) const = 0;
    };

    class SimpleAcquisitionScheme : public IAcquisitionScheme
    {
    public:
        SimpleAcquisitionScheme(
            size_t ms2ScansPerCycle = 9,
            double startPrecursorMz = 500.0,
            double endPrecursorMz = 900.0,
            double startProductMz = 400.0,
            double endProductMz = 1200.0,
            int randomSeed = 0);
        size_t numScans() const override;
        size_t numPrecursors() const override;
        boost::shared_ptr<IScanEvent> scan(size_t scanNum) const override;
    private:
        std::vector<IScanEvent> _scans;
    };

    class IAnalyte
    {
    public:
        virtual ~IAnalyte() = 0 {}
        virtual double precursorMz() const = 0;
        virtual boost::shared_ptr<vector<double>> fragmentMzs() const = 0;
        virtual boost::shared_ptr<vector<double>> fragmentRelIntensities() const = 0;
    };

    class SimpleAnalyte : public IAnalyte
    {
    public:
        SimpleAnalyte(
            int randomSeed = 0,
            double startPrecursorMz = 400.0,
            double endPrecursorMz = 900.0,
            double startFragmentMz = 200.0,
            double endFragmentMz = 1200.0);
        virtual ~SimpleAnalyte() {}
        double precursorMz() const override;
        boost::shared_ptr<vector<double>> fragmentMzs() const override;
        boost::shared_ptr<vector<double>> fragmentRelIntensities() const override;
    private:
        double _precursorMz;
        boost::shared_ptr<vector<double>> _fragmentMzs;
        boost::shared_ptr<vector<double>> _fragmentRelintensities;
        int _randomSeed;
    };

    class IElutionScheme
    {
    public:
        virtual ~IElutionScheme() = 0 {}
        virtual vector<pair<size_t, double>> indexedAnalyteIntensity(double time) const = 0;
        virtual boost::shared_ptr<IAnalyte> analyte(size_t index) const = 0;
    };

    class RegularSineElutionScheme : public IElutionScheme
    {
    public:
        RegularSineElutionScheme();
        virtual ~RegularSineElutionScheme() {};
        vector<pair<size_t, double>> indexedAnalyteIntensity(double time) const override;
        boost::shared_ptr<IAnalyte> analyte(size_t index) const override;

        void setSigma(double sigmaInTime);
        void setPeriod(double timeBetweenPeaks);
        void setSinePeriod(double timeBetweenLongOscillations);

    private:
        double intensity(double time) const;
        /// Index for the eluting peak. This can be used to determine which analyte is eluting.
        size_t peakIndex(double time) const;

        mutable std::map<size_t, boost::shared_ptr<IAnalyte>> _analyteCache;
        double _sigmaInTime;
        double _timeBetweenPeaks;
        double _timeBetweenLongOscillations;
    };

    class SimulatedSpectrum
    {
    public:
        SimulatedSpectrum();
        int msLevel() const;
        double precursorMz() const;
        boost::shared_ptr<vector<double>> mzs() const;
        boost::shared_ptr<vector<double>> intensities() const;
        boost::shared_ptr<IScanEvent> scan() const;
    private:
        boost::shared_ptr<IScanEvent> _scan;
    };

    class SimulatedMassSpec
    {
    public:
        SimulatedMassSpec();
        boost::shared_ptr<SimulatedSpectrum> nextScan();
        void initialize(const IAcquisitionScheme& acquisitionScheme, const IElutionScheme& elutionScheme);

    private:
        boost::shared_ptr<IElutionScheme> _elutionScheme;
        boost::shared_ptr<IAcquisitionScheme> _acquisitionScheme;
        double _minRunDuration;
        double _scanRate;
        size_t _currentScanNum;
    };

    void initializeMSDDemux(msdata::MSData& msd, SimulatedDemuxParams params = SimulatedDemuxParams());

} // namespace test
} // namespace analysis
} // namespace pwiz
#endif // _DEMUXTESTDATA_HPP