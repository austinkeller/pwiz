//
// SpectrumList_Sonar.cpp
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

#include "SpectrumList_Sonar.hpp"
#include "pwiz/analysis/demux/DemuxHelpers.hpp"
#include <boost/format.hpp>
//#include "pwiz/utility/bindings/CLI/common/ParamTypes.hpp"

namespace pwiz
{
namespace analysis
{
    using namespace std;

    class SpectrumList_Sonar::Impl
    {
    public:

        Impl(const msdata::SpectrumListPtr& inner, const Params& p, msdata::DataProcessingPtr dp);
        msdata::SpectrumPtr spectrum(size_t index, bool getBinaryData = false) const;
        msdata::SpectrumPtr spectrum(size_t index, msdata::DetailLevel detailLevel) const;
        size_t size() const;
        const msdata::SpectrumIdentity& spectrumIdentity(size_t index) const;

    private:
        struct SpectrumModifiers
        {
            int msLevel = 1;
            bool ignore = false;
        };

        SpectrumModifiers GetSpectrumModifiers(const msdata::Spectrum& spectrum) const;
        void ReadScheme();
        void CalculatePrecursorRange();
        void MapPrecursors();

        msdata::SpectrumListPtr sl_;
        size_t cycleSize_;

        /// Cache of precursors used in the cycle
        std::vector<msdata::Precursor> precursors_;

        /// Maps from spectrum index to precursor index in precursors cache
        std::map<size_t, size_t> precursorIndexMap_;
        double sonarMzRangeLower_;
        double sonarMzRangeUpper_;
        double sonarWindowSize_;
        std::vector<SpectrumModifiers> spectrumModifiers_;
    };

    SpectrumList_Sonar::Impl::Impl(const msdata::SpectrumListPtr& inner, const Params& p, msdata::DataProcessingPtr dp) :
        cycleSize_(0),
        sonarWindowSize_(p.windowSize),
        sonarMzRangeLower_(p.startMass),
        sonarMzRangeUpper_(p.endMass)
    {
        sl_ = inner;

        // TODO Sanity check that this data is actually Waters data

        // Sanity check input parameters
        if (!(sonarMzRangeLower_ >= 0.0))
            throw runtime_error((boost::format("Sonar range lower m/z must be non-negative. Value: %1%")
            % sonarMzRangeLower_).str());
        if (!(sonarMzRangeUpper_ >= 0.0))
            throw runtime_error((boost::format("Sonar range upper m/z must be non-negative. Value: %1%")
            % sonarMzRangeUpper_).str());
        if (!(sonarMzRangeUpper_ > sonarMzRangeLower_))
            throw runtime_error((boost::format("Sonar upper m/z must be greater than the lower m/z."
            " Lower: %1%, Upper: %2%") % sonarMzRangeLower_ % sonarMzRangeUpper_).str());
        if (!(sonarWindowSize_ > 0.0))
            throw runtime_error((boost::format("Sonar window size must be positive. Value: %1%")
            % sonarWindowSize_).str());

        ReadScheme();
        CalculatePrecursorRange();
        MapPrecursors();
    }

    PWIZ_API_DECL size_t SpectrumList_Sonar::Impl::size() const
    {
        return sl_->size();
    }

    PWIZ_API_DECL const msdata::SpectrumIdentity& SpectrumList_Sonar::Impl::spectrumIdentity(size_t index) const
    {
        return sl_->spectrumIdentity(index);
    }

    PWIZ_API_DECL msdata::SpectrumPtr SpectrumList_Sonar::Impl::spectrum(size_t index, msdata::DetailLevel detailLevel) const
    {
        // TODO: add ability to deal with non-binary-data requests
        return spectrum(index, true);
    }

    inline void replaceCVParam(vector<data::CVParam>& cvParams, cv::CVID oldCVID, cv::CVID newCVID)
    {
        vector<data::CVParam>::iterator it = find_if(cvParams.begin(), cvParams.end(), data::CVParamIs(oldCVID));
        if (it != cvParams.end())
        {
            it->cvid = newCVID;
        }
        else
        {
            cvParams.push_back(data::CVParam(newCVID));
        }
    }

    PWIZ_API_DECL msdata::SpectrumPtr SpectrumList_Sonar::Impl::spectrum(size_t index, bool getBinaryData) const
    {
        auto spectrum = sl_->spectrum(index, getBinaryData);
        auto mods = GetSpectrumModifiers(*spectrum);
        if (!mods.ignore)
        {
            for (int i = 0; i < spectrum->cvParams.size(); ++i)
            {
                auto& cvParam = spectrum->cvParams[i];
                if (cvParam.cvid == cv::MS_ms_level)
                    cvParam.value = to_string(mods.msLevel);
                if (mods.msLevel == 2 && cvParam.cvid == cv::MS_MS1_spectrum)
                    cvParam.cvid = cv::MS_MSn_spectrum;
            }

            // Label this as a profile spectrum
            replaceCVParam(spectrum->cvParams, cv::MS_centroid_spectrum, cv::MS_profile_spectrum);

            if (mods.msLevel == 2)
                spectrum->precursors.push_back(precursors_.at(precursorIndexMap_.at(index)));
        }

        return spectrum;
    }

    SpectrumList_Sonar::Impl::SpectrumModifiers SpectrumList_Sonar::Impl::GetSpectrumModifiers(const msdata::Spectrum& spectrum) const
    {
        SpectrumModifiers returnSpectrumModifiers = SpectrumModifiers();
        string function;
        if (!TryGetScanIDToken(spectrum, "function", function))
            throw runtime_error("Expected spectrum identity to be formatted as a Waters identity");
        auto functionValue = stoi(function);
        int msLevel;
        if (!TryGetMSLevel(spectrum, msLevel))
            throw runtime_error("Could not retrive mse level while interpreting Waters data");

        // Uses logic from https://svn.code.sf.net/p/proteowizard/code/trunk/pwiz/pwiz_tools/Skyline/Model/Results/SpectrumFilter.cs
        // in function UpdateMseLevel()
        switch (functionValue)
        {
        case 1:
            for (auto& param : spectrum.cvParams)
            {
                if (param.cvid == cv::MS_ms_level)
                {
                    //assert(stoi(param.value) == 1); // sanity check, we expect this to already be set to 1
                    returnSpectrumModifiers.msLevel = 1;
                }
            }
            break;
        case 2:
            for (auto& param : spectrum.cvParams)
            {
                if (param.cvid == cv::MS_ms_level)
                {
                    //assert(stoi(param.value) == 1); // sanity check, we expect this to already be set to 1
                    returnSpectrumModifiers.msLevel = 2;
                }
            }
            break;
        case 3:
        {
            returnSpectrumModifiers.ignore = true;
            break;
        }
        default:
            throw runtime_error("Unknown function value in Waters spectrum identity");
        }

        return returnSpectrumModifiers;
    }

    void SpectrumList_Sonar::Impl::ReadScheme()
    {
        // Read through a set of spectra and infer the cycle size from all spectra and infer the number of bins and window size
        int currentMSType = 0;
        int lastMSType = 0;
        size_t lastCycleSize = 0;
        size_t currentCycleSize = 0;
        int fullCyclesObserved = 0;
        for (size_t i = 0; i < sl_->size(); ++i)
        {
            auto currentSpectrum = sl_->spectrum(i);

            // cache the modifications to this spectrum
            auto mods = GetSpectrumModifiers(*currentSpectrum);

            if (mods.ignore)
            {
                currentMSType = 3;
            }
            else if (mods.msLevel == 1)
            {
                currentMSType = 1;
            }
            else if (mods.msLevel == 2)
            {
                currentMSType = 2;
            }
            else
            {
                throw runtime_error("Unknown acquisition method");
            }

            if (lastMSType != currentMSType)
            {
                // Changed spectrum type so the cycle just ended
                if (lastCycleSize != currentCycleSize)
                {
                    // new cycle size was observed, restart count
                    lastCycleSize = currentCycleSize;
                    fullCyclesObserved = 1;
                }
                else
                {
                    ++fullCyclesObserved;
                }

                // restart cycle count
                currentCycleSize = 0;

                // Update MSType of this new cycle
                lastMSType = currentMSType;
            }

            // increment spectrum count for this cycle
            ++currentCycleSize;

            if (fullCyclesObserved >= 4)
            {
                cycleSize_ = lastCycleSize;
            }
        }
    }

    void SpectrumList_Sonar::Impl::CalculatePrecursorRange()
    {
        assert(cycleSize_ > 0);
        assert(sonarMzRangeLower_ > 0.0);
        assert(sonarMzRangeUpper_ > 0.0);
        assert(sonarWindowSize_ > 0.0);

        // Divide the mz range into the number of bins observed in the cycle
        double sonarMzRange = sonarMzRangeUpper_ - sonarMzRangeLower_;
        double sonarOffset = sonarMzRange / (cycleSize_ - 1);
        for (size_t i = 0; i < cycleSize_; ++i)
        {
            double center = sonarOffset * i + sonarMzRangeLower_;
            //double lower = center - sonarWindowSize_ / 2.0;
            //double upper = center + sonarWindowSize_ / 2.0;
            auto p = msdata::Precursor(center);
            p.isolationWindow.set(cv::MS_isolation_window_target_m_z, to_string(center), cv::MS_m_z);
            p.isolationWindow.set(cv::MS_isolation_window_lower_offset, to_string(sonarWindowSize_ / 2.0), cv::MS_m_z);
            p.isolationWindow.set(cv::MS_isolation_window_upper_offset, to_string(sonarWindowSize_ / 2.0), cv::MS_m_z);
            precursors_.push_back(p);
        }
    }

    void SpectrumList_Sonar::Impl::MapPrecursors()
    {
        // Read through all spectra and map the indices of ms2 spectra to their corresponding precursors.
        // This provides quick lookup when demultiplexed spectra are requested with spectrum()
        bool firstMS1Found = false;
        bool firstMS2Found = false;
        size_t cycleSize = 0;
        spectrumModifiers_.clear();
        for (size_t i = 0; i < sl_->size(); ++i)
        {
            auto currentSpectrum = sl_->spectrum(i);
            // cache the modifications to this spectrum
            spectrumModifiers_.push_back(GetSpectrumModifiers(*currentSpectrum));
            if (spectrumModifiers_.back().ignore)
                continue; // don't count lockmass scan as part of the cycle
            auto msLevel = spectrumModifiers_.back().msLevel;
            if (msLevel == 1)
            {
                if (firstMS1Found && firstMS2Found)
                {
                    // Reached the end of the cycle
                    assert(cycleSize == cycleSize_);
                    cycleSize = 0;
                    firstMS2Found = false;
                }
                // found the first spectrum in a cycle
                firstMS1Found = true;
            }
            else if (msLevel == 2 && firstMS1Found)
            {
                firstMS2Found = true;
                precursorIndexMap_.insert(pair<size_t, size_t>(i, cycleSize));
                ++cycleSize;
            }
        }
    }
    SpectrumList_Sonar::SpectrumList_Sonar(const msdata::SpectrumListPtr& inner, const Params& p) : SpectrumListWrapper(inner), impl_(new Impl(inner, p, dp_)) {}
    SpectrumList_Sonar::~SpectrumList_Sonar() {}
    msdata::SpectrumPtr SpectrumList_Sonar::spectrum(size_t index, bool getBinaryData) const { return impl_->spectrum(index, getBinaryData); }
    msdata::SpectrumPtr SpectrumList_Sonar::spectrum(size_t index, msdata::DetailLevel detailLevel) const { return impl_->spectrum(index, detailLevel); }
    size_t SpectrumList_Sonar::size() const { return impl_->size(); }
    const msdata::SpectrumIdentity& SpectrumList_Sonar::spectrumIdentity(size_t index) const { return impl_->spectrumIdentity(index); }
} // namespace analysis
} // namespace pwiz
