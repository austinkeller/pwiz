//
// SpectrumList_Sonar.hpp
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

#ifndef _SPECTRUMLIST_SONAR_HPP
#define _SPECTRUMLIST_SONAR_HPP
#include "pwiz/data/msdata/SpectrumListWrapper.hpp"
#include <boost/smart_ptr/scoped_ptr.hpp>

namespace pwiz
{
namespace analysis
{
    /// A class for reinterpreting Waters data as MS2 data. This functionality would ideally be integrated with
    /// the SpectrumList_Waters class so that conversion from raw to mzML always adds the appropriate attributes.
    class PWIZ_API_DECL SpectrumList_Sonar :
        public msdata::SpectrumListWrapper
    {
    public:

        /// User-defined options
        struct Params
        {
            Params() :
            windowSize(0.0),
            startMass(0.0),
            endMass(0.0)
            {}

            /// Size of each MS2 window in m/z
            double windowSize;
            /// Start mass of isolation scan range
            double startMass;
            /// End mass of isolation scan range
            double endMass;
        };

        SpectrumList_Sonar(const msdata::SpectrumListPtr& inner, const Params& p);
        virtual ~SpectrumList_Sonar();

        /// \name SpectrumList Interface
        ///@{

        msdata::SpectrumPtr spectrum(size_t index, bool getBinaryData = false) const;
        msdata::SpectrumPtr spectrum(size_t index, msdata::DetailLevel detailLevel) const;
        size_t size() const;
        const msdata::SpectrumIdentity& spectrumIdentity(size_t index) const;
        ///@}

    private:
        class Impl;
        boost::scoped_ptr<Impl> impl_;
    };

} // namespace analysis
} // namespace pwiz
#endif // _SPECTRUMLIST_SONAR_HPP
