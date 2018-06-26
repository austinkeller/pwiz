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

#ifndef _SPECTRUMLISTTHREADSAFE_HPP_
#define _SPECTRUMLISTTHREADSAFE_HPP_


#include "pwiz/utility/misc/Export.hpp"
#include "MSData.hpp"
#include "SpectrumListWrapper.hpp"
#include <boost/thread.hpp>


namespace pwiz {
namespace msdata {

/// rate limits the requests of identical spectra
class PWIZ_API_DECL SpectrumListThreadsafe : public SpectrumListWrapper
{
public:

    SpectrumListThreadsafe(const SpectrumListPtr& inner);

    /// returns the requested spectrum which may or may not be cached depending on
    /// the current cache mode
    virtual SpectrumPtr spectrum(size_t index, bool getBinaryData = false) const;

protected:
    //mutable std::vector<boost::mutex> indexedMutex_;
    mutable boost::mutex mutex_;

private:
    SpectrumListThreadsafe(SpectrumListThreadsafe&);
    SpectrumListThreadsafe& operator=(SpectrumListThreadsafe&);
};


} // namespace msdata
} // namespace pwiz


#endif // _SPECTRUMLISTTHREADSAFE_HPP_
