//
// $Id$
//
//
// Original author: Jarrett Egertson <jegertso .@. uw.edu>
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

#define PWIZ_SOURCE

#include "DemuxSolver.hpp"
#include "nnls.h"
#ifdef _PROFILE_PERFORMANCE
#include <iostream>
#endif

namespace pwiz{
namespace analysis{
    typedef NNLS<DemuxTypes::MatrixType> NNLSType;
    void NNLSSolver::Solve(const MatrixPtr& masks, const MatrixPtr& signal, MatrixPtr& solution) const
    {
        NNLSType solver(*masks, numIters_, eps_);
        int numCols = static_cast<int>(signal->cols()); // OpenMP 2.0 only allows signed index variables in for loops

#ifdef _PROFILE_PERFORMANCE
#ifdef _PROFILE_SOLVER_CONVERGENCE

        if (iterationCaching_)
            numItersUsed_.resize(numCols, 0);

#endif
#endif

#pragma omp parallel for firstprivate(solver) schedule(dynamic)
        for (int fragIndex = 0; fragIndex < numCols; ++fragIndex)
        {
            solver.solve(signal->col(fragIndex));
            solution->col(fragIndex).noalias() = solver.x();

#ifdef _PROFILE_PERFORMANCE
#ifdef _PROFILE_SOLVER_CONVERGENCE

            if (iterationCaching_)
                numItersUsed_[fragIndex] = solver.numLS();

#endif
#endif

        }


    }

#ifdef _PROFILE_PERFORMANCE
#ifdef _PROFILE_SOLVER_CONVERGENCE

    void NNLSSolver::setIterationCaching(bool on)
    {
        iterationCaching_ = on;
    }

    void NNLSSolver::writeIterations(std::ostream& out)
    {
        out << "DemuxIterations:";
        for (size_t i = 0; i < numItersUsed_.size(); ++i)
        {
            out << numItersUsed_.at(i);
            if (i + 1 != numItersUsed_.size())
                out << ",";
            else
                out << "\n";
        }
    }

#endif
#endif

} // namespace analysis
} // namespace pwiz