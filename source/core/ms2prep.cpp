/*
 * Copyright (C) 2021  Muhammad Haseeb, and Fahad Saeed
 * Florida International University, Miami, FL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include "ms2prep.hpp"

//
// TODO: insert instrumentation for the work performed in this file
// Specially measure the compute to overhead (comm and I/O) ratio.
//

// external query filenames
extern vector<string_t> queryfiles;

// extern params
extern gParams params;

namespace hcp
{

namespace mpi
{

//
// FUNCTION: getPartitionSize (templated over numbers)
//
template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
T getPartitionSize(T isize)
{
    T loc_size = isize / params.nodes;

    if (params.myid < (isize % params.nodes))
        loc_size++;

    return loc_size;
}

//
// FUNCTION: barrier (wrapper for MPI_Barrier)
//
status_t barrier()
{
    return MPI_Barrier(MPI_COMM_WORLD);
}

//
// TODO: still to implement
//
status_t allgather()
{
    return SLM_SUCCESS;
}

} // namespace mpi

namespace ms2
{

//
// FUNCTION: synchronize
//
status_t synchronize()
{
    status_t status = SLM_SUCCESS;

    // synchronize
    status = hcp::mpi::barrier();

    // write index to file
    MSQuery::write_index();

    return status;
}

//
// FUNCTION: get_instance
// 
MSQuery **& get_instance()
{
    static MSQuery** ptrs = new MSQuery*[queryfiles.size()];
    return ptrs;
}

//
// FUNCTION: initialize
// 
status_t initialize(lwqueue<MSQuery *>** qfPtrs, int_t& nBatches, int_t& dssize)
{
    status_t status = SLM_SUCCESS;
    bool_t qfindex = false;

    int_t nfiles = queryfiles.size();

    // get the ptrs instance
    MSQuery **ptrs = get_instance();

    if (ptrs == nullptr)
        status = ERR_INVLD_PTR;

    if (status == SLM_SUCCESS)
    {
        /* Initialize the queue with already created nfiles */
        *qfPtrs = new lwqueue<MSQuery*>(nfiles, false);

        // initialize the ptrs instances
        for (auto lptr = ptrs; lptr < ptrs + nfiles; lptr++)
            *lptr = new MSQuery;

        // get local partition size
        auto pfiles = hcp::mpi::getPartitionSize(nfiles);

        if (params.nodes > 1)
            MSQuery::init_index();

        // if pfiles > 0
        if (pfiles) 
        {
            // fill the vector with locally processed MS2 file indices
            std::vector<int_t> ms2local(pfiles);

            // first element is params.myid
            ms2local[0] = params.myid;

            // rest of the files in cyclic order
            std::generate(std::begin(ms2local) + 1, std::end(ms2local), 
                          [n=params.myid] () mutable { return n += params.nodes; });

#ifdef USE_OMP
#pragma omp parallel for schedule (dynamic, 1)
#endif/* _OPENMP */
            for (auto fid = 0; fid < pfiles; fid++)
            {
                auto loc_fid = ms2local[fid];
                ptrs[loc_fid]->initialize(&queryfiles[loc_fid], loc_fid);

                // archive the index variables
                ptrs[loc_fid]->archive(loc_fid);
            }
        }

        //
        // Synchronize superstep 2
        //
        if (params.nodes > 1)
        {
            status = hcp::ms2::synchronize();
            info_t *findex = new info_t[nfiles];

            MSQuery::read_index(findex, nfiles);

            // copy global data from the index
            for (auto fid = 0; fid < nfiles; fid ++)
            {
                // use the index information to
                // initialize the remaining index
                if (!ptrs[fid]->isinit())
                {
                    ptrs[fid]->Info() = findex[fid];
                    ptrs[fid]->vinitialize(&queryfiles[fid], fid);
                }
            }
            delete[] findex;
        }

        // -------------------------------------------------------------------------------------- //

        //
        // create/read the MS2 index
        //

        // Push zeroth as is
        (*qfPtrs)->push(ptrs[0]);
        dssize += ptrs[0]->getQAcount();

        // Update batch numbers
        for (auto fid = 1; fid < nfiles; fid++)
        {
            ptrs[fid]->Curr_chunk() = ptrs[fid - 1]->Curr_chunk() + ptrs[fid - 1]->Nqchunks();
            (*qfPtrs)->push(ptrs[fid]);
            dssize += ptrs[fid]->getQAcount();
        }

        // Compute the total number of batches in the dataset
        nBatches = ptrs[nfiles-1]->Curr_chunk() + ptrs[nfiles-1]->Nqchunks();

        if (params.myid == 0)
            std::cout << "\nDataset Size = " << dssize << std::endl << std::endl;
    }
    
    return status;
}

//
// FUNCTION: deinitialize
//
void deinitialize()
{
    MSQuery **ptrs = get_instance();

    /* Delete ptrs */
    if (ptrs != nullptr)
    {
        delete[] ptrs;
        ptrs = nullptr;
    }
}

} // namespace ms2
} // namespace hcp
