#include "lbe.h"

#if defined (USE_TIMEMORY)
#include "timemory/timemory.hpp"

// shorthand
using namespace tim::component;

//
//--------------------------------------------------------------------------------------//
//
TIMEMORY_DECLARE_COMPONENT(inst_per_cycle)
TIMEMORY_STATISTICS_TYPE(component::inst_per_cycle, double)
TIMEMORY_DEFINE_CONCRETE_TRAIT(is_available, component::inst_per_cycle, true_type)
//
//--------------------------------------------------------------------------------------//
//
namespace tim
{
namespace component
{
struct inst_per_cycle : public base<inst_per_cycle, std::array<long long, 2>>
{
    using value_type = std::array<long long, 2>;
    using hw_t       = tim::component::papi_tuple<PAPI_TOT_INS, PAPI_TOT_CYC>;

    static std::string label() { return "inst_per_cycle"; }
    static std::string description() { return "number of instructions per cycle"; }
    static void        thread_init(storage_type*) { hw_t::thread_init(nullptr); }

    void start()
    {
        m_hw.start();
        value = m_hw.get_value();
    }
    void stop()
    {
        m_hw.stop();
        value = m_hw.get_value();
        accum = m_hw.get_accum();
    }
    double get() const { return (accum[1] > 0.0) ? (accum[0] / (1.0 * accum[1])) : 0.0; }

private:
    hw_t m_hw;
};
}  // namespace component
}  // namespace tim

struct hicops_tag
{};

// custom user bundle
using bundle_t = user_bundle<10, hicops_tag>;

// create tuple aliases
using time_tuple_t = tim::auto_tuple<wall_clock, cpu_util, thread_cpu_util, inst_per_cycle, bundle_t>;
using mem_tuple_t = tim::auto_tuple<peak_rss, virtual_memory, written_bytes, read_bytes, num_io_in, num_io_out>;
using wall_tuple_t = tim::auto_tuple<wall_clock>;

#endif // USE_TIMEMORY

#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */

#if defined (USE_TIMEMORY)
#define MARK(mark)
#define ELAPSED(es, m1, m2)
#define ELAPSED_SECONDS(es)
#define PRINT_ELAPSED(es)
#else
#define MARK(mark)                mark = chrono::system_clock::now()
#define ELAPSED(es, m1, m2)       es = m2 - m1
#define ELAPSED_SECONDS(es)       es.count()
#define PRINT_ELAPSED(es)         std::cout << "Elapsed Time: " << es.count() << "s" << std::endl << std::endl
#endif