#ifndef CARE_CORRECT_CPU_HPP
#define CARE_CORRECT_CPU_HPP

#include <config.hpp>

#include <serializedobjectstorage.hpp>
#include <options.hpp>
#include <readlibraryio.hpp>
#include <cpureadstorage.hpp>
#include <cpuminhasher.hpp>

namespace care{
namespace cpu{



	SerializedObjectStorage correct_cpu(
		const ProgramOptions& programOptions,
		CpuMinhasher& minhasher,
		CpuReadStorage& readStorage
	);

	


}
}

#endif
