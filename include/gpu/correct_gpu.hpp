#ifndef CARE_CORRECT_GPU_HPP
#define CARE_CORRECT_GPU_HPP

#include <gpu/gpuminhasher.cuh>
#include <gpu/gpureadstorage.cuh>
#include <gpu/forest_gpu.cuh>

#include <config.hpp>
#include <serializedobjectstorage.hpp>

#include <options.hpp>

#include <mutex>
#include <memory>
#include <vector>

namespace care {
namespace gpu {


SerializedObjectStorage correct_gpu(
	const ProgramOptions& programOptions,
	GpuMinhasher& minhasher,
	GpuReadStorage& readStorage,
	const std::vector<GpuForest>& anchorForests,
	const std::vector<GpuForest>& candidateForests
);

}
}

#endif
