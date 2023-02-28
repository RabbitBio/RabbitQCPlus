#ifndef GLOBAL_CUDA_STREAM_POOL
#define GLOBAL_CUDA_STREAM_POOL

#include <rmm/cuda_device.hpp>
#include <rmm/cuda_stream_pool.hpp>

#include <map>
#include <mutex>

namespace streampool {

namespace detail {

inline rmm::cuda_stream_pool* initial_pool()
{
  static rmm::cuda_stream_pool pool{};
  return &pool;
}

inline std::mutex& map_lock()
{
  static std::mutex map_lock;
  return map_lock;
}

inline auto& get_map()
{
  static std::map<rmm::cuda_device_id::value_type, rmm::cuda_stream_pool*> device_id_to_pool;
  return device_id_to_pool;
}

}  // namespace detail

inline rmm::cuda_stream_pool* get_per_device_pool(rmm::cuda_device_id device_id)
{
  std::lock_guard<std::mutex> lock{detail::map_lock()};
  auto& map = detail::get_map();
  // If a resource was never set for `id`, set to the initial resource
  auto const found = map.find(device_id.value());
  return (found == map.end()) ? (map[device_id.value()] = detail::initial_pool())
                              : found->second;
}

inline rmm::cuda_stream_pool* set_per_device_pool(rmm::cuda_device_id device_id,
                                                       rmm::cuda_stream_pool* new_pool)
{
  std::lock_guard<std::mutex> lock{detail::map_lock()};
  auto& map          = detail::get_map();
  auto const old_itr = map.find(device_id.value());
  // If a resource didn't previously exist for `id`, return pointer to initial_resource
  auto* old_pool           = (old_itr == map.end()) ? detail::initial_pool() : old_itr->second;
  map[device_id.value()] = (new_pool == nullptr) ? detail::initial_pool() : new_pool;
  return old_pool;
}

inline rmm::cuda_stream_pool* get_current_device_pool()
{
  return get_per_device_pool(rmm::detail::current_device());
}

inline rmm::cuda_stream_pool* set_current_device_pool(rmm::cuda_stream_pool* new_pool)
{
  return set_per_device_pool(rmm::detail::current_device(), new_pool);
}

}

#endif
