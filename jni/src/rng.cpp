#include "rialcpp/JVMManager.hpp"
#include "rialcpp/rng.hpp"

namespace rial
{

float gen_float12 (std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_float12 (rnd); }
float gen_float01 (std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_float01 (rnd); }
float gen_gaussian(std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_gaussian(rnd); }

} // rial
