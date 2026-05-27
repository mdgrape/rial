#include "rialcpp/JVMManager.hpp"
#include "rialcpp/rng.hpp"

namespace rial
{

float gen_float12 (std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_float12 (rnd); }
float gen_float01 (std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_float01 (rnd); }
float gen_gaussian(std::uint64_t rnd) { return detail::JVMManager::get_instance().call_gen_gaussian(rnd); }

std::uint32_t squares32(std::uint64_t count, std::uint64_t key)
{
    return detail::JVMManager::get_instance().call_squares32(count, key);
}

float irwin_hall(const std::array<std::uint64_t, 12>& rnds)
{
    return detail::JVMManager::get_instance().call_irwin_hall(rnds.data(), rnds.size());
}

} // rial
