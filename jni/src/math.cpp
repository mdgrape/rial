#include "rialcpp/JVMManager.hpp"
#include "rialcpp/math.hpp"

namespace rial
{

float sqrt       (float x)          {return detail::JVMManager::get_instance().call_sqrt       (x);}
float invsqrt    (float x)          {return detail::JVMManager::get_instance().call_invsqrt    (x);}
float reciprocal (float x)          {return detail::JVMManager::get_instance().call_reciprocal (x);}
float sin        (float x)          {return detail::JVMManager::get_instance().call_sin        (x);}
float cos        (float x)          {return detail::JVMManager::get_instance().call_cos        (x);}
float acos       (float x)          {return detail::JVMManager::get_instance().call_acos       (x);}
float atan2      (float y, float x) {return detail::JVMManager::get_instance().call_atan2      (y, x);}
float exp        (float x)          {return detail::JVMManager::get_instance().call_exp        (x);}
float log        (float x)          {return detail::JVMManager::get_instance().call_log        (x);}

} // rial
