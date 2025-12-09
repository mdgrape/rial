#include "rial/JVMManager.hpp"
#include "rial/math.hpp"

namespace rial
{

float sqrt       (float x)          {return detail::JVMManager::get_instance().call_sqrt       (x);}
float invsqrt    (float x)          {return detail::JVMManager::get_instance().call_invsqrt    (x);}
float reciprocal (float x)          {return detail::JVMManager::get_instance().call_reciprocal (x);}
float sin        (float x)          {return detail::JVMManager::get_instance().call_sin        (x);}
float cos        (float x)          {return detail::JVMManager::get_instance().call_cos        (x);}
float acosPhase1 (float x)          {return detail::JVMManager::get_instance().call_acosPhase1 (x);}
float acosPhase2 (float x)          {return detail::JVMManager::get_instance().call_acosPhase2 (x);}
float acos       (float x)          {return detail::JVMManager::get_instance().call_acos       (x);}
float atan2Phase1(float y, float x) {return detail::JVMManager::get_instance().call_atan2Phase1(y, x);}
float atan2Phase2(float x)          {return detail::JVMManager::get_instance().call_atan2Phase2(x);}
float atan2      (float y, float x) {return detail::JVMManager::get_instance().call_atan2      (y, x);}
float exp        (float x)          {return detail::JVMManager::get_instance().call_exp        (x);}
float log        (float x)          {return detail::JVMManager::get_instance().call_log        (x);}

} // rial
