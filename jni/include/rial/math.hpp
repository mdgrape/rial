#ifndef RIAL_JNI_MATH_HPP
#define RIAL_JNI_MATH_HPP

namespace rial
{

float sqrt       (float x)         ;
float invsqrt    (float x)         ;
float reciprocal (float x)         ;
float sin        (float x)         ;
float cos        (float x)         ;
float acosPhase1 (float x)         ;
float acosPhase2 (float x)         ;
float acos       (float x)         ;
float atan2Phase1(float y, float x);
float atan2Phase2(float x)         ;
float atan2      (float y, float x);
float exp        (float x)         ;
float log        (float x)         ;

} // rial
#endif//RIAL_JNI_MATH_HPP
