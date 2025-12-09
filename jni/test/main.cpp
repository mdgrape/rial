#include <rialcpp/math.hpp>

#include <iomanip>
#include <iostream>

#include <cmath>

template<typename Fref, typename Frial>
void check(std::string name, float from, float to,
        Fref f_ref, Frial f_rial)
{
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << name << std::endl;
    std::cout << "x libm rial diff" << std::endl;
    for(std::size_t i=0; i<100; ++i)
    {
        const float x = 0.01f * (from * i + to * (100 - i));
        const auto libm = f_ref(x);
        const auto rial = f_rial(x);
        const auto diff = std::abs(libm - rial);

        std::cout << std::setw(15) << std::setprecision(10) << x    << " "
                  << std::setw(15) << std::setprecision(10) << libm << " "
                  << std::setw(15) << std::setprecision(10) << rial << " "
                  << diff << std::endl;
    }
}

int main()
{
    constexpr float pi = 3.14159265;

    check("sqrt", 0.0, 4.0,
        [](float x) {return std::sqrt(x); },
        [](float x) {return rial::sqrt(x);}
    );

    check("invsqrt", 0.5, 4.0,
        [](float x) {return 1.0f / std::sqrt(x); },
        [](float x) {return rial::invsqrt(x);}
    );

    check("reciprocal", 0.5, 4.0,
        [](float x) {return 1.0f / x; },
        [](float x) {return rial::reciprocal(x);}
    );

    check("sin", -pi, pi,
        [](float x) {return std::sin(x); },
        [](float x) {return rial::sin(x);}
    );

    check("cos", -pi, pi,
        [](float x) {return std::cos(x); },
        [](float x) {return rial::cos(x);}
    );

    check("acos", -1.0, 1.0,
        [](float x) {return std::acos(x); },
        [](float x) {return rial::acos(x);}
    );


    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "atan2" << std::endl;
    std::cout << "x y libm rial diff" << std::endl;
    for(std::size_t i=0; i<100; ++i)
    {
    for(std::size_t j=0; j<100; ++j)
    {
        const float x = 0.01f * (-1.0f * i + 1.0f * (100 - i));
        const float y = 0.01f * (-1.0f * j + 1.0f * (100 - j));
        const auto libm = std::atan2(y, x);
        const auto rial = rial::atan2(y, x);
        const auto diff = std::abs(libm - rial);

        std::cout << std::setw(15) << std::setprecision(10) << x    << " "
                  << std::setw(15) << std::setprecision(10) << y    << " "
                  << std::setw(15) << std::setprecision(10) << libm << " "
                  << std::setw(15) << std::setprecision(10) << rial << " "
                  << diff << std::endl;
    }
    }

    check("exp", 0.1, 3.1,
        [](float x) {return std::exp(x); },
        [](float x) {return rial::exp(x);}
    );

    check("log", 0.1, 3.1,
        [](float x) {return std::log(x); },
        [](float x) {return rial::log(x);}
    );

    return 0;
}
