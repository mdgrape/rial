# rialcpp

It binds the software simulation of `rial` to C++.

Both the floating point format and the approximation parameters are
configurable, and the result depends on them, so there is no plain
`rial::sqrt(x)`: you build an object for the configuration you want and call
the math functions on it. The approximation tables are generated once, when the
object is constructed, and reused by every call. Objects with different
configurations are independent, so mixing precisions in one program is fine.

```cpp
#include <rialcpp/math.hpp>

// the standard FP32 setup (== MathFuncConfig.standardFP32 on the scala side)
const rial::MathSim fp32;

// a BFloat16 setup with a plain lookup table (nOrder = 0)
const rial::MathSim bf16(rial::RealSpec::bfloat16(), rial::MathFuncConfig::bf16());

// ... or spell the table parameters out to match a specific MathFunctions module
const rial::MathSim custom(rial::RealSpec::float32(),
    rial::MathFuncConfig{/*nOrder=*/1, /*addrW=*/9, /*inBitW=*/14, /*outBitW=*/26});

std::cout << fp32.sqrt(2.0f) << " " << bf16.sqrt(2.0f) << std::endl;
```

`rial::RealSpec` (`exW`, `exBias`, `manW`) mirrors `rial.arith.RealSpec`, and
`rial::MathFuncConfig` (`nOrder`, `addrW`, `inBitW`, `outBitW`) mirrors
`rial.math.MathFuncConfig`, so a `MathSim` built from the same numbers as a
hardware `MathFunctions` module reproduces its results bit-exactly. Presets:
`RealSpec::float32() / float16() / bfloat16()`, `MathFuncConfig::fp32() / bf16()`.

`rial::MathSim` provides:

- `float sqrt      (float x)         ;`
- `float invsqrt   (float x)         ;`
- `float reciprocal(float x)         ;`
- `float sin       (float x)         ;`
- `float cos       (float x)         ;`
- `float acos      (float x)         ;`
- `float atan2     (float y, float x);`
- `float exp       (float x)         ;`
- `float log       (float x)         ;`

## Random number generators

The generators that turn raw random bits into a float depend on the floating
point format too, so they live on `rial::RngSim`, configured the same way
(`rial::RngConfig` carries how many random bits each generator consumes, the
order of the probit approximation, ...).

```cpp
#include <rialcpp/rng.hpp>

const rial::RngSim rng; // the standard FP32 setup
```

- `float rial::RngSim::gen_float12 (std::uint64_t rnd)                       ;` // uniform [1, 2)
- `float rial::RngSim::gen_float01 (std::uint64_t rnd)                       ;` // uniform [0, 1)
- `float rial::RngSim::gen_gaussian(std::uint64_t rnd)                       ;` // N(0, 1)
- `float rial::RngSim::irwin_hall  (const std::array<std::uint64_t, 12>& rnd);` // approx-gaussian

The counter-based generators produce raw bits and never touch the floating point
format, so they stay free functions:

- `std::uint32_t rial::squares32(std::uint64_t count, std::uint64_t key);`
- `std::array<std::uint32_t, 4> rial::threefry4x32(const std::array<std::uint32_t, 4>& ctr, const std::array<std::uint32_t, 4>& key, int rounds = 20);`

`squares32` has no round count because the hardware (`Squares32`) is a fixed
4-round pipeline.

## Notes

Constructing a `MathSim` generates every table on the scala side, which is not
cheap: keep the instance around instead of rebuilding it per call. Neither class
is thread-safe (all calls go through a single JNI environment), and both report
JVM / scala-side errors -- including an invalid configuration, which surfaces
when the object is constructed -- as `std::runtime_error`.

## Usage

First, make `rial.jar`.

```console
$ cd ../
$ sbt assembly
$ cd -
```

Then, build this repo and install it to your favorite location.

installation prefix MUST BE given via `CMAKE_INSTALL_PREFIX` if you use a directory other than the default.
If you use `cmake --install ./build --prefix=/path/to/somewhere`, it does not work.

```console
$ cmake -B ./build -DCMAKE_INSTALL_PREFIX=/home/${you}/.local
$ cmake --build ./build
$ cmake --install ./build
```

Then you can use `find_package(rialcpp)`. see `test/`.
