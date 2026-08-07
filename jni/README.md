# rialcpp

It binds the software simulation of `rial` to C++.

The floating point format, the approximation parameters and the set of
supported functions are all configurable, and the result depends on all three,
so there is no plain `rial::sqrt(x)`: you build an object for the configuration
you want and call the math functions on it. The approximation tables are
generated once, when the object is constructed, and reused by every call.
Objects with different configurations are independent, so mixing precisions in
one program is fine.

```cpp
#include <rialcpp/math.hpp>

// the standard FP32 setup (== MathFuncConfig.standardFP32 on the scala side)
const rial::MathSim fp32;

// a BFloat16 setup with a plain lookup table (nOrder = 0)
const rial::MathSim bf16(rial::RealSpec::bfloat16(), rial::MathFuncConfig::bf16());

// FP32, but supporting sqrt/invsqrt/reciprocal only
const rial::MathSim simple(rial::RealSpec::float32(),
    rial::MathFuncConfig::fp32(rial::Func::simple));

// ... or spell everything out to match a specific MathFunctions module
const rial::MathSim custom(rial::RealSpec::float32(), rial::MathFuncConfig{
    rial::Func::sqrt | rial::Func::exp,
    /*nOrder=*/1, /*addrW=*/9, /*inBitW=*/14, /*outBitW=*/26});

std::cout << fp32.sqrt(2.0f) << " " << bf16.sqrt(2.0f) << std::endl;
```

`rial::RealSpec` (`exW`, `exBias`, `manW`) mirrors `rial.arith.RealSpec`, and
`rial::MathFuncConfig` (`funcs`, `nOrder`, `addrW`, `inBitW`, `outBitW`, plus
the ScaleMixtureGaussian sigmas) mirrors `rial.math.MathFuncConfig`, so a
`MathSim` built from the same numbers as a hardware `MathFunctions` module
reproduces its results bit-exactly. Presets: `RealSpec::float32() / float16() /
bfloat16()`, `MathFuncConfig::fp32() / bf16() / fp32_all() / bf16_all()`.

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
- `float sigmoid   (float x)         ;`
- `float softplus  (float x)         ;`
- `float scale_mixture_gaussian(float x);`

### The function set matters

`MathFuncConfig::funcs` is a bit set of `rial::Func` (`Func::sqrt | Func::exp`,
or the presets `Func::simple / standard / all`, mirroring the ones in the
`MathFuncConfig` companion object). It is not just an access-control switch:
the tables of the supported functions share one merged set of coefficient and
intermediate bit widths, the same way the hardware shares one polynomial
evaluator, so **a narrower function set can change the results**. In FP32, sqrt
from a sqrt-only instance disagrees with the `standard` one on about 0.5% of
inputs, by a ULP. (Not every subset changes something: `Func::simple` already
contains reciprocal, whose table dominates the merge in FP32, so its sqrt is
identical to `standard`'s.)

Functions pull in what they need -- `Func::acos` also gives you `sqrt`, since
its first phase reuses the sqrt table, and `Func::atan2` gives you
`reciprocal` -- exactly as in the hardware. `config().funcs` reports the
effective set after those additions, and `supports(Func)` queries it. Calling a
function outside the set throws `std::runtime_error`.

### Formats a `float` cannot carry

The `float` interface above is fine as long as the emulated format fits in an
IEEE-754 binary32 (FP32, BF16, FP16, ...). For anything wider there is a
parallel interface working on the raw bit pattern of `spec()`, as a
`std::uint64_t`: `sqrt_bits`, `invsqrt_bits`, ... `atan2_bits(y, x)`, and
`gen_float12_bits` / ... / `irwin_hall_bits` on `RngSim`. Bits above the spec
width are ignored.

`to_bits` / `to_double` / `to_float` convert in and out using the instance's
own spec (`rialcpp/real.hpp` has the same as free functions taking a
`RealSpec`, for when you have a bit pattern but no instance).

```cpp
const rial::MathSim fp64(rial::RealSpec{11, 0x3FF, 52},
    rial::MathFuncConfig{rial::Func::sqrt, /*nOrder=*/2, /*addrW=*/8,
                         /*inBitW=*/44, /*outBitW=*/55});

std::cout << fp64.to_double(fp64.sqrt_bits(fp64.to_bits(2.0))) << std::endl;
```

The two interfaces agree wherever both apply. Note that `double` itself tops
out at 52 mantissa bits, so a still wider format can only be inspected as bits
-- though a table-based approximation is not what you would reach for at that
precision anyway.

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

Then build.

```console
$ cmake -B ./build
$ cmake --build ./build
$ ctest --test-dir ./build
```

To install, pick a prefix whichever way you like.

```console
$ cmake --install ./build --prefix /path/to/dir
```

Then you can use `find_package(rialcpp)`.

### path to `rial.jar`

`<directory holding librialcpp>/../share/rialcpp/rial.jar` is automatically loaded.

Set `RIALCPP_JAR` to a path to override the lookup.
