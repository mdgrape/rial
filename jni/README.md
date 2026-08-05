# rialcpp

it binds the following software simulation functions to c++.

- `float rial::sqrt      (float x)         ;`
- `float rial::invsqrt   (float x)         ;`
- `float rial::reciprocal(float x)         ;`
- `float rial::sin       (float x)         ;`
- `float rial::cos       (float x)         ;`
- `float rial::acos      (float x)         ;`
- `float rial::atan2     (float y, float x);`
- `float rial::exp       (float x)         ;`
- `float rial::log       (float x)         ;`

It also binds the random-number generators (input is raw random bits, standard
FP32 config):

- `float         rial::gen_float12 (std::uint64_t rnd)                       ;` // uniform [1, 2)
- `float         rial::gen_float01 (std::uint64_t rnd)                       ;` // uniform [0, 1)
- `float         rial::gen_gaussian(std::uint64_t rnd)                       ;` // N(0, 1)
- `std::uint32_t rial::squares32   (std::uint64_t count, std::uint64_t key)  ;` // uniform u32
- `float         rial::irwin_hall  (const std::array<std::uint64_t, 12>& rnd);` // approx-gaussian
- `std::array<std::uint32_t, 4> rial::threefry4x32(const std::array<std::uint32_t, 4>& ctr, const std::array<std::uint32_t, 4>& key);` // Threefry-4x32-20

Currently, the configurations (e.g. table addrW, extraBits) are fixed.
So the result might be different from customized MathFunctions instance.

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
