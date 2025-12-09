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
