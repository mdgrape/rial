#
.PHONY: all libcrial clean

all : libcrial

libcrial:
	-mkdir build
	cd build; cmake ../clib -DCMAKE_INSTALL_PREFIX=../ ; make ; make install

clean:
	-/bin/rm -rf build
	-/bin/rm -rf include
	-/bin/rm -rf lib
	-/bin/rm -rf test_run_dir

