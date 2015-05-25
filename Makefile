all:
	cd build; make install VERBOSE=1;

clean:
	rm -rf build

cmake-release: clean
	mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Release ..

cmake-debug: clean
	mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Debug ..