LIBS := -lSDL2

debug:
	g++ -march=native -std=c++17 -O0 -g -fopenmp -flax-vector-conversions src/main.cpp -o voxel_game $(LIBS)
profile:
	g++ -march=native -std=c++17 -O3 -g -fopenmp -flax-vector-conversions src/main.cpp -o voxel_game $(LIBS)
release:
	g++ -march=native -std=c++17 -O3 -fopenmp -flax-vector-conversions src/main.cpp -o voxel_game $(LIBS)
release_sse:
	g++ -march=westmere -std=c++17 -O3 -fopenmp -flax-vector-conversions src/main.cpp -o voxel_game $(LIBS)
clean:
	rm voxel_game
