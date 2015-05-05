all:
	make release
	#make debug
debug:
	#g++ -O0 -g -pg -std=c++11 -Wall -Wextra -pedantic -fopenmp kmeans.cpp -o kmeans
	g++ -O0 -g -pg -std=c++11 -Wall -Wextra -pedantic kmeans.cpp -o kmeans
	
release:
	#g++ -O3 -std=c++11 -ffast-math -Wall -Wextra -pedantic -fopenmp -march=native kmeans.cpp -o kmeans
	g++ -O2 -std=c++11 -mtune=native -Wall -Wextra -pedantic kmeans.cpp -o kmeans

clean:
	rm -f ./kmeans
