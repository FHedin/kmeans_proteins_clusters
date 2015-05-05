all:
	make release
debug:
	g++ -O0 -g -pg -std=c++11 -Wall -Wextra -pedantic kmeans.cpp -o kmeans
release:
	g++ -O3 -std=c++11 -ffast-math -Wall -Wextra -pedantic -march=native kmeans.cpp -o kmeans

clean:
	rm -f ./kmeans
