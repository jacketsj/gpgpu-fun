#include <iostream>
using std::cout;
#include <bitset>

#include "algorithm.h"
// #include "fast_bitset.h"
#include "graph.h"
#include "standard_game.h"

#define WIDTH 8
#define HEIGHT 8
#define n 5
const array<size_t, n> lengths = {5, 4, 3, 2, 2};
// #define WIDTH 10
// #define HEIGHT 10
// #define n 5
// const array<size_t, n> lengths = {5, 4, 3, 3, 2};

#define MAX_PLACEMENTS WIDTH*(HEIGHT - 1) + HEIGHT*(WIDTH - 1)
#define fast_set fast_bitset::bitset<MAX_PLACEMENTS>
#define fast_set std::bitset<MAX_PLACEMENTS>
#define HIT_TYPE short
#define COUNT_TYPE long long

int main() {
	standard_game<n, WIDTH, HEIGHT, HIT_TYPE, fast_set> sg(lengths);
	iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set> is;
	result_state<n, MAX_PLACEMENTS, COUNT_TYPE> rs;
	sg.vps.init_is(is);
	COUNT_TYPE total_configurations =
			unroll_algorithm<n, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
											 n>::place_ship(sg.vps, is, rs);
	cout << "total configurations: " << total_configurations << '\n';
	sg.print_output(rs);
}
