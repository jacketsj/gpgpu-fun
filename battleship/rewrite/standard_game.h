#pragma once

#include <array>
using std::array;

#include "graph.h"

// TODO: What info is in each of these again
struct query_miss {};

struct query_hit {};

struct query_sink {};

// The standard game is one that involves ships of 1xk size on a WIDTH * HEIGHT
// sized board. k may vary per-ship, but k>=2 is required.
tempalte<size_t N, size_t WIDTH, size_t HEIGHT> struct standard_game {
	array<size_t, N> lengths;
	constexpr MAX_PLACEMENTS = WIDTH * (HEIGHT - 1) + HEIGHT * (WIDTH - 1);
	valid_placement_subgraph<N, MAX_PLACEMENTS> vps;
	static placement_graph build_ps() {
		placement_graph pg;
		// TODO: Generate pg properly, add in template params
		return pg;
	}
	standard_game(array<size_t, N>&& _lengths)
			: lengths(_lengths), vps(build_ps()) {}
	// TODO: get_overlapping_placements helper function (given a square)
	// TODO: square struct probably
	void update(query_miss q) {
		// TODO: Call remove_vertex a bunch
	}
	void update(query_hit q) {
		// TODO: Call remove_vertex a bunch
		// TODO: Call vps.increment_hit_count() a bunch of times
	}
	void update(query_sink q) {
		// TODO: Call remove_vertex a bunch
	}
};
