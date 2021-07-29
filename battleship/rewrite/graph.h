#pragma once

#include <array>
#include <cassert>
#include <iostream>
#include <utility>
using std::array;
using std::cout;
using std::endl;
using std::pair;

// N is the total number of ships
// MAX_PLACEMENTS is an upper bound on the number of placements for a ship
// In most variants of battleship, WIDTH*(HEIGHT - 1) + HEIGHT*(WIDTH - 1) can
// be used for MAX_PLACEMENTS
template <size_t N, size_t MAX_PLACEMENTS> struct placement_graph {
	// for each ship, store the number of valid placements
	// each should be upper bounded by MAX_PLACEMENTS
	array<size_t, N> placement_counts;
	// for each pair of ships p<q and placements u and v of each respectively,
	// store whether there is an edge between u and v in the placement graph.
	// (i.e. TRUE if the two placements do not intersect, else FALSE)
	// the argument order should be adj[p][u][q][v]
	// For q<=p, the value will always be false
	array<array<array<array<bool, MAX_PLACEMENTS>, N>, MAX_PLACEMENTS>, N> adj;
	// constructor initializes adj to all true
	placement_graph() {
		placement_counts.fill(0);
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < MAX_PLACEMENTS; ++u) {
				for (size_t q = 0; q <= p; ++q)
					for (size_t v = 0; v < MAX_PLACEMENTS; ++v)
						adj[p][u][q][v] = false;
				for (size_t q = p + 1; q < N; ++q)
					for (size_t v = 0; v < MAX_PLACEMENTS; ++v)
						adj[p][u][q][v] = true;
			}
	}
	void set_placement_count(size_t ship, size_t count) {
		assert(ship <= N);
		assert(count <= MAX_PLACEMENTS);
		placement_counts[ship] = count;
		for (size_t u = count + 1; u < MAX_PLACEMENTS; ++u)
			for (size_t q = 0; q < N; ++q)
				for (size_t v = 0; v < MAX_PLACEMENTS; ++v) {
					// Set all out-of-bounds placements to not have edges by default
					adj[ship][u][q][v] = false;
					adj[q][v][ship][u] = false;
				}
	}
	// remove an edge from (ship p, placement u) to (ship q, placement v)
	void remove_edge(size_t p, size_t u, size_t q, size_t v) {
		adj[p][u][q][v] = false;
	}
};

// Data structure to store a set of 'active' vertices within a valid placement
// subgraph, used for clique iteration
template <size_t N, size_t MAX_PLACEMENTS, typename HIT_TYPE, typename fast_set>
struct iteration_state {
	// size_t current_depth = 0;
	// For clique iteration at depth k and ship p>=k, vertex_subsets[k][p] stores
	// the subset of placements of p which could correspond to vertices,
	// as well as the total number of hits covered by all placements of ships q<k
	// current_depth is the maximum value of k for which this is valid
	array<pair<array<fast_set, N>, HIT_TYPE>, N + 1> vertex_subsets;
	// Destructively get the next placement from
	// vertex_subsets[current_depth].first
	/*
	template <size_t current_depth> size_t get_next_placement() {
		return vertex_subsets[current_depth]
				.first[current_depth]
				.bitscan_destructive_any(MAX_PLACEMENTS);
	}
	*/
	template <size_t current_depth> bool is_valid_placement(size_t u) {
		return vertex_subsets[current_depth].first[current_depth].test(u);
	}
	template <size_t current_depth> HIT_TYPE counted_hits() {
		return vertex_subsets[current_depth].second;
	}
	iteration_state() {
		for (size_t p = 0; p <= N; ++p)
			vertex_subsets[p].second = 0;
	}
};
// Data structure used to store total appearance counts of each vertex within a
// valid placement subgraph, used for clique iteration
template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE>
struct result_state {
	array<array<COUNT_TYPE, MAX_PLACEMENTS>, N> appearance_counts;
	void add_result(const size_t& p, const size_t& u, const COUNT_TYPE& count) {
		appearance_counts[p][u] += count;
	}
	result_state() {
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < MAX_PLACEMENTS; ++u)
				appearance_counts[p][u] = 0;
	}
	void sum_from(const result_state& rs) {
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < MAX_PLACEMENTS; ++u)
				appearance_counts[p][u] += rs.appearance_counts[p][u];
	}
};

// Note that unlike placement_graph, the valid_placement_graph will be
// transferred to GPU memory when run on a GPU
// HIT_TYPE should be some relatively small integer type (e.g. unsigned char)
// capable of storing a value up to size WIDTH * HEIGHT
// fast_set should be some sort of bitset, std::bitset<MAX_PLACEMENTS> is fine
// but possibly suboptimal due to lack of speedy bitscan operation
template <size_t N, size_t MAX_PLACEMENTS, typename HIT_TYPE, typename fast_set>
struct valid_placement_subgraph {
	// for each ship, store the number of valid placements
	// each should be upper bounded by MAX_PLACEMENTS
	array<size_t, N> placement_counts;
	// for each pair of ships p<q and a placement u of ship p,
	// store all the placements v of ship q such that there is an edge between the
	// two corresponding vertices in a bitset
	array<array<array<fast_set, N>, MAX_PLACEMENTS>, N> adj;
	array<array<HIT_TYPE, MAX_PLACEMENTS>, N> hit_counts;
	HIT_TYPE total_hit_count = 0;
	valid_placement_subgraph(const placement_graph<N, MAX_PLACEMENTS>& pg)
			: placement_counts(pg.placement_counts) {
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < MAX_PLACEMENTS; ++u)
				for (size_t q = 0; q < N; ++q)
					for (size_t v = 0; v < MAX_PLACEMENTS; ++v) {
						// adj[p][u][q].assign(v, pg.adj[p][u][q][v]);
						adj[p][u][q].set(v, pg.adj[p][u][q][v]);
					}
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < MAX_PLACEMENTS; ++u)
				hit_counts[p][u] = 0;
	}
	// remove a vertex for ship p placement u
	void remove_vertex(size_t p, size_t u) {
		assert(p <= N);
		assert(u <= placement_counts[p]);
		for (size_t q = 0; q < N; ++q)
			for (size_t v = 0; v < MAX_PLACEMENTS; ++v) {
				adj[p][u][q][v] = false;
				adj[q][v][p][u] = false;
			}
	}
	void set_hit_count(size_t p, size_t u, HIT_TYPE count) {
		assert(p <= N);
		assert(u <= placement_counts[p]);
		hit_counts[p][u] = count;
	}
	void increment_hit_count(size_t p, size_t u) {
		assert(p <= N);
		assert(u <= placement_counts[p]);
		hit_counts[p][u]++;
	}
	void increment_total_hit_count() { total_hit_count++; }
	void set_total_hit_count(HIT_TYPE count) { total_hit_count = count; }
	void
	init_is(iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is) const {
		for (size_t depth = 0; depth <= N; ++depth)
			for (size_t p = 0; p < N; ++p) {
				for (size_t u = 0; u < placement_counts[p]; ++u)
					is.vertex_subsets[depth].first[p].set(u);
				for (size_t u = placement_counts[p]; u < MAX_PLACEMENTS; ++u)
					is.vertex_subsets[depth].first[p].reset(u);
			}
	}
	// decrease current depth of pr
	template <size_t current_depth>
	void decrease_depth(
			iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is) const {
		// is.current_depth--;
	}
	// increase current depth of pr, using a placement u of ship pr.current_depth
	template <size_t current_depth>
	void
	increase_depth(iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
								 size_t u) const {
		constexpr size_t p = current_depth;
		is.vertex_subsets[p + 1].second =
				is.vertex_subsets[p].second + hit_counts[p][u];
		for (size_t q = p + 1; q < N; ++q) {
			is.vertex_subsets[p + 1].first[q] =
					is.vertex_subsets[p].first[q] & adj[p][u][q];
		}
		// is.current_depth++;
	}
};
