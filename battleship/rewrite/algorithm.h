#pragma once

#include "graph.h"

#include <iostream>
#include <vector>
using std::cerr;
using std::cout;
using std::endl;
using std::vector;

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set, size_t i>
struct unroll_algorithm {
	static COUNT_TYPE
	place_ship(const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE,
																						fast_set>& vps,
						 iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
						 result_state<N, MAX_PLACEMENTS, COUNT_TYPE>& rs) {
		size_t result = 0;
		size_t u;
		while ((u = is.get_next_placement()) != MAX_PLACEMENTS) {
			vps.increase_depth(is, u);
			COUNT_TYPE current_result =
					unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
													 i - 1>::place_ship(vps, is, rs);
			vps.decrease_depth(is);
			rs.add_result(is.current_depth, u, current_result);
			result += current_result;
		}
		return result;
	}
};

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set>
struct unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set, 0u> {
	static COUNT_TYPE
	place_ship(const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE,
																						fast_set>& vps,
						 iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
						 result_state<N, MAX_PLACEMENTS, COUNT_TYPE>& rs) {
		// 0 if not enough hits, 1 if all hits accounted for
		// not using an if statement to decrease branching (might not matter much)
		if (vps.total_hit_count > 0)
			return is.counted_hits() / vps.total_hit_count;
		return 1;
	}
};

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set, size_t i, size_t MAX_DEPTH>
struct pre_unroll_algorithm {
	static void place_ship(
			const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>&
					vps,
			iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
			vector<iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>>& is_todo) {
		if (MAX_DEPTH == is.current_depth) {
			is_todo.push_back(is);
			return;
		}
		size_t u;
		while ((u = is.get_next_placement()) != MAX_PLACEMENTS) {
			vps.increase_depth(is, u);
			pre_unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
													 i - 1, MAX_DEPTH>::place_ship(vps, is, is_todo);
			vps.decrease_depth(is);
		}
	}
};

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set, size_t MAX_DEPTH>
struct pre_unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
														0u, MAX_DEPTH> {
	static void place_ship(
			const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>&
					vps,
			iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
			vector<iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>>& is_todo) {
		cerr << "Should not have reached here (bottom of pre_unroll)" << endl;
	}
};

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set, size_t i, size_t MAX_DEPTH>
struct post_unroll_algorithm {
	static COUNT_TYPE
	place_ship(const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE,
																						fast_set>& vps,
						 iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
						 result_state<N, MAX_PLACEMENTS, COUNT_TYPE>& rs,
						 vector<result_state<N, MAX_PLACEMENTS, COUNT_TYPE>>& rs_acc,
						 size_t& rs_acc_iter, vector<COUNT_TYPE>& count_acc,
						 size_t& count_acc_iter) {
		if (MAX_DEPTH == is.current_depth) {
			rs.sum_from(rs_acc[rs_acc_iter++]);
			return count_acc[count_acc_iter++];
		}
		size_t result = 0;
		size_t u;
		while ((u = is.get_next_placement()) != MAX_PLACEMENTS) {
			vps.increase_depth(is, u);
			COUNT_TYPE current_result =
					post_unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE,
																fast_set, i - 1,
																MAX_DEPTH>::place_ship(vps, is, rs, rs_acc,
																											 rs_acc_iter, count_acc,
																											 count_acc_iter);
			vps.decrease_depth(is);
			rs.add_result(is.current_depth, u, current_result);
			result += current_result;
		}
		return result;
	}
};

template <size_t N, size_t MAX_PLACEMENTS, typename COUNT_TYPE,
					typename HIT_TYPE, typename fast_set, size_t MAX_DEPTH>
struct post_unroll_algorithm<N, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
														 0u, MAX_DEPTH> {
	static COUNT_TYPE
	place_ship(const valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE,
																						fast_set>& vps,
						 iteration_state<N, MAX_PLACEMENTS, HIT_TYPE, fast_set>& is,
						 result_state<N, MAX_PLACEMENTS, COUNT_TYPE>& rs,
						 vector<result_state<N, MAX_PLACEMENTS, COUNT_TYPE>>& rs_acc,
						 size_t& rs_acc_iter, vector<COUNT_TYPE>& count_acc,
						 size_t& count_acc_iter) {
		cerr << "Should not have reached here (bottom of post_unroll)" << endl;
		return 0;
	}
};
