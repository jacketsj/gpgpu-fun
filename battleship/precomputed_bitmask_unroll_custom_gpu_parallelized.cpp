#include <CL/sycl.hpp>
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

// CONSTANTS
// board dimensions
#define WIDTH 5
#define HEIGHT 5
/*
#define WIDTH 10
#define HEIGHT 10
//boats
const vector<int> lengths = {2,3,3,4,5};
*/
// const vector<int> lengths = {5,4,3,3,2};
const vector<int> lengths = {5, 4, 3, 3, 2};
#define n 5

typedef unsigned long long ull;
// typedef uint64_t ull;

namespace fast_bitset {
template <ll T> struct bitset {
#define sz ll(sizeof(ull) * 8)
#define N (T + sz - 1) / sz // ceiling of T/sz
	ull vals[N];
	bitset() {
		for (int i = 0; i < N; ++i)
			vals[i] = 0;
	}
	void operator&=(const bitset& other) {
		for (int i = 0; i < N; ++i)
			vals[i] &= other.vals[i];
	}
	// any bit is set
	bool any() const {
		for (int i = 0; i < N; ++i)
			if (vals[i] != 0)
				return true;
		return false;
	}
	bool get(ll i) const {
		int k = i / sz;
		i %= sz;
		return (vals[k] & (ll(1) << i)) != 0;
	}
	void set(ll i) {
		int k = i / sz;
		i %= sz;
		vals[k] |= (ll(1) << i);
	}
	void reset(ll i) {
		int k = i / sz;
		i %= sz;
		vals[k] &= (~(ll(1) << i));
	}
	void assign(int i, bool b) {
		if (b)
			set(i);
		else
			reset(i);
	}
	int bitscan_destructive_any() {
		for (int i = 0; i < N; ++i) {
			if (vals[i] != 0) {
				ull t = vals[i] & -vals[i];
				int j = __builtin_ctzll(vals[i]);
				vals[i] ^= t;
				return sz * i + j;
			}
		}
		return -1;
	}
	// do a destructive bitscan
	// returns -1 if nothing found
	int bitscan_destructive() {
		for (int i = 0; i < N; ++i) {
			if (vals[i] != 0) {
				// compute the index
				// int ret = sz*i+(__builtin_ffsll(vals[i])-1);
				// __builtin_ctzll counts trailing 0s
				ull t = vals[i] & -vals[i]; // lowestOneBit
				int j = __builtin_ctzll(vals[i]);
				// destroy it
				vals[i] ^= t;
				return sz * i + j;
				// vals[i] &= ~(1<<j);
				// the fenwick tree code doesn't work (?)
			}
		}
	}
	int ctz(ll v) const {
		for (int i = 0; i < sz; ++i)
			if ((v & (ll(1) << i)) != 0)
				return i;
	}
	// non-destructive bitscan
	int bitscan() {
		for (int i = 0; i < N; ++i)
			if (vals[i] != 0)
				// return sz*i+ctz(vals[i]);
				return sz * i + __builtin_ctzll(vals[i]);
	}
#undef sz
#undef N
};
} // namespace fast_bitset

typedef fast_bitset::bitset<WIDTH*(HEIGHT - 1) + HEIGHT*(WIDTH - 1)> pos_set;

typedef vector<ll> grid_t;

grid_t create_grid() { return grid_t(WIDTH * HEIGHT); }

// produce the number of unique horizontal or vertical positions that a
// top-left coordinate can take on in the corresponding orientation
int horizontal_state_class_count(int length) { return (WIDTH - length + 1); }
int vertical_state_class_count(int length) { return (HEIGHT - length + 1); }
// produce number of horizontal states for a single ship on a clean board
int horizontal_state_count(int length) {
	return horizontal_state_class_count(length) * HEIGHT;
}
int vertical_state_count(int length) {
	return vertical_state_class_count(length) * WIDTH;
}
// produce the number of states for a ship on a clean board
int state_count(int length) {
	return horizontal_state_count(length) + vertical_state_count(length);
}

// (x,y,dx,dy)
typedef tuple<int, int, int, int> dpos_t;

// given a state index and a ship length, produce (x,y,dx,dy)
dpos_t get_state_dpos(int state_index, int length) {
	// define the canonical ordering to be left to right,
	// top to bottom, horizontal to vertical
	// all are defined by their top-left coordinate

	// number of states in a given row for the horizontal orientation
	int horizontal_width = horizontal_state_class_count(length);

	// number of states with horizontal orientation
	int num_in_horizontal = horizontal_width * HEIGHT;
	if (state_index < num_in_horizontal) {
		// horizontal orientation
		return make_tuple(state_index % horizontal_width,
											state_index / horizontal_width, 1, 0);
	}
	// vertical orientation
	state_index %= num_in_horizontal;
	return make_tuple(state_index % WIDTH, state_index / WIDTH, 0, 1);
}

// given a state index and a ship length, draw it onto a grid with a multiplier,
// and return false if a conflict is found, but only after finishing the drawing
bool draw_state(int state_index, int length, ll mult, grid_t& grid) {
	int x, y, dx, dy;
	tie(x, y, dx, dy) = get_state_dpos(state_index, length);
	bool ret = true;
	for (int l = 0; l < length; ++l) {
		grid[x + y * WIDTH] += mult;

		if (grid[x + y * WIDTH] >
				1) // this grid cell is covered by two ships, so the drawing is invalid
			ret = false;

		x += dx;
		y += dy;
	}
	return ret;
}

// returns a list of state indeces that do not conflict with misses
vector<int> find_valid_individual_states(int length, grid_t& misses) {
	vector<int> res; // final result
	// check validity of each state index
	for (int state_index = 0; state_index < state_count(length); ++state_index) {
		// draw and find validity
		if (draw_state(state_index, length, 1, misses))
			res.push_back(state_index);
		// undraw
		draw_state(state_index, length, -1, misses);
	}
	return res;
}

// returns a list of position indeces that do not conflict with misses for each
// ship
vector<vector<int>> find_valid_states(grid_t& misses) {
	vector<vector<int>> res(n);
	for (int i = 0; i < n; ++i)
		res[i] = find_valid_individual_states(lengths[i], misses);
	return res;
}

bool compatible_pair(int length_a, int state_a, int length_b, int state_b,
										 grid_t& test_grid) {
	bool valid = draw_state(state_a, length_a, 1, test_grid);
	valid = draw_state(state_b, length_b, 1, test_grid) && valid;
	draw_state(state_a, length_a, -1, test_grid);
	draw_state(state_b, length_b, -1, test_grid);
	return valid;
}

pos_set all_compatible_pairs(int length_a, int state_a, int length_b,
														 vector<int>& states_b, grid_t& test_grid) {
	pos_set res;
	for (size_t i = 0; i < states_b.size(); ++i) {
		int state_b = states_b[i];
		if (compatible_pair(length_a, state_a, length_b, state_b, test_grid))
			res.set(i);
	}
	return res;
}

vector<vector<vector<pos_set>>>
find_all_pos_sets(vector<vector<int>>& valid_states) {
	// create a single test grid to avoid reallocating memory needlessly
	grid_t test_grid = create_grid();
	// indexed by ship1, position of ship1, ship2
	vector<vector<vector<pos_set>>> res(n);
	for (int i = 0; i < n; ++i) {
		res[i].resize(valid_states[i].size());
		for (size_t state_i = 0; state_i < valid_states[i].size();
				 ++state_i) // for each state of ship i
			res[i][state_i].resize(n);
	}

	for (int i = 0; i < n; ++i) // for each ship
		for (size_t state_i = 0; state_i < valid_states[i].size();
				 ++state_i) // for each state of ship i
			for (int j = i + 1; j < n;
					 ++j) // for each ship that will be seen later than i
				res[i][state_i][j] =
						all_compatible_pairs(lengths[i], valid_states[i][state_i],
																 lengths[j], valid_states[j], test_grid);
	return res;
}

typedef vector<pos_set> place_ship_params;

template <unsigned n_minus_ship_index, unsigned depth> struct unroll_pre {
	static void
	place_ship_pre(const vector<vector<vector<pos_set>>>& validity_masks,
								 vector<pos_set>& currently_valid,
								 vector<place_ship_params>& params_list,
								 const vector<int>& num_valid_states) {
		int ship_index = n - n_minus_ship_index;

		// Save the information about currently_valid that we will need when we
		// return from recursive calls
		vector<pos_set> edits(n - ship_index - 1);
		for (int j = ship_index + 1; j < n; ++j)
			edits[j - ship_index - 1] = currently_valid[j];

		// count of the number of valid placements of all remaining ships
		ll count = 0;

		// iterate over all the ship states that are still valid
		pos_set current = currently_valid[ship_index];
		// current.copy(currently_valid[ship_index]);
		// while (current.any())
		int state_index;
		while ((state_index = current.bitscan_destructive_any()) != -1) {
			// int state_index = current.bitscan_destructive();
			// update legal states for remaining ships
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] &= validity_masks[ship_index][state_index][j];
			// recurse on remaining ships
			unroll_pre<n_minus_ship_index - 1, depth - 1>::place_ship_pre(
					validity_masks, currently_valid, params_list, num_valid_states);
			// count += sub_result;
			// state_frequency[ship_index][state_index] += sub_result;
			// set currently_valid values back
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] = edits[j - ship_index - 1];
		}
		return;
	}
};

template <unsigned n_minus_ship_index, unsigned depth> struct unroll_post {
	static ll
	place_ship_post(const vector<vector<vector<pos_set>>>& validity_masks,
									vector<pos_set>& currently_valid,
									vector<vector<int>>& state_frequency,
									const vector<int>& num_valid_states,
									const vector<ll>& results, size_t& param_index) {
		int ship_index = n - n_minus_ship_index;

		// Save the information about currently_valid that we will need when we
		// return from recursive calls
		vector<pos_set> edits(n - ship_index - 1);
		for (int j = ship_index + 1; j < n; ++j)
			edits[j - ship_index - 1] = currently_valid[j];

		// count of the number of valid placements of all remaining ships
		ll count = 0;

		// iterate over all the ship states that are still valid
		pos_set current = currently_valid[ship_index];
		// current.copy(currently_valid[ship_index]);
		// while (current.any())
		int state_index;
		while ((state_index = current.bitscan_destructive_any()) != -1) {
			// int state_index = current.bitscan_destructive();
			// update legal states for remaining ships
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] &= validity_masks[ship_index][state_index][j];
			// recurse on remaining ships
			ll sub_result =
					unroll_post<n_minus_ship_index - 1, depth - 1>::place_ship_post(
							validity_masks, currently_valid, state_frequency,
							num_valid_states, results, param_index);
			// record counts
			count += sub_result;
			state_frequency[ship_index][state_index] += sub_result;
			// set currently_valid values back
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] = edits[j - ship_index - 1];
		}
		return count;
	}
	static ll place_ship_post(
			const vector<vector<vector<pos_set>>>& validity_masks,
			vector<pos_set>& currently_valid, vector<vector<int>>& state_frequency,
			const vector<int>& num_valid_states, const vector<ll>& results) {
		size_t param_index = 0;
		return place_ship_post(validity_masks, currently_valid, state_frequency,
													 num_valid_states, results, param_index);
	}
};

template <unsigned n_minus_ship_index>
struct unroll_post<n_minus_ship_index, 0u> {
	static ll
	place_ship_post(const vector<vector<vector<pos_set>>>& validity_masks,
									vector<pos_set>& currently_valid,
									vector<vector<int>>& state_frequency,
									const vector<int>& num_valid_states,
									const vector<ll>& results, size_t& param_index) {
		return results[param_index++];
	}
	static ll place_ship_post(
			const vector<vector<vector<pos_set>>>& validity_masks,
			vector<pos_set>& currently_valid, vector<vector<int>>& state_frequency,
			const vector<int>& num_valid_states, const vector<ll>& results) {
		return results[0];
	}
};

template <unsigned n_minus_ship_index>
struct unroll_pre<n_minus_ship_index, 0u> {
	static void
	place_ship_pre(const vector<vector<vector<pos_set>>>& validity_masks,
								 vector<pos_set>& currently_valid,
								 vector<place_ship_params>& params_list,
								 const vector<int>& num_valid_states) {
		params_list.push_back(currently_valid);
	}
};

template <unsigned n_minus_ship_index, typename validity_masks_t,
					typename validity_masks_offsets_t, typename currently_valid_t,
					typename state_frequency_t, typename num_valid_states_t>
struct unroll_gpu {
	// state offset is starting index of states for current ship in
	// state_frequency
	// edit_offset is starting index of edits to store in edit array
	static ll place_ship_gpu(validity_masks_t validity_masks,
													 validity_masks_offsets_t validity_masks_offsets,
													 currently_valid_t currently_valid,
													 state_frequency_t state_frequency,
													 num_valid_states_t num_valid_states,
													 size_t state_offset) {
		int ship_index = n - n_minus_ship_index; // const int?

		// Save the information about currently_valid that we will need when we
		// return from recursive calls

		// count of the number of valid placements of all remaining ships
		ll count = 0;

		pos_set edits[n - ship_index - 1];
		for (int j = ship_index + 1; j < n; ++j)
			edits[j - ship_index - 1] = currently_valid[j];

		// iterate over all the ship states that are still valid
		pos_set current = currently_valid[ship_index];
		int state_index;
		while ((state_index = current.bitscan_destructive_any()) != -1) {
			// int state_index = current.bitscan_destructive();
			// update legal states for remaining ships
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] &=
						validity_masks[validity_masks_offsets[ship_index] +
													 state_index * n + j]; // TODO is this right?
			// recurse on remaining ships
			ll sub_result = unroll_gpu<
					n_minus_ship_index - 1, validity_masks_t, validity_masks_offsets_t,
					currently_valid_t, state_frequency_t,
					num_valid_states_t>::place_ship_gpu(validity_masks,
																							validity_masks_offsets,
																							currently_valid, state_frequency,
																							num_valid_states,
																							state_offset +
																									num_valid_states[ship_index]);
			// record counts
			count += sub_result;
			state_frequency[state_offset + state_index] += sub_result;
			// set currently_valid values back
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] = edits[j - ship_index - 1];
		}
		return count;
	}
};

template <typename validity_masks_t, typename validity_masks_offsets_t,
					typename currently_valid_t, typename state_frequency_t,
					typename num_valid_states_t>
struct unroll_gpu<0u, validity_masks_t, validity_masks_offsets_t,
									currently_valid_t, state_frequency_t, num_valid_states_t> {
	static ll place_ship_gpu(validity_masks_t validity_masks,
													 validity_masks_offsets_t validity_masks_offsets,
													 currently_valid_t currently_valid,
													 state_frequency_t state_frequency,
													 num_valid_states_t num_valid_states,
													 size_t state_offset) {
		// base case: all ships placed successfully
		return 1;
	}
};

template <unsigned n_minus_ship_index> struct unroll {
	static ll place_ship(const vector<vector<vector<pos_set>>>& validity_masks,
											 vector<pos_set>& currently_valid,
											 vector<vector<int>>& state_frequency,
											 const vector<int>& num_valid_states)
	//, ll &total_successful)
	{
		int ship_index = n - n_minus_ship_index;

		// Save the information about currently_valid that we will need when we
		// return from recursive calls
		vector<pos_set> edits(n - ship_index - 1);
		for (int j = ship_index + 1; j < n; ++j)
			edits[j - ship_index - 1] = currently_valid[j];

		// count of the number of valid placements of all remaining ships
		ll count = 0;

		// iterate over all the ship states that are still valid
		pos_set current = currently_valid[ship_index];
		// current.copy(currently_valid[ship_index]);
		// while (current.any())
		int state_index;
		while ((state_index = current.bitscan_destructive_any()) != -1) {
			// int state_index = current.bitscan_destructive();
			// update legal states for remaining ships
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] &= validity_masks[ship_index][state_index][j];
			// recurse on remaining ships
			ll sub_result = unroll<n_minus_ship_index - 1>::place_ship(
					validity_masks, currently_valid, state_frequency, num_valid_states);
			// record counts
			count += sub_result;
			state_frequency[ship_index][state_index] += sub_result;
			// set currently_valid values back
			for (int j = ship_index + 1; j < n; ++j)
				currently_valid[j] = edits[j - ship_index - 1];
		}
		return count;
	}
};

template <> struct unroll<0u> {
	static ll place_ship(const vector<vector<vector<pos_set>>>& validity_masks,
											 vector<pos_set>& currently_valid,
											 vector<vector<int>>& state_frequency,
											 const vector<int>& num_valid_states)
	//, ll &total_successful)
	{
		// base case: all ships placed successfully
		// if (++total_successful % 100000000 == 0)
		//	cout << "Cumulative successful: " << total_successful << endl;
		return 1;
	}
};

// print the final grid
void print_grid(const grid_t& a) {
	for (int y = 0; y < HEIGHT; ++y)
		for (int x = 0; x < WIDTH; ++x)
			cout << a[x + y * WIDTH] << "\t\n"[x == WIDTH - 1];
}

// print the final grid of probabilities
void print_grid_chance(const grid_t& a, const ll& added_count) {
	cout << fixed << setprecision(3); // show 3 decimals for chance
	double count = added_count;
	for (int y = 0; y < HEIGHT; ++y)
		for (int x = 0; x < WIDTH; ++x)
			cout << double(a[x + y * WIDTH]) / count << "\t\n"[x == WIDTH - 1];
}

// list of hits not currently supported, need to derive new algorithm
// count the number of occurrences of each spot on the grid
// do this by counting the frequency of each position for each ship
void count_occurrences(grid_t& misses) {
	vector<vector<int>> valid_states = find_valid_states(misses);
	// count the number of valid states, to be used when placing ships
	vector<int> num_valid_states(n);
	for (int i = 0; i < n; ++i)
		num_valid_states[i] = valid_states[i].size();

	// validity masks: given a ship, a state, and a future ship, provides the mask
	// for that ship
	vector<vector<vector<pos_set>>> validity_masks =
			find_all_pos_sets(valid_states);

	vector<pos_set> currently_valid(
			n); // for each ship, keep track of which states are still valid
	for (int i = 0; i < n; ++i)
		for (int state_i = 0; state_i < num_valid_states[i]; ++state_i)
			currently_valid[i].set(state_i);

	vector<vector<int>> state_frequency(
			n); // for each ship, keep the frequency of each of its states
	for (int i = 0; i < n; ++i) {
		state_frequency[i].resize(num_valid_states[i]);
		for (auto& v :
				 state_frequency[i]) // TODO does this help get rid of garbage data??
			v = 0;
	}

	// ll total_successful = 0;
	// call recursive ship placement routine to iterate through all valid
	// placements

	// procedure plan:
	// call unroll_pre::place_ship_pre, which creates a list of
	// 'currently_valid' params for unroll::place_ship parallel_for on
	// unroll::place_ship starting at ship index of recursion level + 1, store
	// each result in a results vector (and also store stuff in state_frequency)
	// call unroll_post::place_ship_post with the results vector

	// state_frequency is a 2d non-uniform vector by default (intentionally,
	// don't want false sharing) need to unroll it before sending to gpu, then
	// unroll it again
	// can be reversed by using num_valid_states vector
	vector<int> state_frequency_unrolled;
	for (auto& subvec : state_frequency)
		for (auto& elem : subvec)
			state_frequency_unrolled.push_back(elem);

	// currently_valid_sets is the result of a few levels of recursion
	constexpr int PRE_DEPTH = 1;
	vector<place_ship_params> currently_valid_sets;
	vector<pos_set> currently_valid_dupe = currently_valid;
	unroll_pre<n, PRE_DEPTH>::place_ship_pre(validity_masks, currently_valid_dupe,
																					 currently_valid_sets,
																					 num_valid_states);
	vector<pos_set> currently_valid_sets_unrolled;
	for (auto& v : currently_valid_sets)
		for (auto& ps : v)
			currently_valid_sets_unrolled.push_back(ps);

	// done TODO: turn unrolled back into not unrolled currently_valid_sets after
	// completion
	// done TODO: use parallel_for on each contiguous subarray of
	// currently_valid_sets_unrolled modulo n
	// done TODO: define PRE_DEPTH, and modify unroll<n> to be unroll<n-PRE_DEPTH>
	// done TODO: define subproblem_results for each entry in currently_valid_sets
	vector<ll> subproblem_results(currently_valid_sets.size(), 0);

	// TODO: subproblem_results.size() copies of state_frequency_unrolled
	vector<ll> state_frequency_unrolled_arr(
			state_frequency_unrolled.size() * subproblem_results.size(), 0);

	ll total_states = 0;

	{
		cl::sycl::default_selector device_selector;
		cl::sycl::queue queue(device_selector);
		std::cout << "Running on "
							<< queue.get_device().get_info<cl::sycl::info::device::name>()
							<< endl;

		// subproblem_results
		cl::sycl::buffer<long long, 1> subproblem_results_sycl(
				subproblem_results.data(),
				cl::sycl::range<1>(subproblem_results.size()));
		queue.submit([&](cl::sycl::handler& cgh) {
			auto subproblem_results_acc =
					subproblem_results_sycl
							.get_access<cl::sycl::access::mode::read_write>(cgh);
			// validity masks
			vector<pos_set> validity_masks_unrolled;
			vector<int> validity_masks_offsets; // where each first ship index starts
			for (auto& sub1 : validity_masks) {
				validity_masks_offsets.push_back(validity_masks_unrolled.size());
				for (auto& sub2 : sub1)
					for (auto& sub3 : sub2)
						validity_masks_unrolled.push_back(sub3);
			}
			cl::sycl::buffer<pos_set> validity_masks_unrolled_sycl(
					validity_masks_unrolled.data(),
					cl::sycl::range<1>(validity_masks_unrolled.size()));
			cl::sycl::buffer<int> validity_masks_offsets_sycl(
					validity_masks_offsets.data(),
					cl::sycl::range<1>(validity_masks_offsets.size()));
			auto validity_masks_unrolled_acc =
					validity_masks_unrolled_sycl.get_access<cl::sycl::access::mode::read>(
							cgh);
			auto validity_masks_offsets_acc =
					validity_masks_offsets_sycl.get_access<cl::sycl::access::mode::read>(
							cgh);
			// currently_valid sets
			cl::sycl::buffer<pos_set> currently_valid_sets_sycl(
					currently_valid_sets_unrolled.data(),
					cl::sycl::range<1>(currently_valid_sets_unrolled.size()));
			auto currently_valid_sets_acc =
					currently_valid_sets_sycl
							.get_access<cl::sycl::access::mode::read_write>(cgh);
			// state_frequency_unrolled copies
			cl::sycl::buffer<ll> state_frequency_unrolled_arr_sycl(
					state_frequency_unrolled_arr.data(),
					cl::sycl::range<1>(state_frequency_unrolled_arr.size()));
			auto state_frequency_unrolled_arr_acc =
					state_frequency_unrolled_arr_sycl
							.get_access<cl::sycl::access::mode::read_write>(cgh);
			/*
			cl::sycl::buffer<ll> state_frequency_unrolled_sycl(
					state_frequency_unrolled.data(),
					cl::sycl::range<1>(state_frequency_unrolled.size()));
			auto state_frequency_unrolled_acc =
					state_frequency_unrolled_sycl
							.get_access<cl::sycl::access::mode::read_write>(cgh);
			*/
			// total_states
			cl::sycl::buffer<long long, 1> total_states_sycl(&total_states,
																											 cl::sycl::range<1>(1));
			auto total_states_acc =
					total_states_sycl.get_access<cl::sycl::access::mode::read_write>(cgh);
			// num_valid_states
			cl::sycl::buffer<int> num_valid_states_sycl(
					num_valid_states.data(), cl::sycl::range<1>(num_valid_states.size()));
			auto num_valid_states_acc =
					num_valid_states_sycl.get_access<cl::sycl::access::mode::read>(cgh);
			// TODO: convert this to parallel_for
			cgh.parallel_for<class gpu_place_ship>(
					cl::sycl::range<1>{subproblem_results.size()},
					[=](cl::sycl::id<1> item_id) {
						// TODO: fix false sharing with currently_valid_sets_acc
						pos_set currently_valid[n];
						for (int i = 0; i < n; ++i)
							currently_valid[i] = currently_valid_sets_acc[item_id * n + i];
						auto currently_valid_offset_acc =
								&currently_valid_sets_acc[item_id * n];
						auto state_frequency_unrolled_acc =
								&state_frequency_unrolled_arr_acc[item_id];
						subproblem_results_acc[item_id] =
								unroll_gpu<n - PRE_DEPTH, decltype(validity_masks_unrolled_acc),
													 decltype(validity_masks_offsets_acc),
													 decltype(currently_valid_offset_acc),
													 decltype(state_frequency_unrolled_acc),
													 decltype(num_valid_states_acc)>::
										place_ship_gpu(
												validity_masks_unrolled_acc, validity_masks_offsets_acc,
												currently_valid_offset_acc,
												state_frequency_unrolled_acc, num_valid_states_acc, 0);
					});
		});
	}

	// since buffers (and the queue?) are destroyed, queue should wait until all
	// operations using them have finished before continuing

	// put currently_valid_sets_unrolled back into currently_valid_sets
	int cvs_i = 0;
	for (auto& s : currently_valid_sets)
		for (auto& v : s)
			v = currently_valid_sets_unrolled[cvs_i++];
	// run place_ship_post
	auto currently_valid_dupe_post = currently_valid;
	total_states = unroll_post<n, PRE_DEPTH>::place_ship_post(
			validity_masks, currently_valid_dupe_post, state_frequency,
			num_valid_states, subproblem_results);

	// reverse unrolling of state_frequency now that gpu computation is done
	int sfu_i = 0;
	for (auto& copy : state_frequency_unrolled_arr)
		state_frequency_unrolled[sfu_i++] += copy;
	int sf_i = 0;
	for (auto& subvec : state_frequency)
		for (auto& elem : subvec)
			elem += state_frequency_unrolled[sf_i++];

	grid_t frequencies = create_grid();
	for (int i = 0; i < n; ++i)
		for (int state_index = 0; state_index < num_valid_states[i]; ++state_index)
			draw_state(state_index, lengths[i], state_frequency[i][state_index],
								 frequencies);

	cout << "Total states: " << total_states << endl;
	// cout << "Total successful (this should be the same number): " <<
	// total_successful << endl;
	print_grid(frequencies);
	print_grid_chance(frequencies, total_states);
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(NULL);

	// no misses for now, empty grid
	grid_t misses = create_grid();
	count_occurrences(misses);
}
