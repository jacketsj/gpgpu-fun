#include <CL/sycl.hpp>

#include <algorithm>
#include <iostream>
#include <vector>
using std::cout;
using std::vector;
#include <bitset>

#include "algorithm.h"
// #include "fast_bitset.h"
#include "graph.h"
#include "standard_game.h"

// #define WIDTH 3
// #define HEIGHT 3
// #define n 2
// const array<size_t, n> lengths = {3, 2};
#define WIDTH 10
#define HEIGHT 10
#define n 5
const array<size_t, n> lengths = {5, 4, 3, 3, 2};

#define MAX_PLACEMENTS WIDTH*(HEIGHT - 1) + HEIGHT*(WIDTH - 1)
// #define fast_set fast_bitset::bitset<MAX_PLACEMENTS>
#define fast_set std::bitset<MAX_PLACEMENTS>
#define HIT_TYPE short
#define COUNT_TYPE long long

#define MAX_DEPTH 2

void process_batch(
		size_t offset_l, size_t offset_r,
		standard_game<n, WIDTH, HEIGHT, HIT_TYPE, fast_set>& sg,
		vector<iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set>>& is_todo,
		vector<result_state<n, MAX_PLACEMENTS, COUNT_TYPE>>& rs_acc,
		vector<COUNT_TYPE>& count_acc) {
	cl::sycl::default_selector device_selector;
	cl::sycl::queue queue(device_selector);
	// std::cout << "Running on " <<
	// queue.get_device().get_info<cl::sycl::info::device::name>() << endl;
	// getting the total number of compute units
	auto num_groups =
			queue.get_device().get_info<cl::sycl::info::device::max_compute_units>();
	// getting the maximum work group size per thread
	auto work_group_size =
			queue.get_device()
					.get_info<cl::sycl::info::device::max_work_group_size>();
	// building the best number of global thread
	auto total_threads = num_groups * work_group_size;

	// cout << "Total threads available via GPU: " << num_groups << '*' <<
	// work_group_size << '=' << total_threads << endl;
	cl::sycl::buffer<iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set>>
			is_todo_buffer(is_todo.data() + offset_l,
										 cl::sycl::range<1>(offset_r - offset_l));
	cl::sycl::buffer<result_state<n, MAX_PLACEMENTS, COUNT_TYPE>> rs_acc_buffer(
			rs_acc.data() + offset_l, cl::sycl::range<1>(offset_r - offset_l));
	cl::sycl::buffer<COUNT_TYPE> count_acc_buffer(
			count_acc.data() + offset_l, cl::sycl::range<1>(offset_r - offset_l));

	vector<valid_placement_subgraph<n, MAX_PLACEMENTS, HIT_TYPE, fast_set>>
			vps_vec({sg.vps});
	cl::sycl::buffer<
			valid_placement_subgraph<n, MAX_PLACEMENTS, HIT_TYPE, fast_set>>
			vps_buffer(vps_vec.data(), cl::sycl::range<1>(1));

	queue.submit([&](cl::sycl::handler& cgh) {
		auto is_todo_acc =
				is_todo_buffer.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto rs_acc_acc =
				rs_acc_buffer.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto count_acc_acc =
				count_acc_buffer.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto vps_acc =
				vps_buffer.get_access<cl::sycl::access::mode::read_write>(cgh);
		cgh.parallel_for<class gpu_place_ship>(
				cl::sycl::range<1>{offset_r - offset_l},
				[=](cl::sycl::item<1> item_id) {
					count_acc_acc[item_id] =
							unroll_algorithm<n, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE,
															 fast_set,
															 n - MAX_DEPTH>::place_ship(vps_acc[0],
																													is_todo_acc[item_id],
																													rs_acc_acc[item_id]);
				});
	});
}

int main() {
	standard_game<n, WIDTH, HEIGHT, HIT_TYPE, fast_set> sg(lengths);
	iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set> is;
	sg.vps.init_is(is);

	auto run = [=](iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set> is,
								 standard_game<n, WIDTH, HEIGHT, HIT_TYPE, fast_set> sg) {
		iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set> is_pre = is;
		vector<iteration_state<n, MAX_PLACEMENTS, HIT_TYPE, fast_set>> is_todo;
		pre_unroll_algorithm<n, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set, n,
												 MAX_DEPTH>::place_ship(sg.vps, is_pre, is_todo);
		vector<result_state<n, MAX_PLACEMENTS, COUNT_TYPE>> rs_acc(is_todo.size());
		vector<COUNT_TYPE> count_acc(is_todo.size());

		size_t BATCHES = 1;
		for (size_t batch = 0; batch < BATCHES; ++batch) {
			size_t batch_l = batch * is_todo.size() / BATCHES;
			size_t batch_r =
					std::min((batch + 1) * is_todo.size() / BATCHES, is_todo.size());
			// printf("Batch %zu: [%zu to %zu)\n", batch, batch_l, batch_r);
			process_batch(batch_l, batch_r, sg, is_todo, rs_acc, count_acc);
		}

		result_state<n, MAX_PLACEMENTS, COUNT_TYPE> rs;
		size_t rs_acc_iter = 0, count_acc_iter = 0;
		COUNT_TYPE total_configurations =
				post_unroll_algorithm<n, MAX_PLACEMENTS, COUNT_TYPE, HIT_TYPE, fast_set,
															n, MAX_DEPTH>::place_ship(sg.vps, is, rs, rs_acc,
																												rs_acc_iter, count_acc,
																												count_acc_iter);

		// cout << "Total configurations: " << total_configurations << '\n';
		sg.print_output(rs, total_configurations, true, true);
	};
	array<array<char, HEIGHT>, WIDTH> board;
	for (size_t y = 0; y < HEIGHT; ++y)
		for (size_t x = 0; x < WIDTH; ++x) {
			char c;
			do {
				std::cin >> c;
			} while (c != 'o' && c != 'x' && c != '.' &&
							 (c < '0' || c > char(n + '0')));
			board[x][y] = c;
		}
	for (size_t x = 0; x < WIDTH; ++x)
		for (size_t y = 0; y < HEIGHT; ++y) {
			if (board[x][y] == 'o')
				sg.update(query_hit(square(x, y)));
			else if (board[x][y] == 'x')
				sg.update(query_miss(square(x, y)));
		}
	for (size_t x = 0; x < WIDTH; ++x)
		for (size_t y = 0; y < HEIGHT; ++y) {
			if (board[x][y] >= '0' && board[x][y] <= '9')
				sg.update(query_sink(square(x, y), size_t(board[x][y] - '0')));
		}
	run(is, sg);
}
