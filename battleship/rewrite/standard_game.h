#pragma once

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
using std::array;
using std::cout;
using std::fixed;
using std::setprecision;
using std::vector;

#include "graph.h"

struct square {
	size_t x, y;
	square(size_t _x, size_t _y) : x(_x), y(_y) {}
	bool operator==(const square& oth) const { return x == oth.x && y == oth.y; }
};

struct query_miss {
	square loc;
	query_miss(const square& _loc) : loc(_loc) {}
};

struct query_hit {
	square loc;
	query_hit(const square& _loc) : loc(_loc) {}
};

struct query_sink {
	square loc;
	size_t ship;
	query_sink(const square& _loc, const size_t& _ship)
			: loc(_loc), ship(_ship) {}
};

// The standard game is one that involves ships of 1xk size on a WIDTH * HEIGHT
// sized board. k may vary per-ship, but k>=2 is required.
template <size_t N, size_t WIDTH, size_t HEIGHT, typename HIT_TYPE,
					typename fast_set>
struct standard_game {
	array<size_t, N> lengths;
	static constexpr size_t MAX_PLACEMENTS =
			WIDTH * (HEIGHT - 1) + HEIGHT * (WIDTH - 1);
	valid_placement_subgraph<N, MAX_PLACEMENTS, HIT_TYPE, fast_set> vps;
	vector<square> hit_squares;
	// mark a square as hit
	void mark_square_hit(const square& s) { hit_squares.push_back(s); }
	bool is_square_hit(const square& s0) {
		for (const auto& s : hit_squares)
			if (s == s0)
				return true;
		return false;
	}
	// get the number of placements that exist for a ship p
	size_t get_num_placements(size_t p) {
		return (WIDTH - lengths[p] + 1) * HEIGHT +
					 (HEIGHT - lengths[p] + 1) * WIDTH;
	}
	// get the squares used by a ship p at placement u
	vector<square> get_used_squares(size_t p, size_t u) {
		// define the canonical ordering of placements to be left to right,
		// top to bottom, horizontal to vertical
		// each starts at their top-left coordinate
		size_t x, y, dx, dy;
		size_t horizontal_placements_per_row = WIDTH - lengths[p] + 1;
		size_t total_horizontal_placements = horizontal_placements_per_row * HEIGHT;
		if (u < total_horizontal_placements) {
			x = u % horizontal_placements_per_row;
			y = u / horizontal_placements_per_row;
			dx = 1;
			dy = 0;
		} else {
			x = (u - total_horizontal_placements) % WIDTH;
			y = (u - total_horizontal_placements) / WIDTH;
			dx = 0;
			dy = 1;
		}
		vector<square> ret;
		for (size_t l = 0; l < lengths[p]; ++l) {
			ret.emplace_back(x, y);
			x += dx;
			y += dy;
		}
		return ret;
	}
	// check if a ship p at placement u uses a square loc
	bool uses_square(size_t p, size_t u, const square& loc) {
		for (auto& s : get_used_squares(p, u))
			if (s == loc)
				return true;
		return false;
	}
	// returns a list of placements of ship p that overlap square loc
	vector<size_t> get_overlapping_placements(size_t p, const square& loc) {
		vector<size_t> ret;
		for (size_t u = 0; u < get_num_placements(p); ++u) {
			if (uses_square(p, u, loc))
				ret.push_back(u);
		}
		return ret;
	}
	// returns a list of placements of ship p that do not overlap square loc
	vector<size_t> get_nonoverlapping_placements(size_t p, const square& loc) {
		vector<size_t> ret;
		for (size_t u = 0; u < get_num_placements(p); ++u) {
			if (!uses_square(p, u, loc))
				ret.push_back(u);
		}
		return ret;
	}
	// does (ship p, placement u) overlap (ship q, placement v)
	bool overlap(size_t p, size_t u, size_t q, size_t v) {
		for (const auto& s : get_used_squares(q, v))
			if (uses_square(p, u, s))
				return true;
		return false;
	}
	// returns a list of placements v of ship q that overlap (ship p, placement u)
	vector<size_t> get_overlapping_placements(size_t p, size_t u, size_t q) {
		vector<size_t> ret;
		for (size_t v = 0; v < get_num_placements(q); ++v) {
			if (overlap(p, u, q, v))
				ret.push_back(v);
		}
		return ret;
	}
	placement_graph<N, MAX_PLACEMENTS> build_ps() {
		placement_graph<N, MAX_PLACEMENTS> pg;
		for (size_t p = 0; p < N; ++p)
			pg.set_placement_count(p, get_num_placements(p));
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < get_num_placements(p); ++u)
				for (size_t q = 0; q < N; ++q)
					for (size_t v : get_overlapping_placements(p, u, q))
						pg.remove_edge(p, u, q, v);
		return pg;
	}
	standard_game(const array<size_t, N>& _lengths)
			: lengths(_lengths), vps(build_ps()) {}
	void update(query_miss qr) {
		for (size_t p = 0; p < N; ++p)
			for (size_t u : get_overlapping_placements(p, qr.loc))
				vps.remove_vertex(p, u);
	}
	void update(query_hit qr) {
		for (size_t p = 0; p < N; ++p)
			for (size_t u : get_overlapping_placements(p, qr.loc))
				vps.increment_hit_count(p, u);
		vps.increment_total_hit_count();
		mark_square_hit(qr.loc);
	}
	void update(query_sink qr) {
		update(query_hit(qr.loc));
		for (size_t p = 0; p < N; ++p)
			if (p != qr.ship)
				for (size_t u : get_overlapping_placements(p, qr.loc))
					vps.remove_vertex(p, u);
		for (size_t u : get_nonoverlapping_placements(qr.ship, qr.loc))
			vps.remove_vertex(qr.ship, u);
		/*
		// Commented out due to different game behaviour in most online variants
		for (size_t u : get_overlapping_placements(qr.ship, qr.loc))
			for (const auto& s : get_used_squares(qr.ship, u))
				if (!is_square_hit(s)) {
					vps.remove_vertex(qr.ship, u);
					break;
				}
		*/
	}
	template <typename COUNT_TYPE>
	void print_output(const result_state<N, MAX_PLACEMENTS, COUNT_TYPE>& rs,
										COUNT_TYPE total_configurations, bool minimal = false,
										bool col = false) {
		array<array<COUNT_TYPE, HEIGHT>, WIDTH> square_counts;
		for (size_t y = 0; y < HEIGHT; ++y)
			for (size_t x = 0; x < WIDTH; ++x)
				square_counts[x][y] = 0;
		for (size_t p = 0; p < N; ++p)
			for (size_t u = 0; u < get_num_placements(p); ++u) {
				for (const auto& s : get_used_squares(p, u))
					square_counts[s.x][s.y] += rs.appearance_counts[p][u];
			}
		if (!minimal)
			for (size_t y = 0; y < HEIGHT; ++y) {
				for (size_t x = 0; x < WIDTH; ++x)
					cout << '\t' << square_counts[x][y];
				cout << '\n';
			}
		if (!minimal)
			cout << fixed << setprecision(3); // show 3 decimal points
		else
			cout << fixed << setprecision(0); // show 0 decimal points
		cout << "   ";
		for (size_t x = 0; x < WIDTH; ++x)
			if (x + 1 < 10)
				cout << "   " << x + 1;
		for (size_t x = 0; x < WIDTH; ++x)
			if (x + 1 >= 10)
				cout << "  " << x + 1;
		cout << '\n';
		cout << "   ";
		for (size_t x = 0; x < WIDTH; ++x)
			cout << "----";
		cout << '\n';
		size_t highest_square_count = 0;
		// TODO: instead of checking if < total_configurations, check if not hit
		// (since sometimes 100% is possible)
		for (size_t y = 0; y < HEIGHT; ++y)
			for (size_t x = 0; x < WIDTH; ++x)
				if (square_counts[x][y] > highest_square_count &&
						square_counts[x][y] < total_configurations)
					highest_square_count = square_counts[x][y];
		for (size_t y = 0; y < HEIGHT; ++y) {
			if (y + 1 < 10)
				cout << y + 1 << " | ";
			else
				cout << y + 1 << "| ";
			for (size_t x = 0; x < WIDTH; ++x) {
				if (!minimal)
					cout << '\t' << double(square_counts[x][y]) / total_configurations;
				else {
					using std::string;
					bool is_highest = highest_square_count == square_counts[x][y];
					int percent = int(
							100 * double(square_counts[x][y]) / total_configurations + 0.5f);
					string percent_str = std::to_string(percent) + " ";
					if (percent < 10)
						percent_str = string(" ") + percent_str;
					if (percent < 100)
						percent_str = string(" ") + percent_str;
					if (col) {
						string col;
						if (percent < 10)
							col = "34";
						else if (percent < 20)
							col = "36";
						else if (percent < 30)
							col = "35";
						else if (percent < 60)
							col = "33";
						else if (percent < 80)
							col = "91";
						else
							col = "31";
						if (is_highest)
							col += ";4";
						cout << string("\x1B[") + col + "m" + percent_str + "\033[0m";
					} else
						cout << percent_str;
				}
			}
			cout << '\n';
		}
	}
};
