#include <iostream>
using std::cout;

#include "algorithm.h"
#include "graph.h"
#include "standard_game.h"

typedef long long ll;
typedef unsigned long long ull;
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
	friend bitset operator&(const bitset& b1, const bitset& b2) {
		bitset ret;
		for (int i = 0; i < N; ++i)
			ret.vals[i] = b1.vals[i] & b2.vals[i];
		return ret;
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
	size_t bitscan_destructive_any(size_t ERROR) {
		for (size_t i = 0; i < N; ++i) {
			if (vals[i] != 0) {
				ull t = vals[i] & -vals[i];
				size_t j = __builtin_ctzll(vals[i]);
				vals[i] ^= t;
				return sz * i + j;
			}
		}
		return ERROR;
	}
	// do a destructive bitscan
	// returns ERROR if nothing found
	size_t bitscan_destructive(size_t ERROR) {
		for (size_t i = 0; i < N; ++i) {
			if (vals[i] != 0) {
				// compute the index
				// int ret = sz*i+(__builtin_ffsll(vals[i])-1);
				// __builtin_ctzll counts trailing 0s
				ull t = vals[i] & -vals[i]; // lowestOneBit
				size_t j = __builtin_ctzll(vals[i]);
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
	size_t bitscan() {
		for (int i = 0; i < N; ++i)
			if (vals[i] != 0)
				// return sz*i+ctz(vals[i]);
				return sz * i + __builtin_ctzll(vals[i]);
	}
#undef sz
#undef N
};
} // namespace fast_bitset

// #define WIDTH 3
// #define HEIGHT 3
// #define n 3
// const array<size_t, n> lengths = {3, 2, 2};
#define WIDTH 10
#define HEIGHT 10
#define n 5
const array<size_t, n> lengths = {5, 4, 3, 3, 2};

#define MAX_PLACEMENTS WIDTH*(HEIGHT - 1) + HEIGHT*(WIDTH - 1)
#define fast_set fast_bitset::bitset<MAX_PLACEMENTS>
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
