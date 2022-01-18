//#include <CL/sycl.hpp>

#include <bits/stdc++.h>
using namespace std;
//#include <algorithm>
//#include <bitset>
//#include <iostream>
//#include <utility>
//#include <vector>
// using std::array;
// using std::cout;
// using std::endl;
// using std::pair;
// using std::string;
// using std::vector;

//#include "algorithm.h"
// #include "fast_bitset.h"
//#include "graph.h"
//#include "standard_game.h"

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

pair<int, std::string> exec(const char* cmd) {
	std::array<char, 128> buffer;
	std::string result;
	FILE* pipe = popen(cmd, "r");
	if (!pipe)
		throw std::runtime_error("popen() failed!");
	try {
		while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
			result += buffer.data();
		}
	} catch (...) {
		int exit_status = WEXITSTATUS(pclose(pipe));
		// pclose(pipe);
		return make_pair(exit_status, "Error reading compiler output");
		// throw;
	}
	int exit_status = WEXITSTATUS(pclose(pipe));
	// pclose(pipe);
	/*
	std::array<char, 128> buffer; // TODO: Should this size be increased?
	std::string result;
	std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
	if (!pipe) {
		throw std::runtime_error("popen() failed!");
	}
	while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
		result += buffer.data();
	}
	*/
	// return result;
	return make_pair(exit_status, result);
}
pair<int, std::string> exec(const string& cmd) { return exec(cmd.c_str()); }

void run(const string& data) {
	string cmd = "echo \"";
	cmd += data;
	cmd += "\" | ./main_gpu_api";
	int code;
	string res;
	tie(code, res) = exec(cmd);
	if (code != 0)
		cout << "error while running\n";
	cout << res;
}

int main() {
	std::string command_line;
	array<array<char, HEIGHT>, WIDTH> board;
	for (size_t x = 0; x < WIDTH; ++x)
		for (size_t y = 0; y < HEIGHT; ++y)
			board[x][y] = '.';
	cout << "> ";
	while (getline(std::cin, command_line)) {
		if (!command_line.empty()) {
			std::stringstream ss(command_line);
			std::string command;
			ss >> command;
			// run and reset
			if (command == "help") {
				cout << "Available commands:\n";
				cout << "\tdraw: draw the current board state" << '\n';
				cout << "\tmiss [x] [y]: Add a miss to the board" << '\n';
				cout << "\thit [x] [y]: Add a hit to the board" << '\n';
				cout << "\tsink [x] [y] [ship]: (BROKEN) Add a sink to the board of "
								"ship [ship]"
						 << '\n';
				cout << "\tremove [x] [y]: (DEPRECATED) remove a specific point's "
								"query response"
						 << '\n';
				cout << "\tclear: clear the board" << std::endl;
				cout << "\trun: run on gpu" << std::endl;
				cout << "\tquit: quit" << std::endl;
				cout << "Ship sizes:\n";
				for (size_t i = 0; i < n; ++i)
					cout << "\tlength(" << i << ") = " << (lengths[i] + '0') << std::endl;
			} else if (command == "clear") {
				for (size_t x = 0; x < WIDTH; ++x)
					for (size_t y = 0; y < HEIGHT; ++y)
						board[x][y] = '.';
			} else if (command == "run") {
				string s;
				for (size_t y = 0; y < HEIGHT; ++y) {
					for (size_t x = 0; x < WIDTH; ++x) {
						s += board[x][y];
						s += ' ';
					}
					// s += '\n';
				}
				run(s);
			} else if (command == "draw") {
				cout << "Current board state:\n";
				cout << "   ";
				for (size_t x = 0; x < WIDTH; ++x)
					cout << "   " << x + 1;
				cout << '\n';
				cout << "   ";
				for (size_t x = 0; x < WIDTH; ++x)
					cout << "----";
				cout << '\n';
				for (size_t y = 0; y < HEIGHT; ++y) {
					if (y + 1 < 10)
						cout << y + 1 << " | ";
					else
						cout << y + 1 << "| ";
					for (size_t x = 0; x < WIDTH; ++x)
						cout << "  " << board[x][y] << " ";
					cout << '\n';
				}
			} else if (command == "remove") {
				size_t x, y;
				try {
					ss >> x >> y;
					board[x][y] = '.';
				} catch (std::exception e) {
					cout << "malformatted command \"remove [x] [y]\"" << endl;
				}
			} else if (command == "hit") {
				size_t x, y;
				try {
					ss >> x >> y;
					--x;
					--y;
					board[x][y] = 'o';
				} catch (std::exception e) {
					cout << "malformatted command \"hit [x] [y]\"" << endl;
				}
			} else if (command == "sink") {
				size_t x, y, ship;
				try {
					ss >> x >> y >> ship;
					--x;
					--y;
					board[x][y] = char(ship) + '0';
				} catch (std::exception e) {
					cout << "malformatted command \"sink [x] [y] [ship]\"" << endl;
				}
			} else if (command == "miss") {
				size_t x, y;
				try {
					ss >> x >> y;
					--x;
					--y;
					board[x][y] = 'x';
				} catch (std::exception e) {
					cout << "malformatted command \"miss [x] [y]\"" << endl;
				}
			} else {
				cout << "unknown command \"" << command << "\"" << endl;
				cout << "use \"help\" to get a list of available commands" << endl;
			}
		}
		cout << "> ";
	}
}
