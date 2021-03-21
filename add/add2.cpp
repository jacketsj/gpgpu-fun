#include <CL/sycl.hpp>
#include <iostream>
#include <vector>

class vector_exp;

int main(int, char**) {
	std::vector<int> a;
	int MAG = 20;
	for (int i = 0; i < MAG; ++i)
		a.emplace_back(2);
	int n = 2;
	// compute a^n % 1e7+3 (without using squaring exp)
	const int MOD = 1e6 + 3;

	std::vector<int> c = a;

	cl::sycl::default_selector device_selector;

	cl::sycl::queue queue(device_selector);
	std::cout << "Running on "
						<< queue.get_device().get_info<cl::sycl::info::device::name>()
						<< "\n";

	// getting the total number of compute units
	auto num_groups =
			queue.get_device().get_info<cl::sycl::info::device::max_compute_units>();
	// getting the maximum work group size per thread
	auto work_group_size =
			queue.get_device()
					.get_info<cl::sycl::info::device::max_work_group_size>();
	// building the best number of global thread
	auto total_threads = num_groups * work_group_size;

	{
		cl::sycl::buffer<int> a_sycl(a.data(), cl::sycl::range<1>(MAG));
		cl::sycl::buffer<int> c_sycl(c.data(), cl::sycl::range<1>(MAG));

		queue.submit([&](cl::sycl::handler& cgh) {
			auto a_acc = a_sycl.get_access<cl::sycl::access::mode::read>(cgh);
			auto c_acc =
					c_sycl.get_access<cl::sycl::access::mode::discard_write>(cgh);

			cgh.parallel_for<class vector_exp>(
					cl::sycl::range<1>{total_threads}, [=](cl::sycl::item<1> item_id) {
						auto id = item_id.get_id(0);
						for (auto i = id; i < c_acc.get_count();
								 i += item_id.get_range()[0]) {
							auto current = a_acc[i];
							auto initial = current;
							for (auto j = 1; j < n; ++j) // start at 1, since we want a n to
																					 // be the power, not n+1
								current = (current * initial) % MOD;
							c_acc[i] = current;
						}
					});

			// cgh.single_task<class vector_addition>(
			// 		[=]() { c_acc[0] = a_acc[0] + b_acc[0]; });
		});
	}
	for (int i = 0; i < 20; ++i)
		std::cout << "a[" << i << "]=" << a[i] << ",c[" << i << "]=" << c[5]
							<< std::endl;
	// std::cout << "  A { " << a.x() << ", " << a.y() << ", " << a.z() << ", "
	//					<< a.w() << " }\n"
	//					<< "+ B { " << b.x() << ", " << b.y() << ", " << b.z() << ", "
	//					<< b.w() << " }\n"
	//					<< "------------------\n"
	//					<< "= C { " << c.x() << ", " << c.y() << ", " << c.z() << ", "
	//					<< c.w() << " }" << std::endl;

	return 0;
}
