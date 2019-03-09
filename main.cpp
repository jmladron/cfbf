//
// CF-BF Cuckoo filter with integrated Bloom Filter by S. Pontarelli and J. Martinez
//

#include "pch.h"
#include "CF.hpp"
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <string>
 

/*
 * test_setup stores the command line arguments
 */

struct test_setup {
	int filter_size;
	int filter_cells_in_bucket;
	int occupancy;
	int bf_threshold;
	int insertion_burst;
	int run_trials;
	int *fingerprint_bits;
	int fingerprint_bits_len;
};


/*
 * test_data stores the test output
 */

struct test_data {
	int filter_items;
	int cf_items;
	int bf_items;
	int bf_ones;
	double occupancy;
	int total_iterations;
	int max_iterations;
	double avg_iterations;
	int filter_errors;
	int query_errors;
};


/*
 * sets the test setup with the command line arguments
 */

test_setup setup_CFBF_test(int argc, char* argv[]) {
	test_setup setup;

	setup.filter_size            = -1;
	setup.filter_cells_in_bucket = -1;
	setup.occupancy              = -1;
	setup.bf_threshold           = -1;
	setup.insertion_burst        = -1;
	setup.run_trials             = -1;
	int f                        = -1;

	for (int i = 1; i < argc; i++) {
		std::stringstream str;
		str << argv[i];

		std::string argument = str.str();

		// argument 's' is used for filter size, default value is s=8192

		if (argument[0] == 's' || argument[0] == 'S') {
			int len = argument.size();
			setup.filter_size = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'c' is used for the number of cells in bucket, default value is c=4

		if (argument[0] == 'c' || argument[0] == 'C') {
			int len = argument.size();
			setup.filter_cells_in_bucket = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'o' is used for occupancy, default value is o=95

		if (argument[0] == 'o' || argument[0] == 'O') {
			int len = argument.size();
			setup.occupancy = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 't' is used for bf_threshold, default value is t=1000

		if (argument[0] == 't' || argument[0] == 'T') {
			int len = argument.size();
			setup.bf_threshold = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'b' is used for the size of the insertion burst, default value is 1% of the filter cells

		if (argument[0] == 'b' || argument[0] == 'B') {
			int len = argument.size();
			setup.insertion_burst = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'r' is used for runs (trials), default value is r=10000

		if (argument[0] == 'r' || argument[0] == 'R') {
			int len = argument.size();
			setup.run_trials = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'f' is used for fingerprint_bits, default value is f = {12}

		if (argument[0] == 'f' || argument[0] == 'F') {
			int len = argument.size();
			f = atoi(argument.substr(2, len - 1).c_str());
		}
	}

	// set default values for filter size, cells in bucket, occupancy, bf_thershold, runs (trials) and fingerprint_bits

	if (setup.filter_size == -1)
		setup.filter_size = 8192;

	if (setup.filter_cells_in_bucket == -1)
		setup.filter_cells_in_bucket = 4;

	if (setup.occupancy == -1)
		setup.occupancy = 95;

	if (setup.bf_threshold == -1)
		setup.bf_threshold = 1000;

	if (setup.insertion_burst == -1)
		setup.insertion_burst = (setup.filter_size * setup.filter_cells_in_bucket) * (0.01);

	if (setup.run_trials == -1)
		setup.run_trials = 10000;

	if (f == -1) {
		setup.fingerprint_bits_len = 1;
		setup.fingerprint_bits = new int[setup.fingerprint_bits_len];
		setup.fingerprint_bits[0] = 12;
	}
	else {
		setup.fingerprint_bits_len = 1;
		setup.fingerprint_bits = new int[setup.fingerprint_bits_len];
		setup.fingerprint_bits[0] = f;
	}

	return setup;
}


/*
 * shows the test setup on the console
 */

void display_CFBF_test_parameters(test_setup setup) {
	std::cout << "CF-BF Test \n\n";

	std::cout << "Filter size      s = " << setup.filter_size << "\n";
	std::cout << "Cells in bucket  c = " << setup.filter_cells_in_bucket << "\n";
	std::cout << "Occupancy        o = " << setup.occupancy << "%\n";
	std::cout << "BF threshold     t = " << setup.bf_threshold << "\n";
	std::cout << "Insertion burst  b = " << setup.insertion_burst << "\n";
	std::cout << "Runs (trials)    r = " << setup.run_trials << "\n";

	std::cout << "Fingerprint bits f = {";

	for (int i = 0; i < setup.fingerprint_bits_len; i++)
		std::cout << " " << setup.fingerprint_bits[i] << ((i == setup.fingerprint_bits_len - 1) ? " " : ",");

	std::cout << "} \n";
}


/*
 * test: filter initialization, filter warm-up and FPR calculation
 */

test_data test(int total_keys, int filter_size, int filter_cells_in_bucket, int bf_threshold, int insertion_burst, int fingerprint_bits) {
	test_data res;

	CF filter = CF(filter_size, filter_cells_in_bucket, fingerprint_bits);

	int rand_rule;
	uint64_t key;
	const uint64_t base = 17179869184;
	const uint64_t prefix_offset = 8589934592;
	
	int filter_errors = 0;
	int query_errors  = 0;

	// filter initialization using random keys

	for (int i = 1; i <= total_keys; i++) {

		rand_rule = rand() % 1000000 + 1;

		key = (uint64_t)i + (uint64_t)(rand_rule * base) + (uint64_t)prefix_offset;

		// the key is inserted in the CF filter using bf_threshold 1000 to avoid insertions in the BF

		int insert_attemps = -1;

		insert_attemps = filter.insert(key, 1000);

		if (insert_attemps == -1) {
			filter_errors++;
		}
		else {
			if (!filter.query(key)) {
				query_errors++;
			}
		}
	}

	// filter warm-up deleting and inserting filter_size * filter_cells_in_bucket keys

	int offset = total_keys + 1;

	filter_errors = 0;
	query_errors  = 0;

	int items;

	for (int i = 1; i <= (filter_size * filter_cells_in_bucket); i++) {

		// a random key is removed 

		filter.random_remove();

		// a new random key is calculated

		rand_rule = rand() % 1000000 + 1;

		key = (uint64_t)offset + (uint64_t)i + (uint64_t)(rand_rule * base) + (uint64_t)prefix_offset;

		// the key is inserted in the CF filter using bf_threshold 1000 to avoid insertions in the BF

		int insert_attemps = -1;

		insert_attemps = filter.insert(key, 1000);

		if (insert_attemps == -1) {
			filter_errors++;
		}
		else {
			if (!filter.query(key)) {
				query_errors++;
			}
		}
	}

	// dynamic test of the filter inserting 1% of the keys, max iterations and avg. iterations are calculated

	offset = total_keys + 1 + (filter_size * filter_cells_in_bucket);

	int total_iterations = 0;
	int max_iterations   = 0;
	filter_errors        = 0;
	query_errors         = 0;

	for (int i = 1; i <= insertion_burst; i++) {

		rand_rule = rand() % 1000000 + 1;

		key = (uint64_t)offset + (uint64_t)i + (uint64_t)(rand_rule * base) + (uint64_t)prefix_offset;

		// the key is inserted in the filter using original value for bf_threshold

		int insert_attemps = -1;

		insert_attemps = filter.insert(key, bf_threshold);

		if (insert_attemps == -1) {
			filter_errors++;
		}
		else {
			total_iterations = total_iterations + insert_attemps;

			if (max_iterations < insert_attemps)
				max_iterations = insert_attemps;

			if (!filter.query(key)) {
				query_errors++;
			}
		}
	}

	res.cf_items = filter.get_cf_items();
	res.bf_items = filter.get_bf_items();
	res.bf_ones  = filter.get_bf_sets();
	res.filter_items = res.cf_items + res.bf_items;
	res.occupancy = (double) (res.cf_items) / (double) (filter_size * filter_cells_in_bucket);
	res.total_iterations = total_iterations;
	res.max_iterations = max_iterations;
	res.avg_iterations = (double) total_iterations / (double) insertion_burst;
	res.filter_errors = filter_errors;
	res.query_errors = query_errors;

	return res;
}


/*
 * runs the test multiple times for each fingerprint size
 */

void CFBF_run_test(int total_keys, test_setup setup) {
	test_data res;
	int f;

	for (int i = 0; i < setup.fingerprint_bits_len; i++) {
		f = setup.fingerprint_bits[i];

		int max_iterations    = 0;
		int max_bf_items      = 0;
		double avg_iterations = 0.0;
		double avg_bf_items   = 0.0;

		std::cout << "\nTest with s=" << setup.filter_size << ", c=" << setup.filter_cells_in_bucket << ", o=" << setup.occupancy << "%, t=" << setup.bf_threshold << ", r=" << setup.run_trials << ", f=" << f << "\n";

		string data_file;

		std::stringstream ss;
		ss << "s" << setup.filter_size << "_c" << setup.filter_cells_in_bucket << "_o" << setup.occupancy << "_t" << setup.bf_threshold << "_f" << f << "_b" << setup.insertion_burst << "_r" << setup.run_trials << ".txt";
		data_file = ss.str();

		ofstream test_output(data_file);

		test_output << "%keys " << total_keys << ", insertion burst " << setup.insertion_burst << "\n";
		test_output << "%\n";
		test_output << "%items\tcf\tbf\t1s\titer\t max\tavg\terr\tquery\tocc\n";

		for (int r = 1; r <= setup.run_trials; r++) {
			res = test(total_keys, setup.filter_size, setup.filter_cells_in_bucket, setup.bf_threshold, setup.insertion_burst, f);

			if (res.max_iterations > max_iterations)
				max_iterations = res.max_iterations;

			if (res.bf_items > max_bf_items)
				max_bf_items = res.bf_items;

			avg_iterations = avg_iterations + res.avg_iterations;
			avg_bf_items = avg_bf_items + res.bf_items;

			test_output << res.filter_items << "\t" << res.cf_items << "\t" << res.bf_items << "\t" << res.bf_ones << "\t" << res.total_iterations << "\t " << res.max_iterations << "\t" << res.avg_iterations << "\t" << res.filter_errors << "\t" << res.query_errors << "\t" << res.occupancy << "\n";
		}

		avg_iterations = avg_iterations / (double)setup.run_trials;
		avg_bf_items = avg_bf_items / (double)setup.run_trials;

		test_output << "%\n";
		test_output << "%Max iterations " << max_iterations << "\n";
		test_output << "%Avg iterations " << avg_iterations << "\n";
		test_output << "%Max bf items   " << max_bf_items << "\n";
		test_output << "%Avg bf items   " << avg_bf_items << "\n";

		test_output.close();
	}
}


int main(int argc, char* argv[]) {
	test_setup setup;

	// command line arguments for CFBF test program
	//
	// s: filter size, default value is 8192
	// c: filter cells in bucket, default value is 4
	// o: filter occupancy, default value is 95%
	// t: bf threshold, default value is 1000
	// b: insertion burst, default value is 1% of filter size x filter cells in bucket 
	// r: runs (trials), default value is 1000
	// f: fingerprint_bits, default value is f = {12}
	//
	// example: CFBF.exe o=94 t=10 b=491 r=5000 f=12

	setup = setup_CFBF_test(argc, argv);

	// show CF-BF test parameters

	display_CFBF_test_parameters(setup);

	// run CFBF test

	int total_keys = (setup.filter_size * setup.filter_cells_in_bucket) * (double)(setup.occupancy) / 100.0;

	CFBF_run_test(total_keys, setup);

	return 0;

}
