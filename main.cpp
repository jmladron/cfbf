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
	int test_mode;
	int filter_size;
	int filter_cells_in_bucket;
	float occupancy;
	int bf_threshold;
	int *bf_keys;
	int bf_keys_len;
	int insertion_burst;
	int run_trials;
	int *fingerprint_bits;
	int fingerprint_bits_len;
};


/*
 * test1_data stores the output of test1: CFBF filter test with insertion burst
 */

struct test1_data {
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
 * test2_data stores the output of test2: CFBF filter FPR when the BF stores p keys
 */

struct test2_data {
	int filter_items;
	int cf_items;
	int bf_items;
	int bf_ones;
	double occupancy;
	double FPR;
};


/*
 * sets the test setup with the command line arguments
 */

test_setup setup_CFBF_test(int argc, char* argv[]) {
	test_setup setup;

	setup.test_mode              = -1;
	setup.filter_size            = -1;
	setup.filter_cells_in_bucket = -1;
	setup.occupancy              = -1;
	setup.bf_threshold           = -1;
	int p                        = -1;
	setup.insertion_burst        = -1;
	setup.run_trials             = -1;
	int f                        = -1;

	for (int i = 1; i < argc; i++) {
		std::stringstream str;
		str << argv[i];

		std::string argument = str.str();

		// argument 'm' is used for test mode, default value is m=1 (test insertion burst), m=2 (test insertion using bf)

		if (argument[0] == 'm' || argument[0] == 'M') {
			int len = argument.size();
			setup.test_mode = atoi(argument.substr(2, len - 1).c_str());
		}

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
			setup.occupancy = std::stof(argument.substr(2, len - 1).c_str());
		}

		// argument 't' is used for bf_threshold, default value is t=1000

		if (argument[0] == 't' || argument[0] == 'T') {
			int len = argument.size();
			setup.bf_threshold = atoi(argument.substr(2, len - 1).c_str());
		}

		// argument 'p' is used for keys in BF, default value is p = {40, 81, 122, 163, 204, 245, 286, 327, 368, 409} 

		if (argument[0] == 'p' || argument[0] == 'P') {
			int len = argument.size();
			p = atoi(argument.substr(2, len - 1).c_str());
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

	// set default values for test mode, filter size, cells in bucket, occupancy, bf_threshold, runs (trials) and fingerprint_bits

	if (setup.test_mode == -1)
		setup.test_mode = 1;

	if (setup.filter_size == -1)
		setup.filter_size = 8192;

	if (setup.filter_cells_in_bucket == -1)
		setup.filter_cells_in_bucket = 4;

	if (setup.occupancy == -1)
		setup.occupancy = 95;

	if (setup.bf_threshold == -1)
		setup.bf_threshold = 1000;

	if (p == -1) {
		setup.bf_keys_len = 10;
		setup.bf_keys = new int[setup.bf_keys_len];
		setup.bf_keys[0] = 40;  
		setup.bf_keys[1] = 81;
		setup.bf_keys[2] = 122;
		setup.bf_keys[3] = 163;
		setup.bf_keys[4] = 204;
		setup.bf_keys[5] = 245;
		setup.bf_keys[6] = 286;
		setup.bf_keys[7] = 327;
		setup.bf_keys[8] = 368;
		setup.bf_keys[9] = 409;
	}
	else {
		setup.bf_keys_len = 1;
		setup.bf_keys = new int[setup.bf_keys_len];
		setup.bf_keys[0] = p;
	}
 
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
	if (setup.test_mode == 1)
		std::cout << "CF-BF filter test with insertion burst\n\n";
	else
		std::cout << "CF-BF filter FPR having p keys in the BF\n\n";

	std::cout << "Test mode        m = " << setup.test_mode << "\n";
	std::cout << "Filter size      s = " << setup.filter_size << "\n";
	std::cout << "Cells in bucket  c = " << setup.filter_cells_in_bucket << "\n";
	std::cout << "Occupancy        o = " << setup.occupancy << "%\n";
	std::cout << "BF threshold     t = " << setup.bf_threshold << "\n";

	if (setup.test_mode == 1) {
		std::cout << "Insertion burst  b = " << setup.insertion_burst << "\n";
	}
	else {
		std::cout << "Keys in BF       p = {";

		for (int i = 0; i < setup.bf_keys_len; i++)
			std::cout << " " << setup.bf_keys[i] << ((i == setup.bf_keys_len - 1) ? " " : ",");

		std::cout << "} \n";
	}

	std::cout << "Runs (trials)    r = " << setup.run_trials << "\n";

	std::cout << "Fingerprint bits f = {";

	for (int i = 0; i < setup.fingerprint_bits_len; i++)
		std::cout << " " << setup.fingerprint_bits[i] << ((i == setup.fingerprint_bits_len - 1) ? " " : ",");

	std::cout << "} \n";
}


/*
 * test: filter initialization, filter warm-up and FPR calculation
 */

test1_data test1(int total_keys, int filter_size, int filter_cells_in_bucket, int bf_threshold, int insertion_burst, int fingerprint_bits) {
	test1_data res;

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
 * test2: filter initialization, filter warm-up, BF insertion and FPR calculation
 */

test2_data test2(int total_keys, int filter_size, int filter_cells_in_bucket, int bf_threshold, int bf_keys, int fingerprint_bits) {
	test2_data res;

	CF filter = CF(filter_size, filter_cells_in_bucket, fingerprint_bits);

	int rand_rule;
	uint64_t key;
	const uint64_t base = 17179869184;
	const uint64_t prefix_offset = 8589934592;

	int filter_errors = 0;
	int query_errors = 0;

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
	query_errors = 0;

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

	// insertion in the Bloom Filter a maximum of bf_keys 

	offset = total_keys + 1 + (filter_size * filter_cells_in_bucket) + 1;

	int k = 1;

	while (k <= bf_keys) {

		rand_rule = rand() % 1000000 + 1;

		key = (uint64_t)offset + (uint64_t)k + (uint64_t)(rand_rule * base) + (uint64_t)prefix_offset;

		// when the bf_threshold is negative (-1), the key is inserted in the Bloom Filter

		int insert_attemps = -1;

		insert_attemps = filter.insert(key, -1);

		if (insert_attemps != -1) {
			if (!filter.query(key)) {
				query_errors++;
			}
			else {
				k++;
			}
		}

	}

	// FPR calculation querying 1M non existing keys in CFBF

	offset = total_keys + 1 + (filter_size * filter_cells_in_bucket) + bf_keys + 1;

	int max_queries = 1000000;
	int positive_queries = 0;

	for (int i = 1; i <= max_queries; i++) {

		rand_rule = rand() % 1000000 + 1;

		key = (uint64_t)offset + (uint64_t)i + (uint64_t)(rand_rule * base) + (uint64_t)prefix_offset;

		// the key is queryied in the filter
		
		if (filter.query(key)) {
			positive_queries++;
		}

	}

	res.cf_items = filter.get_cf_items();
	res.bf_items = filter.get_bf_items();
	res.bf_ones = filter.get_bf_sets();
	res.filter_items = res.cf_items + res.bf_items;
	res.occupancy = (double)(res.cf_items) / (double)(filter_size * filter_cells_in_bucket);
	res.FPR = (double) positive_queries / (double)max_queries;

	return res;
}


/*
 * runs the test multiple times for each fingerprint size
 */

void CFBF_run_test1(int total_keys, test_setup setup) {
	test1_data res;
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
		ss << "m1_s" << setup.filter_size << "_c" << setup.filter_cells_in_bucket << "_o" << setup.occupancy << "_t" << setup.bf_threshold << "_f" << f << "_b" << setup.insertion_burst << "_r" << setup.run_trials << ".txt";
		data_file = ss.str();

		ofstream test_output(data_file);

		test_output << "%keys " << total_keys << ", insertion burst " << setup.insertion_burst << "\n";
		test_output << "%\n";
		test_output << "%items\tcf\tbf\t1s\titer\t max\tavg\terr\tquery\tocc\n";

		for (int r = 1; r <= setup.run_trials; r++) {
			res = test1(total_keys, setup.filter_size, setup.filter_cells_in_bucket, setup.bf_threshold, setup.insertion_burst, f);

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

void CFBF_run_test2(int total_keys, test_setup setup) {
	test2_data res;
	int f, p;

	for (int i = 0; i < setup.fingerprint_bits_len; i++) {
		f = setup.fingerprint_bits[i];

		for (int k = 0; k < setup.bf_keys_len; k++) {
			p = setup.bf_keys[k];

			std::cout << "\nTest with s=" << setup.filter_size << ", c=" << setup.filter_cells_in_bucket << ", o=" << setup.occupancy << "%, t=" << setup.bf_threshold << ", p=" << p << ", r=" << setup.run_trials << ", f=" << f << "\n";

			string data_file;

			std::stringstream ss;
			ss << "m2_s" << setup.filter_size << "_c" << setup.filter_cells_in_bucket << "_o" << setup.occupancy << "_t" << setup.bf_threshold << "_p" << p << "_f" << f << "_r" << setup.run_trials << ".txt";
			data_file = ss.str();

			ofstream test_output(data_file);

			test_output << "%total keys " << total_keys << ", keys in BF " << p << "\n";
			test_output << "%\n";
			test_output << "%items\tcf\tbf\t1s\tocc\tFPR\n";

			for (int r = 1; r <= setup.run_trials; r++) {
				res = test2(total_keys, setup.filter_size, setup.filter_cells_in_bucket, setup.bf_threshold, p, f);

				test_output << res.filter_items << "\t" << res.cf_items << "\t" << res.bf_items << "\t" << res.bf_ones << "\t" << res.occupancy << "\t" << res.FPR << "\n";
			}

			test_output.close();
		}
	}
}

int main(int argc, char* argv[]) {
	test_setup setup;

	// command line arguments for CFBF test program
	//
	// m: test mode, 1:CFBF filter test with insertion burst, 2:CFBF filter FPR having p keys in the BF, default value is 1 
	// s: filter size, default value is 8192
	// c: filter cells in bucket, default value is 4
	// o: filter occupancy, default value is 95%
	// t: bf threshold, default value is 1000
	// b: insertion burst, default value is 1% of filter size x filter cells in bucket 
	// p: keys in BF, default value is p = {40, 81, 122, 163, 204, 245, 286, 327, 368, 409} 
	// r: runs (trials), default value is 1000
	// f: fingerprint_bits, default value is f = {12}
	//
	// example: CFBF.exe o=94 t=10 b=491 p=40 r=5000 f=12

	setup = setup_CFBF_test(argc, argv);

	// show CF-BF test parameters

	display_CFBF_test_parameters(setup);

	// run CFBF test

	int total_keys = (setup.filter_size * setup.filter_cells_in_bucket) * (double)(setup.occupancy) / 100.0;

	if (setup.test_mode == 1)
		CFBF_run_test1(total_keys, setup);
	else
		CFBF_run_test2(total_keys, setup);

	return 0;

}
