//
// CF.cpp : Cuckoo Filter by S. Pontarelli and CFBF by J. Martinez 
//

#include "pch.h"
#include <iostream>
#include "CF.hpp"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <functional>
#include <cstring>
#include <smmintrin.h>


/*
 * CF (Cuckoo Filter) constructor
 */

CF::CF(int size, int cells, int f)
{
	cf_size  = size;
	cf_cells = cells;
	fp_size  = (1 << f);

	// CF memory uses a table with M rows and c columns, BF flag is for the Bloom filter
	
	table = new int*[cf_size];
	bf = new unsigned char[cf_size];
	
	for (int i = 0; i < cf_size; i++)
		table[i] = new int[cf_cells];
	
	// the CF memory is set to -1, counters are set to 0

	cf_num_items = 0;
	bf_num_items = 0;
	victim_fingerprint = -1;
	victim_pointer = -1;

	for (int i = 0; i < cf_size; i++)
	{
		bf[i] = 0;

		for (int j = 0; j < cf_cells; j++)
			table[i][j] = -1;
	}
}


/*
 * CF (Cuckoo Filter) destructor
 */

CF::~CF()
{
	for (int i = 0; i < cf_size; i++)
		delete[] table[i];
	
	delete[] table;
	delete[] bf;
}


/*
 * insert calculates the fingerprint and calls insert2 function, which performs the insertion in the
 * table or sets the bits used by the Bloom filter
 */ 

int CF::insert(uint64_t key, int bf_threshold)
{
	int fingerprint = hash(key, 1, fp_size);

	int p = hash(key, 2, cf_size);
	p = p % cf_size;

	// bf_threshold == - 1 indicates that the key should be inserted in the Bloom Filter
	// otherwise, the key is inserted in the CFBF filter

	if (bf_threshold == -1)
		return insert_bf(p, fingerprint);
	else
		return insert_cfbf(p, fingerprint, bf_threshold);
}


/*
 * insert_cfbf inserts an element in the table:
 *
 * if there is an empty cell at positions p or p1, the element is inseted in the CF, otherwise
 * if the number of attemps is higher than the BF Threshold, the element is inserted in the BF if bf[p] and bf[p1] are zero
 * when both previous insertion fail, an element previously inserted is pushed to place the new element
 *
 * if the insertion fails insert_cfbf returns -1, otherwise it returns the number of insertion attempts run
 */

int CF::insert_cfbf(uint64_t p, int fingerprint, int bf_threshold)
{
	int p1;
	int j = 0;
	int jj = 0;
	int pushed_fingerprint = -1;

	for (int t = 1; t <= 1001; t++) {

		// the fingerprint may be inserted in 2 buckets given by 2 hash functions
		// each fingerprint may be stored in 2 x cf_cells different cells

		for (jj = 0; jj < cf_cells; jj++) {
			p = p % cf_size;

			if (table[p][jj] == -1) {
				table[p][jj] = fingerprint;

				cf_num_items++;

				return t;
			}

			// if table[p][jj] is already occupied, the alternative p1 position is checked

			p1 = p ^ hash(fingerprint, 2, cf_size);

			p1 = p1 % cf_size;

			if (table[p1][jj] == -1) {
				table[p1][jj] = fingerprint;

				cf_num_items++;

				return t;
			}
		}

		// at this point the 2 x cf_cells buckets are full and the new fingerprint is inserted in the Bloom filter
		// if the number of insertion attempts run is higher than the BF threshold

		// tries to insert on the BF with no increase on the FPR

		if (t > bf_threshold) {
			if ((bf[p] == 1) && (bf[p1] == 1)) {
				bf[p] = 1;
				bf[p1] = 1;

				bf_num_items++;

				return t;
			}
		}

		// tries to insert on the BF with small increase on the FPR

		if (t > bf_threshold + 50) {
			if ((bf[p] == 1) || (bf[p1] == 1)) {
				bf[p] = 1;
				bf[p1] = 1;

				bf_num_items++;

				return t;
			}
		}

		// inserts when we can if we are close to the limit

		if (t > 950) {
			bf[p] = 1;
			bf[p1] = 1;

			bf_num_items++;

			return t;
		}

		// at this point the insertion has failed and the new fingerprint pushes an element previously stored

		j = rand() % 2;
		p = p ^ (j*hash(fingerprint, 2, cf_size));
		p = p % cf_size;

		jj = rand() % cf_cells;

		pushed_fingerprint = table[p][jj];

		table[p][jj] = fingerprint;

		// check if the pushed fingerprint is a fingerprint

		if (pushed_fingerprint == -1) {
			cf_num_items++;
			return t;
		}

		// find a new place for the pushed fingerprint

		fingerprint = pushed_fingerprint;
	}

	// at this point the insertion has failed after 1001 attempts

	victim_pointer = p;
	victim_fingerprint = fingerprint;

	return -1;
}

int CF::insert_bf(uint64_t p, int fingerprint)
{
	int p1;

	p = p % cf_size;
	p1 = (p ^ hash(fingerprint, 2, cf_size)) % cf_size;

	// the key is inserted in the BF only if both positions are free

	if ((bf[p] == 0) && (bf[p1] == 0)) {
		bf[p] = 1;
		bf[p1] = 1;

		bf_num_items++;

		return 1;
	} 
	else {
		return -1;
	}
}

/*
 * random_remove deletes a random element
 */

void CF::random_remove()
{
	int p1, p2;

	do {
		p1 = rand() % cf_size;
		p2 = rand() % cf_cells;

	} while (table[p1][p2] == -1);

	table[p1][p2] = -1;
	cf_num_items--;
}


/*
 * remove deletes an element from the table when the fingerprint matches with one of the two positions given
 * by the hash functions, if the fingerprint is not found, the BF is checked and reset
 */

bool CF::remove(uint64_t key)
{
	int fingerprint = hash(key, 1, fp_size);
	fingerprint = fingerprint % fp_size;
	int p = hash(key, 2, cf_size);
	int *p1 = new int[2];
	bool found = false;

	for (int j = 0; j < 2; j++) {
		p = hash(key, 2, cf_size) ^ (j*hash(fingerprint, 2, cf_size));
		p = p % cf_size;

		// p1 stores both positions given by the hash functions

		p1[j] = p;

		// all the elements of the bucket should be cheched, from table[p][0] to table[p][cf_cells-1]

		for (int jj = 0; jj < cf_cells; jj++) {
			fingerprint = fingerprint % fp_size;

			// if the fingerprint matches table[p][jj], it is removed

			if (table[p][jj] == fingerprint) {
				table[p][jj] = -1;
				cf_num_items--;
				found = true;

				break;
			}
		}

		// if the fingerprint was found and removed, the loop ends

		if (found)
			break;
	}

	if (found)
		return true;
	else {

		// if the fingerprint was not found, the BF is checked and reset

		if ((bf[p1[0]] == 1) && (bf[p1[1]] == 1)) {
			bf[p1[0]] = 0;
			bf[p1[1]] = 0;
			bf_num_items--;

			return true;
		}
	}

	return false;
}


/*
 * query lookups an element in the table and the Bloom filter when necessary
 */

bool CF::query(uint64_t key)
{
	int fingerprint = hash(key, 1, fp_size);
	fingerprint = fingerprint % fp_size;
	int p = hash(key, 2, cf_size);
	int p1[2];

	if ((fingerprint == victim_fingerprint) && (p == victim_pointer)) return true;

	for (int j = 0; j < 2; j++) {
		p = hash(key, 2, cf_size) ^ (j*hash(fingerprint, 2, cf_size));
		p = p % cf_size;

		// p1 stores both positions given by the hash functions

		p1[j] = p;

		// all the elements of the bucket should be cheched, from table[p][0] to table[p][cf_cells-1]

		for (int jj = 0; jj < cf_cells; jj++) {
			fingerprint = fingerprint % fp_size;

			if (table[p][jj] == fingerprint)
				return true;
		}
	}

	if ((bf[p1[0]] == 1) && (bf[p1[1]] == 1))
		return true;

	return false;
}


/*
 * get_fingerprints returns the number of elements stored in the table
 */

int CF::get_fingerprints()
{
	int fingerprints = 0;

	for (int i = 0; i < cf_size; i++)
		for (int j = 0; j < cf_cells; j++)
			if (table[i][j] != -1)
				fingerprints++;

	return fingerprints;
};


/*
* get_bfsets returns the number of 1s stored in the Bloom filter
*/

int CF::get_bf_sets()
{
	int bfsets = 0;

	for (int i = 0; i < cf_size; i++)
		if (bf[i] == 1)
			bfsets++;

	return bfsets;
};


/*
 * hash functions
 */

int CF::hash(uint64_t key, int i, int s)
{
	unsigned int  val = RSHash(key) + i * JSHash(key);
	
	if (i == 2) val = RSHash(key);
	if (i == 1) val = JSHash(key);
	
	return (val % s);
}

// RSHash function

int CF::RSHash(uint64_t key)
{
	int b = 378551;
	int a = 63689;
	int hash = 0;
	int i = 0;
	char k[8];

	memcpy(k, &key, 8);

	for (i = 0; i < 8; i++)
	{
		hash = hash * a + k[i];
		a = a * b;
	}

	return hash;
}

// JSHash function

int CF::JSHash(uint64_t key)
{
	int hash = 1315423911;
	int i = 0;
	char k[8];

	memcpy(k, &key, 8);

	for (i = 0; i < 8; i++)
	{
		hash ^= ((hash << 5) + k[i] + (hash >> 2));
	}

	return hash;
}
