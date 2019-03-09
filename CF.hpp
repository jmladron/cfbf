//
// Cuckoo Filter by S. Pontarelli and CFBF by J. Martinez
//

#include <cstdint>

using namespace std;

class CF {	
	public:

		int **table;               // CF memory
		unsigned char *bf;         // flags for Bloom Filter
		int cf_size;               // size of CF memory
		int cf_cells;			   // number of cells in each bucket
		int fp_size;               // 1 << f
		int cf_num_items;   	   // number of items inserted in CF
		int bf_num_items;          // number of items insserted in BF
		int victim_fingerprint;
		int victim_pointer;

		CF(int size, int cells, int f);
		virtual ~CF();
		int insert(uint64_t key, int bf_thresold);
		void random_remove();
		bool remove(uint64_t key);
		bool query(uint64_t key);
		int get_cf_items() { return cf_num_items; }
		int get_bf_items() { return bf_num_items; }
		int get_fingerprints();
		int get_bf_sets();
		int get_size() { return cf_cells*cf_size; }


	private:

		int insert2(uint64_t p, int fingerprint, int bf_thresold);
        int hash(uint64_t key, int i, int s);
        int RSHash(uint64_t key);
        int JSHash(uint64_t key);

};
