#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include <stdlib.h>
#include <time.h>
#include <iterator>
#include <cmath>

//enumeration of data types for super position values
enum data_type : int{
	type_a,
	type_b,
	type_c,
	type_d,
	type_e,
	type_f,
	type_g,
	type_h,
	type_i,
	type_j,
	end_data_type
};

//info for wave function
struct wave_index{
	data_type collapse_val;
	bool is_collapsed;
	bool is_changed; // try to use later for entropy computation redundancy reduction
	bool is_tracked; // for entropy matrix reduction
	float shan_entropy;
	std::multiset<data_type> superpositions;
};

//class for matrix mapping
class wfc_mapping{
	public:
	wfc_mapping();
	~wfc_mapping();
	//collapse matrix index
	bool collapse(int row, int col);
	
	//propagate the collapse
	bool propagate(int row, int col);
	
	//recalc entropy
	std::pair<int, int> entropy();
	
	//random entry point selection
	std::pair<int,int> start();
	
	//display bitmap
	void display();
	
	protected:
	//propagation helper funcs
	void prop_index(wave_index & adj_index, data_type collapse_val);
	bool prop_unstuck(data_type upd_val); //when collapse corners itself and no adjacent indecies are uncollapsed
	std::pair<int,int> region_ent(); //when collapse count for region reaches max size ie collapse_count == collapse_size, 
	float entropy_formula(int row, int col);
	
	// temp hard code values
	static const int rows = 40; //height
	static const int cols= 40; //width
	int d_count = 2; //count size of duplicate data sections
	int collapse_size = (rows * cols) /  (end_data_type * d_count);
	int max_collapse_count = rows * cols;
	int collapse_count = 0;
	
	//matrix data
	char data[end_data_type] ={'A','B','C','D','E','F','G','H','I','J'};
	data_type bitmap[rows][cols];
	wave_index wave_matrix[rows][cols];
	//add indecies in propagation func call, has indecies which need entropy recalc
	std::vector<std::pair<int,int>> prop_set;
	//list of all indecies with smallest entropy value at back. in future use another var to track when a region finish populating so 2 of the same region not adjacent to each other. add indecies greater than back to begin. can sort them later as needed when vector is size 1
	std::vector<std::vector<std::pair<int, int>>> entropy_list;
	
	//remove a data type from set once the count reaches max size
	std::multiset<data_type> loop_history;
	
};
