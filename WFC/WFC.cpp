#include "WFC.h"

/*
to do:
-refactor so everything not in one file

- change the containter types so a list of uncollapsed indecies is populated as a hash map in the constructor and deleted an index
from the hash made when it is collapsed. Keep a list of adjacent indecies from collapsed indecies for next collapse and entropy
recalc

- change superposition multiset to a hash map and erase values when region finished, repoluation the super position when full
loop and region collapses finishes

- minor bug with collapse count on duolicate data regions > 1. yields uneven dispersal of data types. typically sets of two types share the offset.

-clean up debug code

-try using is changed flag to increase computation efficiency

-change region size so that after filling all regions from all loops there are exactly 1 region worth of tiles remaining to file,
doing this will remove the strange sparcity of the last region fill being so spread out. in this last fill just populate the
remaining uncollapsed tiles with any collapsed adjecent tile use a while loop untill no tiles remain in uncollapsed list

*/

///////////////////////////////////////////
//constructor for wfc_mapping class
wfc_mapping::wfc_mapping(){
	//default superpostions
	std::multiset<data_type> superpos_set;
	//populate superpos set, d_count number of each data_type
	for(int i = type_a; i != end_data_type; i++){
		for(int j = 1; j <= d_count; j++){
			data_type temp = static_cast<data_type>(i);
			superpos_set.insert(temp);
			loop_history.insert(temp);
		}	
	}
	
	//max entropy val, ent = sum(p(x) * log2(1/p(x)), so with 10 options ent = 10 * [(1/10) * log2(1/(1/10))], done in the entropy_formula func called in entropy func
	float start_ent = 1000.f;
	//default values in wave matrix and bitmap
	struct wave_index def_val = {end_data_type, false, false, false, start_ent, superpos_set};
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			wave_matrix[i][j] = {def_val};
			bitmap[i][j] = end_data_type;
		}
	}
}


///////////////////////////////////////////
wfc_mapping::~wfc_mapping(){
	
}

///////////////////////////////////////////
//matrix starting index selection
std::pair<int,int> wfc_mapping::start(){
	auto index = std::make_pair(rand() % rows, rand() % cols);
	return index;
}
///////////////////////////////////////////
//display all data, includes values during generatuon, helpful for debugging
void wfc_mapping::display(){
	int data_counts[end_data_type] = {0};
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++){
			auto val = bitmap[i][j];
			if(val != end_data_type){
				data_counts[val] += 1;	
				std::cout << data[static_cast<int>(val)];
			}
			else
				std::cout << "0"; //null otherwise from end_data_type
		}
		std::cout << std::endl;
	}
	
	// debug
	std::cout << "entropy and superpositions: " << std::endl;
	int level = 0;
	for(auto i : entropy_list){
		std::cout << "matrix level: " << level << std::endl;
		for(auto j : i){
			std::cout << "index: " << j.first << "," << j.second <<":" <<wave_matrix[j.first][j.second].shan_entropy << ", superpositions: ";
			for(auto k : wave_matrix[j.first][j.second].superpositions){
				std::cout << data[k];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		level++;
	}	
	
	//debug
	std::cout << std::endl;
	std::cout << "collapse count remaining: " << max_collapse_count <<std::endl;
	
	std::cout << "data_types remaining: " << std::endl;
	for(auto i: loop_history){
		std::cout << data[static_cast<int>(i)]<< " ";
	}
	std::cout << std::endl;
	
	//data counts debug
	std::cout << "data type counts: " << std::endl;
	for (int i = 0; i < end_data_type; i++){
		std::cout << data[i] << ": "<< data_counts[i] << std::endl;
	}
	std::cout << std::endl;
}