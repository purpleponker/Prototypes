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
//collapse, close index and return true else return false for contradiction
bool wfc_mapping::collapse(int row, int col){
	if(wave_matrix[row][col].is_collapsed)
		return false;
	if(wave_matrix[row][col].superpositions.empty())
		return false;
	//else not empty superposition values or collapsed so collapse index
	int rand_val = rand() % wave_matrix[row][col].superpositions.size();
	auto it = wave_matrix[row][col].superpositions.begin();
	auto val = std::next(it, rand_val);
	
	//after selection from valid superpositions clear the set and change flags
	wave_matrix[row][col].is_collapsed = true;
	wave_matrix[row][col].is_changed = false;
	wave_matrix[row][col].collapse_val = (*val);
	wave_matrix[row][col].superpositions.clear();
	 bitmap[row][col] = (*val);
	 max_collapse_count -= 1;
	 collapse_count += 1;
	
	return true;
}

///////////////////////////////////////////
//propogation on adjacent indecies of collpased index
//returns true if propagation successful, and false if all propagation targets are already collapsed
bool wfc_mapping::propagate(int row, int col){
	if(collapse_size == collapse_count){
		//remove data type from loop history if region size reaches max
		auto found = loop_history.find(wave_matrix[row][col].collapse_val);
		if(found != loop_history.end())
			loop_history.erase(found);
	}	
	
	int adj =0;
	data_type tar_val = wave_matrix[row][col].collapse_val;
	//check bounds on adjacent indecies and add to entropy set if valid
	if(row - 1 >=0 ){//top index
		auto & index = wave_matrix[row-1][col];
		if(!index.is_collapsed){
			prop_index(index, tar_val);
			prop_set.push_back(std::make_pair(row-1, col));
			adj++;
		}
	}
	if(row + 1 < rows){//bottom
		auto & index = wave_matrix[row + 1][col];
		if(!index.is_collapsed){
			prop_index(index, tar_val);
			prop_set.push_back(std::make_pair(row+1, col));
			adj++;
		}
	}
	if(col - 1 >= 0){//left index
		auto & index = wave_matrix[row][col-1];
		if(!index.is_collapsed){
			prop_index(index, tar_val);
			prop_set.push_back(std::make_pair(row, col-1));
			adj++;
		}
	}
	if(col + 1 < cols){//right index
		auto & index = wave_matrix[row][col+1];
		if(!index.is_collapsed){
			prop_index(index, tar_val);
			prop_set.push_back(std::make_pair(row, col+1));
			adj++;
		}
	}
	//return tue if any adjacent indecies valid else return false
	if(adj  > 0)
		return true;
	return prop_unstuck(tar_val);
}

//helper function for propagate
void wfc_mapping::prop_index(wave_index & adj_index, data_type collapse_value){
	//adjust superpostion of index and add to prop set
	adj_index.superpositions.clear();
	adj_index.superpositions.insert(collapse_value);
}

//propagatiom helper when collapse corners itself
bool wfc_mapping::prop_unstuck(data_type upd_val){
	if(entropy_list.empty())
		return false;
	//could find a better optimized solution
	for(auto it : entropy_list.back()){
		wave_matrix[it.first][it.second].superpositions.clear();
		wave_matrix[it.first][it.second].superpositions.insert(upd_val);
	}
	return true;
}

///////////////////////////////////////////
//entropy recalc entropy where balid super positions are changed
std::pair<int, int> wfc_mapping::entropy(){
	//go thru prop_set and recalc entropy for each index, push back values to 2D vector depending on state of 2D vectoe
	
		//return -1,-1 if entropy list is still empty at end. either propagation failed or wfc is complete. else return random min value of entropy set
	if(entropy_list.empty() && max_collapse_count == 0){//catch last collapse on col_count-- each collapse
		return std::make_pair(-1,-1);
	}	
	//if region size reaches max entropy over whole region
	if(collapse_size == collapse_count){
		return region_ent();
	}
	
	while(!prop_set.empty()){// while not empty
		//recalc entropy, entropy = sum(p(x) * log2(1/(p(x))))
			float shan_ent = entropy_formula(prop_set.back().first, prop_set.back().second);
			wave_matrix[prop_set.back().first][prop_set.back().second].shan_entropy = shan_ent;
	
		//entropy list is empty
		if(entropy_list.empty()){
			std::vector<std::pair<int,int>> temp;
			temp.push_back(prop_set.back());
			entropy_list.push_back(temp);
			prop_set.pop_back();
		}
		else{//entropy list not empty
			//case 1: current region aleady has already started to populate the entropy list and its not the first region
			if(entropy_list.size() > 1){
				//case 1.a new min, add a new index at the back, push back new index extending vector size by 1
				if(shan_ent < wave_matrix[entropy_list.back()[0].first][entropy_list.back()[0].second].shan_entropy){
					std::vector<std::pair<int,int>> temp;
					temp.push_back(prop_set.back());
					entropy_list.push_back(temp);
					prop_set.pop_back();	
				}
				//case 1.b same value min, add to back() has current min values for region
				else if (shan_ent == wave_matrix[entropy_list.back()[0].first][entropy_list.back()[0].second].shan_entropy){
					entropy_list.back().push_back(prop_set.back());
					prop_set.pop_back();
				}
				//case 1.c greater than min add to index 1, has current region values
				else{
					entropy_list[1].push_back(prop_set.back());
					prop_set.pop_back();
				}
			}
			//case 2 first region
			else if(end_data_type * d_count == loop_history.size()){
				entropy_list.back().push_back(prop_set.back());
				prop_set.pop_back();
			}
			
			//case 3 only old indecies are stored in entropy list at index 0 from region entropy function reduction
			else{
				std::vector<std::pair<int,int>> temp;
				temp.push_back(prop_set.back());
				entropy_list.push_back(temp);
				prop_set.pop_back();	
			}	
		}
	}

	//check entropylist has valid data	
	if(entropy_list.empty())
		return std::make_pair(-1,-1);
	
	//if entropy list has valid data find.next index for collapse
	std::pair<int,int> last;
	bool valid = false;
	//select rand value from min entropy values
	do{
		if(entropy_list.empty())
			return std::make_pair(-1,-1);
		int ran = rand() % entropy_list.back().size();
		last = entropy_list.back()[ran];
		//swap the rand index with back if size greater than 1 then pop back
		if(entropy_list.back().size() > 1)
			std::iter_swap(entropy_list.back().begin() + ran, entropy_list.back().end()-1);
		//after conditional swap pop back
		entropy_list.back().pop_back();
		if(!wave_matrix[last.first][last.second].is_collapsed){
			valid = true;
		}
		if(entropy_list.back().empty()){
			entropy_list.pop_back();
		}
	}while(valid == false);

	return last;
}

//entropy over region when size of collapses = max count
std::pair <int,int> wfc_mapping::region_ent(){
	//added prop set to entropy list and recalc entropy for whole entropy list
	while(!prop_set.empty()){
		entropy_list.back().push_back(prop_set.back());
		prop_set.pop_back();
	}

	//set all superpositions in index 0 to loop history
	for(int i = entropy_list[0].size()-1; i >= 0; i--){
		auto ind = entropy_list[0][i];
		wave_matrix[ind.first][ind.second].superpositions = loop_history;
		wave_matrix[ind.first][ind.second].shan_entropy = entropy_formula(ind.first, ind.second);
		wave_matrix[ind.first][ind.second].is_tracked= true;	
	}
		
	//reduce 2D array down to size 2, all into index 0, and level index 1 there bit empty for next data region valuea. also set superposotions accodring to loophistory
	for(auto it = entropy_list.end()-1; it != entropy_list.begin(); it-=1){
		while(!(*it).empty()){
			//dont add duplicates to list
			if(wave_matrix[(*it).back().first][(*it).back().second].is_tracked){
				(*it).pop_back();
			}
			
			//reduce down to 1d and marked is_tracked as true
			else{ 
				wave_matrix[(*it).back().first][(*it).back().second].superpositions = loop_history;
				wave_matrix[(*it).back().first][(*it).back().second].shan_entropy = entropy_formula((*it).back().first, (*it).back().second);
				entropy_list[0].push_back((*it).back());
				wave_matrix[(*it).back().first][(*it).back().second].is_tracked = true;
				(*it).pop_back();
			}
		}
		entropy_list.pop_back();		
	}
	//return if list is empty
	if(entropy_list[0].empty())
		return std::make_pair(-1,-1);
		
	int ran =0;
	bool valid = false;
	std::pair<int,int> next_ind;
	//select new index at random
	do{
		ran = rand() % entropy_list[0].size();
		if(entropy_list[0].empty())
			return std::make_pair(-1,-1);	
		next_ind = entropy_list[0][ran];
		//see if index is collapsed
		if(wave_matrix[next_ind.first][next_ind.second].is_collapsed == false){
			valid = true;
		}
		//swap index with back and pop
		else{
			std::iter_swap(entropy_list[0].begin() + ran, entropy_list[0].end()-1);
			entropy_list[0].pop_back();		
		}
	}while(valid == false);
	//swap selected index with back for pop if size greater than 1
	if(entropy_list[0].size() > 1){
		std::iter_swap(entropy_list[0].begin() + ran, entropy_list[0].end()-1);
	}
	entropy_list[0].pop_back();
	collapse_count = 0;				
	return next_ind;
}

//entropy formula calc entropy value for a goven index
float wfc_mapping::entropy_formula(int row, int col){
	auto supers = wave_matrix[row][col].superpositions;
	float len = static_cast<float>(supers.size());
	float prob = 1.f/len;
	float log_prob = std::abs(static_cast<float>(log2(prob)));
	float shan_ent = len * (prob * log_prob);
	//std::cout << len << ", " << prob << ", " << log_prob << std::endl;
	return shan_ent;
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