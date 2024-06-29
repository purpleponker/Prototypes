#include "WFC.h"

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
