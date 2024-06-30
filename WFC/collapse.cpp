#include "WFC.h"
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