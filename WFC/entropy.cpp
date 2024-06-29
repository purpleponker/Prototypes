#include "WFC.h"



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

