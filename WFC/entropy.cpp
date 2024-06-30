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
	
	bool needs_sort = false;
	while(!prop_set.empty()){// while not empty
		//recalc entropy, entropy = sum(p(x) * log2(1/(p(x))))
			float shan_ent = entropy_formula(prop_set.back().first, prop_set.back().second);
			wave_matrix[prop_set.back().first][prop_set.back().second].shan_entropy = shan_ent;
	
		//entropy list is empty
		if(entropy_list.empty()){
			std::pair<int,int> temp;
			temp = prop_set.back();
			entropy_list.push_back(temp);
			prop_set.pop_back();
		}
		else{//entropy list not empty
			//case 1: current region aleady has already started to populate the entropy list and its not the first region
			if(entropy_list.size() > 1){
				//case 1.a new min, add a new index at the back, push back new index extending vector size by 1
				if(shan_ent < wave_matrix[entropy_list.back().first][entropy_list.back().second].shan_entropy){
					std::pair<int,int> temp;
					temp = prop_set.back();
					entropy_list.push_back(temp);
					prop_set.pop_back();	
				}
				//case 1.b same value min, add to back() has current min values for region
				else if (shan_ent == wave_matrix[entropy_list.back().first][entropy_list.back().second].shan_entropy){
					entropy_list.push_back(prop_set.back());
					prop_set.pop_back();
				}
				//case 1.c greater than min add to index 1, has current region values
				else{
					entropy_list.push_back(prop_set.back());
					prop_set.pop_back();
					needs_sort = true;
				}
			}
			//case 2 first region
			else if(end_data_type * d_count == loop_history.size()){
				entropy_list.push_back(prop_set.back());
				prop_set.pop_back();
			}
		}
	}
	if(needs_sort){
		std::sort(entropy_list.begin(), entropy_list.end());
	}

	//check entropylist has valid data	
	if(entropy_list.empty() && uncollapsed.empty() == false)
		return std::make_pair(-1,-1);
	
	//if entropy list has valid data find.next index for collapse
	std::pair<int,int> last;
	bool valid = false;
	//select min value or if more than 1 duplicate min value select rand entry from min entropy values
	do{
		if(entropy_list.empty())
			return std::make_pair(-1,-1);
		int ran = rand() % entropy_list.size();
		last = entropy_list[ran];
		//swap the rand index with back if size greater than 1 then pop back
		if(entropy_list.size() > 1)
			std::iter_swap(entropy_list.begin() + ran, entropy_list.end()-1);
		//after conditional swap pop back
		entropy_list.pop_back();
		if(!wave_matrix[last.first][last.second].is_collapsed){
			valid = true;
		}
		if(entropy_list.empty()){
			entropy_list.pop_back();
		}
	}while(valid == false);                       
	return last;
}

//entropy over region when size of collapses = max count
std::pair <int,int> wfc_mapping::region_ent(){
	std::sort(entropy_list.begin(), entropy_list.end());
	std::pair<int,int> next_reg;
	if(uncollapsed.empty()){
		next_reg.first = -1;
		next_reg.second = -1;
	}
	//get min value entropy or select from pool of duplicate minimum values
	//check if more than 1 element in list, 1 element is size 2 including /null terminating element
	if(entropy_list.size() > 2){
		//pool min value elements into vector and select rand element
		int i = entropy_list.size() - 1;
		int j = entropy_list.size() -2;
		next_reg.first = entropy_list.back().first;
		next_reg.second = entropy_list.back().second;
		std::vector<std::pair<int,int>> next_reg_list;
		next_reg_list.push_back(next_reg);
		while(wave_matrix[entropy_list[i].first][entropy_list[i].second].shan_entropy == wave_matrix[entropy_list[j].first]
		[entropy_list[j].second].shan_entropy){
			next_reg.first=entropy_list[j].first;
			next_reg.second=entropy_list[j].second;
			next_reg_list.push_back(next_reg);
		}
		int rand_sel = rand() % next_reg_list.size();
		next_reg.first = next_reg_list[rand_sel].first;
		next_reg.second = next_reg_list[rand_sel].second;
	}
	else{
		next_reg.first = entropy_list.back().first;
		next_reg.second = entropy_list.back().second;
	}
	return next_reg;
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

