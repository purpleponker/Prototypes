#include "WFC.h"

//main function
int main(int argc, char *argv[]){
	wfc_mapping matrix;
	
	//seed random based on current time
	srand(time(NULL));
	
	//random entry point
	std::pair index = matrix.start();
	
	//bool flag for WFC complete
	bool wfc_valid= true;
	do{
		//collapse index
		bool do_prop = matrix.collapse(index.first, index.second);
		// collapse check
		// break on wfc contradiction
		if(!do_prop){
			wfc_valid = false;
			std::cout  << "WFC contradiction" << std::endl;
			break;
		}
		
		//propogate collapse
		wfc_valid = matrix.propagate(index.first, index.second);
		
		// calculate new entropy vals	
		index = matrix.entropy();
		// repeat untill controdiction or completion
		//change the check from entropy list to the uncollapsed list, either ending proper or on controdiction
		if(index == std::make_pair(-1,-1)){
			std::cout << "collapse selection contradiction" << std::endl;
			wfc_valid = false;
		}
	}while(wfc_valid);
	
	// show map
	matrix.display();
}
