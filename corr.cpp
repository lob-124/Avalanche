#include <vector>
#include <complex>
#include <list>
#include <tuple>
#include <iostream>
#include <fstream>
#include <random>
#include <time.h>
#include <omp.h>
#include "Chain.hxx"
#include "Mapping.hxx"



using namespace std;

vector<tuple<int,int>> get_hoppings(int L,int N,int L_b){
    // ****
    //      Computes & returns a vector of matrix elements for the hoppings between the 
    //          right & left ends of the bath  


	vector<tuple<int,int>> hoppings_R;
	vector<tuple<int,int>> hoppings_L;

	int start_point = ((1 << L_b) - 1) << (L - L_b);
	int right_end = 1 << (L-L_b);
	int left_end = 1 << (L-1);

	for(int i= L-L_b; i > 0; i--){

		int init_state = start_point;

        //Construct the inital site, by placing the remaining particles as far left as possible (while still 
        //	permitting the desired hopping)
        int to_place = N - L_b;
        int start_position = L - L_b;
        for(int j=start_position; j > 0; j--){
        	if (j==i){ continue; }
         	else{
         		start_position = j;
         		init_state += (1 << (j-1));
         		to_place--;
         		if(to_place == 0){ break; }
         	}
        }


        //Construct the final state hopping from the right end
        int final_state_R = init_state - right_end + (1 << (i-1));

        //Construct the final state hopping from the left end
        int final_state_L = init_state - left_end + (1 << (i-1));

        //Append the (initial,final) pairs to the result
        hoppings_R.push_back(make_tuple(init_state,final_state_R));
        hoppings_L.push_back(make_tuple(init_state,final_state_L));

        }

        //Concatenate the results and return
        hoppings_R.insert(hoppings_R.end(),hoppings_L.begin(),hoppings_L.end());
        return hoppings_R;

}


vector<tuple<int,int>> get_indices(vector<tuple<int,int>> hoppings, vector<int> basis){

	vector<tuple<int,int>> res;

	//Walk through the array of hoppings, finding the index in the basis of each hopping
	for(auto tup: hoppings){
		int first_index = find(basis.begin(), basis.end(), get<0>(tup)) - basis.begin();
		int second_index = find(basis.begin(), basis.end(), get<1>(tup)) - basis.begin();

		res.push_back(make_tuple(first_index,second_index));
	}


	return res;
}


double find_g(Chain c, double epsilon, double g_step, double max_g, double tol = 1e-8){
/*
 *    Finds the critical g for a given energy density
 */
    for(double curr_g = g_step; curr_g <= max_g; curr_g += g_step){
        c.new_g(curr_g);
        vector<complex<double>> spectrum = c.spectrum();
        double E_min = real(spectrum[0]);
        double E_max = real(spectrum[c.dim-1]);
        double best_dist = 1.0;
        bool is_complex = false;
        for(complex<double> E: spectrum){
            double curr_eps = (real(E)-E_min)/(E_max - E_min);
            double curr_dist = abs(epsilon - curr_eps);            
            if(curr_dist < best_dist){            
                best_dist = curr_dist;                
                if(abs(imag(E)) > tol){ is_complex = true; }                   
                else{ is_complex = false; }
            }
        }
        if(is_complex){ return curr_g; }
    }
 
    return 0.0;
}




int main(int argc, char* argv[]){

if (argc == 1){
    	cout << "Usage: L N L_b W U g <realizations> <iterations> <threads> <path> <tol (op)>" << endl;
    	return 0;
	}

	//Read in command line parameters
	int L = stoi(argv[1]);
	int N = stoi(argv[2]);
	int L_b = stoi(argv[3]);
	double W = stod(argv[4]);
	double U = stod(argv[5]);
	double g = stod(argv[6]);
        int realizations = stoi(argv[7]);
	int iterations = stoi(argv[8]);
	int num_threads = stoi(argv[9]);
	string path = string(argv[10]);

	double tol = 1e-6;
	if(argc == 12){
		tol = stod(argv[11]);
	}

        //Parameters for the critical g search
        double max_g = 2.0;
        double delta_g = .01;
        double eps = 0.5;

        //Initialize the seeds
        srand(time(NULL));
        vector<unsigned int> seeds;
        for(int i=0; i < num_threads; i++){ seeds.push_back((unsigned int) rand()); }

	//Initialize the chains
	vector<Chain> chains;
	for(int i=0; i < num_threads; i++){ chains.push_back(Chain(L,N,W,U,g,seeds[i])); }

	//Open the output file and write some parameters to it 
	ofstream outfile(path,ios::out|ios::binary);
	outfile.write((char*)&realizations,sizeof(int));
	outfile.write((char*)&num_threads,sizeof(int));
	for(unsigned int seed: seeds){ outfile.write((char*)&seed, sizeof(unsigned int)); }


	//Mask for the "bath"
	int bath_mask = ((1 << L_b) - 1) << (L-L_b);
	//Basis in the fock space
	vector<int> basis = chains[0].basis;

	//Vector to store the data & record the thread numbers (so that we can compute the critical g's later for the same
	//	disorder realizations)
	vector<double> largest_amps(realizations);
	//vector<double> largest_amps_g(realizations);
        vector<vector<double>> data(realizations);
        //vector<vector<double>> data_g(realizations);
	vector<int> thread_numbers(realizations);
        //vector<double> g_cs(realizations);

	//Get the indices indicating the matrix elements we want to record
        vector<tuple<int,int>> hoppings = get_hoppings(L,N,L_b);
        vector<tuple<int,int>> indices = get_indices(hoppings,basis);


	//Set the number of threads
	omp_set_dynamic(0);
	omp_set_num_threads(num_threads); 


	//**** Main for loop ****//
	#pragma omp parallel for 
	for(int i=0; i < realizations; i++){
	    //Get and record the thread number
	    int thread_num = omp_get_thread_num();
	    thread_numbers[i] = thread_num;

	    //Extract the Hamiltonian
	    cx_mat H_prime = conv_to<cx_mat>::from(chains[thread_num].H);

	    //Peform the iterated displacement transformations
	    tuple<vector<int>,vector<int>,double> largest;
            for(int j=0; j < iterations; j++){
    	        largest = find_largest_bath(H_prime,basis,bath_mask,L);
    	        if(get<2>(largest) < tol){ break; }

    	        H_prime = remove_term(H_prime,get<0>(largest),get<1>(largest),basis);
    	        filter_zeroes(H_prime);
    	    } 

           //Record the largest amplitude (to track convergence)
           largest_amps[i] = get<2>(largest);
    	    
           //Record the matrix elements
    	   for(auto pair: indices){ data[i].push_back(abs(H_prime(get<0>(pair), get<1>(pair)))); }
 

           //****
           //     Now, find the ciritcal g_c at the given energy density, and repeat for that g
           //****
           // double g_c = find_g(chains[thread_num],eps,delta_g,max_g);
            //chains[thread_num].new_g(g_c);
            //g_cs[i] = g_c;

            //Extract the Hamiltonian
	    //H_prime = conv_to<cx_mat>::from(chains[thread_num].H);

	    //Peform the iterated displacement transformations
            //for(int j=0; j < iterations; j++){
    	      //  largest = find_largest_bath(H_prime,basis,bath_mask,L);
    	        //if(get<2>(largest) < tol){ break; }

    	        //H_prime = remove_term(H_prime,get<0>(largest),get<1>(largest),basis);
    	        //filter_zeroes(H_prime);
    	    //} 

           //Record the largest amplitude (to track convergence)
           //largest_amps_g[i] = get<2>(largest);
    	    
           //Record the matrix elements
    	    //for(auto pair: indices){ data_g[i].push_back(abs(H_prime(get<0>(pair), get<1>(pair)))); }
 



    	    //Generate a new disorder realization
     	    //chains[thread_num].new_g(0.0);
            chains[thread_num].new_disorder();

	}

        //Write out the largest amplitude (to track convergence)
        for(double amp: largest_amps){ outfile.write((char*)&amp, sizeof(amp)); }
        //for(double amp: largest_amps_g){ outfile.write((char*)&amp, sizeof(amp)); }


	//Write out the thread numbers (so that we can compute the critical g's later for the same
	//	disorder realizations)
	for(int thread_num: thread_numbers){ outfile.write((char*)&thread_num, sizeof(int)); }

        //Write out the critical g_cs
        //for(double g_c: g_cs){ outfile.write((char*)&g_c, sizeof(double)); }
 
	//Write out the amplitudes
	for(vector<double> amplitudes: data ){
		for(double d: amplitudes){ outfile.write((char*)&d, sizeof(double)); }
	}
	//for(vector<double> amplitudes: data_g ){
	//	for(double d: amplitudes){ outfile.write((char*)&d, sizeof(double)); }
	//}


	outfile.close();
}

