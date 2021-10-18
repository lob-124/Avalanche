#include <vector>
#include <complex>
#include <list>
#include <numeric>
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
    //      Computes & returns a vector of matrix elements for the hoppings from the 
    //          right end of the bath  


	vector<tuple<int,int>> hoppings_R;
//	vector<tuple<int,int>> hoppings_L;

	int start_point = ((1 << L_b) - 1) << (L - L_b);
	int right_end = 1 << (L-L_b);
//	int left_end = 1 << (L-1);

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
//        int final_state_L = init_state - left_end + (1 << (i-1));

        //Append the (initial,final) pairs to the result
        hoppings_R.push_back(make_tuple(init_state,final_state_R));
//        hoppings_L.push_back(make_tuple(init_state,final_state_L));

        }

        //Concatenate the results and return
//        hoppings_R.insert(hoppings_R.end(),hoppings_L.begin(),hoppings_L.end());
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


vector<int> find_closest(vector<complex<double>> spectrum, vector<double> targets){
    /**
        Finds the indices of elements of spectrum whose real part is closest to each value in targets.
        We assume that both spectrum and target are sorted in ascending order (of real part, in the case of spectrum)
     **/

    int len = targets.size();
    vector<int> closest_elems(len);
    vector<double> distances(len,INFINITY); 

    //Walk through the spectrum, comparing the distance from the real part of each element
    //    to the targets
    int found_up_to = 0;
    for(int i=0; i < spectrum.size(); i++){
        //Walk through the targets, recording new elements as we get closer to each
        for(int j=found_up_to; j < len; j++){
                double dist = abs(real(spectrum[i]) - targets[j]);
            if (dist <= distances[j]){    //Update while we're getting closer
                closest_elems[j] = i;
                distances[j] = dist;
            }
            else{ found_up_to = j; }    //If the distance is now increasing, we've found the closest elements for this and preceding targets
        }

        if(found_up_to == len-1){ return closest_elems; }
    }

    return closest_elems;
}


vector<double> critical_g_vals(Chain c, vector<double> epsilons, double g_step, double g_max, double tol=1e-8){
    vector<double> results(epsilons.size(),0.0);                    //Vector of g_c's at which energy densities coalesce
    list<int> indices;                                              //List of indices of energy densities remaining to check
    for(int i=0; i<epsilons.size(); i++){ indices.push_back(i); }

    for(double curr_g = g_step; curr_g <= g_max; curr_g += g_step){
        //Find the spectrum for the current g
        c.new_g(curr_g);
        vector<complex<double>> spectrum = c.spectrum();
        
        //Find the energies in the spectra closest to the (remaining) target energy densities
        double E_min =  real(spectrum.front()); double E_max = real(spectrum.back());
        vector<double> targets;
        for(int i: indices){ targets.push_back(E_min + epsilons[i]*(E_max-E_min)); }
        vector<int> closest_vals_indices = find_closest(spectrum,targets);

        //For each of these energies that is complex, record the current g value, and remove the 
        //  index of the energy from the list of indices remaining to be found
        vector<int> indices_copy(indices.begin(),indices.end());
        for(int i=0; i<closest_vals_indices.size(); i++){
            if(abs(imag(spectrum[closest_vals_indices[i]])) > tol){
                results[indices_copy[i]] = curr_g;
                indices.remove(indices_copy[i]);
            }
    }
    
        if(indices.empty()){ return results; } 
    }
    
    return results;
}




int main(int argc, char* argv[]){

if (argc == 1){
    	cout << "Usage: L N L_b W U <eps_start> <eps_stop> <eps_step> <realizations> <iterations> <threads> <path> <tol (op)>" << endl;
    	return 0;
	}

	//Read in command line parameters
	int L = stoi(argv[1]);
	int N = stoi(argv[2]);
	int L_b = stoi(argv[3]);
	double W = stod(argv[4]);
	double U = stod(argv[5]);
        double eps_start = stod(argv[6]);
        double eps_stop = stod(argv[7]);
        double eps_step = stod(argv[8]);
        int realizations = stoi(argv[9]);
	int iterations = stoi(argv[10]);
	int num_threads = stoi(argv[11]);
	string path = string(argv[12]);

	double tol = 1e-6;
	if(argc == 14){
		tol = stod(argv[13]);
	}

        //Parameters for the critical g search
        double max_g = 2.0;
        double delta_g = .01;
        vector<double> epsilons;
        for(double eps = eps_start; eps <= eps_stop + eps_step/2; eps += eps_step){ epsilons.push_back(eps); } 

        //Initialize the seeds
        vector<int> seeds;
        srand(time(NULL));
        for(int i=0; i < num_threads; i++){ seeds.push_back((unsigned int) rand()); }
        
	//Initialize the open chains
	vector<Open_Chain> chains;
	for(int i=0; i < num_threads; i++){ chains.push_back(Open_Chain(L,N,W,U,seeds[i])); }
        
        //Initialize the rings
        vector<Chain> rings;
        for(int i=0; i < num_threads; i++){ rings.push_back(Chain(L,N,W,U,0,seeds[i])); }




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
        vector<vector<double>> hopping_amplitudes(realizations);
        vector<vector<double>> open_spectra(realizations);
        //vector<vector<double>> data_g(realizations);
	vector<int> thread_numbers(realizations);
        vector<vector<double>> g_cs(realizations);

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

     //        //Set the disorder distributions to be the same
     //        //chains[thread_num].set_disorder(rings[thread_num].on_site);

     //        //Extract the spectrum of the open chain
            open_spectra[i] = chains[thread_num].spectrum();

	    // //Extract the Hamiltonian
	    cx_mat H_prime = conv_to<cx_mat>::from(chains[thread_num].H);

	    // //Peform the iterated displacement transformations
	    tuple<vector<int>,vector<int>,double> largest;
            double prev_largest = 1.0;
            for(int j=0; j < iterations; j++){
    	        largest = find_largest_bath(H_prime,basis,bath_mask,L);
    	        double largest_amp = get<2>(largest);
                if((largest_amp < tol) || ((largest_amp == prev_largest) && (largest_amp < 1.0))){ break; }

                prev_largest = largest_amp;
    	        H_prime = remove_term(H_prime,get<0>(largest),get<1>(largest),basis);
    	        filter_zeroes(H_prime);
    	    } 

           //Record the largest amplitude (to track convergence)
           largest_amps[i] = get<2>(largest);
    	    
           //Record the matrix elements
    	    for(auto pair: indices){ hopping_amplitudes[i].push_back(abs(H_prime(get<0>(pair), get<1>(pair)))); }
 

           //****
           //     Now, find the ciritcal g_c at the given energy density, and repeat for that g
           //****
           vector<double> g_c_vals = critical_g_vals(rings[thread_num],epsilons,delta_g,max_g);
           //chains[thread_num].new_g(g_c);
           g_cs[i] = g_c_vals;

           /*
            //Extract the Hamiltonian
	    H_prime = conv_to<cx_mat>::from(chains[thread_num].H);

	    //Peform the iterated displacement transformations
            for(int j=0; j < iterations; j++){
    	        largest = find_largest_bath(H_prime,basis,bath_mask,L);
    	        if(get<2>(largest) < tol){ break; }

    	        H_prime = remove_term(H_prime,get<0>(largest),get<1>(largest),basis);
    	        filter_zeroes(H_prime);
    	    } 

           //Record the largest amplitude (to track convergence)
           largest_amps_g[i] = get<2>(largest);
    	    
           //Record the matrix elements
    	    for(auto pair: indices){ data_g[i].push_back(abs(H_prime(get<0>(pair), get<1>(pair)))); }
 

            */

    	    //Generate a new disorder realization
     	   // chains[thread_num].new_g(0.0);
            chains[thread_num].new_disorder();
            rings[thread_num].new_disorder();

	}

        //Write out the largest amplitude (to track convergence)
        for(double amp: largest_amps){ outfile.write((char*)&amp, sizeof(amp)); }
        //for(double amp: largest_amps_g){ outfile.write((char*)&amp, sizeof(amp)); }


	//Write out the thread numbers (so that we can compute the critical g's later for the same
	//	disorder realizations)
	for(int thread_num: thread_numbers){ outfile.write((char*)&thread_num, sizeof(int)); }

        //Write out the critical g_cs
        for(vector<double> _gcs: g_cs){ 
           for(double g_c: _gcs) {
                outfile.write((char*)&g_c, sizeof(double)); 
           }
        }

        //Write out the spectra (NOT INCLUDED PREVIOUSLY)
        for(vector<double> _spectrum: open_spectra){
            for(double d: _spectrum){ outfile.write((char*)&d, sizeof(double)); }
        }

	//Write out the amplitudes
	for(vector<double> amplitudes: hopping_amplitudes ){
	    for(double d: amplitudes){ outfile.write((char*)&d, sizeof(double)); }
	}
	//for(vector<double> amplitudes: data_g ){
	//	for(double d: amplitudes){ outfile.write((char*)&d, sizeof(double)); }
	//}


	outfile.close();
}

