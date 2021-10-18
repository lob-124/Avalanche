#ifndef MAPPING
#define MAPPING
 
#include <complex>
#include <vector>
#include <tuple>
#include <algorithm>
#include <armadillo>
//#include <iostream>
//#include <fstream>
#include "generate-matrices.hxx"

using namespace std;
using namespace arma;


vector<tuple<vector<int>,vector<int>>> gen_hop_rec(vector<int>sites, vector<int> sites_1 , vector<int> sites_2, int num_bath_sites, int num_hopping, int start_index, vector<int> hopping_distances){
	// ****
	//		Recursively generate all possible hopping pairs for the given list of hopping distances.
	//		Returns a vector of tuples, each tuple containing two vectors, each containing the sites being hopped to/from
	// ****

	//The recursive step
	//Essentially, we add (to,from) site pairs (to to sites_1 & from to sites_2), with to & from separated by
	//	the current hopping distance (hopping_distances.back()) and make a recursive call with the new sites arrays, one
	//	fewer particle hopping, and the next hopping distance
	if (num_hopping > 0){
		int L_loc = sites.size() - num_bath_sites;
		vector<tuple<vector<int>,vector<int>>> res;	

		//Loop through the sites
		for(int i=start_index; i < sites.size(); i++){
			int site_1 = sites[i]; 
			int site_2 = sites[(i + hopping_distances.back() - num_bath_sites) % L_loc + num_bath_sites];

			//Check that these sites aren't already included in either of the sites arrays
			if ( (find(sites_1.begin(),sites_1.end(),site_1) == sites_1.end()) && (find(sites_2.begin(),sites_2.end(),site_1) == sites_2.end()) ){
				if ( (find(sites_1.begin(),sites_1.end(),site_2) == sites_1.end()) && (find(sites_2.begin(),sites_2.end(),site_2) == sites_2.end()) ){
					//Add the sites to the site arrays and make the recursive call
					vector<int> sites_1_new(sites_1); vector<int> sites_2_new(sites_2);
					sites_1_new.push_back(site_1); sites_2_new.push_back(site_2);
					vector<int> distances_new(hopping_distances);
					distances_new.pop_back();					
					vector<tuple<vector<int>,vector<int>>> temp = gen_hop_rec(sites, sites_1_new, sites_2_new, num_bath_sites, num_hopping-1,i+1,distances_new);
					res.insert(res.end(),temp.begin(),temp.end());
				}
			} 
		}
		return res;

	}
	//Base case - there are no more particles hopping, so just retrun the sites arrays
	else { return { make_tuple(sites_1,sites_2) }; }

}


tuple<int,int,double> find_largest(cx_mat H){
	// ****
	//		Finds the largest (by absolute value) off-diagonal matrix element
	//			of the given matrix
	// ****

	double max_amplitude = 0.0;
	int first_index; int second_index;

	//Loop through all off-diagonal elements and keep track of the largest (by absolute
	//	value) one
	for(int i=0; i < H.n_rows; i++){
		for(int j=0; j < H.n_cols; j++){
			if ( i==j){ continue; }

			double amplitude = abs(H(i,j));

			//if(isnan(amplitude)) { cout << "WOAH  DUDE" << endl; return make_tuple(-1,-1,0.0);}
			if (amplitude > max_amplitude){
				max_amplitude = amplitude;
				first_index = i;
				second_index = j;
			}
		}
	}

	return make_tuple(first_index,second_index,max_amplitude);
}


tuple<vector<int>,vector<int>> get_sites(int basis_elem_1, int basis_elem_2, int L){
	// ****
	//		Given two basis elements corresponding to a matrix element H_{ij},
	//			determines the sites the hopping is between
	// ****

	//Extract masks for the sites involved in the hopping
	int all_sites = basis_elem_1 ^ basis_elem_2;
	int sites_mask_1 = basis_elem_1 & all_sites;
	int sites_mask_2 = basis_elem_2 & all_sites;

	//Split up the site masks into lists of individual sites
	//This is essentially just a decomposition into powers of two
	vector<int> sites_1; vector<int> sites_2;
	for(int i = (1 << (L-1)); i >= 1; i = i >> 1){
		if(sites_mask_1 & i){ sites_1.push_back(i); }
		else if(sites_mask_2 & i){ sites_2.push_back(i); }
	}

	return make_tuple(sites_1,sites_2);

}



tuple<vector<int>,vector<int>,double> find_largest_bath(cx_mat H, vector<int> basis, int bath_mask, int L){
	// ****
	//		Finds the largest (by absolute value) off-diagonal matrix element
	//			of the given matrix pertaining to hopping between non-bath sites. 
	//		Returns the sites involved in said hopping
	// ****

	double max_amplitude = 0.0;
	vector<int> sites1; vector<int> sites2;

	//Loop through all off-diagonal elements and keep track of the largest (by absolute
	//	value) one
	for(int i=0; i < H.n_rows; i++){
		for(int j=0; j < H.n_cols; j++){
			if ( i==j){ continue; }

			double amplitude = abs(H(i,j));

			if (amplitude > max_amplitude){
				//If the current matrix element has a smaller amplitude, check if it 
				// involves any bath sites
				tuple<vector<int>,vector<int>> sites = get_sites(basis[i],basis[j],L);
				vector<int> sites1_temp = get<0>(sites);
				vector<int> sites2_temp = get<1>(sites);
				bool is_not_bath_hopping = true;
				for(int k=0; k < sites1_temp.size(); k++){
					if((sites1_temp[k] & bath_mask) || (sites2_temp[k] & bath_mask)){
						is_not_bath_hopping = false;
						break;
					}
				}

				//If it doesn't, update the largest element & element
				if(is_not_bath_hopping){
					max_amplitude = amplitude;
					sites1 = sites1_temp;
					sites2 = sites2_temp;
				}
			}
		}
	}

	return make_tuple(sites1,sites2,max_amplitude);
}


cx_mat remove_term(cx_mat H, vector<int> sites_1, vector<int> sites_2, vector<int> basis){
	// ****
	//		Removes an off-diagonal term of the form
	//					(\prod_i c_i^{\dag})(\prod_j c_j) + \alpha(\prod_j c_j^{\dag})(\prod_i c_i)
	//			from the given Hamiltonian
	// ****

	//Bit masks to identify the hopping sites
	int all_sites_mask = 0;
	int sites_1_mask = 0;
	int sites_2_mask = 0;
	for(int i=0; i < sites_1.size(); i++){ 
		all_sites_mask += sites_1[i] + sites_2[i]; 
		sites_2_mask += sites_2[i];
		sites_1_mask += sites_1[i];
	}


	//alpha, delta V, and t are all diagonal in the subspace spanned by the sites involved in the hopping. 
	//	Hence, we need to enumerate all the possible states in this subspace to construct the alpha, delta V, and t
	//	matrices

	//Generate a list of all the sites involved in the hopping
	vector<int> all_hopping_sites(sites_1);
	all_hopping_sites.insert(all_hopping_sites.end(),sites_2.begin(),sites_2.end());
	sort(all_hopping_sites.begin(),all_hopping_sites.end());

	//Generate the basis for the subspace spanned by the appropriate number of particles distributed among these sites
	vector<int> subspace_basis = generate_combos(all_hopping_sites,0,-1,all_hopping_sites.size()/2,all_hopping_sites.size()); 


	//****
	//		Compute alpha - the "ratio" between right & left hopping terms
	//		Note that in the non-hermitian interacting case, this may be an operator!
	//****
	cx_vec alpha(basis.size(),fill::ones);	
	for (int i=0; i < basis.size(); i++){
		int basis_elem_1 = basis[i];
		complex<double> mat_elem = 0.0;
		for (int j = i+1; j < basis.size(); j++){
			int basis_elem_1 = basis[i];
			int basis_elem_2 = basis[j];

			//Check if these basis elements include the desired hopping, and extract alpha if they do
			if (((basis_elem_1 ^ basis_elem_2) == all_sites_mask) && ((basis_elem_2 & all_sites_mask) == sites_2_mask)) { 
				mat_elem = H(j,i)/H(i,j);
				if(::isnan(real(mat_elem)) || ::isnan(imag(mat_elem))){ mat_elem = 1.0;}
				break; 
			}
			else if (((basis_elem_1 ^ basis_elem_2) == all_sites_mask) && ((basis_elem_1 & all_sites_mask) == sites_2_mask)) {
				mat_elem = H(i,j)/H(j,i);
				if(::isnan(real(mat_elem)) || ::isnan(imag(mat_elem))){ mat_elem = 1.0;}
				break;
			}		
		}  
		if(mat_elem == complex<double>(0.0,0.0)){ continue; }	
	
		//Fill in the entire subspace spanned by arbitrary occupation of the sites
		//	involved in the hopping, while keeping the occupation of the "other sites"
		//	fixed
		int other_sites = basis_elem_1 & (~all_sites_mask);
		for(int elem : subspace_basis){
			int index = find(basis.begin(), basis.end(), other_sites + elem) - basis.begin();
			alpha(index) = mat_elem;
		}	

	}

	
	//****
	//        Construct delta V
	//****
	cx_vec delta_V(basis.size(),fill::zeros);
	for (int j=0; j < basis.size(); j++){
		int basis_elem = basis[j];
		complex<double> mat_elem;

		//Check if we are either in the initial or final state of the hopping
		if ((basis_elem & all_sites_mask) == sites_1_mask){ mat_elem = H(j,j); }	//Sites 1 occupied, sites 2 unoccupied
		else if ((basis_elem & all_sites_mask) == sites_2_mask){ mat_elem = -H(j,j); }	//Sites 2 occupied, sites 1 unoccupied
		else { continue; }

		//Fill in the entire subspace spanned by arbitrary occupation of the sites
		//	involved in the hopping, while keeping the occupation of the "other sites"
		//	fixed
		int other_sites = basis_elem & (~all_sites_mask);
		for(int elem : subspace_basis){
			int index = find(basis.begin(), basis.end(), other_sites + elem) - basis.begin();
			delta_V(index) += mat_elem;
		}
	} //end delta V loop


	//****
	//       Construct t
	//****
	cx_vec t(basis.size(),fill::zeros);
	for (int i=0; i < basis.size(); i++){
		for (int j = i+1; j < basis.size(); j++){
			int basis_elem_1 = basis[i];
			int basis_elem_2 = basis[j];
			complex<double> mat_elem;

			//Check if these basis elements include the desired hopping, and extract the diagonal part if they do
			if ((basis_elem_1 ^ basis_elem_2) == all_sites_mask){
				if((basis_elem_2 & all_sites_mask) == sites_2_mask) { mat_elem = H(i,j); }	//H(i,j) ~ c_k^{\dag}c_l
				else if((basis_elem_1 & all_sites_mask) == sites_2_mask) { mat_elem = H(j,i); }	//H(i,j) ~ \alpha c_l^{\dag} c_k
				else { continue; }
			}
			else { continue; }

			//Now we need to count the number of particles between the hopping sites for each of the hoppings, to get 
			//	the signs of the matrix elements right.
			int num = 0;
			for(int i=0; i < sites_1.size(); i++){
				if(sites_1[i] > sites_2[i]){		//site 2 right of site 1
					for(int curr_site = (sites_1[i] >> 1); curr_site > sites_2[i]; curr_site = curr_site >> 1){ if(basis[i] &  curr_site) { num ++; } }
				}      
				else{		//site 2 left of site 1
					for(int curr_site = (sites_1[i] << 1); curr_site < sites_2[i]; curr_site = curr_site << 1){ if(basis[i] &  curr_site) { num ++; } } 
				}
			} 
			double sign = pow(-1.0,num);
			mat_elem *= sign;  

			//Fill in the entire subspace spanned by arbitrary occupation of the sites
			//	involved in the hopping, while keeping the occupation of the "other sites"
			//	fixed
			int other_sites = basis_elem_1 & (~all_sites_mask);
			for (int elem : subspace_basis){
				int index = find(basis.begin(), basis.end(), other_sites + elem) - basis.begin();
				t(index) = mat_elem;
			}
			break;
		} 
	}//end t loop


	//****
	//      Construct the sine/cosine terms appearing in the displacement transforms
	//****
	cx_vec cos_lambda_m1(basis.size()); cx_vec sin_lambda(basis.size());
	for(int i=0; i<basis.size(); i++){
		//Compute \cos(2\lambda) & \sin(2\lambda)
		complex<double> cos_2lambda = delta_V(i)/sqrt(pow(delta_V(i),2) + alpha(i)*pow(2.0*t(i),2));
		complex<double> sin_2lambda = 2.0*sqrt(alpha(i))*t(i)/sqrt(pow(delta_V(i),2) + alpha(i)*pow(2.0*t(i),2));
		
		//If t_{ij} = \delta_V = 0.0, then set set \sin(\lambda) = 0.0, \cos(\lambda) = I
		if (::isnan(abs(cos_2lambda))) { cos_lambda_m1(i) = 0.0; sin_lambda(i) = 0.0; continue; }

		//Otherwise, use half-angle formulae to compute \sin(\lambda) & \cos(\lambda) - I
		cos_lambda_m1(i) = sqrt((1.0+cos_2lambda)/2.0)-1.0;
		sin_lambda(i) = sqrt((1.0-cos_2lambda)/2.0);

		//Check the sign of sin(\lambda) by computing \tan(2\lambda) and comparing against the expected value
		//	of 2\sqrt{\alpha} t_{ij}/\Delta_{ij} V
		complex<double> tan_2lambda = 2.0*sin_lambda(i)*(cos_lambda_m1(i)+1.0)/(pow(cos_lambda_m1(i)+1.0,2.0) - pow(sin_lambda(i),2.0));
		complex<double> tan_2lambda_expected = 2.0*sqrt(alpha(i))*t(i)/delta_V(i);
		if(real(tan_2lambda/tan_2lambda_expected) < 0.0){ sin_lambda(i) *= -1.0; }
	}


	//****
	//        Construct the projector P and the hopping matrix J
	//****
	vec P(basis.size(),fill::zeros);
	cx_mat J(basis.size(),basis.size(),fill::zeros);
	for(int i=0; i < basis.size(); i++){
		//Check if this state has the right occupations for a hopping to occurr
		if( ((basis[i] & all_sites_mask) == sites_1_mask) || ((basis[i] & all_sites_mask) == sites_2_mask) ){ P(i) = 1;}

		//Loop through all other basis sites, determine if there's a hopping between basis[i] & it, and compute
		//	the matrix element
		for(int j=i+1; j < basis.size(); j++){
			//Check if there's a hopping between basis[i] & basis[j]
			//If there is, fill in the matrix elements accordingly. If there's not, continue to the next basis element
			if (((basis[i] ^ basis[j]) == all_sites_mask) && ((basis[j] & all_sites_mask) == sites_2_mask)) { 
				J(i,j) = complex<double>(1.0,0.0); 
				J(j,i) = -alpha(i);
			}
			else if (((basis[i] ^ basis[j]) == all_sites_mask) && ((basis[i] & all_sites_mask) == sites_2_mask)) { 
				J(i,j) = -alpha(i);
				J(j,i) = complex<double>(1.0,0.0);
			}
			else continue;

			//Now we need to count the number of particles between the hopping sites for each of the hoppings, to get 
			//	the signs of the matrix elements right.
			int num = 0;
			for(int i=0; i < sites_1.size(); i++){
				if(sites_1[i] > sites_2[i]){		//site 2 right of site 1
					for(int curr_site = (sites_1[i] >> 1); curr_site > sites_2[i]; curr_site = curr_site >> 1){ if(basis[i] &  curr_site) { num ++; } }
				}      
				else{		//site 2 left of site 1
					for(int curr_site = (sites_1[i] << 1); curr_site < sites_2[i]; curr_site = curr_site << 1){ if(basis[i] &  curr_site) { num ++; } } 
				}
			} 

			//Multiply by the sign factor
			double sign = pow(-1.0,num);
			J(i,j) *= sign;  
			J(j,i) *= sign;
			break; 			//Each basis state hops to exactly one other basis state, so once we've found it we break
		}
	} //end P & J loop


	//****
	//        Construct the transformation matrices & perform the transform
	//****
	
	//We first construct the operator sin(lambda)/\sqrt{\alpha}
	cx_vec sin_over_sqrt_alpha(basis.size());
	for(int i=0; i < basis.size(); i++){ 
		if(abs(alpha(i)) != 0.0){ sin_over_sqrt_alpha(i) = sin_lambda(i)/sqrt(alpha(i)); }
		else{ sin_over_sqrt_alpha(i) = 0.0; }
	}
	
	//Construct & apply the transformation matrices
	cx_mat S_plus = eye(basis.size(),basis.size()) + diagmat(cos_lambda_m1)*diagmat(P) + diagmat(sin_over_sqrt_alpha)*J;
	cx_mat S_minus = eye(basis.size(),basis.size()) + diagmat(cos_lambda_m1)*diagmat(P) - diagmat(sin_over_sqrt_alpha)*J;
	cx_mat result = S_plus*H*S_minus;
	
	return result;

} //end remove_term


vector<tuple<int,int>> expected_zeroes(vector<int> sites_1, vector<int> sites_2, vector<int> basis){
	// ****
	//		Returns a list of the (upper triangular) indices expcted to be (near) zero after removing 
	//			hopping between the given sets of sites
	// ****

	//Bit masks to identify the hopping sites
	int all_sites_mask = 0;
	int sites_1_mask = 0;
	int sites_2_mask = 0;
	for(int i=0; i < sites_1.size(); i++){ 
		all_sites_mask += sites_1[i] + sites_2[i]; 
		sites_2_mask += sites_2[i];
		sites_1_mask += sites_1[i];
	}

	//Loop through the pairs of basis elements and add if they represent the desired hopping
	vector<tuple<int,int>> res;
	for(int i=0; i < basis.size(); i++){
		for(int j=i+1; j < basis.size(); j++){
			if((basis[i] ^ basis[j]) == all_sites_mask){
				if (((basis[i] & sites_1_mask) == sites_1_mask) || ((basis[i] & sites_2_mask) == sites_2_mask)){	
					res.push_back(make_tuple(i,j));	
				}
			}
		}
	}

	return res;
} // end expected_zeroes


void filter_zeroes(cx_mat &H, double tol=1e-10){
	//****
	//		Zero out any matrix elements whose absolute value is less than tol
	//****
	for(int i = 0; i < H.n_rows; i++){
		for(int j = 0; j < H.n_cols; j++){
			if(abs(H(i,j)) < tol){ H(i,j) = 0.0; H(j,i) = 0.0; }
		}
	}
} //end filter_zeroes


#endif
