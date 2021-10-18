#ifndef GENERATE_MATRICES
#define GENERATE_MATRICES

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include <iostream>

using namespace std;


/* ***
 *
 *    GENERATE THE BASIS
 *
 *** */
vector<int> generate_combos(vector<int> sites,int site_val,int site,int n,int L){

if (n > 0){
    vector<int> res;
    for(int i=site+1; i < L+1-n; i++){
        vector<int> combos = generate_combos(sites,sites[i],i,n-1,L);
        for(vector<int>::iterator combo = combos.begin(); combo != combos.end(); ++combo){
            res.push_back(site_val + *combo); 
        }
    }
    return res;
}
else{
    vector<int> res = {site_val};
    return res;
}
} //end generate_combos


vector<int> generate_basis(int L){
  // ****
  //      Generate the full basis
  // ****
  vector<int> basis;
  for(int i=0; i < pow(2,L)-1; i++){ basis.push_back(i); }
 
  return basis;
}


/* ***
 *
 *    GENERATE HOPPING TERMS   
 *
 * ***/
vector<double> generate_hopping(vector<int> sites, vector<int> basis, vector<tuple<double,double>> amplitudes, double alpha, double g, int dim, int L, int N, int L_bath){
    vector<double> hopping(dim*dim,0);
    int left_bath_end = sites[0];
    int right_bath_end = sites[L_bath-1];

    //Loop over all possible non-diagonal matrix elements
    for (int i=0; i<dim; i++){
    	for(int j=0 ;j<dim;j++){
        	
        	//Check if the two basis states are related by a hopping 
           	for(int k=0; k<L-L_bath; k++){
           		
                //****
                //    Check for a hopping to/from the right end of the bath
                //****
               	if((basis[i] ^ basis[j]) == (right_bath_end + sites[L_bath+k])){
               		
               		  //Check how many particles are between the two hopping sites & compute sign
               		  //	from anticommuting past the other particle's creation ops
               		  int num = 0;
               		  for(int l =L_bath; l<L_bath + k; l++){ if(sites[l] & basis[i]){ num++; } }
               		  int sign = pow(-1,num);	

               		  //Check whether i->j is left or right hopping and act accordingly
               		  if(basis[i] & right_bath_end){	//left hopping
                        hopping[dim*i+j] = sign*get<1>(amplitudes[k])*pow(alpha,k+1)*exp((k+1)*g);
               		  }
               		else{	//right hopping
               			hopping[dim*i+j] = sign*get<1>(amplitudes[k])*pow(alpha,k+1)*exp(-(k+1)*g);
               		}

               		break;
               	}
                //****
                //    Check for a hopping to/from the left end of the bath
                //****
                else if((basis[i] ^ basis[j]) == (left_bath_end + sites[L_bath+k])){
                  
                    //Check how many particles are between the two hopping sites & compute sign
                    //  from anticommuting past the other particle's creation ops
                    int num = 0;
                    for(int l =1; l<L_bath + k; l++){ if(sites[l] & basis[i]){ num++; } }
                    int sign = pow(-1,num);  

                    //Check whether this is left or right hopping and act accordingly
                    if(basis[j] & left_bath_end){  //left hopping
                        hopping[dim*i+j] = sign*get<0>(amplitudes[k])*pow(alpha,L-L_bath-k)*exp((L-L_bath-k)*g);
                    }
                    else{ //right hopping
                        hopping[dim*i+j] = sign*get<0>(amplitudes[k])*pow(alpha,L-L_bath-k)*exp(-(L-L_bath-k)*g);
                    }

                    break;
                }
           	} //end loop over sites
   		}    
	}
    return hopping;

} //end generate_hopping




/* ***
 *
 *    GENERATE LEFT HOPPING TERMS   
 *
 * ***/
 
//****
//    No random amplitudes
//****
vector<int> generate_left_hopping(vector<int> sites, vector<int> basis, int dim, int L, int N, int L_b){
    vector<int> left_hopping(dim*dim,0);
    int left_hopping_masks [L-1][2];
    for(int i=0; i < L-L_b; i++){
        left_hopping_masks[i][0] = sites[i+L_b];
        left_hopping_masks[i][1] = sites[i+L_b+1];
    }
    int left_end_mask [2] = {sites[L-1] , sites[L_b]};    //Note the change sites[0] -> sites[L_b] here to "skip over" the bath at site 0

    for (int i=0; i<dim; i++){
    for(int j=0;j<dim;j++){
        if (i == j){ continue; }
        if(((basis[i] ^ basis[j]) == (left_end_mask[1] + left_end_mask[0])) && (basis[i] & left_end_mask[0])){
            left_hopping[dim*i+j] = pow(-1,N-1);
        }
        else{
            for(int k=0; k<L-1; k++){
                if(((basis[i] ^ basis[j]) == (left_hopping_masks[k][1] + left_hopping_masks[k][0])) && (basis[i] & left_hopping_masks[k][0])){
                    left_hopping[dim*i+j] = 1;
                    break;   
                }
            }
        }
    }
    }
    return left_hopping; 
} //end generate_left_hopping



//****
//    Random amplitudes (requires L_b > 0)
//****
vector<double> generate_left_hopping(vector<int> sites, vector<int> basis, vector<double> amplitudes, int dim, int L, int N, int L_b){
    vector<double> left_hopping(dim*dim,0);
    int left_hopping_masks [L-L_b][2];
    for(int i=0; i < L-L_b; i++){  //NB: i starts at -1 to incorporate bath into ring
        left_hopping_masks[i][0] = sites[i+L_b-1];
        left_hopping_masks[i][1] = sites[i+L_b];
    }
    int left_end_mask [2] = {sites[L-1] , sites[0]};



    for (int i=0; i<dim; i++){
    for(int j=0;j<dim;j++){
        if (i == j){ continue; }
        if(((basis[i] ^ basis[j]) == (left_end_mask[1] + left_end_mask[0])) && (basis[i] & left_end_mask[0])){
            if(N){ left_hopping[dim*i+j] = amplitudes.back()*pow(-1,N-1); }   //fixed particle number sector
            else{     //arbitrary particle number
              int num = 0;
              for(int k=1; k < L-1; k++){
                if(sites[k] & basis[i]){ num++; }
              }
              left_hopping[dim*i+j] = amplitudes.back()*pow(-1,num);
            }
        }
        else{
            for(int k=0; k<L-L_b; k++){
                if(((basis[i] ^ basis[j]) == (left_hopping_masks[k][1] + left_hopping_masks[k][0])) && (basis[i] & left_hopping_masks[k][0])){
                    left_hopping[dim*i+j] = amplitudes[k];
                    break;   
                }
            }
        }
    }
    }
    return left_hopping; 
} //end generate_left_hopping


//****
//    No bath
//****
vector<double> generate_left_hopping(vector<int> sites, vector<int> basis, vector<double> amplitudes, int dim, int L, int N){
    vector<double> left_hopping(dim*dim,0);
    int left_hopping_masks [L-1][2];
    for(int i=0; i < L-1; i++){
        left_hopping_masks[i][0] = sites[i];
        left_hopping_masks[i][1] = sites[i+1];
    }
    int left_end_mask [2] = {sites[L-1] , sites[0]};    

    for (int i=0; i<dim; i++){
    for(int j=0;j<dim;j++){
        if (i == j){ continue; }
        if(((basis[i] ^ basis[j]) == (left_end_mask[1] + left_end_mask[0])) && (basis[i] & left_end_mask[0])){
            left_hopping[dim*i+j] = amplitudes.back()*pow(-1,N-1);
        }
        else{
            for(int k=0; k<L-1; k++){
                if(((basis[i] ^ basis[j]) == (left_hopping_masks[k][1] + left_hopping_masks[k][0])) && (basis[i] & left_hopping_masks[k][0])){
                    left_hopping[dim*i+j] = amplitudes[k];
                    break;   
                }
            }
        }
    }
    }
    return left_hopping; 
} //end generate_left_hopping


//****
//    No bath (open chain)
//****
vector<double> generate_left_hopping_open(vector<int> sites, vector<int> basis, vector<double> amplitudes, int dim, int L, int N){
    vector<double> left_hopping(dim*dim,0);
    int left_hopping_masks [L-1][2];
    for(int i=0; i < L-1; i++){
        left_hopping_masks[i][0] = sites[i];
        left_hopping_masks[i][1] = sites[i+1];
    }

    for (int i=0; i<dim; i++){
        for(int j=0;j<dim;j++){
            if (i == j){ continue; }
            for(int k=0; k<L-1; k++){
                if(((basis[i] ^ basis[j]) == (left_hopping_masks[k][1] + left_hopping_masks[k][0])) && (basis[i] & left_hopping_masks[k][0])){
                    left_hopping[dim*i+j] = amplitudes[k];
                    break;   
                }
            }
        }
    } 
    return left_hopping; 
} //end generate_left_hopping_open




/* ***
 *
 *    GENERATE DIAGONAL TERMS
 *
 * ***/
vector<double> generate_diagonal(vector<int> sites, vector<int> basis, vector<double> on_site, double U, int L, int L_b,int dim){
    vector<double> diagonal(dim,0.0);
    for (int i =0; i<dim; i++){
        //Check for an interaction between rightmost bath site & first non-bath site
        if((basis[i] & sites[L_b-1]) && (basis[i] & sites[L_b])){ diagonal[i] += U; }
        
        //Check the rest of the diagonal terms
        for (int j=0; j<L-L_b; j++){
            if (basis[i] & sites[j+L_b]){
                diagonal[i] += on_site[j];
                if (basis[i] & sites[(j+L_b+1)%L]){ diagonal[i] += U; }
            }
        }
    } 
    return diagonal;
} //end generate_diagonal



/* ***
 *
 *    GENERATE DIAGONAL TERMS (no bath)
 *
 * ***/
vector<double> generate_diagonal(vector<int> sites, vector<int> basis, vector<double> on_site, double U, int L, int dim){
    vector<double> diagonal(dim,0.0);
    for (int i =0; i<dim; i++){
        for (int j=0; j<L; j++){    //Note the change j=0 -> j=1 to ignore the bath site
            if (basis[i] & sites[j]){
                diagonal[i] += on_site[j];
                if(j == L-1){
                  if(basis[i] & sites[0]){ diagonal[i] += U; }  //Note change to skip bath site with PBC
                }
                else if (basis[i] & sites[j+1]){ diagonal[i] += U; }
            }
        }
    } 
    return diagonal;
} //end generate_diagonal


/* ***
 *
 *    GENERATE DIAGONAL TERMS (open)
 *
 * ***/
vector<double> generate_diagonal_open(vector<int> sites, vector<int> basis, vector<double> on_site, double U, int L, int dim){
    vector<double> diagonal(dim,0.0);
    for (int i =0; i<dim; i++){
        for (int j=0; j<L; j++){   // Loop through basis elements and then sites 
            if (basis[i] & sites[j]){
                diagonal[i] += on_site[j];
                if(j != L-1){
                  if(basis[i] & sites[j+1]){ diagonal[i] += U; }  //Check for interaction
                }
            }
        }
    } 
    return diagonal;
} //end generate_diagonal






#endif
