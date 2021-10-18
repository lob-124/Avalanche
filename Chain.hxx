#ifndef CHAIN
#define CHAIN

#include <vector>
#include <tuple>
#include <algorithm>
#include <complex>
#include <random>
#include <time.h>
#include <armadillo>
#include "generate-matrices.hxx"
#include <stdlib.h>


#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;


bool compare(complex<double> i, complex<double> j){
    // ****
    //      Comparator functions for complex double's. Compares by real part
    //****
    if(real(i) == real(j))
        return imag(i) < imag(j);
    else
        return real(i) < real(j);
}


bool compare_pair(tuple<int,complex<double>> first, tuple<int,complex<double>> second){
    // ****
    //      Comparator function for <integer, complex double> pairs. Compares by real part
    //          of the associated complex numbers
    //****
    complex<double> i = get<1>(first); complex<double> j = get<1>(second);
    if(real(i) == real(j))
        return imag(i) < imag(j);
    else
        return real(i) < real(j);
}


// /****
//         Class Definition
// ****/
// class Bath_Chain{
// /****
//         Bath_Chain has the bath OUTSIDE of the chain, connected only to one site
// ****/
// private:
// mt19937 engine;
// uniform_real_distribution<double> on_site_dist;
// uniform_real_distribution<double> hopping_dist;
// normal_distribution<double> bath_dist;
// int bath_mask;

// public:
// int L; int N; int L_b; double W; double U; double g; double beta; unsigned int seed; 
// int dim = 1; int bath_dim = 1;
// vector<int> sites; vector<int> basis; 
// vector<double> on_site; vector<double> hopping_amplitudes;

// Mat<double> H;
// vector<double> hopping; vector<double> diagonal; 
// Mat<double> H_bath;

// Bath_Chain(int l, int n, int l_b, double w, double u, double G, double Beta, unsigned int s){
//     L = l; N = n; L_b = l_b; W = w; U = u; g = G; beta = Beta; seed = s;
    
//     //Seed RNG
//     if(seed){
//         //cout << "uhh yup" <<  endl;
//         engine = mt19937(seed);
//     }
//     else{
//         //cout << "uhh nope" <<  endl;
//         srand(time(NULL));
//         seed = (unsigned int) rand();
//         //cout << seed << endl;
//         engine = mt19937(seed);
//     }

//     //Initialize distributions for on-site potentials, hopping amplitudes, bath
//     on_site_dist = uniform_real_distribution<double>(-W,W);
//     hopping_dist = uniform_real_distribution<double>(0.0,1.0);
//     bath_dist = normal_distribution<double>(0.0,1.0);


//     //Compute dimensionality
//     for(int i=L; i > N; i--){dim *= i;}
//     for(int i=2; i <= L-N; i++){dim /= i;}

//     //Compute bath dimensionality
//     bath_dim = bath_dim << L_b;
//     bath_mask = bath_dim - 1;

//     //Compute sites array - an array of integers which, when written in binary,
//     //  have a single 1 in position i, representing a particle at site i
//     sites.resize(L);
//     for(int i=0; i < L; i++){sites[i] = 1 << (L-i-1);}

//     //Compute basis 
//     basis = generate_combos(sites, 0, -1,N,L);

//     //Create the Hamiltonian 
//     H.zeros(dim,dim);  
//     hopping.resize(dim*dim);
//     diagonal.resize(dim);

//     //Initialize random on-site potential
//     on_site.resize(L-L_b);
//     for (int i=0; i<L-L_b; i++){ on_site[i] = on_site_dist(engine);}

//     //Initialize random power-law hopping coefficients
//     hopping_amplitudes.resize(L-L_b);
//     for (int i=0; i<L-L_b; i++){ hopping_amplitudes[i] = hopping_dist(engine); cout << hopping_amplitudes[i] << endl;}

//     //Generate matrices for hopping & diagonal (on-site + interactions) terms
//     hopping = generate_hopping(sites,basis,dim,L,N,L_b,g,hopping_amplitudes); 
//     diagonal = generate_diagonal(sites,basis,on_site,U,L,L_b,dim); 

//     //Initialize bath Hamiltonian
//     H_bath.zeros(bath_dim, bath_dim);
//     for(int i=0; i<bath_dim; i++){
//         for (int j = i+1; j < bath_dim; j++){
//             H_bath(i,j) = beta*bath_dist(engine);
//             H_bath(j,i) = H_bath(i,j);
//         }

//         H_bath(i,i) = 1.4142135623730951*beta*bath_dist(engine);    //GOE has  H_ii ~ N(0,2) = \sqrt{2}*N(0,1)
//     }

//     //Copy the matrix elements into the Hamiltonian
//     for(int i=0; i<dim; i++){
//         for(int j=i+1; j<dim; j++){
//             //Hopping
//             H(i,j) += hopping[dim*i + j];
//             H(j,i) += hopping[dim*j+i];

//             //Copy in the bath Hamiltonian terms
//             if((~bath_mask & basis[i]) == (~bath_mask & basis[j])){
//                 int bath_index_1 = (bath_mask & basis[i]) >> (L - L_b);
//                 int bath_index_2 = (bath_mask & basis[j]) >> (L - L_b);
//                 H(i,j) += H_bath(bath_index_1,bath_index_2);
//                 H(j,i) += H_bath(bath_index_1,bath_index_2);
//             }
//         }

//         //Add in the diagonal terms
//         H(i,i) += diagonal[i];
//         int bath_index = (bath_mask & basis[i]) >> (L-L_b);
//         H(i,i) += H_bath(bath_index,bath_index);
//     } 
    
// } //End constructor


// void set_g(double new_g){ 
//     // ****
//     //      Set a new value of the tilt g for this Chain & update
//     //           the Hamiltonian
//     // ****
//     this->g = new_g; 
//     this->gen_hopping();
// }


// void set_W(double new_W){
//     // ****
//     //      Set a new value of the disorder width W for this Chain 
//     //          & update the Hamiltonian
//     // ****
//     this->W = new_W;
//     this->on_site_dist = uniform_real_distribution<double>(-new_W,new_W);
//     this->new_disorder();
// }


// void set_U(double new_U){
//     // ****
//     //      Set a new value of the interaction strength U for this Chain 
//     //          & update the Hamiltonian
//     // ****
//     this->U = new_U;
//     this->gen_diagonal();
// }


// void set_beta(double new_beta){
//     // ****
//     //      Set a new value of the bath "strength" beta for this Chain 
//     //          & update the Hamiltonian
//     // ****
//     this->beta = new_beta;
//     this->new_bath();
// }


// void new_disorder(){
//     // ****
//     //      Generate a new disorder realization & update the Hamiltonian
//     // ****
//     for (int i=0; i<this->L-this->L_b; i++){ this->on_site[i] = this->on_site_dist(this->engine);}
//     this->gen_diagonal();   
// }

// void new_amplitudes(){
//     // ****
//     //      Generate new random hopping amplitudes & update the Hamiltonian
//     // ****
//     for (int i=0; i<this->L - this->L_b; i++){ this->hopping_amplitudes[i] = this->hopping_dist(this->engine);}
//     this->gen_hopping();   
// }


// void new_bath(){
//     // ****
//     //      Generate a new random bath Hamiltonian & update the Hamiltonian
//     // ****
//     this->H_bath.zeros(this->bath_dim, this->bath_dim);
//     for(int i=0; i<this->bath_dim; i++){
//         for (int j = i+1; j < this->bath_dim; j++){
//             this->H_bath(i,j) = this->beta*this->bath_dist(this->engine);
//             this->H_bath(j,i) = this->H_bath(i,j);
//         }

//         this->H_bath(i,i) = 1.4142135623730951*this->beta*this->bath_dist(this->engine);
//     }
//     this->gen_bath();
// }


// void gen_diagonal(){
//     // ****
//     //      Compute & update the diagonal terms (on-site potential + interactions) for this 
//     //          Chain's Hamiltonian
//     // ****
//     this->diagonal = generate_diagonal(this->sites,this->basis,this->on_site,this->U,this->L,this->L_b, this->dim);
//     for(int i=0; i<dim; i++){ 
//         this->H(i,i) = this->diagonal[i]; 

//         //Copy in the bath Hamiltonian terms
//         int bath_index = (this->bath_mask & this->basis[i]) >> (this->L - this->L_b);
//         this->H(i,i) += this->H_bath(bath_index,bath_index);
//     }
// }


// void gen_hopping(){ 
//     // ****
//     //      Compute & update the hopping terms for this Chain's Hamiltonian
//     // ****   
//     this->hopping = generate_hopping(this->sites,this->basis, this->dim, this->L, this->N, this->L_b, this->g, this->hopping_amplitudes); 
//     for(int i=0; i<this->dim; i++){
//         for(int j=i+1; j<this->dim; j++){
//             this->H(i,j) = this->hopping[this->dim*i + j];
//             this->H(j,i) = this->hopping[this->dim*j + i];

//             //Copy in the bath Hamiltonian terms
//             if((~this->bath_mask & this->basis[i]) == (~this->bath_mask & this->basis[j])){
//                 int bath_index_1 = (this->bath_mask & this->basis[i]) >> (this->L - this->L_b);
//                 int bath_index_2 = (this->bath_mask & this->basis[j]) >> (this->L - this->L_b);
//                 this->H(i,j) += this->H_bath(bath_index_1,bath_index_2);
//                 this->H(j,i) += this->H_bath(bath_index_1,bath_index_2);
//             }
//         }
//     }
// }


// void gen_bath(){
//     //****
//     //     Compute & update the bath terms for this Chain's Hamiltonian
//     //**** 
//     for(int i=0; i<this->dim; i++){
//         for(int j=i+1; j<this->dim; j++){
//             this->H(i,j) = this->hopping[this->dim*i + j];
//             this->H(j,i) = this->hopping[this->dim*j + i];

//             //Copy in the bath Hamiltonian terms
//             if((~this->bath_mask & this->basis[i]) == (~this->bath_mask & this->basis[j])){
//                 int bath_index_1 = (this->bath_mask & this->basis[i]) >> (this->L - this->L_b);
//                 int bath_index_2 = (this->bath_mask & this->basis[j]) >> (this->L - this->L_b);
//                 this->H(i,j) += this->H_bath(bath_index_1,bath_index_2);
//                 this->H(j,i) += this->H_bath(bath_index_1,bath_index_2);
//             }
//         }
    
//         this->H(i,i) = this->diagonal[i];
//         int bath_index = (this->bath_mask & this->basis[i]) >> (this->L - this->L_b);
//         this->H(i,i) += this->H_bath(bath_index,bath_index);
//     }
// }




// vector<complex<double>> spectrum(){
//     // ****
//     //      Compute & return the full spectrum, sorted by ascending real part 
//     // ****
//     cx_vec eigvals = eig_gen(this->H);
//     vector<complex<double>> spectrum = conv_to<vector<complex<double>>>::from(eigvals);
//     sort(spectrum.begin(), spectrum.end(), compare);
    
//     return spectrum;
// } //end spectrum()


// tuple<vector<complex<double>>, vector<vector<complex<double>>>> states(){
//     // ****
//     //      Compute & return the full spectrum and eigenvectors, sorted by ascending real part of 
//     //          the eigenvalues
//     // ****

//     //Compute the eigenvalues and eigenvectors
//     cx_vec eigvals; cx_mat eigvecs;
//     eig_gen(eigvals,eigvecs,this->H);

//     //Sort the spectrum (note we pair each energy with its original index so we 
//     //	are effectively argsort-ing)
//     vector<tuple<int,complex<double>>> pairs(this->dim);
//     for(int i=0; i<dim; i++){
//     	pairs[i] = make_tuple(i,eigvals(i));
//     }
//     sort(pairs.begin(), pairs.end(), compare_pair);
//     vector<complex<double>> spectrum(this->dim);
//     vector<vector<complex<double>>> states(this->dim);
//     for(int i=0; i<dim; i++){
//     	spectrum[i] = get<1>(pairs[i]);
//     	states[i] = conv_to<vector<complex<double>>>::from(eigvecs.col(get<0>(pairs[i])));
//     }
    
//     //Return the sorted spectrum and sorted eigenvectors
//     return make_tuple(spectrum,states);
// } //end states()


// }; //END BATH_CHAIN CLASS



/****
        Class Definition
****/
class Chain{
private:
mt19937 engine;
uniform_real_distribution<double> dist;
uniform_real_distribution<double> hopping_dist;

public:
int L; int N; double W; double U; double g; unsigned int seed; 
int dim = 1; 
vector<int> sites; vector<int> basis; vector<double> on_site; vector<double> hopping_amps; 
Mat<double> H;
vector<double> left_hopping; vector<double> diagonal;

Chain(int l, int n, double w, double u, double G, unsigned int s){
    L = l; N = n; W = w; U = u; g = G; seed = s;
    
    //Seed RNG
    if(seed){
        engine = mt19937(seed);
    }
    else{
        srand(time(NULL));
        seed = (unsigned int) rand();
        engine = mt19937(seed);
    }

    dist = uniform_real_distribution<double>(-W,W);
    hopping_dist = uniform_real_distribution<double>(0.0,1.0);

    //Compute dimensionality
    if(N==1){ dim = L; }
    else{
        for(int i=L; i > N; i--){ dim *= i; }
        for(int i=2; i <= L-N; i++){ dim /= i; }
    }

    //Compute sites array
    sites.resize(L);
    for(int i=0; i < L; i++){sites[i] = 1 << (L-i-1);}
    

    //Compute basis 
    basis = generate_combos(sites, 0, -1,N,L);

    /* ***
     *    Create the Hamiltonian
     * ***/  
    H.zeros(dim,dim);   
    left_hopping.resize(dim*dim);
    diagonal.resize(dim);

    //Initialize random on-site potential
    on_site.resize(L);
    for (int i=0; i<L; i++){ on_site[i] = dist(engine);}

    //Initialize random hopping amplitudes
    hopping_amps.resize(L);
    for(int i=0; i<L; i++){ hopping_amps[i] = 1.0; }//hopping_dist(engine); }


    left_hopping = generate_left_hopping(sites,basis,hopping_amps,dim,L,N); 
    diagonal = generate_diagonal(sites,basis,on_site,U,L,dim);   

    /* ****
     *     Copy the matrix elements into the Hamiltonian
     * ****/
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            H(i,j) += exp(g)*left_hopping[dim*i + j];
            H(j,i) += exp(-g)*left_hopping[dim*i + j];
        }
    } 
    
    for(int i=0; i<dim; i++){ H(i,i) = diagonal[i]; }
   
} //End constructor


void new_g(double new_g){ 
    this->g = new_g; 
    this->make_H();

} //end new_g()


void new_W(double new_W){
    this->W = new_W;
    this->dist = uniform_real_distribution<double>(-new_W,new_W);
    this->new_disorder();
} //end new_W()


void new_U(double new_U){
    this->U = new_U;
    this->new_disorder();
} //end new_U()


void new_disorder(){
    for (int i=0; i<L; i++){ on_site[i] = dist(engine);}
    this->diagonal = generate_diagonal(this->sites,this->basis,this->on_site,this->U,this->L, this->dim);
    for(int i=0; i<dim; i++){ this->H(i,i) = this->diagonal[i]; }

}// end new_disorder()


void make_H(){

    /* ****
     *     Copy the matrix elements into the Hamiltonian
     * ****/
    this->H.zeros(this->dim,this->dim); //=  Eigen::MatrixXd::Zero(this->dim,this->dim);  
    for(int i=0; i<this->dim; i++){
        for(int j=0; j<this->dim; j++){
            this->H(i,j) += exp(this->g)*this->left_hopping[this->dim*i + j];
            this->H(j,i) += exp(-this->g)*this->left_hopping[this->dim*i + j];
        }
    } 
    
    for(int i=0; i<this->dim; i++){ this->H(i,i) = this->diagonal[i]; }

} // end make_H()


vector<complex<double>> spectrum(){
    
    cx_vec eigvals = eig_gen(this->H);
    vector<complex<double>> spectrum = conv_to<vector<complex<double>>>::from(eigvals);//(this->dim);
    sort(spectrum.begin(), spectrum.end(), compare);
    
    return spectrum;
} //end spectrum()


tuple<vector<complex<double>>, vector<vector<complex<double>>>> states(){
    //Compute the eigenvalues and eigenvectors
    cx_vec eigvals; cx_mat eigvecs;
    eig_gen(eigvals,eigvecs,this->H);

    //Sort the spectrum (note we pair each energy with its original index so we 
    //  are effectively argsort-ing)
    vector<tuple<int,complex<double>>> pairs(this->dim);
    for(int i=0; i<dim; i++){
        pairs[i] = make_tuple(i,eigvals(i));
    }
    sort(pairs.begin(), pairs.end(), compare_pair);
    vector<complex<double>> spectrum(this->dim);
    vector<vector<complex<double>>> states(this->dim);
    for(int i=0; i<dim; i++){
        spectrum[i] = get<1>(pairs[i]);
        states[i] = conv_to<vector<complex<double>>>::from(eigvecs.col(get<0>(pairs[i])));
    }
    
    //Return the sorted spectrum and sorted eigenvectors
    return make_tuple(spectrum,states);
} //end states()



}; //END CHAIN CLASS


/****
        Class Definition
****/
class Open_Chain{
private:
mt19937 engine;
uniform_real_distribution<double> dist;
uniform_real_distribution<double> hopping_dist;

public:
int L; int N; double W; double U; unsigned int seed; 
int dim = 1; 
vector<int> sites; vector<int> basis; vector<double> on_site; vector<double> hopping_amps; 
Mat<double> H;
vector<double> left_hopping; vector<double> diagonal;

Open_Chain(int l, int n, double w, double u, unsigned int s){
    L = l; N = n; W = w; U = u; seed = s;
    
    //Seed RNG
    if(seed){
        engine = mt19937(seed);
    }
    else{
        srand(time(NULL));
        seed = (unsigned int) rand();
        engine = mt19937(seed);
    }

    dist = uniform_real_distribution<double>(-W,W);
    hopping_dist = uniform_real_distribution<double>(0.0,1.0);

    //Compute dimensionality
    if(N==1){ dim = L; }
    else{
        for(int i=L; i > N; i--){ dim *= i; }
        for(int i=2; i <= L-N; i++){ dim /= i; }
    }

    //Compute sites array
    sites.resize(L);
    for(int i=0; i < L; i++){sites[i] = 1 << (L-i-1);}
    

    //Compute basis 
    basis = generate_combos(sites, 0, -1,N,L);

    /* ***
     *    Create the Hamiltonian
     * ***/  
    H.zeros(dim,dim);   
    left_hopping.resize(dim*dim);
    diagonal.resize(dim);

    //Initialize random on-site potential
    on_site.resize(L);
    for (int i=0; i<L; i++){ on_site[i] = dist(engine);}

    //Initialize random hopping amplitudes
    hopping_amps.resize(L);
    for(int i=0; i<L; i++){ hopping_amps[i] = 1.0; }//hopping_dist(engine); }


    left_hopping = generate_left_hopping_open(sites,basis,hopping_amps,dim,L,N); 
    diagonal = generate_diagonal_open(sites,basis,on_site,U,L,dim);   

    /* ****
     *     Copy the matrix elements into the Hamiltonian
     * ****/
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            H(i,j) += left_hopping[dim*i + j];
            H(j,i) += left_hopping[dim*i + j];
        }
    } 
    
    for(int i=0; i<dim; i++){ H(i,i) = diagonal[i]; }
   
} //End constructor

void set_W(double new_W){
    this->W = new_W;
    this->dist = uniform_real_distribution<double>(-new_W,new_W);
    this->new_disorder();
} //end new_W()


void new_U(double new_U){
    this->U = new_U;
    this->new_disorder();
} //end new_U()


void new_disorder(){
    for (int i=0; i<L; i++){ on_site[i] = dist(engine);}
    this->diagonal = generate_diagonal_open(this->sites,this->basis,this->on_site,this->U,this->L, this->dim);
    for(int i=0; i<dim; i++){ this->H(i,i) = this->diagonal[i]; }

}// end new_disorder()


void set_disorder(vector<double> realization){
    for(int i=0; i<L; i++){ this->on_site[i] = realization[i]; }
    this->diagonal = generate_diagonal_open(this->sites,this->basis,this->on_site,this->U,this->L,this->dim);
    for(int i=0; i<dim; i++){ this->H(i,i) = this->diagonal[i]; }
} //end set_disorder()


void make_H(){

    /* ****
     *     Copy the matrix elements into the Hamiltonian
     * ****/
    this->H.zeros(this->dim,this->dim); //=  Eigen::MatrixXd::Zero(this->dim,this->dim);  
    for(int i=0; i<this->dim; i++){
        for(int j=0; j<this->dim; j++){
            this->H(i,j) += this->left_hopping[this->dim*i + j];
            this->H(j,i) += this->left_hopping[this->dim*i + j];
        }
    } 
    
    for(int i=0; i<this->dim; i++){ this->H(i,i) = this->diagonal[i]; }

} // end make_H()


vector<double> spectrum(){
    //Compute the spectrum (note it will be sorted already) 
    vec eigvals = eig_sym(this->H);
    vector<double> spectrum = conv_to<vector<double>>::from(eigvals);
    
    return spectrum;
} //end spectrum()


tuple<vector<double>, vector<vector<double>>> states(){
    //Compute the eigenvalues and eigenvectors
    vec eigvals; mat eigvecs;
    eig_sym(eigvals,eigvecs,this->H);

    //Convert to the return types 
    vector<double> spectrum(this->dim);
    spectrum = conv_to<vector<double>>::from(eigvals);
    vector<vector<double>> states(this->dim);
    states = conv_to<vector<vector<double>>>::from(eigvecs);
    
    //Return the sorted spectrum and sorted eigenvectors
    return make_tuple(spectrum,states);
} //end states()



}; //END OPEN_CHAIN CLASS




#endif
