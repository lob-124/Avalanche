compile-SC: Spin-Chain.hxx generate-data.cpp
	g++ -std=c++11 -O3 generate-data.cpp -o generate-data -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-C: Chain.hxx generate-data-chain-p.cpp
	g++ -std=c++11 -O3 -fopenmp generate-data-chain-p.cpp -o generate-data-chain-p -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-chain: Chain.hxx generate-data-p.cpp
	g++ -std=c++11 -O3 -fopenmp generate-data-p.cpp -o generate-data-p -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-convergence: Bath-Chain.hxx Mapping.hxx generate-matrices.hxx convergence.cpp
	g++ -std=c++11 -O3 -fopenmp convergence.cpp -o convergence -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-corr: Chain.hxx Mapping.hxx generate-matrices.hxx corr.cpp 
	g++ -std=c++11 -O3 -fopenmp corr.cpp -o corr -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-corr-open: Chain.hxx Mapping.hxx generate-matrices.hxx corr-open.cpp 
	g++ -std=c++11 -O3 -fopenmp corr-open.cpp -o corr-open -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
ccompile-correlations: Chain.hxx Mapping.hxx generate-matrices.hxx correlations.cpp 
	g++ -std=c++11 -O3 -fopenmp correlations.cpp -o correlations -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-corr-conv: Chain.hxx Mapping.hxx generate-matrices.hxx corr-conv.cpp
	g++ -std=c++11 -O3 -fopenmp corr-conv.cpp -o corr-conv -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-gs: Bath-Chain.hxx generate-matrices.hxx critical-gs.cpp
	g++ -std=c++11 -O3 -fopenmp critical-gs.cpp -o critical-gs -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack
compile-corr-open-gs: Chain.hxx Mapping.hxx generate-matrices.hxx corr-open-gs.cpp
	g++ -std=c++11 -O3 -fopenmp corr-open-gs.cpp -o corr-open-gs -I ../../armadillo-10.1.2/include -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

