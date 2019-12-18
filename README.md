# Nonlinear Markov Clustering By Minimum Curvilinear Sparse Similarity

Steps to apply sucessfuly MC_MCL in MATLAB (tested using cygwin in a windows 10 environment).

1. for windows users: Install cygwin, a collection of GNU and open source tools which provide functionality like in linux.
2. Install C/C++ compilers in the cygwin/linux environment (required for step 3).
3. Install MCL following the instructions as in https://micans.org/mcl/
4. MC_MCL should be ready to use. Its I/O are:

Inputs:

x:  matrix data where in the rows are the samples and in the
columns the features.

C: number of clusters that the user is trying to find.

dist:  number that represent the distance applied on "x" to construct the graph
      [1] apply the MC matrix created with euclidean distance
      [2] apply the MC matrix created with correlation distance
      [3] apply the MC matrix created with correlation distance
default: 1

cygwin_path: cygwin path to run MCL (folder where cygwin was installed)

factor: the option for the MC distance:
      [1] nothing
      [2] sqrt
      [3] log(1+x)
default: 1

cutoff:  threshold in the similarity of the samples used to compute
the graph. The higher the cutoff it is, the more number of components the graph may have.
default: the program calculate for itselfs the higher cutoff where the
graph contains 1 component.

max_nc: maximum number of components of the graph for choosing the cutoff
default : 1

mod: perform modality:
         [1] choosing cutoff automatically, matrix-dataset goes from 0
         until 1-cutoff.
         [2] sparcifying network using cutoff selected by the user
default: 1

Outputs:

comm: vector of size of the lengths of the rows of x containing the associated clusters to the samples.

nc: number of components of the network introuced to MCL.

