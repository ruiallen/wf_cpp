This is a automated version of power.f, grave.f and medoc.f, written in C++.
Takes two state (N.L.M) and performs the following operations:
  1. Calcualte potential and separation constant.
  2. calulate semi-analytical wave function expansion.
  3. calculate the coupling matrix given the previous results.
  4. (TODO) Radial Coupling and ETF correction terms based on MEDOC integral.

To compile:
g++ -std=c++11 one_elec_waveFunction.cpp -o someNameforExecutible
It is important to use standard 11 and above. And then
./someNameforExecutible
 
