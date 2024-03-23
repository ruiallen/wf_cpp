## Program overvies
This is a automated version of power.f, grave.f and medoc.f, written in C++.

  Takes two state (N.L.M) and performs the following operations:
  1. Calcualte potential and separation constant.
  2. calulate semi-analytical wave function expansion.
  3. calculate the coupling matrix given the previous results.

<details>
<summary>Todo list</summary>


Radial Coupling and ETF correction terms based on MEDOC integral.

</details>

## Inputfiles:
The only input file required is input.txt, and inside it contains detailed instruction.<br />
## how to run and compile
To compile: <br />
<br />
g++ -std=c++11 one_elec_waveFunction.cpp -o someNameforExecutible <br />
<br />
It is important to use standard 11 and above. And then <br />
<br />
./someNameforExecutible
 
