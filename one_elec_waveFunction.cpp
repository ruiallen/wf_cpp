// one_elec_waveFunction.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include<cmath>
using namespace std;

void CORR(double ZA, double ZB, double N, int L, int M, double &NS, int K, int ICEN) {
    //Correltates UA and SA states using analytic formula of J. Power (Phil. Trans. Royal Society of London 1973)
    if (ZA == ZB) { //Homonuclear OEDM
        ICEN = 0;
        K = (L - M) / 2;
        NS = N - L + M + K;
    }
    else { //Heteronuclear OEDM
        int NU = ((1 + N) * ZB - ZA * (M + N - 1 - L)) / (ZA + ZB);
        if (NU < 0) {
            NU = 0;
        }
        int   NUTEST = (N * ZB - ZA * (M + N - 1 - L)) / (ZA + ZB);
        if (NU == 0 || NUTEST == NU) { //Dissociates with the e- localised on centre A
            ICEN = 1;
            K = L - M - NU;
            NS = N - NU;
        }
        else {//Disccociates with the e- localised on centre B
            ICEN = 2;
            K = NU - 1;
            NS = N - L + M + K;
        }
    }
    return;
}


void WRLG1(double Z, double ZP, int N, int D, double M, double R, double &W, double &W0, double &W1, double &W2, double &W3, double &W4, double &W5, double &W6) {
    //ASYMPTOTIC LARGE R EXPANSION OF W IN POWERS OF 1/R
    int DD = D*D;
    double ratio = ZP / Z;
    int NN = N * N;
    double MM = M * M;
    W0 = -0.5 * Z * Z / NN;
    W1 = -ZP;

    W2 = 1.5 * N * D * ratio;

    W3 = -NN * (6 * DD + 1 - NN) * ratio / (2 * Z);

    W4 = N * NN * D * (109.0 * DD - 39.0 * NN - 9.0 * MM + 59.0) * ratio/16.0 + NN*NN*(3.0 * DD + 9 * MM - 17 * NN -19) * pow(ratio,2.0) / 16.0;

    W5 = NN * NN * (594 * NN * DD - 1065 * DD * DD - 1230 * DD + 234 * MM * DD - 9 * MM * MM - 33 * NN * NN
        + 18 * MM * (1 + NN) - 105 + 138 * NN) * ratio / 64.0 + N * NN * NN * D * (111 * NN - 21 * DD - 63 * MM
            + 189) * pow(ratio, 2.0) / 16;

    W6 = N * NN * NN * D * (2727 * DD * DD - 2076 * NN * DD + 5544 * DD - 1056 * DD * MM + 93 * MM * MM
        + 273 * NN * NN - 78 * NN * MM - 450 * MM + 1533 - 1470 * NN) * ratio / 64 + NN * NN * NN * (207 *
            DD * DD - 1044 * NN * DD - 2436 * DD + 576 * MM * DD + 42 * NN - 371 + 162 * MM - 42 * NN * MM
            + 89 * NN * NN - 15 * MM * MM) * pow(ratio, 2.0) / 32.0 + NN * NN * NN * N * D * (69 * NN - 3 * DD + 117
                + 33 * MM) * pow(ratio, 3.0) / 32;

    W4 = W4 / pow(Z, 2.0);
    W5 = W5 / pow(Z, 3.0);
    W6 = W6 / pow(Z, 4.0);

    double T = 1 / R;
    W = W0 + T * (W1 + T * (W2 + T * (W3 + T * (W4 + T * (W5 + T * W6)))));
    return;


}
void WRLG2(double R, double W, double W0, double W1, double W2, double W3, double W4, double W5, double W6) {
    double T = 1 / R;
    W = W0 + T * (W1 + T * (W2 + T * (W3 + T * (W4 + T * (W5 + T * W6)))));
    return;

}

double CRLG1(int N, int L, int M, double P, double RZ) {
	double NU = 2 * (N - L - 1) + M + 1;
	double MM = M * M;
	double NN = NU * NU;
	double C1 = -0.5 * (NN + 1 - MM);
	double C2 = NU * (NN + 1 - MM);
	double C3 = MM - 1 - 3 * NN;
	double C4 = 6 * NN * MM - (1 - MM) * (1 - MM) - 5 * NN * NN - 10 * NN;
	double C5 = NU * (20 * (NN + 1) - 12 * MM);
	double C6 = 8 * (MM - 1) - 24 * NN;
	double C7 = NU * (NN * (33 * NN + 114 - 46 * MM) + 37 + MM * (13 * MM - 50));
	double C8 = NN * (138 * MM - 342 - 165 * NN) - 37 + MM * (50 - 13 * MM);
	double C9 = NU * (284 * NN + 292 - 156 * MM);
	double C10 = NN * (MM * (230 - 39 * MM + 100 * NN) - 239 - NN * (340 + 63 * NN)) - 14 + MM * (30 - MM * (18 - 2 * MM));
	double C11 = NU * (NN * (1360 + 378 * NN - 400 * MM) + 478 + MM * (78 * MM - 460));
	double C12 = NN * (630 * MM - 1810 - 845 * NN) - 209 + MM * (250 - 41 * MM);
	double C13 = NU * (860 * NN + 900 - 460 * MM);
	double C14 = NU * (NN * (5221 + MM * (465 * MM - 939 * NN - 3750) + NN * (4139 + 527 * NN)) + 1009 - MM * (1591 - MM * (635 - 53 * MM)));
	double C15 = NN * (MM * (11250 - 1395 * MM + 4695 * NN) - 15663 - NN * (20695 + NN * 3689)) - 1009 + MM * (1591 - MM * (635 - 53 * MM));
	double C16 = NU * (14072 + NN * (37640 - 9520 * MM + 10128 * NN) - MM * (11640 - 1440 * MM));
	double C17 = NN * (9780 * MM - 30140 - 13750 * NN) - 3630 + MM * (4140 - 510 * MM);
	double C18 = NU * (9520 * NN + 10080 - 5040 * MM);

	double S = RZ / (2 * P);
	return 2 * P * (S - NU) + NU * S + C1 + (C2 + S * (C3 + 2 * NU * S)) / (8 * P) +
		(C4 + S * (C5 + S * (C6 + 8 * NU * S))) / (64 * P * P) +
		(C7 + S * (C8 + S * (C9 + S * (8 * C6 + 40 * NU * S)))) / (512 * P * P * P) +
		(C10 + S * (C11 + S * (C12 + S * (C13 + S * (16 * C6 + 56 * NU * S))))) / (1024 * P * P * P * P) +
		(C14 + S * (C15 + S * (C16 + S * (C17 + S * (C18 + S * (128 * C6 + 336 * NU * S)))))) / (8192 * P * P * P * P * P);
}


void EXTRAP(vector<double> & X, vector<double> F, int N) {
    // given values of X for i = 0,1,...N-1 and F for i= 0,1,...,N-2 
    int NM2 = N - 2;
    vector<double> DELTA(999);
    for (int i = 0; i < NM2; ++i) {
        //Assume the functions are approximately liner in X so the divided differences will be nearly constant.
        //First extrapolate to get an estimate of the divided difference corresponding to i = N-1
        DELTA[i] = (F[i + 1] - F[i]) / (X[i + 1] - X[i]);
    }
    double SUM = DELTA[NM2];
    int J = NM2;
    while (J > 0) {
        --J;
        if (abs(DELTA[J]) > 0.2 * abs(DELTA[J + 1])) {
            break;
        }
        SUM += DELTA[J];
    }
    F[N] = F[N - 1] + (X[N] - X[N - 1]) * SUM;


}



int main()
{    
    //define constants, As of now not sure their functionality
    vector<int> IACDFT{ 11,11,10,10 };
    vector<int> IACMAX{ 13,13,10,10 };
    vector<int> IACMIN{ 10,10,10,10 };
    const int IIN{ 2 }, IOUT{ 3 }, IHOLD{1}, IXTRAP{ 10 }, ITMAX{ 50 };
    // ZA, ZB: nuclear charges;
    // N,L,M: UNITED atom quantum number
    int QU{ 1 };
    int N{ 1 }, L{ 0 }, M{ 0 };
    int lStart = (L - M != 2 * (L - M) / 2) ? M + 1 : M;

    // some subroutines require ZA >= Zb
    double ZA = 1.0, ZB = QU;
    int IGO = (ZA != ZB) ? 1 : 2;
    double rStart{ 0.0 };
    int Z = ZA + ZB;
    int nPts = 100;
    vector <double> PP(999),CSEP(999), RR(999), ICVGY(999), ICVGX(999);

    double rIncr = 0.3;
    for (int i = 1; i < nPts+1; ++i) {
        RR[i] = i* rStart;
    }
    int NOPTS{ 1 };
    int ISTART = NOPTS + 1;
    NOPTS += nPts;
    double NS{ 0 }, KS{ 0 }, ICEN{ 0 };
    CORR(ZA, ZB, N, L, M, NS, KS, ICEN); 


    //Now we have found the dissociation products of the OEDM state

    double SIGN = pow(-1.0, (ICEN - 1));
    double ZZ = ZA * (2 - ICEN) + ZB * (ICEN - 1);
    double ZP = ZA * (ICEN - 1) + ZB * (2 - ICEN);
    int iDelta = N - L - 1 - KS;

    double rJump = 2 * NS * NS * (1 + 2 * sqrt(abs(ZP / ZZ))) / abs(ZZ);
    double W{ 0 }, W0{ 0 }, W1{ 0 }, W2{ 0 }, W3{ 0 }, W4{ 0 }, W5{ 0 }, W6{0};
    WRLG1(ZZ, ZP, NS, iDelta, M, 1.0, W, W0 , W1, W2, W3, W4, W5, W6);

    W1;
    // if RJ =0, use RJUMP, else RJ
    double RJ = 1.0;
    rJump = (RJ == 0) ? rJump : RJ;

    

    if (rStart == 0) {
        PP[1] = 0;
        CSEP[1] = L * (L + 1);
        ICVGX[1] = 0;
        ICVGY[1] = 0;
        ISTART = 2;
    }
    else {
        ISTART = 1;
        WRLG2(1.0, W, W0, W1, W2, W3, W4, W5, W6);
        PP[1] = rStart * sqrt(-W / 2);
        CSEP[1] = CRLG1(N, L, M, PP[0], rStart * Z);
    }
    CSEP[1];
    //Begin loop over R values
    int NR = NOPTS - ISTART + 1;
    for (int i = ISTART; i <= NOPTS; ++i) {
        double R = RR[i];
        //IGO=1,2, OR 3 RESPECTIVELY FOR HETERONUCLEAR SMALL R, 
        //HOMONUCLEAR OR HETERONUCLEAR LARGE R EXPANSION FOR Y(ETA).
        if (IGO == 3) {
            IGO = 1;
        }
        if (IGO == 1 && R >= rJump) {
            IGO = 3;
        }
        double RZ = R * Z;
        double RZDIF = R * (ZA - ZB) * SIGN;
        if (i>1){
            //extrapolote for P and C at current R
            //USE DATA FROM THE IXTRAP(= 10) PREVIOUS POINTS, IF AVAILABLE
            nPts = min(i,IXTRAP);
            int J = i + 1 - nPts;
            EXTRAP(RR, PP, nPts);
        }
        PP;
        //evaluate the second term
        PP[1] = 0.5 * R * Z / N;
        double ZK = RZDIF * 0.5 / PP[1];
        double DC{ 0 };
        if (L == 0) {
            DC = 2 * (1 - ZK * ZK) / 3 * pow(PP[2], 2.0);
        }
        else {// L!=0
            DC = 2 * (((L + 1)*(L+1) - M * M) * ((L + 1)*(L+1) - ZK * ZK) / ((L + 1) * (L + L + 3))
                - (L * L - M * M) * (L * L - ZK * ZK) / (L * (L + L - 1))) / (L + L + 1) * pow(PP[1], 2);

        }
        CSEP[1] = CSEP[0] + DC;

        //goto 7, this is where 7 starts
        double P = PP[i];
        double C = CSEP[i];

    }
    PP;
    
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
