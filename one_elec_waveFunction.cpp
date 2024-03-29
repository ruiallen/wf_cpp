// one_elec_waveFunction.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include<cmath>
#include<tuple>
#include <fstream>
#include <string>
using namespace std;

void CORR(double ZA, double ZB, double N, int L, int M, double &NS, int &K, int &ICEN) {
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
void WRLG2(double R, double &W, double W0, double W1, double W2, double W3, double W4, double W5, double W6) {
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


void EXTRAP(vector<double> & X, vector<double> &F, int N) {
    // given values of X for i = 0,1,...N-1 and F for i= 0,1,...,N-2 
    int NM2 = N - 2;
    vector<double> DELTA(999);
    for (int i = 1; i <= NM2; ++i) {
        //Assume the functions are approximately liner in X so the divided differences will be nearly constant.
        //First extrapolate to get an estimate of the divided difference corresponding to i = N-1
        DELTA[i] = (F[i + 1] - F[i]) / (X[i + 1] - X[i]);
    }
    double SUM = DELTA[NM2];
    int J = NM2;
    J = J - 1;
    while (J != 0) {
        for (int i = 1; i <= J; ++i) {
            DELTA[i] = DELTA[i + 1] - DELTA[i];
        }
        J--;
    }//J = 0, goto 4
    //4
    double curSum = DELTA[NM2];
    J = NM2;
    J = J - 1;
	while (J!=0 and not (abs(DELTA[J])>0.2*abs(DELTA[J+1]))) {
        curSum += DELTA[J];
		J--;
	}//J = 0, goto 6
    F[N] = F[N - 1] + (X[N] - X[N - 1]) * curSum;
    return;
}


void CTDFRN(int IGO, double& P, double& C, double RZ, double RZDIF, int M, double LCHAIN, double LSTART, double ACCY, double CUTOFF, double& F, double &DFDC, double &DFDP,int &ICHK) {

    bool TEST = false;
    const double SCALE{ 1e-7 };
    ICHK = 0;
    //ignore all outputs in the original Fortran subroutine
    double ACY = ACCY;
    int MM = M * M;
    double P2 = 2 * P * P;
    int JINCR = 1;
    //....................................................//
    //....................................................//
    //....................................................//
    //....................................................//
    double AO = 1;
    double DAODC = 0;
    double DAODP = 0;
    double DADC = 1;
    int JGO; // JGO is used as a label



    // initialize vars used in later cases
    int J{ 0 };
    double RZDIF2{ 0 };
    double A{ 0 };
    double DADP{ 0 };
    double P4{ 0 };
    int JJ{ 0 };
    double ZK{ 0 };
    double S1{ 0 };
    double S2{ 0 };
    double DS1{ 0 };
    double SIGMA{ 0 };
    double DSIGMA{ 0 };
    double DS2{ 0 };
    double ASAVE{ 0 };
    double AJ{ 0 };
    double BJ{ 0 };
    double DAJ{ 0 };
    double DBJ{ 0 };
    double BO = 0;
    double DBODC = 0;
    double DBODP = 0;
    double B = 1;
    double DBDC = 0;
    double DBDP = 0;
    double RATOLD = A;
    double DIFOLD = 1;
    double AC = 0;
    double BC = 0;
    double AP = 0;
    double BP = 0;
    double ACP = 0;
    double BSAVE{ 0 };
    double RATNEW{ 0 };
    double DENOM{ 0 };
    double DIFNEW{ 0 };
    double RATVO{ 0 };
    switch (IGO) {
    case 1:
        J = M;
        RZDIF2 = RZDIF * RZDIF;
        A = C - MM - M;
        DADP = 0;
        JGO = 6;
        goto label40;
    case 2:
        J = LSTART;
        P4 = P2 * P2 / 4;
        JJ = J * J;
        ASAVE = (MM + JJ + J - 1.0) / (4 * (JJ + J) - 3);
        A = C - JJ - J - P2 * ASAVE;
        DADP = -4 * P * ASAVE;
        //HOmonuclear Y(ETA) requires incr = 2
        JINCR = 2;
        JGO = 7;
        goto label40;
    case 3:
        J = 0;
        ZK = 0.5 * RZDIF / P;
        S1 = C - MM - M - 2 * P * (M + 1 - ZK);
        S2 = -4 * P - M - M - 1;
        DS1 = -2 * (M + 1);
        A = S1;
        DADP = DS1;
        JGO = 8;
        goto label40;
    case 4:
        J = 0;
        SIGMA = 0.5 * RZ / P - M - 1;
        S1 = C - RZ - (M + 1) * (M + SIGMA - P - P);
		S2 = 4 * P - 2 * SIGMA;
		DSIGMA = -0.5 * RZ / pow(P, 2);
		DS1 = (M + 1) * (2 - DSIGMA);
		DS2 = 4 - 2 * DSIGMA;
		A = S1;
		DADP = DS1;
        JGO = 9;
        goto label40;;
    }//all cases go to label 40
label5:
    J = J + JINCR;
    ICHK = ICHK + 1;
    if (ICHK == 10000) { cout << "Errors" << endl; return; }
    JJ = J * J;
    switch (JGO) {
    case 6:
        ASAVE = (JJ - MM) / (4 * JJ - 1);
        AJ = (2 * P2 * JJ - RZDIF2) * ASAVE;
        BJ = C - JJ - J;
        DAJ = 8 * P * JJ * ASAVE;
        DBJ = 0;
        goto label10;

    case 7:
        AJ = -P4 * (JJ - J - J + 1 - MM) * (JJ - MM) / ((4 * JJ - 8 * J + 3) * (4 * JJ - 1));
        ASAVE = (MM + JJ + J - 1.0) / (4 * (JJ + J) - 3);
        BJ = C - JJ - J - P2 * ASAVE;
        DAJ = 4 * AJ / P;
        DBJ = -4 * P * ASAVE;
        goto label10;
    case 8:
        ASAVE = -4 * J * (J + M);
        DAJ = ASAVE * (J + M);
        AJ = P * ASAVE * (J + M - ZK);
        BJ = -JJ + J * S2 + S1;
        DBJ = DS1 - 4 * J;
        goto label10;
    case 9:
        ASAVE = -J * (J + M);
        AJ = ASAVE * (J - 1 - SIGMA) * (J - 1 - SIGMA - M);
        BJ = S1 + J * S2 + 2 * JJ;
        DAJ = ASAVE * (M + 2 * (1 - J + SIGMA)) * DSIGMA;
        DBJ = DS1 + J * DS2;
        goto label10;
    }
label10:
    if (TEST) { goto label12; }
    ASAVE = DADC;
    DADC = A + BJ * DADC + AJ * DAODC;
    DAODC = ASAVE;
    ASAVE = DADP;
    DADP = BJ * DADP + AJ * DAODP + DBJ * A + DAJ * AO;
    DAODP = ASAVE;
    ASAVE = A;
    A = BJ * A + AJ * AO;
    AO = ASAVE;
label40:
    if (J < LCHAIN) { goto label5; }
    TEST = true;
    A = A / AO;
    DADP = (DADP - A * DAODP) / AO;
    DADC = (DADC - A * DAODC) / AO;
    AO = 1;
    DAODC = 0;
    DAODP = 0;
    BO = 0;
    DBODC = 0;
    DBODP = 0;
    B = 1;
    DBDC = 0;
    DBDP = 0;
    RATOLD = A;
    DIFOLD = 1;
    goto label5;
label12://start backward evaluation
    AC = DADC;
    DADC = A + BJ * DADC + AJ * DAODC;
    BC = DBDC;
    DBDC = B + BJ * DBDC + AJ * DBODC;
    AP = DADP;
    DADP = BJ * DADP + AJ * DAODP + DBJ * A + DAJ * AO;
    BP = DBDP;
    DBDP = BJ * DBDP + AJ * DBODP + DBJ * B + DAJ * BO;
    ASAVE = A;
    A = BJ * A + AJ * AO;
    BSAVE = B;
    B = BJ * B + AJ * BO;
    RATNEW = A / B;
    DENOM = abs(RATNEW) + SCALE; //avoid possible divide by zero
       
    DIFNEW = abs(RATNEW - RATOLD) / DENOM;
    if (DIFNEW > 0.1) { goto label13; }
       
    if ((DIFNEW + DIFOLD) <= ACY) {goto label15;}
    ACP = CUTOFF / DENOM;
    if (ACP > 1e-2) { ACP = 1e-4; }
    if (ACY < ACP) { ACY = ACP; }
label13:
    DIFOLD = DIFNEW;
    RATVO = RATOLD;
    RATOLD = RATNEW;
    if (abs(B) > 1e5) { goto label14; }
    AO = ASAVE;
    BO = BSAVE;
    DAODC = AC;
    DAODP = AP;
    DBODC = BC;
    DBODP = BP;
    goto label5;
label14://Rescale terms to avoid exponent overflow
    A = A * SCALE;
    B = B * SCALE;
    AO = ASAVE * SCALE;
    BO = BSAVE * SCALE;
    DADC = DADC * SCALE;
    DBDC = DBDC * SCALE;
    DAODC = AC * SCALE;
    DBODC = BC * SCALE;
    DADP = DADP * SCALE;
    DBDP = DBDP * SCALE;
    DAODP = AP * SCALE;
    DBODP = BP * SCALE;
    goto label5;

label15:
    F = RATNEW;
    
    if (abs(F - RATOLD) < 1e-12 * abs(F)) { goto label16; }
        
    F = RATVO - (RATOLD - RATVO) * (RATOLD - RATVO) / (F - RATOLD - RATOLD + RATVO);
 label16:   
    DFDC = (DADC - F * DBDC) / B;
    DFDP = (DADP - F * DBDP) / B;
    return;
}

void CHAIN(int& FCHAIN, int& GCHAIN, int IGO, double &RZ, double &RZDIF, int N, int L, int M, double KS, double &P, double &C) {
    //Function that determines the values of the continued fraction indices at which we chain
    int J{ 0 };
    double SIGMA{ 0 };
    double A{ 0 };
    double B{ 0 };
    double ALPHA{ 0 };
    double BETA{ 0 };
    J = N - L - 1;
    SIGMA = 0.5 * RZ / P - M - 1;
    A = 2 * P * SIGMA + (M + 1) * (M + SIGMA) - C;
    B = 2 * (SIGMA - P - P);
label1:
    J = J + 1;
    ALPHA = (J - 1 - SIGMA) * (J - 1 - SIGMA - M);
    BETA = -2 * J * J + J * B + A;
    if (abs(ALPHA) > abs(BETA)) { goto label1; }
    GCHAIN = J - 1;
    if (IGO != 3) { goto label3; }
    J = KS;
    //heteronuclear large R expansion for Y(ETA)
    A = C - 2 * P * (M + 1) + RZDIF - M * (M + 1);
label2:
    J = J + 1;
    ALPHA = 2 * P * (J + M) - RZDIF;
    BETA = A - 4 * P * J - J * (J + 1 + M + M);
    if (abs(ALPHA) > abs(BETA)) { goto label2;}
    FCHAIN = J - 1;
    return;
label3:
    J = L;
    if (IGO == 1) { goto label5; }
    //homonuclear expansion for Y(ETA)
label4:
    J = J + 2;
    ALPHA = P * P * (J - 1 - M) * (J - M) / ((J + J - 3) * (J + J - 1));
    BETA = C - J * (J + 1) - 2 * P * P * (M * M + J * (J + 1) - 1) / ((J + J + 3) * (J + J - 1));
    if (abs(ALPHA) > abs(BETA)) { goto label4; }
    FCHAIN = J - 2;
    return;
    //Heternucler Small R expansion for Y(ETA)
label5:
    J = J + 1;
    int TEST = J * J;
    if (TEST < C) { goto label5; }
    FCHAIN = J - 1;
    return;
}

//to do pack main() into function wf

void wave_function(int QU, double RINCR, int N, int L, int M, int NR, vector<double>& RR, vector<double>& PP, vector<double>& CSEP, vector<double>& GraveP, vector<double>& GraveC)
{
    //define constants, As of now not sure their functionality
    vector<int> IACDFT{ 0,11,11,10,10 };
    vector<int> IACMAX{ 0,13,13,10,10 };
    vector<int> IACMIN{ 0,10,10,10,10 };
    vector<double> IAC(5);
    vector<double> AC(5);
    double* PAC = &AC[1];
    double* CAC = &AC[2];
    double* ACCIN = &AC[3];
    double* ACCOUT = &AC[4];
    int J = 0;
    double F{ 0 };
    double G{ 0 };
    const int IIN{ 2 }, IOUT{ 3 }, IHOLD{ 1 }, IXTRAP{ 10 }, ITMAX{ 50 };
    int FCHAIN{ 0 };
    int GCHAIN{ 0 };
    int lStart{ 0 };
    double DFDC{ 0 };
    double DFDP{ 0 };
    int NCVIN{ 0 };
    double DGDC{ 0 };
    double DGDP{ 0 };
    int NCVOUT{ 0 };
    int KOUNT = 0;
    int FCUT = 0;
    int GCUT = 0;
    double P = 0;
    double C = 0;
    double DENOM = 0;
    double DELTAP = 0;
    double DELTAC = 0;
    int ISTART = 0;
    double RZ = 0;
    double RZDIF = 0;
    double DC = 0;
    double ZK = 0;
    double NS{ 0 };
    int KS{ 0 }, ICEN{ 0 };
    int NOPTS{ 1 };
    int IRISK = 0;
    int nPts = NR;

    vector <double> ICVGY(999), ICVGX(999);
    // ZA, ZB: nuclear charges;
    // N,L,M: UNITED atom quantum number


    // some subroutines require ZA >= Zb
    double ZA = 1.0, ZB = QU;
    // IGO = 1 , 2, 3 corresponds to Heteronuclear small R, Homo/Heteronuclear larger R for Y(ETA) 
    int IGO = (ZA != ZB) ? 1 : 2;
    lStart = M;
    if ((L - M) != 2 * int((L - M) / 2)) lStart = M + 1;

    double rStart{ 0.0 };
    int Z = ZA + ZB;
    for (int i = 1; i <= 4; i++) {
        if (IAC[i] < IACMIN[i] || IAC[i] > IACMAX[i]) {
            IAC[i] = IACDFT[i];
        }
        AC[i] = 1.0 / pow(10, IAC[i]);
    }

    
    // the above should be in a input file in the future
    
    for (int i = 1; i < nPts+1; ++i) {
        RR[i] = (i-1)* RINCR;
    }
    
    ISTART = NOPTS + 1;
    NOPTS += nPts;
    
    CORR(ZA, ZB, N, L, M, NS, KS, ICEN); 


    //Now we have found the dissociation products of the OEDM state

    double SIGN = pow(-1.0, (ICEN - 1));
    double ZZ = ZA * (2 - ICEN) + ZB * (ICEN - 1);
    double ZP = ZA * (ICEN - 1) + ZB * (2 - ICEN);
    int iDelta = N - L - 1 - KS;

    double rJump = 2 * NS * NS * (1 + 2 * sqrt(abs(ZP / ZZ))) / abs(ZZ);
    double W{ 0 }, W0{ 0 }, W1{ 0 }, W2{ 0 }, W3{ 0 }, W4{ 0 }, W5{ 0 }, W6{0};
    WRLG1(ZZ, ZP, NS, iDelta, M, 1.0, W, W0 , W1, W2, W3, W4, W5, W6);

    // if RJ =0, use RJUMP, else RJ
    double RJ = 1.0;
    rJump = (RJ == 0) ? rJump : RJ;

    
    // Now assume rstart = 0. Later work on the case for rstart != 0
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

    //Begin loop over R values
    // we have already initialized i = 0 case. 
    //int NR = NOPTS - ISTART + 1;
    for (int i = ISTART; i <= NR; ++i) {
        double R = RR[i];
        //IGO=1,2, OR 3 RESPECTIVELY FOR HETERONUCLEAR SMALL R, 
        //HOMONUCLEAR OR HETERONUCLEAR LARGE R EXPANSION FOR Y(ETA).
        if (IGO == 3) {
            IGO = 1;
        }
        if (IGO == 1 && R >= rJump) {
            IGO = 3;
        }
        RZ = R * Z;
        RZDIF = R * (ZA - ZB) * SIGN;
        if (i > 2) { goto label6; }
        if (i == 1) { goto label7; }
        PP[2] = 0.5 * R * Z / N;
        ZK = RZDIF * 0.5 / PP[2];
        if (L == 0) {
            DC = 2 * (1 - ZK * ZK) / 3 * pow(PP[2], 2.0);
        }
        else {// L!=0
            DC = 2 * (((L + 1) * (L + 1) - M * M) * ((L + 1) * (L + 1) - ZK * ZK) / ((L + 1) * (L + L + 3))
                - (L * L - M * M) * (L * L - ZK * ZK) / (L * (L + L - 1))) / (L + L + 1) * pow(PP[2], 2);

        }
        CSEP[2] = CSEP[1] + DC;
        goto label7;
        //extrapolate for P and C at current R
        //Use data from the previous i points if available.
label6:
        nPts = min(i,IXTRAP);
        
        EXTRAP(RR, PP, i);
        EXTRAP(RR, CSEP, i);
label7:
        P = PP[i];
        C = CSEP[i];
        
        CHAIN(FCHAIN, GCHAIN, IGO, RZ, RZDIF, N, L, M, KS, P, C);
        KOUNT = 0;
        FCUT = 0;
        GCUT = 0;
label8: 
        //cout <<i<<' '<< P << ' ' << C << " " << RZ << " " << RZDIF << ' ' << endl;
        CTDFRN(IGO, P, C, RZ, RZDIF, M, FCHAIN, lStart, *ACCIN, FCUT, F, DFDC, DFDP,  NCVIN);

        CTDFRN(4, P, C, RZ, RZDIF, M, GCHAIN, lStart , *ACCOUT, GCUT, G, DGDC, DGDP, NCVOUT);
        //BUild an array of the convergents for transfer to WFNC
        ICVGY[i] = NCVIN;
        ICVGX[i] = NCVOUT;
        DENOM = DFDC * DGDP - DFDP * DGDC;
        DELTAP = (F * DGDC - G * DFDC) / DENOM;
        DELTAC = (G * DFDP - F * DGDP) / DENOM;

        
        P +=DELTAP;
        C +=DELTAC;
        if (abs(DELTAP) / (1. + P) < *PAC and abs(DELTAC) / (1. + abs(C)) < *CAC) { goto label9; }
        KOUNT += 1;
        if (KOUNT > ITMAX) { 
            cout << "max iteration reached at R index: " << i; 
            return; }
        FCUT = 0.5 * abs(DENOM) * max(*PAC * (1 + P) / abs(DGDC), *CAC * (1 + abs(C)) / abs(DGDP));
            
        GCUT = 0.5 * abs(DENOM) * max(*PAC * (1 + P) / abs(DFDC), *CAC * (1 + abs(C)) / abs(DFDP));
        if (IRISK == 0) { goto label8; }
        if (abs(F) * (*ACCIN / FCUT) < 1.0 and abs(G) * 100 * (*ACCOUT / GCUT) < 1.0) { goto label9; }
        goto label8;
label9:
        CSEP[i] = C;
        C = -C + P * P;
        //GraveC and GraveP are for later Grave and MEDOC functions to calcualte couplings
        GraveC[i] = C;
        GraveP[i] = P;
        PP[i] = P;
    }
    //iteration over R is now complete, store energies in PP
    //PP[i] = P;
    //CSEP[i] = C;
    PP[1] = -0.5 * pow((Z / N),2);
    for (int i = ISTART; i < NOPTS; i++) {
        PP[i] = -2 * pow((PP[i] / RR[i]), 2);
    }
    //todo: print to file
    //suspicious, PP[1] is set to an arbitray number for Rstart = 0;
    PP[1] = 1e50;
    vector<double> POT;
    POT.push_back(1e50); //infinity
    int ZAB = ZA * ZB;
    for (int i = 1; i < NR; i++) {
        POT.push_back(PP[i] + ZAB / RR[i]);
    }
    string potname1;
    potname1 = "pot" + to_string(N) + to_string(L) + to_string(M);
    //ignore potential at R = 0
    std::ofstream outFilePot1(potname1);
    if (outFilePot1.is_open()) {
        for (int i = ISTART; i < NR; ++i) {
            outFilePot1 << RR[i] << "\t" << POT[i] << endl;
        }

        outFilePot1.close(); // Close the file
    }
    else {
        std::cerr << "Unable to open file";
    }
}





struct CommonBlockTRAP {
    //used to transfer quantities defining the wave fucntion to different subroutines in GRAVE
    double QU, R,  IVER, E, PE, A, NG, NF;
    int ME, N, L;
    std::vector<double> XG = vector<double>(41);
    std::vector<double> XF = vector<double>(41);

};


struct CommonBlockINIT {
    double SIGMA, EM, CEM, PE2, RQU, ME2, MEC;
};


struct CommonBlockMAT {
   
    vector<double> A = vector<double>(201);
    vector<double> AM = vector<double>(201);
    vector<double> INT = vector<double>(41);
};

vector<vector<double>> array2MatrixColumnWise(const vector<double>& vec, int n) {
    int N = vec.size()-1;
    // Calculate the number of columns. If N is not perfectly divisible by n, add 1 more column to fit all elements.
    int m = (N % n == 0) ? (N / n) : (N / n + 1);

    vector<vector<double>> matrix(n, vector<double>(m, 0)); // Initialize matrix with zeros

    for (int i = 1; i <= N; ++i) {
        // Calculate the row and column index based on column-wise filling
        int row = (i-1) % n;
        int col = (i-1) / n;
       
        matrix[row][col] = vec[i];
    }

    return matrix;
}

vector<double> flattenMatrixColumnWise(const vector<vector<double>>& matrix) {
    if (matrix.empty()) return {};

    int n = matrix.size(); // Number of rows
    int m = matrix[0].size(); // Assuming all rows have the same number of columns
    vector<double> flattened;
//placeholder for 0-indexed array
    flattened.push_back(0);
    for (int col = 0; col < m; ++col) {
        for (int row = 0; row < n; ++row) {
            flattened.push_back(matrix[row][col]);
        }
    }

    return flattened;
}

//NN, 0, M2, 0, 1, A, AM, INT, BG
void BANSOL(int N, int M1, int M2, int IE, int IR, vector<double>& AA, vector<double>& MM, vector<double>& INT, vector<double>& BB){
    vector<vector<double>> A = array2MatrixColumnWise(AA, N);
    vector<vector<double>> M = array2MatrixColumnWise(MM, N);
    vector<vector<double>> B = array2MatrixColumnWise(BB, N);
    double X;
    int W,I, IK,KK;
    int L = M1;
    if (IE != 0) { goto label31; }
    for (int K = 1;K<=N;K++){
        I = int(K);
        if (I == K) { goto label11; }
            
        X = B[K-1][0];
        B[K-1][0] = B[I-1][0];
        B[I-1][0] = X;
    label11:
        if (L < N) { L = L + 1; }
        I = K + 1;
    label30:
        if (I > L) { break; }
        IK = I - K;
        X = M[K-1][IK];
        B[I-1][0] = B[I-1][0] - X * B[K-1][0];
        I = I + 1;
        goto label30;
    }
label31:


    L = -M1;
    for (int II = 1; II <= N; II++) {
        I = N - II + 1;
        X = B[I-1][0];
        W = I + M1;
        int K = 1 - M1;
    label62:
        if (K > L) { goto label61; }
        IK = K + M1 + 1;
        KK = K + W;
        X = X - A[I-1][IK-1]* B[KK-1][0];
        K = K + 1;
        goto label62;
    label61:
        B[I-1][0] = X / A[I-1][0];
        if (L < M2) { L += 1; }
    }
    AA = flattenMatrixColumnWise(A);
    BB = flattenMatrixColumnWise(B);
    MM = flattenMatrixColumnWise(M);



    return;
}


void BANDET(int N, int M1, int M2, int IE, double LAM, double MAC, vector<double>& A, double& D1, int& ID2, vector<double>& M, vector<double>& INT, bool &FAIL) {
    /* PROGRAM BANDET1 (Handbook for automatic computation Vol II, I/6)
    * Calculates the determinant of a Band Matrix and put the matrix in a form suitable for use by BANSOL
    * N: dimension of the matrix
    * M1: number of subdiagonal lines
    * M2: number of Superdiagonal lines
    * IE=1 if the determinant is requred and 0 otherwise
    * LAM: Calculations are performed on A-A-LAM*I instead of A
    * MAC: smallest number s.t. 1.0+MAC!=1
    * A: Matrix of dimension N*(M1+M2+1) 
    * D1, ID2: give the value of the determinant as D1*2^ID2
    * M: Real working matrix of the size N*(M1+M2+1)
    * INT: Working array if dimension N
    * FAIL: if BANSOL fails, FAIL=True
    */


    //probably no need to convert A to a column vector
    double DABS,  X,  NORM;
    vector<vector<double>> AA = array2MatrixColumnWise(A,N);
    vector<vector<double>> MM = array2MatrixColumnWise(M,N);
    FAIL = false;
    int JJ, JP;
    int MSUP = M1+M2 + 1;
    NORM = 0.0;
    int J;
    for (int i = 1; i <= N; ++i) {
        AA[i-1][0] = AA[i-1][0] - LAM;
    }
    int M0 = M1 + 1;
    int L = M1;
    int I = 1;
label42:
    if (I > M1) { goto label41; }
    J = 1 - I;
label38:
    while (J <= M2) {
        JJ = J - L + M0;
        JP = J + M0;
        AA[I-1][JJ - 1] = AA[I-1][JP - 1];
        J += 1;
    }
label39:
    L -= 1;
    J = M2 - L;
label48:
    if (J > M2) { goto label49; }
    JJ = J + M0;
    AA[I-1][JJ-1] = 0.0;
    J = J + 1;
    goto label48;
label49:
    I += 1;
    goto label42;
label41:
    D1 = 1.0;
    ID2 = 0;
    L = M1;
    for (int K = 1; K <= N; ++K) {
        X = AA[K-1][0];
        I = K;
        if (L < N) { L += 1; }
        J = K + 1;
    label52:
        if (J > L) { goto label51; }
        if (abs(AA[J-1][0]) <= abs(X)) { goto label53; }
        X = AA[J-1][0];
        I = J;
    label53:
        J += 1;
        goto label52;
    label51:
        K = I;
        D1 = D1 * X;
        if (X != 0) { goto label54; }
        ID2 = 0;
        if (IE == 1) { goto label55; }
        FAIL = true;
    label55:
        AA[K-1][0] = NORM * MAC;
    label54:
        if (D1 == 0) { goto label56; }
    label58:
        if (abs(D1) < 1.0) { goto label57; }
        ID2 = ID2 + 4;
        D1 = D1 * 0.0625;
        goto label58;
    label57:
        if (abs(D1) < 1.0) { goto label56; }
        ID2 = ID2 - 4;
        D1 = D1 * 16.0;
        goto label57;
    label56:
        if (I == K) { goto label61; }
        D1 = -D1;
        X = AA[K-1][0];
        AA[K-1][0] = AA[I-1][0];
        AA[I-1][0] = X;
    label61:
        I = K + 1;
    label63:
        if (I > L) { goto labelend; }
        JJ = I - K;
        M[K-1] = AA[I-1][0] / AA[K-1][0];
        X = M[K-1];
        J = 2;
    label68:
        if (J > MSUP) { goto label69; }
        JJ = J - 1;
        AA[I-1][0] = AA[I-1][0] - X * AA[K-1][0];
        J += 1;
        goto label68;
    label69:
        AA[I-1][0] = 0;
        I += 1;
        goto label63;

    }

labelend:
    A = flattenMatrixColumnWise(AA);
    M = flattenMatrixColumnWise(MM);
    return;
}


void DGE(int N, double AG, double PE, vector<double> &BG, CommonBlockINIT &INIT, CommonBlockMAT &MAT) {
    /*************************************************************************
    Calculates the tridiagonal matrix G and solves the system of linear equations
    which gives the coefficient of the expansion of the outer wavefunction
    *************************************************************************/
    double AA, D1{0}, HI;
    const double ZERO = 0.0;
    const double MACH = 2.3e-16; // MACH is machine dependent, can make it smaller in the future
    bool FAIL;
    int NN = N - 1;
    int KS, NS, M2;
    int J = 0;
    int ID2 = 0;


    for (int i = 1; i <= NN; i++) {
        J += 1;
        AA = double(i - 1) - INIT.SIGMA;
        MAT.A[J] = AA * (AA - INIT.EM);
    }
    AA = INIT.SIGMA - 2.0 * PE;
    if (NN == 1) { goto label21; }
    KS = NN - 1;
    for (int i = 1; i <= KS; ++i) {
        J += 1;
        HI = double(i);
        MAT.A[J] = 2.0 * HI * (AA - HI) + INIT.CEM - AG;
    }
    J += 1;
    MAT.A[J] = 0;
    if (NN == 2) { goto label31; }
    KS = NN - 2;
    for (int i = 1; i <= KS; i++) {
        J = J + 1;
        HI = double(i + 1);
        MAT.A[J] = HI * (HI + INIT.EM);
    }
    J += 1;
    MAT.A[J] = 0;
    J += 1;
    MAT.A[J] = 0;
    M2 = 2;
label32:
    BANDET(NN, 0, M2, 0, ZERO, MACH, MAT.A, D1, ID2, MAT.AM, MAT.INT, FAIL=false);
    if (FAIL) { goto label99; }
    if (NN == 2) { goto label51; }
    NS = NN - 2;
    for (int i = 0; i <= NS; i++) {
        BG[i] = 0;
    }
label51:
    NS = NN - 1;
    HI = double(NS + 1);
    BG[NS] = -HI * (HI + INIT.EM);
    HI = double(NN);
    BG[NN] = -2.0 * HI * (AA - HI) - INIT.CEM + AG;
    BANSOL(NN, 0, M2, 0, 1, MAT.A, MAT.AM, MAT.INT, BG);
    BG[N] = 1.0;
    goto labelreturn;
label31:
    M2 = 1;
    goto label32;
label21:
    MAT.A[2] = 2.0 * (AA - 1) + INIT.CEM - AG;
    BG[1] = -MAT.A[2] / MAT.A[1];
    BG[2] = 1.0;
    goto labelreturn;
label99:
    cout << "ERROR. QUIT" << endl;
labelreturn:
    return;
}

void DFE(int N, double AF, double& PE, vector<double> BF, CommonBlockINIT &INIT, CommonBlockTRAP &TRAP, CommonBlockMAT &MAT) {
    //Calculates the pentadiagonal matrix F and solves the systems of linear eqns which give the coeff
    //of the expansion of the inner wave-function

    //for asymmetric system only
    double D1,  I2, IME,JJ,JNU,JDEN;
    const double MACH = 2.3e-16;
    bool FAIL=false;
    int NN, J, KS, M2, NS,ID2;
    const double ZERO = 0.0;
    NN = N - 1;
    if (NN == 1) { goto label11; }
    J = 1;
    MAT.A[1] = 0;
    for (int I = 2; I <= NN; I++) {
        J += 1;
        I2 = I + I + INIT.ME2 - 1;
        MAT.A[J] = INIT.PE2 * double(I * (I - 1)) / double(I2 * (I2 - 2));
    }
    for (int I = 1; I <= NN; I++) {
        J += 1;
        MAT.A[J] = INIT.RQU * float(I) / float(INIT.ME2 + I + I - 1);
    }
    KS = NN - 1;
    for (int I = 1; I <= KS; I++) {
        J = J + 1;
        IME = I + TRAP.ME;
        JJ = IME * (IME + 1);
        JNU = JJ + JJ - INIT.MEC;
        I2 = I + I + INIT.ME2 - 1;
        JDEN = I2 * (I2 + 4);
        MAT.A[J] = -JJ + INIT.PE2 * JNU / JDEN - AF;
    }
    J = J + 1;
    MAT.A[J] = 0.0;
    if(NN == 2) goto label31;
    KS = NN - 2;
    for (int I = 1; I <= KS; I++) {
        J = J + 1;
        JNU = INIT.ME2 + I + 1;
        MAT.A[J] = INIT.RQU * JNU / (JNU + I + 2);
    }
    J = J + 1;
    MAT.A[J] = 0.0;
    J = J + 1;
    MAT.A[J] = 0.0;
    if (NN == 3) { goto label41; };
    KS = NN - 3;
    for (int I = 1; I <= KS; I++) {
        J += 1;
        JNU = INIT.ME2 + I + 1;
        JDEN = JNU + I + 2;
        MAT.A[J] = INIT.PE2 * (JNU * (JNU + 1)) / (JDEN * (JDEN + 2));

    }
    J = J + 1;
    MAT.A[J] = 0.0;
    J = J + 1;
    MAT.A[J] = 0.0;
    J = J + 1;
    MAT.A[J] = 0.0;
    M2 = 3;
    NS = NN - 3;
    for (int i = 1; i <= NS; i++) {
        BF[i] = 0.0;
    }
label42:
    NS = NN - 2;
    JNU = INIT.ME2 + NS + 1;
    JDEN = JNU + NS + 2;
    BF[NS] = -INIT.PE2 * (JNU * (JNU + 1)) / (JDEN * (JDEN + 2));
label32:
    NS = NN - 1;
    JNU = INIT.ME2 + NS + 1;
    BF[NS] = -INIT.RQU * JNU / (JNU + NS + 2);
    IME = NN + TRAP.ME;
    JJ = IME * (IME + 1);
    I2 = INIT.ME2 + NN + NN - 1;
    JNU = JJ + JJ - INIT.MEC;
    JDEN = (I2 + 4) * I2;
    BF[NN] = (JJ) - INIT.PE2 * (JNU) / (JDEN) + AF;
    BANDET(NN, 1, M2, 0, ZERO, MACH, MAT.A, D1, ID2, MAT.AM, MAT.INT, FAIL);
    if (FAIL) { cout << "FAIL in subroutine DGE" << endl; }
    BANSOL(NN, 1, M2, 0, 1, MAT.A, MAT.AM, MAT.INT, BF);
    BF[N] = 1.0;
    goto labelend;
label41:
    M2 = 2;
    goto label42;
label31:
    goto label32;
label11:
    MAT.A[1] = INIT.RQU /double(INIT.ME2 + 1);
    JJ = (TRAP.ME + 1) * (TRAP.ME + 2);
    JNU = JJ + JJ - INIT.MEC;
    JDEN = (INIT.ME2 + 5) * (INIT.ME2 + 1);
    MAT.A[2] = -double(JJ) + INIT.PE2 * JNU / JDEN - AF;
    BF[1] = -MAT.A[2] / MAT.A[1];
    BF[2] = 1.0;
labelend:
    return;
}

void DFSYM(int N,double AF,double PE, vector<double> &BF, int NPAR, CommonBlockINIT &INIT, CommonBlockTRAP &TRAP, CommonBlockMAT &MAT) {
    //Calcualtes the traditional Matrix F and solves the system of linear eqns that gives the
    //coeff of the inner wave function
    //For QU = 1 (symmetrical system) only.
    double D1,JNU,JDEN;
    bool FAIL;
    const double MACH = 2.3e-16;
    const double ZERO = 0.0;
    int I,II,I2,NN, NS, J, JJ,KS,IME,M2,ID2;
    NN = N - 1;
    if (NN == 1) { goto label11; }
    J = 0;
    for (int ii = 1; ii <= NN; ii++) {
        I = 2 * ii + NPAR;
        J += 1;
        I2 = I + I + INIT.ME2 - 1;
        MAT.A[J] = INIT.PE2 * (I * (I - 1)) / (I2 * (I2 - 2));
    }
    KS = NN - 1;
    for (int II = 1; II <= KS; II++) {
        J = J + 1;
        I = 2 * II + NPAR;
        IME = I + TRAP.ME;
        JJ = IME * (IME + 1);
        JNU = JJ + JJ - INIT.MEC;
        I2 = I + I + INIT.ME2 - 1;
        JDEN = I2 * (I2 + 4);
        MAT.A[J] = -double(JJ) + INIT.PE2 * JNU / JDEN - AF;
    }
    J += 1;
    MAT.A[J] = 0.0;
    if (NN == 2) { goto label31; }
    KS = NN - 2;
    for (int II = 1; II <= KS; ++II) {
        I = 2 * II + NPAR;
        J = J + 1;
        JNU = INIT.ME2 + I + 1;
        JDEN = JNU + I + 2;
        MAT.A[J] = INIT.PE2 * (JNU * (JNU + 1)) / (JDEN * (JDEN + 2));
    }
    J += 1;
    MAT.A[J] = 0;
    J += 1;
    MAT.A[J] = 0;
    M2 = 2;
label32:
    BANDET(NN, 0, M2, 0, ZERO, MACH, MAT.A, D1, ID2, MAT.AM, MAT.INT, FAIL);
    if (FAIL) { cout << "FAIL inside function DFE" << endl; }
    if (NN == 2) { goto label51; }
    KS = NN - 2;
    for (int i = 1; i <= KS; i++) {
        BF[i] = 0.0;
    }
label51:
    KS = NN - 1;
    I = 2 * KS + NPAR;
    JNU = INIT.ME2 + I + 1;
    JDEN = JNU + I + 2;
    BF[KS] = INIT.PE2 * (JNU * (JNU + 1)) / (JDEN * (JDEN + 2));
    KS = NN;
    I = KS + KS + NPAR;
    IME = I + TRAP.ME;
    JJ = IME * (IME + 1);
    JNU = JJ + JJ - INIT.MEC;
    I2 = I + I + INIT.ME2 - 1;
    JDEN = I2 * (I2 + 4);
    BF[KS] = -double(JJ) + INIT.PE2 * JNU / JDEN - AF;
    BANSOL(NN, 0, M2, 0, 1, MAT.A, MAT.AM, MAT.INT, BF);
label12:
    BF[N] = 1.0;
    for (int I = 2; I <= N; ++I) {
        II = N - I + 2;

        KS = II + II + NPAR - 1;
        BF[KS] = BF[II];
        KS = KS - 1;
        BF[KS] = 0.0;
    }
    BF[NPAR + 1] = BF[1];
    if (NPAR == 1) { BF[1] = 0.0; }
    goto labelend;
label31:
    M2 = 1;
    goto label32;
label11:
    I = 2 + NPAR;
    I2 = I + I + INIT.ME2 - 1;
    MAT.A[1] = INIT.PE2 * double(I * (I - 1)) / double(I2 * (I2 - 2));
    IME = I + TRAP.ME;
    JJ = IME * (IME + 1);
    JNU = JJ + JJ - INIT.MEC;
    JDEN = I2 * (I2 + 4);
    MAT.A[2] = -double(JJ) + INIT.PE2 * JNU / JDEN - AF;
    BF[1] = -MAT.A[2] / MAT.A[1];
    goto label12;
labelend:
    return;
}


tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> GRAVE(const double QU, int N, int L, int ME, int NR, vector<double> RR, vector<double> GraveP, vector<double> GraveC) {
    double AG, DABS, XF1;
    int NG, NF, MYINDEX;
    int JWRIT, NGMAX, NFMAX, NGIN, NFIN, NSTATE, XGMAX, XFMAX, NFPMAX;
    bool NGTEST, NFTEST;
    int NPAR, NCPAR;
    int JPAR;
    int NFF;
    vector<double> SWITCH(999);
    vector<double> AOUT(88);
    // &QU = AOUT, not sure why placeholder for future review.
    //output array for MEDOC
    vector<vector<double>> params;
    vector<vector<double>> XGRes;
    vector<vector<double>> XFRes;

    const int JAOUT = 88;


    CommonBlockTRAP TRAP = {};
    TRAP.QU = QU;
    TRAP.N = N;
    TRAP.L = L;
    TRAP.ME = ME;
    CommonBlockINIT INIT = {};
    CommonBlockMAT MAT = {};
    //declare reference objects for simplicity
    double& SIGMA = INIT.SIGMA;
    double& EM = INIT.EM;
    double& CEM = INIT.CEM;
    double& PE2 = INIT.PE2;
    double& RQU = INIT.RQU;
    double& ME2 = INIT.ME2;
    double& MEC = INIT.MEC;

    //reference for TRAP
    /*
    
    
    int& N = TRAP.N;
    int& L = TRAP.L;
    int& ME = TRAP.ME;
    double& E = TRAP.E;
    
    */
    double& PE = TRAP.PE;
    double& A = TRAP.A;
    double& R = TRAP.R;
    double& E = TRAP.E;
    const int IVER = 2;
    //for now, ignore all reads in and set parameters directly. Need to be more dynamic in the future. 
    /****************************************************************************************************************************** 
    ******************************************************************************************************************************
    ******************************************************************************************************************************
    JWRIT: The Coeffs of the semi-analytic expansions of the WF is pringted if JWRIT!=0.

    XGMAX: THE TRUNCATION OF THE SEMI-ANALYTIC EXPANSION OF THE OUTER WAVE-FUNCTION IS SUCH THAT THE RATIO OF THE FIRST COEFFICIENT
    TO THE LAST IN THE EXPANSION IS GREATER THAN XGMAX.

    XFMAX: SAME AS XFMAX FOR THE INNER WAVE-FUNCTION.

    NGMAX: MAXIMUM NUMBER OF TERMS IN THE EXPANSION OF THE OUTER WAVE-FUNCTION.

    NFMAX: SAME AS NGMAX FOR INNER WAVE-FUNCTION.

    NGIN: INITIAL GUESS FOR THE NUMBER OF TERMS IN THE EXPANSION OF THE OUTER WAVE-FUNCTION.
    NFIN: SAME AS NFIN FOR INNER WAVER-FUNCTION.

    NSTATE: NUMBER OF STATES FOR WHICH THE CALCULATIONS ARE TO BE PERFORMED.
     ******************************************************************************************************************************
    ******************************************************************************************************************************
    *******************************************************************************************************************************
    *******************************************************************************************************************************/

    JWRIT = 1;
    XGMAX = 10000;
    XFMAX = 10000;
    NGMAX = 80;
    NFMAX = 80;
    NGIN = 7;
    NFIN = 8;
    NSTATE = 1;

    NG = NGIN;
    NF = NFIN;

    INIT.EM = float(ME);
    JPAR = 22;
    NFPMAX = NFMAX;
    if (abs(QU - 1.0) >= 2e-5) { goto label101; }
    NPAR = 0;
    JPAR = 220;
    if (pow(-1.0, L + ME) < 0) { NPAR = 1; }
    NCPAR = 1 - NPAR;
    NFPMAX = 2 * ((NFMAX + NCPAR) / 2) - NCPAR;
label101:
    INIT.ME2 = 2*ME;
    INIT.MEC = 2 * ME * ME + 1;
    //start looping on internuclear distance
    // i = 2 because GraveP, GraceC and RR start with the second index
    for (int i = 1; i < NR; i++) {
        TRAP.R = RR[i];
        TRAP.PE = GraveP[i];
        TRAP.A = GraveC[i];
        INIT.RQU = TRAP.R * (QU - 1.0);
        INIT.PE2 = PE * PE;
        AG = -TRAP.A;
        E = -2 * PE * PE / (TRAP.R * TRAP.R);
        NGTEST = true;
        NFTEST = true;
        SIGMA = TRAP.R * (1.0 + QU) / (2.0 * TRAP.PE) - 1.0 - EM;
        CEM = EM * (EM + SIGMA + 1.0) + SIGMA * (1.0 + 2.0 * PE) - INIT.PE2;
        DGE(NG, AG, TRAP.PE, TRAP.XG, INIT, MAT);
    labelJPAR:

        switch (JPAR) {
             case 22:
                DFE(NF, TRAP.A, TRAP.PE, TRAP.XF, INIT,  TRAP, MAT);
            case 220:
                if (NF < 4) { NF = 4;}
                NFF = (NF + NCPAR) / 2;
                NF = 2 * NFF - NCPAR;
                DFSYM(NFF, A, PE, TRAP.XF, NPAR, INIT, TRAP, MAT);       
        }
        //check if enough terms in series to achieve required precision
    label20:
        if (abs(TRAP.XG[1]) > XGMAX) { goto label11; }
        NGTEST = false;
        if (NG < NGMAX) { goto label19; }
        goto label11;
    label19:
        XF1 = TRAP.XF[1];
        if (abs(XF1) < 1e-10) { XF1 = TRAP.XF[2]; }
        if (abs(XF1) > XFMAX) { goto label12; }
        NGTEST = false;
        if (NF < NFPMAX) { goto label16; }
        goto label12;
    label11:
        XF1 = TRAP.XF[1];
        if(abs(XF1)<1e-10) XF1 = TRAP.XF[2];
        if(abs(XF1)>XFMAX) goto label13;
        NFTEST = false;
        if (NF < NFPMAX) {
            goto label17;
        };
        goto label13;
    label12:
        NG += 2;
        if (NG > NGMAX) { NG = NGMAX; }
        DGE(NG, AG, PE, TRAP.XG, INIT, MAT);
        goto label20;
    label16:
        NG += 2;
        if (NG > NGMAX) { NG = NGMAX; }
        NF += 2;
        if (NF > NFPMAX) { NF = NFPMAX; }
        DGE(NG, AG, PE, TRAP.XG, INIT, MAT);
        goto labelJPAR;
    label17:
        NF += 2;
        if (NF > NFPMAX) { NF = NFPMAX; }
        goto labelJPAR;
    label13:
        //printing begins
        //cout << E << " " << PE << " " << A << " " << NG << " " << NF << endl;
        //cout << NF<<" XF "<<TRAP.XF[1] << " " << TRAP.XF[2] << " " << TRAP.XF[3] << " " << TRAP.XF[4] << " " << TRAP.XF[5] << endl;
        //cout <<NG<<" XG "<<TRAP.XG[1] << " " << TRAP.XG[2] << " " << TRAP.XG[3] << " " << TRAP.XG[4] << " " << TRAP.XG[5] << endl;
        //output results
        vector<double> curParam = { E,PE,A,double(NG),double(NF) };
        params.push_back(curParam);
        
        XGRes.push_back(TRAP.XG);
        XFRes.push_back(TRAP.XF);

        if (NGTEST and NG != 2) { NG -= 1; }
        if (NFTEST and NF != 2) { NF -= 1; }
    }
    cout << "GRAVE complete successfully" << endl;
    return make_tuple(params,XGRes,XFRes);
}

struct Gauss
{ //can be expanded with more terms in the future
    int NXI = 20;
    //the first term in XI and OMEGA are place holder for simplicity, because of 0-indexed and 1-indexed. 
    vector<double> XI{ 0,0.07053989,0.372126818,0.916582102,1.707306531,2.749199255,4.048925314,5.615174971,7.459017454,9.59439287,12.03880255,14.81429344,17.94889552,21.47878824,25.45170279,
            29.93255463,35.01343424,40.83305706,47.61999405,55.81079575,66.52441653 };
    vector<double> OMEGA{0, 0.181080062,0.422556768,0.666909547,0.915352373,1.169539707,1.431354986,1.702981138,1.987015891,2.286635781,2.605834728,2.949783734,3.325395782,
            3.742255471,4.21423671,4.762518461,5.421726044,6.254012357,7.387314389,9.151328731,12.89338865 };

};

void NORDIF(int ME, double PE, double SIGMA, int NG, vector<double>&XG, int NF, vector<double>&XF, double& FNOR, double R) {
    //calcualte the norm of a given wavefunciton
    //in GRAVE, the last term in semi-analytical expansion is set equal to 1. Here the first non-zero term is equal to 1
    //The norm is calculated for this new form of expansion.
    double ALAM, ALAMXI, ALAM1, ALAM2, ALPHA, BETA, DEXP, DLOG,  FONC, PP,SAM, SG, SIM, SPG, SUM, TROC, TRUC, XF2, XN, Y;
    double IS, IM, IME, NFP;

    Gauss Gauss;
    int NXI = Gauss.NXI;
    vector<double> XI = Gauss.XI;
    vector<double> OMEGA = Gauss.OMEGA;
    XN = XG[1];
    for (int K = 1; K <= NG;K++) {
        XG[K] = XG[K] / XN;
    }
    XN = XF[1];
    if (XN == 0) { XN = XF[2]; }
    for (int k = 1; k <= NF; k++) {
        XF[k] = XF[k] / XN;
    }
    //integration over mu in spheridol prolate coordinate
    double IALPHA = 1;
    int ME2 = 2 * ME;
    for (int K = 1; K <= ME2; K++) {
        IALPHA = IALPHA * K;
    }
    SAM = 0;
    SUM = 0;
    ALPHA = IALPHA;
    for (int IT = 1; IT <= NF; IT++) {
        IS = IT - 1;
        XF2 = XF[IT] * XF[IT];
        IM = ME2 + IS;
        IME = IM + IS + 1;
        SUM = SUM + ALPHA * 2.0 * XF2 / IME;
        TROC = ((IT) * (IM + 1)) / (IME + 2);
        TROC = TROC + (IS * IM) / (IME - 2);
        TROC = TROC / (IME * IME);
        SAM = SAM + ALPHA * 2.0 * XF2 * TROC;
        ALPHA = ALPHA * (IM + 1) /IT;
    }
    BETA = (IALPHA * (ME2 + 1) * (ME2 + 2));
    SIM = 0;
    if (NF <= 2) { goto label21; }
    NFP = NF - 2;
    for (int IT = 1; IT <= NFP; IT++) {
        IS = IT - 1;
        IME = ME2 + IS + IS + 1;                                                  
        SIM = SIM + 4.0 * XF[IT] * XF[IT + 2] * BETA / (IME * (IME + 2) * (IME + 4));
        BETA = BETA * (ME2 + IS + 3) /double(IT);
    }
label21:
    SAM = SAM + SIM;
    //integration over Lambda coord
    SG = 0.0;
    SPG = 0.0;
    PP = 2.0 * PE;
    for (int N = 1; N <= Gauss.NXI; N++) {
        ALAM = XI[N] / PP + 1.0;
        ALAM1 = ALAM - 1.0;
        ALAM2 = ALAM + 1.0;
        ALAMXI = ALAM1 / ALAM2;
        FONC = exp(-PE * ALAM + SIGMA * log(ALAM2));
        //(FONC = (lam+1)^sig*exp(-p*lam))
        FONC = FONC * FONC;
        TRUC = ALAM1 * ALAM2;
        //TRUC = (lam+1)*(lam-1)

        TRUC = pow(TRUC, ME);
        //^m/2*^m/2 = ^m
        FONC = FONC * TRUC;
        Y = 0;
        for (int JN = 1; JN <= NG; JN++) {
            Y += XG[JN] * pow(ALAMXI, JN - 1);
        }
        FONC = FONC * Y * Y * OMEGA[N];
        SG += FONC;
        SPG += FONC * ALAM * ALAM;
    }
    SG = SG / PP;
    SPG = SPG / PP;
    FNOR = R * R * R * (SPG * SUM - SG * SAM) / 8.0;
    FNOR = sqrt(1 / FNOR);
    return;
}

void ROTDIS(double& CROT, double& GAMROT, Gauss Gaussian, vector<double> XFP, vector<double> XGP, int NFP, int NGP, vector<double> XF, vector<double> XG, int NF, int NG, double R, double PE, double PEP, double SIGMA, double SIGMAP, double DELTA, int L, int ME, double OMEZUT) {
    double  AI, ALAM, ALAMXI, ALAM1, ALAM2, ALPHA,  BI,  CI, COEF, COEFF, DEXP, DLOG,  DY, FISC, FONC, GDEU,GPDEU, GPUN, 
        GTROI, GUN, SAM, SG, SIM, SJM, SKM, SP, SPG, SPI, SPJ, SPK, SPN, SUMP1, SUMP2, SUM1, SUM2,TRIC, TRUC,  Y, YP, SPK2;
    double SUM{ 0 }, SUMP{ 0 };
    int IBMAX, IN,NMAX,IS,IS2,MIM,IM,IME,IT,NPMAX,NN;
    auto A = [](int I, int I2, int ME, int ME2) -> double {
        return -2.0*double(I*(ME2+I+2.))/double((ME2+I2-1.)*(ME2+I2+1.)*(ME2+I2+3.));
        };

    auto B = [](int I, int I2, int ME, int ME2) -> double {
        return 2.0*double((ME2+I+2)*(-4*ME*ME-4*ME*(I+2)-2*I-3))/double((ME2+I2+3)*(ME2+I2+3)*(ME2+I2+5)*(ME2+I2+1));
        };

    auto C = [](int I, int I2, int ME, int ME2) -> double {
        return 2*double((ME2+I+3)*(ME2+I+2))/double((ME2+I2+7)*(ME2+I2+3)*(ME2+I2+5));
        };
    int ME2 = ME * 2;
    int IALPHA = 1;
    /*
    Start Iegrating over
    *  
     B and B''
    *   
    */
    for (int k = 1; k <= ME2; k++) {
        IALPHA = IALPHA * k;
    }
    ALPHA = IALPHA;
    GDEU = -2.0 * ALPHA * double((ME + 2) * (ME2 + 1)) / double (ME2 + 3);
    GPDEU = -2.0 * ALPHA * double(ME2 + 1) / double(ME2 + 3);                   
    SP = 0.0;                  
    IBMAX = (NFP - 1) / 2 + 1;
    for (int IB = 1; IB <= IBMAX; IB++) {
        IN = 2 * IB - 1;
        SP = SP + XFP[IN];
    }
    SUM2 = XF[2] * GDEU * SP;
    SUMP2 = XF[2] * GPDEU * SP;                                              
    ALPHA = ALPHA * double(ME2 + 1);                                    
    SUM1 = 0.0; 
    SUMP1 = 0.0; 
    NMAX = min(NFP, NF - 1);
    if (NMAX < 2) { goto label21; }
    for (int IT = 2; IT <= NMAX; IT++) {
        IS = IT - 1;                                                           
        IS2 = 2 * IS;
        MIM = ME + IS;
        IM = MIM + ME + 1;
        IME = IM + IS;
        GUN = ALPHA * 2.0 * double((MIM - 1) * IS) / double(IME - 2);
        GPUN = -ALPHA * double(IS2) / double(IME - 2);
        GDEU = -ALPHA * 2.0 * double((MIM + 2) * IM) / double(IME + 2);
        GPDEU = -2.0 * ALPHA * double(IM) / double(IME + 2);
        SP = 0.0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP = SP + XFP[IN];
        }
        SUM1 = SUM1 + XF[IT - 1] * GUN * SP;
        SUMP1 = SUMP1 + XF[IT - 1] * GPUN * SP;
        SUM2 = SUM2 + XF[IT + 1] * GDEU * SP;
        SUMP2 = SUMP2 + XF[IT + 1] * GPDEU * SP;
        ALPHA = ALPHA * double(IM) / double(IT);
    }
label21:
    if (NFP - NF < 0) {
        goto label24;
    }
    else if (NFP - NF == 0) {
        goto label27;
    }
    else {
        goto label26;
    }
label27:
    IT = NF;                                                             
    IS = IT - 1;
    IS2 = 2 * IS;
    GUN = ALPHA * 2.0 * double((ME + IS - 1) * IS) / double(ME2 + IS2 - 1);
    GPUN = -ALPHA * double(IS2) / double(ME2 + IS2 - 1);
    SP = XFP[IT];
    SUM1 = SUM1 + XF[IT - 1] * GUN * SP;
    SUMP1 = SUMP1 + XF[IT - 1] * GPUN * SP;
    goto label24;
label26:
    IT = NF;
    IS = IT - 1;                                                        
    IS2 = 2 * IS;                                                 
    GUN = ALPHA * 2.0 * double((ME + IS - 1) * IS) / double(ME2 + IS2 - 1);
    GPUN = -ALPHA * double(IS2) / double(ME2 + IS2 - 1);
    SP = 0.0;
    IBMAX = (NFP - IT) / 2 + 1;
    for (int IB = 1; IB <= IBMAX; IB++) {
        IN = 2 * IB - 2 + IT;
        SP += XFP[IN];
    }
    SUM1 = SUM1 + XF[IT - 1] * GUN * SP;
    SUMP1 = SUMP1 + XF[IT - 1] * GPUN * SP;
    ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
    IT = NF + 1;
    IS = IT - 1;
    IS2 = 2 * IS;
    GUN = ALPHA * 2.0 * double((ME + IS - 1) * IS) / double(ME2 + IS2 - 1);
    GPUN = -ALPHA * double(IS2) / double(ME2 + IS2 - 1);
    SP = 0.0;
    IBMAX = (NFP - IT) / 2 + 1;
    for (int IB = 1; IB <= IBMAX; IB++) {
        IN = 2 * IB - 2 + IT;
        SP += XFP[IN];
    }
    SUM1 += XF[IT - 1] * GUN * SP;
    SUMP += SUMP1 + SUMP2;
label24:
    SUM = SUM1 + SUM2;
    SUMP = SUMP1 + SUMP2;
    /*
        Start Iegrating over
        *
        B'
        *
    */
    ALPHA = double(IALPHA);
    if(NF<4) XF[4] = 0.0;
    if(NF<2) XF[2] = 0.0;                                            
    if(NFP<2) XFP[2] = 0.0;                                         
    COEF = double(ME2 + 1);                                               
    SAM = COEF * XFP[1] * (B(0, 0,ME,ME2) * XF[2] + C(0, 0,ME,ME2) * XF[4]) * ALPHA; 
    ALPHA = ALPHA * double(ME2 + 1);                                        
    NMAX = min(NFP, NF - 3);                                             
    if(NMAX<2) goto label59;      
    for (int IT = 2; IT <= NMAX; IT++) {
        IS = IT - 1;
        MIM = ME + IS;                                                         
        IM = ME + MIM + 1;                                                      
        IME = IM + IS;                                                         
        AI = -2.0 * double(IS * (IM + 1)) / double((IME - 2) * IME * (IME + 2));            
        BI = 2.0 * double((IM + 1) * (-4 * ME * ME - 4 * ME * (IS + 2) - 2 * IS - 3)) / double((IME + 2) * (IME + 2) * (IME + 4) * IME);
        CI = 2.0 * double((IM + 2) * (IM + 1)) / double((IME + 6) * (IME + 4) * (IME + 2));
        COEF = double(IM);
        SAM = SAM + COEF * XFP[IT] * (AI * XF[IT - 1] + BI * XF[IT + 1] + CI * XF[IT + 3]) * ALPHA;
        ALPHA = ALPHA * double(IM) / double(IT);
    }
label59:
    if (NF-3 < 0) {
        goto label64;
    }
    else if (NF - 3 == 0) {
        goto label67;
    }
    else {
        goto label68;
    }
label68:
    if (NMAX == NFP) goto label61;
    if (NFP - NF+1 < 0) {
        goto label62;
    }
    else if (NFP - NF + 1 == 0) {
        goto label63;
    }
    else {
        goto label63;
    }
label62:
    IT = NF - 2;                                                          
    IS = IT - 1;                                                           
    IS2 = 2 * IS;                                                        
    COEF = double(ME2 + IS + 1);                                             
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2,ME,ME2) * XF[IT - 1] + B(IS, IS2,ME,ME2) * XF[IT + 1]);
    goto label61;
label63:
    IT = NF - 2;                                                           
    IS = IT - 1;                                                     
    IS2 = 2 * IS;                                                          
    COEF = double(ME2 + IS + 1);  
    //cout << A(IS, IS2, ME, ME2)<<" B "<<B(IS, IS2, ME, ME2) << endl;
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2,ME,ME2) * XF[IT - 1] + B(IS, IS2,ME,ME2) * XF[IT + 1]);     
    ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
label67:
    IT = NF - 1;
    IS = IT - 1;
    IS2 = 2 * IS;
    COEFF = double(ME2 + IS + 1);
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2, ME, ME2) * XF[IT - 1] + B(IS, IS2, ME, ME2) * XF[IT + 1]);
    ALPHA = ALPHA *double(ME2 + IS + 1) / double(IS + 1);
    
label64:
    if (NFP-NF < 0) {
        goto label61;
    }
    else if (NFP - NF == 0) {
        goto label65;
    }
    else {
        goto label66;
    }
label65:
    IT = NF;
    IS = IT - 1;                                                         
    IS2 = 2 * IS;                                                         
    COEF = double(ME2 + IS + 1);                                            
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2,ME,ME2) * XF[IT - 1]);
    
    goto label61;
label66:
    IT = NF;                                                             
    IS = IT - 1;
    IS2 = 2 * IS;
    COEF = double(ME2 + IS + 1);
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2,ME,ME2) * XF[IT - 1]);
    ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
    IT = NF + 1;
    IS = IT - 1;
    IS2 = 2 * IS;
    COEF = double(ME2 + IS + 1);
    SAM = SAM + COEF * XFP[IT] * (A(IS, IS2,ME,ME2) * XF[IT - 1]);
label61:
    /*
        Start Iegrating over
        *
                B*
        *
    */
    SIM = 0.0;                                                         
    ALPHA = double(IALPHA * (ME2 + 1) * (ME2 + 2));                            
    NMAX = min(NFP, NF - 2);
    if (NMAX < 1) { goto label12; }
    for (int IT = 1; IT <= NMAX; IT++) {
        IS = IT - 1;
        IME = ME2 + IS + IS + 1;
        COEFF = 2.0 / double(IME + 2);
        SIM = SIM + COEFF * XFP[IT] * ALPHA * (XF[IT + 2] / double(IME + 4) - XF[IT] / double(IME));
        ALPHA = ALPHA * double(ME2 + IS + 3) / double(IT);
    }
label12:
    if (NFP - NF +1 < 0) {
        goto label231;
    }
    else if (NFP - NF + 1 == 0) {
        goto label232;
    }
    else {
        goto label233;
    }
label232:
    IT = NF - 1;
    IS = IT - 1;
    IS2 = 2 * IS;
    COEFF = 2.0 / double(ME2 + IS2 + 3);
    SIM = SIM - COEFF * XFP[IT] * ALPHA * XF[IT] / double(ME2 + IS2 + 1);
    goto label231;
label233:
    IT = NF - 1;                                                        
    IS = IT - 1;                                                           
    IS2 = 2 * IS;
    COEFF = 2.0 / double(ME2 + IS2 + 3);
    SIM = SIM - COEFF * XFP[IT] * ALPHA * XF[IT] / double(ME2 + IS2 + 1);
    ALPHA = ALPHA * double(ME2 + IS + 3) / double(IS + 1);
    IT = NF;
    IS = IT - 1;
    IS2 = 2 * IS;
    COEFF = 2.0 / double(ME2 + IS2 + 3);
    SIM = SIM - COEFF * XFP[IT] * ALPHA * XF[IT] / double(ME2 + IS2 + 1);

label231:
    //start integratoin over B''*
    SJM = 0.0;
    ALPHA = double(ALPHA);
    NMAX = min(NFP, NF);
    for (IT = 1; IT <= NMAX; IT++) {
        IBMAX = (NFP - IT) / 2 + 1;
        SP = 0;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP = SP + XFP[IN];
        }
        SP = SP + XFP[IN];
        SJM = SJM - SP * 2.0 * ALPHA * XF[IT];
        ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
    }
    //integration over B'*
    SKM = 0.0;
    ALPHA = double(IALPHA);
    if (NF < 4) { XF[4] = 0.0; }
    if (NF < 3) { XF[3] = 0.0; }
    for (IT = 1; IT <= 2; IT++) {
        IS = IT - 1;
        IS2 = 2 * IS;
        COEFF = -2.0 * ALPHA;
        GUN = (double((ME2 + IS) * (ME + IS + 1) * IS) / double(ME2 + IS2 - 1)
            - double((ME + IS) * (ME2 + IS + 1) * (IS + 1)) / double(ME2 + IS2 + 3)) / double(ME2 + IS2 + 1);
        GDEU = double((ME2 + IS + 2) * (ME + IS + 3) * (ME2 + IS + 1)) / double((ME2 + IS2 + 3) * (ME2 + IS2 + 5));
        SP = 0.0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP += XFP[IN];
        }
        SP = 0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP = SP + XFP[IN];
        }
        SKM += COEFF * SP * (XF[IT] * GUN + XF[IT + 2] * GDEU);
        ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
    }
    NMAX = min(NFP, NF - 2);
    if (NMAX < 3) { goto label125; }
    for (int IT = 3; IT <= NMAX; IT++) {
        IS = IT - 1;
        COEFF = -2.0 * ALPHA;
        MIM = IS + ME;
        IM = MIM + ME + 1;
        IME = IM + IS;
        GUN = (double((IM - 1) * (MIM + 1) * IS) / double(IME - 2)
            - double(MIM * IM * IT) / double(IME + 2)) / double(IME);
        GDEU = double((IM + 1) * (MIM + 3) * IM) / double((IME + 2) * (IME + 4));
        GTROI = -double((MIM - 2) * (IS - 1) * IS) / double((IME - 2) * (IME - 4));
        SP = 0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP = SP + XFP[IN];
        }
        SKM = SKM + COEFF * SP * (XF[IT] * GUN + XF[IT + 2] * GDEU + XF[IT - 2] * GTROI);
        ALPHA = ALPHA * double(IM) / double(IT);
    }
label125:
    if (NMAX == NFP) goto label130;
    NMAX += 1;
    if (NMAX < 3) { NMAX = 3; }
    NPMAX = min(NFP, NF);
    if (NPMAX < NMAX) { goto label132; }
    for (int IT = NMAX; IT <= NPMAX; IT++) {
        IS = IT - 1;
        IS2 = 2 * IS;
        COEFF = -2.0 * ALPHA;
        GUN = (double((ME2 + IS) * (ME + IS + 1) * IS) / double(ME2 + IS2 - 1)
            - double((ME + IS) * (ME2 + IS + 1) * (IS + 1)) / double(ME2 + IS2 + 3)) / double(ME2
                + IS2 + 1);
        GTROI = -double((ME + IS - 2) * (IS - 1) * IS) / double((ME2 + IS2 - 1) * (ME2 + IS2 - 3));
        SP = 0.0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP += XFP[IN];
        }
        SKM += COEFF * SP * (XF[IT] * GUN + XF[IT - 2] * GTROI);
        ALPHA *= double(ME2 + IS + 1) / double(IS + 1);
    }
label132:
    if (NPMAX == NFP) { goto label130; }
    NMAX = min(NFP, NF + 2);
    NPMAX = NPMAX + 1;
    if (NMAX < NPMAX) { goto label130; }
    for (int IT = NPMAX; IT <= NMAX; IT++) {
        IS = IT - 1;
        IS2 = 2 * IS;
        COEFF = -2 * ALPHA;
        GTROI = -double((ME + IS - 2) * (IS - 1) * IS) / double((ME2 + IS2 - 1) * (ME2 + IS2 - 3));
        SP = 0;
        IBMAX = (NFP - IT) / 2 + 1;
        for (int IB = 1; IB <= IBMAX; IB++) {
            IN = 2 * IB - 2 + IT;
            SP = SP + XFP[IN];
        }
        SKM = SKM + COEFF * SP * XF[IT - 2] * GTROI;
        ALPHA = ALPHA * double(ME2 + IS + 1) / double(IS + 1);
    }
label130:
    //start integration LAMBDA
    double SIGPE = PE + PEP;
    SG = 0.0;                                                           
    SPG = 0.0;
    SPI = 0.0;
    SPJ = 0.0;
    SPK = 0.0;
    SPN = 0.0;
    SPK2 = 0.0;
    for (int N = 1; N <= Gaussian.NXI; N++) {
        ALAM = Gaussian.XI[N] / SIGPE + 1;
        ALAM1 = ALAM - 1.0;
        ALAM2 = ALAM + 1.0;
        ALAMXI = ALAM1 / ALAM2;
        FONC = exp(-SIGPE * ALAM + (SIGMA + SIGMAP) * log(ALAM2));
        TRUC = ALAM1 * ALAM2;
        TRIC = pow(TRUC,ME);
        TRUC = TRIC * TRUC;
        FISC = FONC * TRIC * Gaussian.OMEGA[N];
        FONC = FONC * TRUC * Gaussian.OMEGA[N];
        YP = XGP[NGP];
        for (int JN = 2; JN <= NGP; JN++) {
            NN = NGP - JN + 1;
            YP = YP * ALAMXI + XGP[NN];
        }
        Y = XG[NG];
        DY = 0.0;
        for (int JN = 2; JN <= NG; JN++) {
            NN = NG - JN + 1;
            Y = Y * ALAMXI + XG[NN];
            DY = DY * ALAMXI + XG[NN + 1] * double(NN);
        }
        DY = 2 * DY;
        DY = DY + Y * (ME * ALAM / ALAMXI + ALAM2 * (SIGMA - PE * ALAM2));
        DY = DY / (ALAM2 * ALAM2);
        FISC = FISC * Y * YP;
        SG = SG + FONC * DY * YP;
        SPG = SPG + FONC * Y * YP * ALAM;
        SPI = SPI + FISC * ALAM;
        SPJ = SPJ + FONC * DY * YP * ALAM;
        SPK = SPK + FONC * Y * YP;
        SPK2 = SPK2 + FONC * Y * YP * ALAM * ALAM;
        SPN = SPN + FISC;
    }
    SG = SG / SIGPE;
    SPG = SPG / SIGPE;
    SPI = SPI / SIGPE;
    SPJ = SPJ / SIGPE;
    SPK = SPK / SIGPE;
    SPN = SPN / SIGPE;
    SPK2 = SPK2 / SIGPE;

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!   NEED SPK and SPK2 values for ETF correction calculation     !!!!*/
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    CROT = +R * R * R * (SG * SAM - SPG * SUM + ME * (SPI * SAM + SUMP * SPG)) / 16;
    
    GAMROT = +R * R * R * (SPJ * SIM - SPK * SKM + ME * (SPK * SJM + SIM * SPN)) / 16.0;
    
    CROT = CROT + OMEZUT * GAMROT;
    //cout << CROT << " " << GAMROT << endl;
    if (ME == 0) {
        CROT = CROT * sqrt(2.0);
        GAMROT = GAMROT * sqrt(2.0);
    }
    return;
}

void MEDOC(double QU, double AO, int N1, int L1, int M1, int N2, int L2, int M2, int NR, vector<double> RR, tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> GRAVERes1, tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> GRAVERes2, vector<double> &couplings,string&name) {
    
    double OMEZUT = 1.0 - 2 * AO;
    int JPERF;
    int NN, L, ME, NNP, LP, MEP;//state quantum number
    double R, RP;
    //state specific results from function GRAVE.
    int NG, NF, NGP, NFP;
    double E, EP, PE, PEP;
    vector<double> XG, XF, XGP, XFP, param1, param2;
    double SIGMA, SIGMAP;
    double FNOR, FNORP;
    double CRAD, GAMRRAD;
    double CROT, GAMROT;
    double DELTA{ 0 };

    //gaussian weights
    Gauss Gaussian;
    couplings.push_back(0.0);//place holder
    //loop through internuclear distances
    for (int JN = 1; JN < NR; JN++) {
        R = RR[JN];
        RP = R;
        //reads in outputs from Grave
        if ( abs(M1-M2)-1 <0 or (abs(M1-M2) == 1 and M1 <M2)) {
            //bra
            NN = N1;
            L = L1;
            ME = M1;
            param1 = get<0>(GRAVERes1)[JN - 1];
            E = param1[0];
            PE = param1[1];
            NG = param1[3];
            NF = param1[4];
            vector<double> tempXG = get<1>(GRAVERes1)[JN - 1];
            vector<double> tempXF = get<2>(GRAVERes1)[JN - 1];
            XG = vector<double>(tempXG.begin(), tempXG.begin() + NG+1);
            XF = vector<double>(tempXF.begin(), tempXF.begin() + NF+1);
            //ket
            NNP = N2;
            LP = L2;
            MEP = M2;
            param2 = get<0>(GRAVERes2)[JN - 1];
            EP = param2[0];
            PEP = param2[1];
            NGP = param2[3];
            NFP = param2[4];
            vector<double> tempXGP = get<1>(GRAVERes2)[JN - 1];
            vector<double> tempXFP = get<2>(GRAVERes2)[JN - 1];
            XGP = vector<double>(tempXGP.begin(), tempXGP.begin() + NGP+1);
            XFP = vector<double>(tempXFP.begin(), tempXFP.begin() + NFP+1);
        }
        else if (abs(M1-M2)-1 ==0 and M1>M2) {
            //bra is the second input set
            NN = N2;
            L = L2;
            ME = M2;
            param2 = get<0>(GRAVERes2)[JN - 1];
            E = param2[0];
            PE = param2[1];
            NG = param2[3];
            NF = param2[4];
            vector<double> tempXG = get<1>(GRAVERes2)[JN - 1];
            vector<double> tempXF = get<2>(GRAVERes2)[JN - 1];
            XG = vector<double>(tempXG.begin(), tempXG.begin() + NG+1);
            XF = vector<double>(tempXF.begin(), tempXF.begin() + NF+1);
            //ket is the first input set
            NNP = N1;
            LP = L1;
            MEP = M1;
            param1 = get<0>(GRAVERes1)[JN - 1];
            EP = param1[0];
            PEP = param1[1];
            NGP = param1[3];
            NFP = param1[4];
            vector<double> tempXGP = get<1>(GRAVERes1)[JN - 1];
            vector<double> tempXFP = get<2>(GRAVERes1)[JN - 1];
            XGP = vector<double>(tempXGP.begin(), tempXGP.begin() + NGP+1);
            XFP = vector<double>(tempXFP.begin(), tempXFP.begin() + NFP+1);
        }
        
        else {
            cout << "Error, check if state quantum numbers are valid " << endl;
        }
        
        //data reads in complete, keep looping over internuclear distance.
        SIGMA = R * (1.0 + QU) / (2.0 * PE) - 1.0 - double(ME);
        SIGMAP = RP * (1.0 + QU) / (2.0 * PEP) - 1.0 - double(MEP);
        if (MEP != ME) {
            //rotational coupling calculation
            NORDIF(ME, PE, SIGMA, NG, XG, NF, XF, FNOR, R);
            NORDIF(MEP, PEP, SIGMAP, NGP, XGP, NFP, XFP, FNORP, R);
            
            ROTDIS(CROT, GAMROT, Gaussian, XFP,XGP,NFP,NGP, XF,XG,NF,NG, R,PE,PEP,SIGMA,SIGMAP,DELTA,L,ME,OMEZUT);
            //cout << FNOR << " " << FNORP << " " << CROT << " " << GAMROT << endl;
            CROT = CROT * FNOR * FNORP;
            couplings.push_back( CROT);
        }
        else {
            //radial coupling calculation
        }

    
    }
    name = "coupling"+to_string(NN) + to_string(L) + to_string(ME) + "-" + to_string(NNP) + to_string(LP) + to_string(MEP);
    return;
}

int main() {
    vector<double> RR(999);
    int NR;
    double RINCR;
    float QU ,AO;
    int N1, L1,M1;
    int N2, L2, M2;

    std::ifstream file("input.txt");
    if (!file) {
        std::cerr << "Unable to open file" << std::endl;
        return 1;
    }
    file >> NR >> RINCR;
    file >> QU >> AO;
    file >> N1 >> L1 >> M1;
    file >> N2 >> L2 >> M2;
    file.close();
    //state1 declaration
    
    vector <double> PP1(999), CSEP1(999);
    vector <double> GraveP1(999), GraveC1(999);
    wave_function(QU,RINCR,N1, L1, M1, NR,RR, PP1, CSEP1,GraveP1,GraveC1);
    //post-processing
    RR.erase(RR.begin());
    PP1.erase(PP1.begin());
    CSEP1.erase(CSEP1.begin());
    GraveP1.erase(GraveP1.begin());
    GraveC1.erase(GraveC1.begin());
    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> GRAVERes1;
    GRAVERes1 = GRAVE(QU, N1, L1, M1, NR, RR, GraveP1, GraveC1);

    

    //state2
    vector<double> RR2(999);
    vector <double> PP2(999), CSEP2(999);
    vector <double> GraveP2(999), GraveC2(999);
    wave_function(QU,RINCR, N2, L2, M2, NR,RR2, PP2, CSEP2, GraveP2, GraveC2);
    tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> GRAVERes2;
    RR2.erase(RR2.begin());
    PP2.erase(PP2.begin());
    CSEP2.erase(CSEP2.begin());
    GraveP2.erase(GraveP2.begin());
    GraveC2.erase(GraveC2.begin());
    GRAVERes2 = GRAVE(QU, N2, L2, M2, NR, RR2, GraveP2, GraveC2);

    vector<double> couplings{};
    string name;
    MEDOC(QU,AO, N1, L1, M1, N2, L2, M2, NR, RR, GRAVERes1, GRAVERes2, couplings,name);
    
    std::ofstream outFileCoup(name);
    if (outFileCoup.is_open()) {
        for (int i = 1; i < NR; ++i) {
        //cout << RR[i] << "\t" << couplings[i]<<endl;
        outFileCoup << RR[i] <<"\t" <<couplings[i]<<endl; 
        }

        outFileCoup.close(); // Close the file
    }
    else {
        std::cerr << "Unable to open file";
    }
    return 0;
}
