#include <math.h>

void initialize_states_default(double *STATES) {
    STATES[0] =  7.02128101897185673e-04;
    STATES[1] =  3.94923428392655786e-03;
    STATES[2] =  1.35538532457244482e-01;
    STATES[3] =  1.03674364292988680e-01;
    STATES[4] =  1.90759804527589089e-01;
    STATES[5] =  1.35640688636079511e-02;
    STATES[6] =  2.14063418881809235e-02;
    STATES[7] =  4.45327242854324807e-03;
    STATES[8] =  1.27856586024588575e-01;
    STATES[9] =  5.69999505293381902e-03;
    STATES[10] =  1.83143535034222225e-02;
    STATES[11] =  2.10808768153058460e-04;
    STATES[12] =  3.25814677291117296e-04;
    STATES[13] =  2.33018340557575125e-04;
    STATES[14] =  3.61396062660070427e+00;
    STATES[15] =  7.88607791910409195e-01;
    STATES[16] =  9.15153381546177336e+00;
    STATES[17] =  9.15182798281732346e+00;
    STATES[18] =  5.02305826642838293e-01;
    STATES[19] =  1.13337536953687845e+00;
    STATES[20] = -7.34336366728778671e+01;
    STATES[21] =  2.16850216379767157e-05;
    STATES[22] =  9.98384427312367095e-01;
    STATES[23] =  4.49572164109603364e-02;
    STATES[24] =  3.28512098597005947e-02;
    STATES[25] = 120.0;
    STATES[26] =  1.31290096227093382e-03;
    STATES[27] =  7.49436760722081534e-03;
    STATES[28] =  9.15199678386256998e+00;
    STATES[29] =  3.93548562883350357e-04;
    STATES[30] =  9.58234428284286399e-01;
    STATES[31] =  3.15482710277587786e-01;
    STATES[32] =  2.48034071360795916e-01;
    STATES[33] =  1.89326933812916480e-02;
    STATES[34] =  3.79829335413739144e-02;
    STATES[35] =  1.01974216400706526e-02;
    STATES[36] =  1.37939236359928058e-03;
    STATES[37] =  9.45874848392074696e-01;
    STATES[38] =  5.01323282772066123e-07;
    STATES[39] =  2.01567245823636694e-06;
    STATES[40] =  8.00819151705148946e-01;
}

void initialize_constants_default(double* CONSTANTS) {
CONSTANTS[0] = 0.024;
CONSTANTS[1] = 0.0171;
CONSTANTS[2] = 0.14;
CONSTANTS[3] = 0.07;
CONSTANTS[4] = 0.14;
CONSTANTS[5] = 0.238;
CONSTANTS[6] = 0.00046;
CONSTANTS[7] = 5.7e-05;
CONSTANTS[8] = 0.03;
CONSTANTS[9] = 1.3;
CONSTANTS[10] = 0.06;
CONSTANTS[11] = 3.2e-05;
CONSTANTS[12] = 0.00333;
CONSTANTS[13] = 34.0;
CONSTANTS[14] = 13.8;
CONSTANTS[15] = 0.0157;
CONSTANTS[16] = 100.0;
CONSTANTS[17] = 100.0;
CONSTANTS[18] = 100.0;
CONSTANTS[19] = 2.37;
CONSTANTS[20] = 0.003;
CONSTANTS[21] = 32.7;
CONSTANTS[22] = 1.0;
CONSTANTS[23] = 0.0;
CONSTANTS[24] = 7.561;
CONSTANTS[25] = 1.65;
CONSTANTS[26] = 0.001;
CONSTANTS[27] = 0.0001;
CONSTANTS[28] =  3.72425607984805052e-12;
CONSTANTS[29] = 1.1e-10;
CONSTANTS[30] =  8.24130542277896849e-13;
CONSTANTS[31] = 96485.0;
CONSTANTS[32] = 65.0;
CONSTANTS[33] = 100.0;
CONSTANTS[34] = 0.0;
CONSTANTS[35] = 0.0;
CONSTANTS[36] =  1.83127823220607955e-14;
CONSTANTS[37] =  1.63862792221979433e-12;
CONSTANTS[38] = 100.0;
CONSTANTS[39] = 10.25;
CONSTANTS[40] =  3.14159265358979312e+00;
CONSTANTS[41] =  6.06430000000000033e-04;
CONSTANTS[42] = 0.11;
CONSTANTS[43] = 1.8;
CONSTANTS[44] = 1.0;
CONSTANTS[45] = 0.0;
CONSTANTS[46] = 0.0;
CONSTANTS[47] = 0.00027;
CONSTANTS[48] = 1.35e-07;
CONSTANTS[49] = 7.5e-09;
CONSTANTS[50] = 0.9;
CONSTANTS[51] = 1.8;
CONSTANTS[52] = 5.4;
CONSTANTS[53] = 140.0;
CONSTANTS[54] = 0.009;
CONSTANTS[55] = 0.0548;
CONSTANTS[56] = 0.1;
CONSTANTS[57] = 0.0525;
CONSTANTS[58] = 0.002;
CONSTANTS[59] = 0.035;
CONSTANTS[60] = 0.0035;
CONSTANTS[61] = 0.01833;
CONSTANTS[62] = 0.045;
CONSTANTS[63] = 23.0;
CONSTANTS[64] = 0.000597;
CONSTANTS[65] = 3.15;
CONSTANTS[66] = 0.000384;
CONSTANTS[67] = 0.00359;
CONSTANTS[68] = 1.3;
CONSTANTS[69] = 12.29;
CONSTANTS[70] = 87.5;
CONSTANTS[71] = 1.57;
CONSTANTS[72] = 0.27;
CONSTANTS[73] = 0.35;
CONSTANTS[74] = 1.26;
CONSTANTS[75] = 1.5;
CONSTANTS[76] = 0.0025;
CONSTANTS[77] = 600.0;
CONSTANTS[78] = 15.0;
CONSTANTS[79] = 150.0;
CONSTANTS[80] = 0.0471;
CONSTANTS[81] = 0.0005;
CONSTANTS[82] = 2.35;
CONSTANTS[83] = 0.165;
CONSTANTS[84] = 8314.0;
CONSTANTS[85] = 310.0;
CONSTANTS[86] = 5.348e-06;
CONSTANTS[87] = 1.7;
CONSTANTS[88] = 15.0;
CONSTANTS[89] = 1.0;
CONSTANTS[90] = 2.6;
CONSTANTS[91] = 0.0053114;
CONSTANTS[92] = 0.45;
CONSTANTS[93] = 1.787;
CONSTANTS[94] = 0.5;
CONSTANTS[95] = 0.005;
CONSTANTS[96] = 0.06;
CONSTANTS[97] = 25.0;
CONSTANTS[98] = -12.5;
CONSTANTS[99] = 5.0;
CONSTANTS[100] = 50.0;
CONSTANTS[101] = 1000.0;
CONSTANTS[102] =  (1.00000+ 0.500000*CONSTANTS[23])*0.0196000;
CONSTANTS[103] =  ( ( CONSTANTS[40]*pow(CONSTANTS[39], 2.00000))*CONSTANTS[38])*1.00000e-15;
CONSTANTS[104] = pow( CONSTANTS[81]*1.00000, 1.60000);
CONSTANTS[105] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[47];
CONSTANTS[106] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[48];
CONSTANTS[107] =  ( (1.00000+ 0.500000*CONSTANTS[23])*(1.00000 -  0.500000*CONSTANTS[34]))*CONSTANTS[49];
CONSTANTS[108] = 1.00000 - CONSTANTS[50];
CONSTANTS[109] = (CONSTANTS[31]/CONSTANTS[84])/CONSTANTS[85];
CONSTANTS[110] =  ( (1.00000+CONSTANTS[34])* pow((CONSTANTS[52]/5.40000), 1.0 / 2))*CONSTANTS[57];
CONSTANTS[111] =  CONSTANTS[59]* pow((CONSTANTS[52]/5.40000), 1.0 / 2);
CONSTANTS[112] =  ((1.00000+CONSTANTS[34])+ 2.00000*CONSTANTS[23])*CONSTANTS[60];
CONSTANTS[113] =  ((1.00000+CONSTANTS[34])+ 2.00000*CONSTANTS[23])*CONSTANTS[60];
CONSTANTS[114] =  ( ( (1.00000 -  0.500000*CONSTANTS[34])*(1.00000+ 2.00000*CONSTANTS[23]))*(1.00000+ 0.200000*CONSTANTS[35]))*CONSTANTS[62];
CONSTANTS[115] =  CONSTANTS[63]*(1.00000 -  0.100000*CONSTANTS[34]);
CONSTANTS[116] =  (1.00000+ 0.400000*CONSTANTS[34])*CONSTANTS[65];
CONSTANTS[117] =  11.0000*(1.00000 -  0.250000*CONSTANTS[23]);
CONSTANTS[118] = (exp(CONSTANTS[53]/67.3000) - 1.00000)/7.00000;
CONSTANTS[119] =  CONSTANTS[76]*CONSTANTS[34];
CONSTANTS[120] =  (1.00000 -  0.700000*CONSTANTS[34])*CONSTANTS[83];
CONSTANTS[121] =  (2.50000 -  1.25000*CONSTANTS[23])*0.000246000;
CONSTANTS[122] =  ((10.0000+ 20.0000*CONSTANTS[34])+ ( 10.0000*CONSTANTS[23])*(1.00000 - CONSTANTS[34]))*1.00000;
CONSTANTS[123] =  ( 0.0539000*0.0100000)*CONSTANTS[103];
CONSTANTS[124] = pow( CONSTANTS[81]*1.00000, 1.60000);
CONSTANTS[125] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[78]/CONSTANTS[79]);
CONSTANTS[126] =  0.650000*CONSTANTS[103];
CONSTANTS[127] = 1.00000 - CONSTANTS[42];
CONSTANTS[128] =  (( 0.00165000*CONSTANTS[126])/CONSTANTS[123])*0.100000;
CONSTANTS[129] =  (( 0.00460000*CONSTANTS[126])/CONSTANTS[123])*0.100000;
CONSTANTS[130] =  0.0200000*CONSTANTS[103];
CONSTANTS[131] =  0.0350000*CONSTANTS[103];
CONSTANTS[132] = (CONSTANTS[85] - 310.000)/10.0000;
CONSTANTS[133] =  (CONSTANTS[126]/CONSTANTS[130])*0.0134000;
CONSTANTS[134] =  (CONSTANTS[126]/CONSTANTS[130])*0.0374000;
CONSTANTS[135] =  (CONSTANTS[126]/CONSTANTS[131])*0.140000;
}


void computeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC) {
RATES[25] = 0;
RATES[0] =  ( CONSTANTS[13]*STATES[11])*(CONSTANTS[0] - STATES[0]) -  CONSTANTS[5]*STATES[0];
RATES[23] =  ( ( 1.70000*STATES[12])*(1.00000 - STATES[23]) -  0.0119000*STATES[23])*1.00000;
RATES[24] =  ( ( 1.70000*STATES[13])*(1.00000 - STATES[24]) -  0.0119000*STATES[24])*1.00000;
RATES[28] =  (CONSTANTS[37]/CONSTANTS[126])*(STATES[17] - STATES[28]);
RATES[1] =  ( CONSTANTS[14]*STATES[11])*((CONSTANTS[4] - STATES[1]) - STATES[2]) -  CONSTANTS[6]*STATES[1];
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[20]+91.0000)/6.10000));
RATES[34] = (ALGEBRAIC[9] - STATES[34])/CONSTANTS[77];
RATES[2] =  ( CONSTANTS[15]*CONSTANTS[22])*((CONSTANTS[4] - STATES[1]) - STATES[2]) -  CONSTANTS[7]*STATES[2];
ALGEBRAIC[0] = 1.00000/(1.00000+exp(- ((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000)/6.00000));
ALGEBRAIC[14] = ( ALGEBRAIC[0]*(1.00000 - exp(- ((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000)/6.00000)))/( 0.0350000*((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000));
RATES[21] = (ALGEBRAIC[0] - STATES[21])/ALGEBRAIC[14];
ALGEBRAIC[1] = 1.00000/(1.00000+exp(((STATES[20]+ 3.00000*CONSTANTS[23])+30.0000)/7.00000))+0.200000/(1.00000+exp(((50.0000 - STATES[20]) -  3.00000*CONSTANTS[23])/20.0000));
ALGEBRAIC[15] = 1.00000/( 0.0197000*exp(- pow( 0.0337000*((STATES[20]+ 3.00000*CONSTANTS[23])+25.0000), 2.00000))+0.0200000);
RATES[22] = (ALGEBRAIC[1] - STATES[22])/ALGEBRAIC[15];
ALGEBRAIC[2] = 1.00000/(1.00000+exp(- (STATES[20]+10.0000)/5.00000));
ALGEBRAIC[16] = ( ((550.000)/(1.00000+exp((- 22.0000 - STATES[20])/9.00000)))*6.00000)/(1.00000+exp((STATES[20]+11.0000)/9.00000))+230.000/(1.00000+exp((STATES[20]+40.0000)/20.0000));
RATES[26] = (ALGEBRAIC[2] - STATES[26])/ALGEBRAIC[16];
ALGEBRAIC[3] = 1.00000/(1.00000+exp(- ((STATES[20]+ 40.0000*CONSTANTS[23])+3.80000)/14.2500));
ALGEBRAIC[17] = 990.100/(1.00000+exp(- ((STATES[20]+ 40.0000*CONSTANTS[23])+2.43600)/14.1200));
RATES[27] = (ALGEBRAIC[3] - STATES[27])/ALGEBRAIC[17];
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[20]+6.00000)/- 8.60000));
ALGEBRAIC[18] = 9.00000/(1.00000+exp((STATES[20]+5.00000)/12.0000))+0.500000;
RATES[29] = (ALGEBRAIC[4] - STATES[29])/ALGEBRAIC[18];
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[20]+7.50000)/10.0000));
ALGEBRAIC[19] = 590.000/(1.00000+exp((STATES[20]+60.0000)/10.0000))+3050.00;
RATES[30] = (ALGEBRAIC[5] - STATES[30])/ALGEBRAIC[19];
ALGEBRAIC[8] = 1.00000/pow(1.00000+exp(- (56.8600+STATES[20])/9.03000), 2.00000);
ALGEBRAIC[22] =  0.129200*exp(- pow((STATES[20]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[20] - 4.82300)/51.1200, 2.00000));
RATES[33] = (ALGEBRAIC[8] - STATES[33])/ALGEBRAIC[22];
ALGEBRAIC[10] = ( 0.320000*(STATES[20]+47.1300))/(1.00000 - exp( - 0.100000*(STATES[20]+47.1300)));
ALGEBRAIC[23] =  0.0800000*exp(- STATES[20]/11.0000);
RATES[35] =  ALGEBRAIC[10]*(1.00000 - STATES[35]) -  ALGEBRAIC[23]*STATES[35];
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[20]+1.00000)/11.0000));
ALGEBRAIC[24] =  3.50000*exp(- pow(STATES[20]/30.0000, 2.00000))+1.50000;
RATES[36] = (ALGEBRAIC[11] - STATES[36])/ALGEBRAIC[24];
ALGEBRAIC[12] = 1.00000/(1.00000+exp((STATES[20]+40.5000)/11.5000));
ALGEBRAIC[25] =  25.6350*exp(- pow((STATES[20]+52.4500)/15.8827, 2.00000))+24.1400;
RATES[37] = (ALGEBRAIC[12] - STATES[37])/ALGEBRAIC[25];
RATES[3] =  ( CONSTANTS[16]*STATES[12])*(CONSTANTS[128] - STATES[3]) -  CONSTANTS[8]*STATES[3];
RATES[4] =  ( CONSTANTS[16]*STATES[13])*(CONSTANTS[133] - STATES[4]) -  CONSTANTS[8]*STATES[4];
ALGEBRAIC[27] = 1.00000/pow(1.00000+exp((STATES[20]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[6] =  (STATES[20]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[20]+80.0000)/6.80000))*1.00000;
ALGEBRAIC[20] =  (STATES[20]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[20]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[20])+ ( 3.10000*pow(10.0000, 5.00000))*exp( 0.348500*STATES[20]))*1.00000;
ALGEBRAIC[30] = 1.00000/(ALGEBRAIC[6]+ALGEBRAIC[20]);
RATES[31] = (ALGEBRAIC[27] - STATES[31])/ALGEBRAIC[30];
ALGEBRAIC[28] = 1.00000/pow(1.00000+exp((STATES[20]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[7] =  (STATES[20]>=- 40.0000 ? 0.00000 :  (( ( ( - 2.54280*pow(10.0000, 4.00000))*exp( 0.244400*STATES[20]) -  ( 6.94800*pow(10.0000, - 6.00000))*exp( - 0.0439100*STATES[20]))*(STATES[20]+37.7800))/(1.00000+exp( 0.311000*(STATES[20]+79.2300))))*1.00000)*1.00000;
ALGEBRAIC[21] =  (STATES[20]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[20]))/(1.00000+exp( - 0.100000*(STATES[20]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[20]))/(1.00000+exp( - 0.137800*(STATES[20]+40.1400))))*1.00000;
ALGEBRAIC[31] = 1.00000/(ALGEBRAIC[7]+ALGEBRAIC[21]);
RATES[32] = (ALGEBRAIC[28] - STATES[32])/ALGEBRAIC[31];
ALGEBRAIC[13] = ((1.00000 - STATES[40]) - STATES[39]) - STATES[38];
ALGEBRAIC[26] = CONSTANTS[88] - (CONSTANTS[88] - CONSTANTS[89])/(1.00000+pow(CONSTANTS[92]/STATES[18], 2.50000));
ALGEBRAIC[29] =  CONSTANTS[94]*ALGEBRAIC[26];
ALGEBRAIC[32] = CONSTANTS[122]/ALGEBRAIC[26];
RATES[38] = ( ( ALGEBRAIC[29]*STATES[12])*STATES[39] -  CONSTANTS[95]*STATES[38]) - ( CONSTANTS[96]*STATES[38] -  ( ALGEBRAIC[32]*pow(STATES[12], 2.00000))*ALGEBRAIC[13]);
RATES[39] = ( ( ALGEBRAIC[32]*pow(STATES[12], 2.00000))*STATES[40] -  CONSTANTS[96]*STATES[39]) - ( ( ALGEBRAIC[29]*STATES[12])*STATES[39] -  CONSTANTS[95]*STATES[38]);
RATES[40] = ( CONSTANTS[95]*ALGEBRAIC[13] -  ( ALGEBRAIC[29]*STATES[12])*STATES[40]) - ( ( ALGEBRAIC[32]*pow(STATES[12], 2.00000))*STATES[40] -  CONSTANTS[96]*STATES[39]);
RATES[5] =  ( CONSTANTS[17]*STATES[12])*(CONSTANTS[129] - STATES[5]) -  CONSTANTS[9]*STATES[5];
RATES[6] =  ( CONSTANTS[17]*STATES[13])*(CONSTANTS[134] - STATES[6]) -  CONSTANTS[9]*STATES[6];
RATES[7] =  ( CONSTANTS[18]*STATES[11])*(CONSTANTS[1] - STATES[7]) -  CONSTANTS[10]*STATES[7];
RATES[8] =  ( CONSTANTS[19]*STATES[11])*((CONSTANTS[2] - STATES[8]) - STATES[9]) -  CONSTANTS[11]*STATES[8];
RATES[9] =  ( CONSTANTS[20]*CONSTANTS[22])*((CONSTANTS[2] - STATES[8]) - STATES[9]) -  CONSTANTS[12]*STATES[9];
RATES[10] =  ( CONSTANTS[21]*STATES[11])*(CONSTANTS[3] - STATES[10]) -  CONSTANTS[102]*STATES[10];
RATES[14] =  ( CONSTANTS[27]*STATES[16])*(CONSTANTS[24] - STATES[14]) -  CONSTANTS[26]*STATES[14];
RATES[15] =  ( CONSTANTS[27]*STATES[17])*(CONSTANTS[25] - STATES[15]) -  CONSTANTS[26]*STATES[15];
RATES[19] =  ( CONSTANTS[33]*STATES[18])*(CONSTANTS[135] - STATES[19]) -  CONSTANTS[32]*STATES[19];
ALGEBRAIC[33] = RATES[5]+RATES[3];
ALGEBRAIC[80] =  ((1.00000/CONSTANTS[109])/2.00000)*log(CONSTANTS[51]/STATES[12]);
ALGEBRAIC[81] =  ( CONSTANTS[42]*CONSTANTS[41])*(STATES[20] - ALGEBRAIC[80]);
ALGEBRAIC[36] =  (( ( ( CONSTANTS[105]*4.00000)*( ( STATES[20]*CONSTANTS[31])*CONSTANTS[109]))*( ( 0.341000*STATES[12])*exp( ( 2.00000*STATES[20])*CONSTANTS[109]) -  0.341000*CONSTANTS[51]))/(exp( ( 2.00000*STATES[20])*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[37] =  ( ( ( ( ( CONSTANTS[50]*ALGEBRAIC[36])*STATES[21])*STATES[22])*((1.00000 - STATES[23])+CONSTANTS[46]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[59] = 1.00000/(1.00000+pow(CONSTANTS[66]/STATES[12], 2.00000));
ALGEBRAIC[61] =  ( exp( ( CONSTANTS[73]*STATES[20])*CONSTANTS[109])*pow(STATES[16], 3.00000))*CONSTANTS[51];
ALGEBRAIC[62] =  ( exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109])*pow(CONSTANTS[53], 3.00000))*STATES[12];
ALGEBRAIC[63] = ((( ( CONSTANTS[67]*pow(CONSTANTS[53], 3.00000))*(1.00000+pow(STATES[16]/CONSTANTS[69], 3.00000))+ ( pow(CONSTANTS[70], 3.00000)*STATES[12])*(1.00000+STATES[12]/CONSTANTS[67]))+ CONSTANTS[68]*pow(STATES[16], 3.00000))+ pow(STATES[16], 3.00000)*CONSTANTS[51])+ pow(CONSTANTS[53], 3.00000)*STATES[12];
ALGEBRAIC[64] = (( ( ( ( CONSTANTS[42]*CONSTANTS[116])*pow(CONSTANTS[71], CONSTANTS[132]))*ALGEBRAIC[59])*(ALGEBRAIC[61] - ALGEBRAIC[62]))/ALGEBRAIC[63])/(1.00000+ CONSTANTS[72]*exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109]));
ALGEBRAIC[74] = pow( STATES[12]*1.00000, 1.60000);
ALGEBRAIC[75] = ( ( ( CONSTANTS[42]*pow(CONSTANTS[82], CONSTANTS[132]))*CONSTANTS[80])*ALGEBRAIC[74])/(ALGEBRAIC[74]+CONSTANTS[104]);
ALGEBRAIC[82] = ((ALGEBRAIC[37]+ALGEBRAIC[81])+ALGEBRAIC[75]) -  2.00000*ALGEBRAIC[64];
ALGEBRAIC[83] =  ( CONSTANTS[97]*STATES[39])*(STATES[18] - STATES[12]);
ALGEBRAIC[85] =  ( (1.00000+ 0.250000*CONSTANTS[34])*(STATES[18] - STATES[12]))*CONSTANTS[86];
RATES[12] = (((( - ALGEBRAIC[82]*CONSTANTS[29])/( ( CONSTANTS[123]*2.00000)*CONSTANTS[31])+ (CONSTANTS[30]/CONSTANTS[123])*(STATES[13] - STATES[12])) - ALGEBRAIC[33])+( ALGEBRAIC[83]*CONSTANTS[131])/CONSTANTS[123])+( ALGEBRAIC[85]*CONSTANTS[126])/CONSTANTS[123];
ALGEBRAIC[35] = (((((RATES[10]+RATES[8])+RATES[9])+RATES[0])+RATES[1])+RATES[2])+RATES[7];
ALGEBRAIC[87] = ( ( pow(CONSTANTS[90], CONSTANTS[132])*CONSTANTS[91])*(pow(STATES[11]/CONSTANTS[121], CONSTANTS[93]) - pow(STATES[18]/CONSTANTS[87], CONSTANTS[93])))/((1.00000+pow(STATES[11]/CONSTANTS[121], CONSTANTS[93]))+pow(STATES[18]/CONSTANTS[87], CONSTANTS[93]));
RATES[11] = (( - ALGEBRAIC[87]*CONSTANTS[131])/CONSTANTS[126] - ALGEBRAIC[35])+ (CONSTANTS[28]/CONSTANTS[126])*(STATES[13] - STATES[11]);
ALGEBRAIC[34] = RATES[6]+RATES[4];
ALGEBRAIC[84] =  ((1.00000/CONSTANTS[109])/2.00000)*log(CONSTANTS[51]/STATES[13]);
ALGEBRAIC[86] =  ( CONSTANTS[127]*CONSTANTS[41])*(STATES[20] - ALGEBRAIC[84]);
ALGEBRAIC[38] =  (( ( ( CONSTANTS[105]*4.00000)*( ( STATES[20]*CONSTANTS[31])*CONSTANTS[109]))*( ( 0.341000*STATES[13])*exp( ( 2.00000*STATES[20])*CONSTANTS[109]) -  0.341000*CONSTANTS[51]))/(exp( ( 2.00000*STATES[20])*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[39] =  ( ( ( ( ( CONSTANTS[108]*ALGEBRAIC[38])*STATES[21])*STATES[22])*((1.00000 - STATES[24])+CONSTANTS[45]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[60] = 1.00000/(1.00000+pow(CONSTANTS[66]/STATES[13], 2.00000));
ALGEBRAIC[65] =  ( exp( ( CONSTANTS[73]*STATES[20])*CONSTANTS[109])*pow(STATES[17], 3.00000))*CONSTANTS[51];
ALGEBRAIC[66] =  ( exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109])*pow(CONSTANTS[53], 3.00000))*STATES[13];
ALGEBRAIC[67] = ((( ( CONSTANTS[67]*pow(CONSTANTS[53], 3.00000))*(1.00000+pow(STATES[17]/CONSTANTS[69], 3.00000))+ ( pow(CONSTANTS[70], 3.00000)*STATES[13])*(1.00000+STATES[13]/CONSTANTS[67]))+ CONSTANTS[68]*pow(STATES[17], 3.00000))+ pow(STATES[17], 3.00000)*CONSTANTS[51])+ pow(CONSTANTS[53], 3.00000)*STATES[13];
ALGEBRAIC[68] = (( ( ( ( CONSTANTS[127]*CONSTANTS[116])*pow(CONSTANTS[71], CONSTANTS[132]))*ALGEBRAIC[60])*(ALGEBRAIC[65] - ALGEBRAIC[66]))/ALGEBRAIC[67])/(1.00000+ CONSTANTS[72]*exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109]));
ALGEBRAIC[76] = pow( STATES[13]*1.00000, 1.60000);
ALGEBRAIC[77] = ( ( ( CONSTANTS[127]*pow(CONSTANTS[82], CONSTANTS[132]))*CONSTANTS[80])*ALGEBRAIC[76])/(ALGEBRAIC[76]+CONSTANTS[124]);
ALGEBRAIC[89] = ((ALGEBRAIC[39]+ALGEBRAIC[86])+ALGEBRAIC[77]) -  2.00000*ALGEBRAIC[68];
RATES[13] = ((( - ALGEBRAIC[89]*CONSTANTS[29])/( ( CONSTANTS[130]*2.00000)*CONSTANTS[31])+ (CONSTANTS[30]/CONSTANTS[130])*(STATES[12] - STATES[13]))+ (CONSTANTS[28]/CONSTANTS[130])*(STATES[11] - STATES[13])) - ALGEBRAIC[34];
RATES[18] = (ALGEBRAIC[87] - (( ALGEBRAIC[85]*CONSTANTS[126])/CONSTANTS[131]+ALGEBRAIC[83])) - RATES[19];
ALGEBRAIC[43] =  (( ( ( ( CONSTANTS[107]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[16])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[53]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[44] =  ( ( ( ( ( CONSTANTS[50]*ALGEBRAIC[43])*STATES[21])*STATES[22])*((1.00000 - STATES[23])+CONSTANTS[46]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[102] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[53]/STATES[16]);
ALGEBRAIC[103] =  ( ( ( ( CONSTANTS[42]*CONSTANTS[115])*pow(STATES[33], 3.00000))*STATES[31])*STATES[32])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[104] =  ( CONSTANTS[42]*CONSTANTS[64])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[69] = 1.00000/((1.00000+ 0.124500*exp( ( - 0.100000*STATES[20])*CONSTANTS[109]))+ ( 0.0365000*CONSTANTS[118])*exp( - STATES[20]*CONSTANTS[109]));
ALGEBRAIC[71] = (( ( ( CONSTANTS[42]*CONSTANTS[74])*ALGEBRAIC[69])*CONSTANTS[52])/(1.00000+pow(CONSTANTS[117]/STATES[16], 4.00000)))/(CONSTANTS[52]+CONSTANTS[75]);
ALGEBRAIC[105] =  ( ( ( CONSTANTS[42]*CONSTANTS[119])*pow(STATES[35], 3.00000))*STATES[34])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[114] = ((((ALGEBRAIC[103]+ALGEBRAIC[104])+ 3.00000*ALGEBRAIC[64])+ 3.00000*ALGEBRAIC[71])+ALGEBRAIC[44])+ALGEBRAIC[105];
RATES[16] = (( - ALGEBRAIC[114]*CONSTANTS[29])/( CONSTANTS[123]*CONSTANTS[31])+ (CONSTANTS[36]/CONSTANTS[123])*(STATES[17] - STATES[16])) - RATES[14];
ALGEBRAIC[45] =  (( ( ( ( CONSTANTS[107]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[17])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[53]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[46] =  ( ( ( ( ( CONSTANTS[108]*ALGEBRAIC[45])*STATES[21])*STATES[22])*((1.00000 - STATES[24])+CONSTANTS[45]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[106] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[53]/STATES[17]);
ALGEBRAIC[107] =  ( ( ( ( CONSTANTS[127]*CONSTANTS[115])*pow(STATES[33], 3.00000))*STATES[31])*STATES[32])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[109] =  ( CONSTANTS[127]*CONSTANTS[64])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[72] = (( ( ( CONSTANTS[127]*CONSTANTS[74])*ALGEBRAIC[69])*CONSTANTS[52])/(1.00000+pow(CONSTANTS[117]/STATES[17], 4.00000)))/(CONSTANTS[52]+CONSTANTS[75]);
ALGEBRAIC[111] =  ( ( ( CONSTANTS[127]*CONSTANTS[119])*pow(STATES[35], 3.00000))*STATES[34])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[115] = ((((ALGEBRAIC[107]+ALGEBRAIC[109])+ 3.00000*ALGEBRAIC[68])+ 3.00000*ALGEBRAIC[72])+ALGEBRAIC[46])+ALGEBRAIC[111];
RATES[17] = ((( - ALGEBRAIC[115]*CONSTANTS[29])/( CONSTANTS[130]*CONSTANTS[31])+ (CONSTANTS[36]/CONSTANTS[130])*(STATES[16] - STATES[17]))+ (CONSTANTS[37]/CONSTANTS[130])*(STATES[28] - STATES[17])) - RATES[15];
ALGEBRAIC[90] = ALGEBRAIC[82]+ALGEBRAIC[89];
ALGEBRAIC[47] =  CONSTANTS[54]*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[49] =  (( CONSTANTS[42]*CONSTANTS[55])/(1.00000+CONSTANTS[56]/STATES[12]))*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[51] =  (( CONSTANTS[127]*CONSTANTS[55])/(1.00000+CONSTANTS[56]/STATES[13]))*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[52] = ALGEBRAIC[49]+ALGEBRAIC[51];
ALGEBRAIC[78] = ALGEBRAIC[52]+ALGEBRAIC[47];
ALGEBRAIC[116] = ALGEBRAIC[114]+ALGEBRAIC[115];
ALGEBRAIC[41] =  (( ( ( ( CONSTANTS[106]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[25])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[52]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[42] =  ( ( ( ( ALGEBRAIC[41]*STATES[21])*STATES[22])*( CONSTANTS[50]*(CONSTANTS[46]+(1.00000 - STATES[23]))+ CONSTANTS[108]*(CONSTANTS[45]+(1.00000 - STATES[24]))))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[91] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[52]/STATES[25]);
ALGEBRAIC[92] =  1.00000*(1.02000/(1.00000+exp( 0.238500*((STATES[20] - ALGEBRAIC[91]) - 59.2150))));
ALGEBRAIC[93] =  1.00000*(( 0.491240*exp( 0.0803200*((STATES[20]+5.47600) - ALGEBRAIC[91]))+exp( 0.0617500*((STATES[20] - ALGEBRAIC[91]) - 594.310)))/(1.00000+exp( - 0.514300*((STATES[20] - ALGEBRAIC[91])+4.75300))));
ALGEBRAIC[94] = ALGEBRAIC[92]/(ALGEBRAIC[92]+ALGEBRAIC[93]);
ALGEBRAIC[95] =  ( CONSTANTS[110]*ALGEBRAIC[94])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[53] = 1.00000/(1.00000+exp(7.48800 - STATES[20]/5.98000));
ALGEBRAIC[96] =  ( ( CONSTANTS[42]*CONSTANTS[58])*ALGEBRAIC[53])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[97] =  ( ( CONSTANTS[127]*CONSTANTS[58])*ALGEBRAIC[53])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[98] = ALGEBRAIC[96]+ALGEBRAIC[97];
ALGEBRAIC[54] = 1.00000/(1.00000+exp((STATES[20]+74.0000)/24.0000));
ALGEBRAIC[99] =  ( ( CONSTANTS[111]*STATES[26])*ALGEBRAIC[54])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[55] =  (1.00000/CONSTANTS[109])*log((CONSTANTS[52]+ CONSTANTS[61]*CONSTANTS[53])/(STATES[25]+ CONSTANTS[61]*STATES[28]));
ALGEBRAIC[56] =  ( ( CONSTANTS[42]*CONSTANTS[112])*pow(STATES[27], 2.00000))*(STATES[20] - ALGEBRAIC[55]);
ALGEBRAIC[57] =  ( ( CONSTANTS[127]*CONSTANTS[113])*pow(STATES[27], 2.00000))*(STATES[20] - ALGEBRAIC[55]);
ALGEBRAIC[58] = ALGEBRAIC[56]+ALGEBRAIC[57];
ALGEBRAIC[100] =  ( ( CONSTANTS[114]*STATES[29])*STATES[30])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[73] = ALGEBRAIC[71]+ALGEBRAIC[72];
ALGEBRAIC[101] =  ( ( CONSTANTS[120]*STATES[36])*STATES[37])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[113] = ((((((ALGEBRAIC[101]+ALGEBRAIC[99])+ALGEBRAIC[58])+ALGEBRAIC[95]) -  2.00000*ALGEBRAIC[73])+ALGEBRAIC[42])+ALGEBRAIC[98])+ALGEBRAIC[100];
ALGEBRAIC[117] = ((ALGEBRAIC[116]+ALGEBRAIC[78])+ALGEBRAIC[90])+ALGEBRAIC[113];
// ALGEBRAIC[118] = ((VOI - CONSTANTS[100]) -  CONSTANTS[101]*floor((VOI - CONSTANTS[100])/CONSTANTS[101])<CONSTANTS[99] ? 1.00000 : 0.00000);
ALGEBRAIC[119] =  ALGEBRAIC[118]*CONSTANTS[98];
RATES[20] = - (ALGEBRAIC[117]+ALGEBRAIC[119]);
}

void computeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[20]+91.0000)/6.10000));
ALGEBRAIC[0] = 1.00000/(1.00000+exp(- ((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000)/6.00000));
ALGEBRAIC[14] = ( ALGEBRAIC[0]*(1.00000 - exp(- ((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000)/6.00000)))/( 0.0350000*((STATES[20]+ 3.00000*CONSTANTS[23])+9.00000));
ALGEBRAIC[1] = 1.00000/(1.00000+exp(((STATES[20]+ 3.00000*CONSTANTS[23])+30.0000)/7.00000))+0.200000/(1.00000+exp(((50.0000 - STATES[20]) -  3.00000*CONSTANTS[23])/20.0000));
ALGEBRAIC[15] = 1.00000/( 0.0197000*exp(- pow( 0.0337000*((STATES[20]+ 3.00000*CONSTANTS[23])+25.0000), 2.00000))+0.0200000);
ALGEBRAIC[2] = 1.00000/(1.00000+exp(- (STATES[20]+10.0000)/5.00000));
ALGEBRAIC[16] = ( ((550.000)/(1.00000+exp((- 22.0000 - STATES[20])/9.00000)))*6.00000)/(1.00000+exp((STATES[20]+11.0000)/9.00000))+230.000/(1.00000+exp((STATES[20]+40.0000)/20.0000));
ALGEBRAIC[3] = 1.00000/(1.00000+exp(- ((STATES[20]+ 40.0000*CONSTANTS[23])+3.80000)/14.2500));
ALGEBRAIC[17] = 990.100/(1.00000+exp(- ((STATES[20]+ 40.0000*CONSTANTS[23])+2.43600)/14.1200));
ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[20]+6.00000)/- 8.60000));
ALGEBRAIC[18] = 9.00000/(1.00000+exp((STATES[20]+5.00000)/12.0000))+0.500000;
ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[20]+7.50000)/10.0000));
ALGEBRAIC[19] = 590.000/(1.00000+exp((STATES[20]+60.0000)/10.0000))+3050.00;
ALGEBRAIC[8] = 1.00000/pow(1.00000+exp(- (56.8600+STATES[20])/9.03000), 2.00000);
ALGEBRAIC[22] =  0.129200*exp(- pow((STATES[20]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[20] - 4.82300)/51.1200, 2.00000));
ALGEBRAIC[10] = ( 0.320000*(STATES[20]+47.1300))/(1.00000 - exp( - 0.100000*(STATES[20]+47.1300)));
ALGEBRAIC[23] =  0.0800000*exp(- STATES[20]/11.0000);
ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (STATES[20]+1.00000)/11.0000));
ALGEBRAIC[24] =  3.50000*exp(- pow(STATES[20]/30.0000, 2.00000))+1.50000;
ALGEBRAIC[12] = 1.00000/(1.00000+exp((STATES[20]+40.5000)/11.5000));
ALGEBRAIC[25] =  25.6350*exp(- pow((STATES[20]+52.4500)/15.8827, 2.00000))+24.1400;
ALGEBRAIC[27] = 1.00000/pow(1.00000+exp((STATES[20]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[6] =  (STATES[20]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[20]+80.0000)/6.80000))*1.00000;
ALGEBRAIC[20] =  (STATES[20]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[20]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[20])+ ( 3.10000*pow(10.0000, 5.00000))*exp( 0.348500*STATES[20]))*1.00000;
ALGEBRAIC[30] = 1.00000/(ALGEBRAIC[6]+ALGEBRAIC[20]);
ALGEBRAIC[28] = 1.00000/pow(1.00000+exp((STATES[20]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[7] =  (STATES[20]>=- 40.0000 ? 0.00000 :  (( ( ( - 2.54280*pow(10.0000, 4.00000))*exp( 0.244400*STATES[20]) -  ( 6.94800*pow(10.0000, - 6.00000))*exp( - 0.0439100*STATES[20]))*(STATES[20]+37.7800))/(1.00000+exp( 0.311000*(STATES[20]+79.2300))))*1.00000)*1.00000;
ALGEBRAIC[21] =  (STATES[20]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[20]))/(1.00000+exp( - 0.100000*(STATES[20]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[20]))/(1.00000+exp( - 0.137800*(STATES[20]+40.1400))))*1.00000;
ALGEBRAIC[31] = 1.00000/(ALGEBRAIC[7]+ALGEBRAIC[21]);
ALGEBRAIC[13] = ((1.00000 - STATES[40]) - STATES[39]) - STATES[38];
ALGEBRAIC[26] = CONSTANTS[88] - (CONSTANTS[88] - CONSTANTS[89])/(1.00000+pow(CONSTANTS[92]/STATES[18], 2.50000));
ALGEBRAIC[29] =  CONSTANTS[94]*ALGEBRAIC[26];
ALGEBRAIC[32] = CONSTANTS[122]/ALGEBRAIC[26];
ALGEBRAIC[33] = RATES[5]+RATES[3];
ALGEBRAIC[80] =  ((1.00000/CONSTANTS[109])/2.00000)*log(CONSTANTS[51]/STATES[12]);
ALGEBRAIC[81] =  ( CONSTANTS[42]*CONSTANTS[41])*(STATES[20] - ALGEBRAIC[80]);
ALGEBRAIC[36] =  (( ( ( CONSTANTS[105]*4.00000)*( ( STATES[20]*CONSTANTS[31])*CONSTANTS[109]))*( ( 0.341000*STATES[12])*exp( ( 2.00000*STATES[20])*CONSTANTS[109]) -  0.341000*CONSTANTS[51]))/(exp( ( 2.00000*STATES[20])*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[37] =  ( ( ( ( ( CONSTANTS[50]*ALGEBRAIC[36])*STATES[21])*STATES[22])*((1.00000 - STATES[23])+CONSTANTS[46]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[59] = 1.00000/(1.00000+pow(CONSTANTS[66]/STATES[12], 2.00000));
ALGEBRAIC[61] =  ( exp( ( CONSTANTS[73]*STATES[20])*CONSTANTS[109])*pow(STATES[16], 3.00000))*CONSTANTS[51];
ALGEBRAIC[62] =  ( exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109])*pow(CONSTANTS[53], 3.00000))*STATES[12];
ALGEBRAIC[63] = ((( ( CONSTANTS[67]*pow(CONSTANTS[53], 3.00000))*(1.00000+pow(STATES[16]/CONSTANTS[69], 3.00000))+ ( pow(CONSTANTS[70], 3.00000)*STATES[12])*(1.00000+STATES[12]/CONSTANTS[67]))+ CONSTANTS[68]*pow(STATES[16], 3.00000))+ pow(STATES[16], 3.00000)*CONSTANTS[51])+ pow(CONSTANTS[53], 3.00000)*STATES[12];
ALGEBRAIC[64] = (( ( ( ( CONSTANTS[42]*CONSTANTS[116])*pow(CONSTANTS[71], CONSTANTS[132]))*ALGEBRAIC[59])*(ALGEBRAIC[61] - ALGEBRAIC[62]))/ALGEBRAIC[63])/(1.00000+ CONSTANTS[72]*exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109]));
ALGEBRAIC[74] = pow( STATES[12]*1.00000, 1.60000);
ALGEBRAIC[75] = ( ( ( CONSTANTS[42]*pow(CONSTANTS[82], CONSTANTS[132]))*CONSTANTS[80])*ALGEBRAIC[74])/(ALGEBRAIC[74]+CONSTANTS[104]);
ALGEBRAIC[82] = ((ALGEBRAIC[37]+ALGEBRAIC[81])+ALGEBRAIC[75]) -  2.00000*ALGEBRAIC[64];
ALGEBRAIC[83] =  ( CONSTANTS[97]*STATES[39])*(STATES[18] - STATES[12]);
ALGEBRAIC[85] =  ( (1.00000+ 0.250000*CONSTANTS[34])*(STATES[18] - STATES[12]))*CONSTANTS[86];
ALGEBRAIC[35] = (((((RATES[10]+RATES[8])+RATES[9])+RATES[0])+RATES[1])+RATES[2])+RATES[7];
ALGEBRAIC[87] = ( ( pow(CONSTANTS[90], CONSTANTS[132])*CONSTANTS[91])*(pow(STATES[11]/CONSTANTS[121], CONSTANTS[93]) - pow(STATES[18]/CONSTANTS[87], CONSTANTS[93])))/((1.00000+pow(STATES[11]/CONSTANTS[121], CONSTANTS[93]))+pow(STATES[18]/CONSTANTS[87], CONSTANTS[93]));
ALGEBRAIC[34] = RATES[6]+RATES[4];
ALGEBRAIC[84] =  ((1.00000/CONSTANTS[109])/2.00000)*log(CONSTANTS[51]/STATES[13]);
ALGEBRAIC[86] =  ( CONSTANTS[127]*CONSTANTS[41])*(STATES[20] - ALGEBRAIC[84]);
ALGEBRAIC[38] =  (( ( ( CONSTANTS[105]*4.00000)*( ( STATES[20]*CONSTANTS[31])*CONSTANTS[109]))*( ( 0.341000*STATES[13])*exp( ( 2.00000*STATES[20])*CONSTANTS[109]) -  0.341000*CONSTANTS[51]))/(exp( ( 2.00000*STATES[20])*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[39] =  ( ( ( ( ( CONSTANTS[108]*ALGEBRAIC[38])*STATES[21])*STATES[22])*((1.00000 - STATES[24])+CONSTANTS[45]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[60] = 1.00000/(1.00000+pow(CONSTANTS[66]/STATES[13], 2.00000));
ALGEBRAIC[65] =  ( exp( ( CONSTANTS[73]*STATES[20])*CONSTANTS[109])*pow(STATES[17], 3.00000))*CONSTANTS[51];
ALGEBRAIC[66] =  ( exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109])*pow(CONSTANTS[53], 3.00000))*STATES[13];
ALGEBRAIC[67] = ((( ( CONSTANTS[67]*pow(CONSTANTS[53], 3.00000))*(1.00000+pow(STATES[17]/CONSTANTS[69], 3.00000))+ ( pow(CONSTANTS[70], 3.00000)*STATES[13])*(1.00000+STATES[13]/CONSTANTS[67]))+ CONSTANTS[68]*pow(STATES[17], 3.00000))+ pow(STATES[17], 3.00000)*CONSTANTS[51])+ pow(CONSTANTS[53], 3.00000)*STATES[13];
ALGEBRAIC[68] = (( ( ( ( CONSTANTS[127]*CONSTANTS[116])*pow(CONSTANTS[71], CONSTANTS[132]))*ALGEBRAIC[60])*(ALGEBRAIC[65] - ALGEBRAIC[66]))/ALGEBRAIC[67])/(1.00000+ CONSTANTS[72]*exp( ( (CONSTANTS[73] - 1.00000)*STATES[20])*CONSTANTS[109]));
ALGEBRAIC[76] = pow( STATES[13]*1.00000, 1.60000);
ALGEBRAIC[77] = ( ( ( CONSTANTS[127]*pow(CONSTANTS[82], CONSTANTS[132]))*CONSTANTS[80])*ALGEBRAIC[76])/(ALGEBRAIC[76]+CONSTANTS[124]);
ALGEBRAIC[89] = ((ALGEBRAIC[39]+ALGEBRAIC[86])+ALGEBRAIC[77]) -  2.00000*ALGEBRAIC[68];
ALGEBRAIC[43] =  (( ( ( ( CONSTANTS[107]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[16])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[53]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[44] =  ( ( ( ( ( CONSTANTS[50]*ALGEBRAIC[43])*STATES[21])*STATES[22])*((1.00000 - STATES[23])+CONSTANTS[46]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[102] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[53]/STATES[16]);
ALGEBRAIC[103] =  ( ( ( ( CONSTANTS[42]*CONSTANTS[115])*pow(STATES[33], 3.00000))*STATES[31])*STATES[32])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[104] =  ( CONSTANTS[42]*CONSTANTS[64])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[69] = 1.00000/((1.00000+ 0.124500*exp( ( - 0.100000*STATES[20])*CONSTANTS[109]))+ ( 0.0365000*CONSTANTS[118])*exp( - STATES[20]*CONSTANTS[109]));
ALGEBRAIC[71] = (( ( ( CONSTANTS[42]*CONSTANTS[74])*ALGEBRAIC[69])*CONSTANTS[52])/(1.00000+pow(CONSTANTS[117]/STATES[16], 4.00000)))/(CONSTANTS[52]+CONSTANTS[75]);
ALGEBRAIC[105] =  ( ( ( CONSTANTS[42]*CONSTANTS[119])*pow(STATES[35], 3.00000))*STATES[34])*(STATES[20] - ALGEBRAIC[102]);
ALGEBRAIC[114] = ((((ALGEBRAIC[103]+ALGEBRAIC[104])+ 3.00000*ALGEBRAIC[64])+ 3.00000*ALGEBRAIC[71])+ALGEBRAIC[44])+ALGEBRAIC[105];
ALGEBRAIC[45] =  (( ( ( ( CONSTANTS[107]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[17])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[53]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[46] =  ( ( ( ( ( CONSTANTS[108]*ALGEBRAIC[45])*STATES[21])*STATES[22])*((1.00000 - STATES[24])+CONSTANTS[45]))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[106] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[53]/STATES[17]);
ALGEBRAIC[107] =  ( ( ( ( CONSTANTS[127]*CONSTANTS[115])*pow(STATES[33], 3.00000))*STATES[31])*STATES[32])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[109] =  ( CONSTANTS[127]*CONSTANTS[64])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[72] = (( ( ( CONSTANTS[127]*CONSTANTS[74])*ALGEBRAIC[69])*CONSTANTS[52])/(1.00000+pow(CONSTANTS[117]/STATES[17], 4.00000)))/(CONSTANTS[52]+CONSTANTS[75]);
ALGEBRAIC[111] =  ( ( ( CONSTANTS[127]*CONSTANTS[119])*pow(STATES[35], 3.00000))*STATES[34])*(STATES[20] - ALGEBRAIC[106]);
ALGEBRAIC[115] = ((((ALGEBRAIC[107]+ALGEBRAIC[109])+ 3.00000*ALGEBRAIC[68])+ 3.00000*ALGEBRAIC[72])+ALGEBRAIC[46])+ALGEBRAIC[111];
ALGEBRAIC[90] = ALGEBRAIC[82]+ALGEBRAIC[89];
ALGEBRAIC[47] =  CONSTANTS[54]*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[49] =  (( CONSTANTS[42]*CONSTANTS[55])/(1.00000+CONSTANTS[56]/STATES[12]))*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[51] =  (( CONSTANTS[127]*CONSTANTS[55])/(1.00000+CONSTANTS[56]/STATES[13]))*(STATES[20] - CONSTANTS[125]);
ALGEBRAIC[52] = ALGEBRAIC[49]+ALGEBRAIC[51];
ALGEBRAIC[78] = ALGEBRAIC[52]+ALGEBRAIC[47];
ALGEBRAIC[116] = ALGEBRAIC[114]+ALGEBRAIC[115];
ALGEBRAIC[41] =  (( ( ( ( CONSTANTS[106]*STATES[20])*CONSTANTS[31])*CONSTANTS[109])*( ( 0.750000*STATES[25])*exp( STATES[20]*CONSTANTS[109]) -  0.750000*CONSTANTS[52]))/(exp( STATES[20]*CONSTANTS[109]) - 1.00000))*CONSTANTS[44];
ALGEBRAIC[42] =  ( ( ( ( ALGEBRAIC[41]*STATES[21])*STATES[22])*( CONSTANTS[50]*(CONSTANTS[46]+(1.00000 - STATES[23]))+ CONSTANTS[108]*(CONSTANTS[45]+(1.00000 - STATES[24]))))*pow(CONSTANTS[43], CONSTANTS[132]))*0.450000;
ALGEBRAIC[91] =  (1.00000/CONSTANTS[109])*log(CONSTANTS[52]/STATES[25]);
ALGEBRAIC[92] =  1.00000*(1.02000/(1.00000+exp( 0.238500*((STATES[20] - ALGEBRAIC[91]) - 59.2150))));
ALGEBRAIC[93] =  1.00000*(( 0.491240*exp( 0.0803200*((STATES[20]+5.47600) - ALGEBRAIC[91]))+exp( 0.0617500*((STATES[20] - ALGEBRAIC[91]) - 594.310)))/(1.00000+exp( - 0.514300*((STATES[20] - ALGEBRAIC[91])+4.75300))));
ALGEBRAIC[94] = ALGEBRAIC[92]/(ALGEBRAIC[92]+ALGEBRAIC[93]);
ALGEBRAIC[95] =  ( CONSTANTS[110]*ALGEBRAIC[94])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[53] = 1.00000/(1.00000+exp(7.48800 - STATES[20]/5.98000));
ALGEBRAIC[96] =  ( ( CONSTANTS[42]*CONSTANTS[58])*ALGEBRAIC[53])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[97] =  ( ( CONSTANTS[127]*CONSTANTS[58])*ALGEBRAIC[53])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[98] = ALGEBRAIC[96]+ALGEBRAIC[97];
ALGEBRAIC[54] = 1.00000/(1.00000+exp((STATES[20]+74.0000)/24.0000));
ALGEBRAIC[99] =  ( ( CONSTANTS[111]*STATES[26])*ALGEBRAIC[54])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[55] =  (1.00000/CONSTANTS[109])*log((CONSTANTS[52]+ CONSTANTS[61]*CONSTANTS[53])/(STATES[25]+ CONSTANTS[61]*STATES[28]));
ALGEBRAIC[56] =  ( ( CONSTANTS[42]*CONSTANTS[112])*pow(STATES[27], 2.00000))*(STATES[20] - ALGEBRAIC[55]);
ALGEBRAIC[57] =  ( ( CONSTANTS[127]*CONSTANTS[113])*pow(STATES[27], 2.00000))*(STATES[20] - ALGEBRAIC[55]);
ALGEBRAIC[58] = ALGEBRAIC[56]+ALGEBRAIC[57];
ALGEBRAIC[100] =  ( ( CONSTANTS[114]*STATES[29])*STATES[30])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[73] = ALGEBRAIC[71]+ALGEBRAIC[72];
ALGEBRAIC[101] =  ( ( CONSTANTS[120]*STATES[36])*STATES[37])*(STATES[20] - ALGEBRAIC[91]);
ALGEBRAIC[113] = ((((((ALGEBRAIC[101]+ALGEBRAIC[99])+ALGEBRAIC[58])+ALGEBRAIC[95]) -  2.00000*ALGEBRAIC[73])+ALGEBRAIC[42])+ALGEBRAIC[98])+ALGEBRAIC[100];
ALGEBRAIC[117] = ((ALGEBRAIC[116]+ALGEBRAIC[78])+ALGEBRAIC[90])+ALGEBRAIC[113];
// ALGEBRAIC[118] = ((VOI - CONSTANTS[100]) -  CONSTANTS[101]*floor((VOI - CONSTANTS[100])/CONSTANTS[101])<CONSTANTS[99] ? 1.00000 : 0.00000);
ALGEBRAIC[119] =  ALGEBRAIC[118]*CONSTANTS[98];
ALGEBRAIC[40] = ALGEBRAIC[37]+ALGEBRAIC[39];
ALGEBRAIC[48] = ALGEBRAIC[44]+ALGEBRAIC[46];
ALGEBRAIC[50] = (ALGEBRAIC[40]+ALGEBRAIC[42])+ALGEBRAIC[48];
ALGEBRAIC[70] = ALGEBRAIC[64]+ALGEBRAIC[68];
ALGEBRAIC[79] = ALGEBRAIC[75]+ALGEBRAIC[77];
ALGEBRAIC[88] = ALGEBRAIC[81]+ALGEBRAIC[86];
ALGEBRAIC[108] = ALGEBRAIC[103]+ALGEBRAIC[107];
ALGEBRAIC[110] = ALGEBRAIC[104]+ALGEBRAIC[109];
ALGEBRAIC[112] = ALGEBRAIC[105]+ALGEBRAIC[111];
}