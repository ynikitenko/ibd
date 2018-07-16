#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "CrossSection.h"

int main() {
    print_dsigma_dcostheta("dcs.txt"); 

    switchCVC = true;
    neutrinoEMax = 4;
    neutrinoERes = 20;
    print_sigma("cs_cvc.txt");
    switchCVC = false;
    print_sigma("cs_no_cvc.txt");

    return 0;
}

