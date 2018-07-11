#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "CrossSection.h"

int main() {
    print_dsigma_dcostheta("dcs.txt"); 
    print_sigma("cs.txt");
    return 0;
}

