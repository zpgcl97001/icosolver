//
// Created by Lin on 2020/7/27.
//

#ifndef ICOSOLVER_COUETTEFLOW_H
#define ICOSOLVER_COUETTEFLOW_H
#include <iostream>
#include "math.h"
#include <vector>
#include <stdio.h>
#include<algorithm>
#include <fstream>

struct Points{
    double x;
    double y;
    int I;
    int J;
    double Q;
    double QPlus;
    double P;
    double F;//nablaF
    double u;
    double v;
    double NablaV_;
};

struct Bluntbody{
    std::vector<int> PW_point;
    std::vector<int> PE_point;
    std::vector<int> PS_point;
    std::vector<int> PN_point;
    std::vector<int> uW_point;
    std::vector<int> uE_point;
    std::vector<int> uS_point;
    std::vector<int> uN_point;
    std::vector<int> vW_point;
    std::vector<int> vE_point;
    std::vector<int> vS_point;
    std::vector<int> vN_point;
    std::vector<int> vInner_point;
    std::vector<int> uInner_point;
    std::vector<int> PInner_point;
};


//Create Staggared grid




#endif //ICOSOLVER_COUETTEFLOW_H
