//
// Created by Lin on 2020/7/21.
//

#include "Matsolver.h"
#include <iostream>

int main(){
    int Totalpoints = 3;
    MatSol matsol;

    double **AMat = new double *[Totalpoints];
    double **BMat = new double *[Totalpoints];
    double **CMat = new double *[Totalpoints];
    for (int i = 0; i < Totalpoints; ++i) {
        AMat[i] = new double[Totalpoints];
        BMat[i] = new double[Totalpoints];
        CMat[i] = new double[Totalpoints];
    };

    for (int i = 0; i <Totalpoints ; ++i) {
        for (int j = 0; j <Totalpoints ; ++j) {
            BMat[i][i] =  1;
        }

    }
AMat[0][0]= 1;
    AMat[0][1]=2;
    AMat[0][2]=3;
    AMat[1][0]=2;
    AMat[1][1]=4;
    AMat[1][2]=7;
    AMat[2][0]=3;
    AMat[2][1]=5;
    AMat[2][2]=3;

    matsol.PLUdecompose(BMat,CMat,AMat,3,0);
    std::cout<<"P-----------------"<<std::endl;
    for (int i = 0; i <Totalpoints ; ++i) {
        for (int j = 0; j <Totalpoints ; ++j) {
            std::cout<<BMat[i][j]<<std::endl;
        }

    }
    std::cout<<"L-------------------"<<std::endl;
    for (int i = 0; i <Totalpoints ; ++i) {
        for (int j = 0; j <Totalpoints ; ++j) {
            std::cout<<CMat[i][j]<<std::endl;
        }

    }

    std::cout<<"U------------------"<<std::endl;
    for (int i = 0; i <Totalpoints ; ++i) {
        for (int j = 0; j <Totalpoints ; ++j) {
            std::cout<<AMat[i][j]<<std::endl;
        }

    }


    double pauss;
    return 0;

}


