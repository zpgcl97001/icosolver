//
// Created by Lin on 2020/7/21.
//

#include "Matsolver.h"
#include <iostream>

int main(){
    int Totalpoints = 3;
    MatSol matsol;

    double **AMat = new double *[Totalpoints];
    for (int i = 0; i < Totalpoints; ++i) {
        AMat[i] = new double[Totalpoints];
    };


AMat[0][0]= 1;
    AMat[0][1]=2;
    AMat[0][2]=3;
    AMat[1][0]=2;
    AMat[1][1]=5;
    AMat[1][2]=7;
    AMat[2][0]=3;
    AMat[2][1]=5;
    AMat[2][2]=3;

    double **LMat = new double*[Totalpoints];
    double **UMat = new double*[Totalpoints];
    double *BMat = new double[Totalpoints];
    double *xMat = new double[Totalpoints];
    BMat[0] = 1;
    BMat[1] = 7;
    BMat[2] = 8;

    for (int i = 0; i < Totalpoints; ++i) {
        LMat[i] = new double[Totalpoints];
        UMat[i] = new double[Totalpoints];
    }

    for (int i = 0; i <Totalpoints ; ++i) {
        xMat[i] = 0;
        for (int j = 0; j <Totalpoints ; ++j) {
            LMat[i][j] =0;
            UMat[i][j] =0;
        }

    }
    std::cout<<"decompose is run "<<std::endl;
    matsol.LUdecomp(AMat,Totalpoints,LMat,UMat);

    for (int i = 0; i <Totalpoints ; ++i) {
        for (int j = 0; j <Totalpoints ; ++j) {
            std::cout<<LMat[i][j]<<"\t"<<UMat[i][j]<<std::endl;
        }

    }

    matsol.SolveLU(LMat,UMat,BMat,xMat,Totalpoints);

//    matsol.PLUdecompose(BMat,CMat,AMat,3,0);
//    std::cout<<"P-----------------"<<std::endl;
//    for (int i = 0; i <Totalpoints ; ++i) {

//
//    }
//    std::cout<<"L-------------------"<<std::endl;
//    for (int i = 0; i <Totalpoints ; ++i) {
//        for (int j = 0; j <Totalpoints ; ++j) {
//            std::cout<<CMat[i][j]<<std::endl;
//        }
//
//    }

//    std::cout<<"U------------------"<<std::endl;
//    for (int i = 0; i <Totalpoints ; ++i) {
//        for (int j = 0; j <Totalpoints ; ++j) {
//            std::cout<<AMat[i][j]<<std::endl;
//        }
//
//    }


    double pauss;
    return 0;

}


