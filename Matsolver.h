//
// Created by Lin on 2020/7/16.
//

#ifndef ICOSOLVER_MATSOLVER_H
#define ICOSOLVER_MATSOLVER_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>

class MatSol{
private:
    double **PMat;
    double **AMat;
    double *bMat;
    double *SolMat;
    int sizen;
    double **MatA;
    double **MatB;
    double **MatC;
    int rowA;
    int colA;
    int rowB;
    int colB;
    int rowC;
    int colC;
public:

void PLUdecompose(double **PMat,double **LMat,double **AMat,int sizen,int iter){

    double **Lninv = new double *[sizen];
    double **NewPMat = new double *[sizen];
    double **NewAMat = new double *[sizen];
    double **NewLMat = new double *[sizen];
    double **AMat_ =new double *[sizen];
    for(int i = 0;i<sizen;++i){
        NewPMat[i] = new double[sizen];
        NewAMat[i] = new double[sizen];
        Lninv[i]   = new double[sizen];
        NewLMat[i] = new double[sizen];
        AMat_[i]   = new double[sizen];
    }

    //copy AMat value
    for (int u = 0; u <sizen ; ++u) {
        for (int v = 0; v < sizen ; ++v) {
            AMat_[u][v] = AMat[u][v];
            if(u==v){
                LMat[u][v] = 1;
            }else {
                LMat[u][v] = 0;
            }
            }

    }



    //End initialize
    for (int i=0;i<sizen-1; ++i) {
        double aii = AMat[i][i];
        std::cout << "Decomposing process is on! Iter=" << iter << "   i=" << i << std::endl;
        if (abs(aii - 0) < 1e-5) {
//            std::cout<<"Aii=0! return"<<std::endl;
//            return;

            iter++;
            for (int j = 0; j < sizen; ++j) {
                for (int k = 0; k < sizen; ++k) {
                    if (j < i) {
                        AMat[j][k] = AMat_[j][k];
                        NewPMat[j][k] = PMat[j][k];
                    } else {
                        int I = i + (j - i + 1) % (sizen - i);
                        AMat[j][k] = AMat_[I][k];
                        NewPMat[j][k] = PMat[I][k];
                    }
                }
            }
            for (int u = 0; u < sizen; ++u) {
                for (int v = 0; v < sizen; ++v) {
                    PMat[u][v] = NewPMat[u][v];
                }
            }
            if (iter < sizen) {
                std::cout << "Iter=" << iter << std::endl;
                MatSol::PLUdecompose(PMat, LMat, AMat, sizen, iter);
            } else {
                std::cout << "Singular Matrix" << std::endl;
                return;
            }
        }else {


            for (int j = 0; j < sizen; ++j) {
                for (int k = 0; k < sizen; ++k) {
                    NewLMat[j][k] = 0;
                    if (j == k) {
                        Lninv[j][k] = 1;
                    } else {
                        Lninv[j][k] = 0;
                    }
                    if (j > i) {
                        Lninv[j][i] = AMat[j][i] / AMat[i][i];
                        NewAMat[j][k] = AMat[j][k] - AMat[j][i] / AMat[i][i] * AMat[i][k];

                    }

                }

            }

            MatSol::MatMul(LMat, Lninv, NewLMat, sizen, sizen, sizen, sizen, sizen, sizen);

            //update LMat
            for (int u = 0; u < sizen; ++u) {
                for (int v = 0; v < sizen; ++v) {
                    LMat[u][v] = NewLMat[u][v];
                    if (u > i) {
                        AMat[u][v] = NewAMat[u][v];
                    }
                }
            }
            std::cout << "Iter round " << i << " finish!" << std::endl;

        }
    }

    //Calculate P L U
    delete[] Lninv;
    delete[] NewAMat;
    delete[] NewPMat;
    delete[] NewLMat;
    delete[] AMat_;

    };

void LUdecomp(double** AMat,int sizen,double **LMat,double **UMat){


    for (int i = 0; i <sizen ; ++i) {
        std::cout<<"LU---Iter"<<i<<std::endl;
        for (int k = i; k < sizen; ++k) {
            double sum = 0;
            for (int j = 0; j <sizen ; ++j) {
                sum += (LMat[i][j] * UMat[j][k]);
            }
            UMat[i][k] = AMat[i][k]-sum;
        }


        for (int k = i; k < sizen; ++k) {
            if (i == k){
                LMat[i][i] = 1;
            } else {
                double sum = 0;
                for (int j = 0; j <  i ; ++j) {
                    sum += (LMat[k][j] * UMat[j][i]);
                }
                LMat[k][i] = (AMat[k][i] - sum)/UMat[i][i];
            }

        }

    }




}

void SolveLU(double** LMat,double** UMat,double *BMat,double* xMat,int sizen){
    double * yMat = new  double[sizen];
    for (int i = 0; i <sizen ; ++i) {
        yMat[i] = 0;
    }

    //sovle Y
    for (int i = 0; i <sizen ; ++i) {
        if(i==0){
            yMat[0] = BMat[0]/LMat[0][0];
        }else {
            yMat[i] += BMat[i];
            for (int j = 0; j < i ; ++j) {
                yMat[i] -=   LMat[i][j] * yMat[j];
            }
        }
    }

    for (int k = 0; k <sizen ; ++k) {
        std::cout<<yMat[k]<<std::endl;
    }
    //Solve X
    std::cout<<"--------------------------"<<std::endl;
    for (int i = sizen-1; i >=0 ; --i) {
        if(i == sizen-1){
            xMat[i] = yMat[i]/UMat[i][i];
        }
        else{
//            std::cout<<"Line 200 is run"<<std::endl;
            xMat[i] += yMat[i]/UMat[i][i];
            for (int j = i+1; j <sizen ; ++j) {
                xMat[i] -= 1/UMat[i][i] * UMat[i][j]*xMat[j];
            }
        }

    }
    for (int k = 0; k <sizen ; ++k) {
        std::cout<<xMat[k]<<std::endl;
    }
}


    void MatMul(double **MatA,double **MatB,double **MatC,int rowA,int colA,int rowB,int colB,int rowC,int colC){

        if(colA != rowB || rowA !=rowC || colB !=colC){
            std::cout<<"MatMul error"<<std::endl;
        }
        for (int i=0;i<rowA;++i){
            for (int j=0;j<colA;++j){
                for (int k = 0; k < colB; ++k) {
                    MatC[i][k] += MatA[i][j]*MatB[j][k];
                }
            }
        }
    }







};








#endif //ICOSOLVER_MATSOLVER_H
