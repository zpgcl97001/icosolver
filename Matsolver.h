//
// Created by Lin on 2020/7/16.
//

#ifndef ICOSOLVER_MATSOLVER_H
#define ICOSOLVER_MATSOLVER_H

#include <math.h>
#include <iostream>

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
    int b;

void PLUdecompose(double **PMat,double **LMat,double **AMat,int sizen,int iter){

    iter++;
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
                LMat[u][v]= 1;
            }else {
                LMat[u][v] = 0;
            }
            }

    }

    //End initialize
    for (int i=0;i<sizen-1; ++i) {
        double aii = AMat[i][i];
        std::cout<<"aii="<<aii<<"      step="<<i<<std::endl;
        if(abs(aii - 0)<1e-5){
            iter++;
            for (int j = 0; j <sizen ; ++j) {
                for (int k = 0; k <sizen ; ++k) {
                    if(j<i){
                       AMat[j][k] = AMat_[j][k];
                       NewPMat[j][k] = PMat[j][k];
                    }else{
                        int I = i+ (j-i+1)%(sizen-i);
                        AMat[j][k] = AMat_[I][k];
                        NewPMat[j][k] = PMat[I][k];
                    }
                    std::cout<<"Ajk"<<AMat[j][k]<<std::endl;
                }
            }
            for (int u = 0; u <sizen ; ++u) {
                for (int v = 0; v <sizen ; ++v) {
                    PMat[u][v] = NewPMat[u][v];
                }
            }
            if(iter<sizen){
                MatSol::PLUdecompose(PMat,LMat,AMat,sizen,iter);
            }else{
                std::cout<<"Singular Matrix"<<std::endl;
                return;
            }
        }

        for (int u = 0; u < sizen; ++u) {
            for (int v = 0; v <sizen ; ++v) {
                if(v==u){
                    Lninv[u][v] = 1;
                }
                else{
                    Lninv[u][v] = 0;
                }

            }
        }
            for(int j =i+1;j<sizen;++j) {
                Lninv[j][i] = AMat[j][i]/ AMat[i][i];
                std::cout<<j<<"   ani/aii="<<AMat[j][i]/ AMat[i][i]<<std::endl;
                for(int k = 0;k<sizen;++k) {
                    NewAMat[j][k] =AMat[j][k] - AMat[j][i] / AMat[i][i] * AMat[i][k];
                }
            }
        for (int u=i+1;u<sizen;++u){
            for (int v=0;v<sizen;++v){
                AMat[u][v]=NewAMat[u][v];
            }
        }
            //intital MatC;
            for (int u=0;u<sizen;++u){
                for (int v=0;v<sizen;++v){
                    NewLMat[u][v] = 0;
                }
            }
            //calculate newL
            MatSol::MatMul(LMat,Lninv,NewLMat,sizen,sizen,sizen,sizen,sizen,sizen);
            std::cout<<"Lninv"<<std::endl;
        for (int u=0;u<sizen;++u){
            for (int v=0;v<sizen;++v){
                std::cout<<Lninv[u][v]<<std::endl ;
            }
        }
            //update LMat
            for (int u=0;u<sizen;++u){
                for (int v=0;v<sizen;++v){
                    LMat[u][v] = NewLMat[u][v];
                }
            }


        }
    //Calculate P L U
    delete Lninv;
    delete NewAMat;
    delete NewPMat;
    delete NewLMat;
    delete AMat_;
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
