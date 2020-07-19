//
// Created by Lin on 2020/7/16.
//

#ifndef ICOSOLVER_MATSOLVER_H
#define ICOSOLVER_MATSOLVER_H



Class MatSol{
private:

public:
void MatMul(double**A,int A_row,int A_col,double**B,int B_row,int B_col){

    if(A_row != B_col){
        std::cout<<"Arow not equal to Bcol"<<std::endl;
        return ;
    }
    static double**  C = new double*[A_row];
    for (int i = 0; i <A_row ; ++i) {
        C[i] = new double[B_col];
    }
    for (int i = 0; i <A_row ; ++i) {
        for (int j = 0; j <A_col ; ++j) {
            for(int k=0;k<B_col;++k){
               C[i][k] += A[i][j]*B[j][k];
            }
        }
    }

    return C;



}

void Gaussian_Elimination(double **PMat,double **AMat,double *bMat,double *SolMat,int size,int iter){

    double **NewPMat = new double *[size];
    for(int i = 0;i<size;++i){
        PMat[i] = new double[size];
    }
    for(int i = 0;i<size;++i){
        NewPMat[i] = new double[size];
    }



    double **NewAMat = new double *[size];
    for(int i =0;i<size;++i) {
        NewAMat[i] = new double[size];
    }
    for (;iter <size -1; ++iter) {
        aii = AMat[iter][iter];
        if(abs(aii-0)>1e-10){
            for(int j =iter+1;j<size;++j) {
                bMat[j] -= AMat[j][iter]/AMat[iter][iter]*bMat[iter];
                for(k = iter+1;k<size;++k) {
                    AMat[j][k] -= AMat[j][iter] / AMat[iter][iter] * AMat[iter][k];
                }
            }
        }else{
            for (int j = 0; j <size ; ++k) {
                for (int k = 0; k <size ; ++k) {
                if(j<iter){
                    NewAMat[j][k] = AMat[j][k];
                    NewPMat[j][k] = PMat[j][k];
                }else{
                    int I = iter+ (j-iter+1)%(size-iter);
                    NewAMat[j][k] = Amat[I][k];
                }
                NewAMat[I][j] = AMat[i][j];
                NewPMat[I][j] = PMat[i][j];
                }
            }
            for (int i = 0; i <size ; ++i) {
                for (int j = 0; j <size ; ++j) {
                    AMat[i][j] = NewAMat[i][j];
                    PMat[i][j] = NewPMat[i][j];
                }
            }
            MatSol::Gaussian_Elimination(PMat,AMat,bMat,SolMat,size,iter);
        }

    }



}





};





#endif //ICOSOLVER_MATSOLVER_H
