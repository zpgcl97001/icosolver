#include <iostream>
#include "math.h"
#include <vector>
#include <stdio.h>
#include "Deriv.h"
#include "icoSolver.h"
#include<algorithm>

int main(int argc, char *argv[]){


//Initalize and Mesh creating

double L = 1.0;//Unit Length
int XL = 10;
int YL = 3;
double Lx = 10 * L;
double Ly = 3 * L;
double Ng = 10;//number of grid per L
double Np = Ng +1 ;//number of points per L
double delx = L/(Ng);
double dely = L/(Ng);
int Nx = Ng * XL + 1;
int Ny = Ng * YL + 1;
double dt = 1e-3;
double Stoptime = 1;
int Totalsteps = 1000;

//Air property
double mu =  17.9e-6;//Pa¡¤s
double rho = 1.293; //kg/m3
double Re = 100;
double U_init;
U_init = Re * mu /(rho * L);
double atm =1.013e+5;



int Totalpoints = Nx*Ny - ( pow((Np-2),2));
int Totalcells  = (Nx-1) * (Ny-1) - pow(Ng,2);
std::cout<<"TotalPoint:"<<Totalpoints<<std::endl;

Bluffbody Square;
Bluffbody Domain;
std::vector<Points> PointDomain;
std::vector<Cells> CellDomain;

PointDomain.push_back(Points());
PointDomain[0].x = 0;
PointDomain[0].y = 0;
std::cout<<"Initialize Point Domain"<<std::endl;
for(int i=1; i<Totalpoints ; ++i){
    PointDomain.push_back(Points());
    if(abs(PointDomain[i-1].x - (Ng*XL*delx)) <1e-10 ){
        PointDomain[i].x = 0;
        PointDomain[i].y = PointDomain[i-1].y + dely;
    }
    else if(abs(PointDomain[i-1].x - 3 * Ng * delx)<1e-10 && PointDomain[i-1].y >= Ng *dely &&  PointDomain[i-1].y <= 2 * Ng *dely  ){
        PointDomain[i].x = PointDomain[i-1].x + Ng * delx;
        PointDomain[i].y = PointDomain[i-1].y;
    }
    else{
        PointDomain[i].x = PointDomain[i-1].x + delx;
        PointDomain[i].y = PointDomain[i-1].y;
    }
//    std::cout<<i<<"---"<<PointDomain[i].x<<","<<PointDomain[i].y<<std::endl;
}

//std::cout<<"Initialize cell domain"<<std::endl;
//int iter_ = 0;
//for(int i=0;i<Totalpoints;++i){
//    if(abs(PointDomain[i].x-PointDomain[i+1].x)< delx+1e-10 && abs(PointDomain[i+1].y -3)>1e-10 && Coor2pid(PointDomain[i+1].x,PointDomain[i+1].y+dely,PointDomain)>0){
//            CellDomain.push_back(Cells());
//            CellDomain[iter_].v1id = i;
//            CellDomain[iter_].v2id = i + 1;
//            CellDomain[iter_].v3id = Coor2pid(PointDomain[i + 1].x, PointDomain[i + 1].y + dely, PointDomain);
//            CellDomain[iter_].v4id = Coor2pid(PointDomain[i].x, PointDomain[i].y + dely, PointDomain);
//            CellDomain[iter_].x = 1 / 4 *
//                                  (PointDomain[CellDomain[iter_].v1id].x + PointDomain[CellDomain[iter_].v2id].x +
//                                   PointDomain[CellDomain[iter_].v3id].x + PointDomain[CellDomain[iter_].v4id].x);
//            CellDomain[iter_].y = 1 / 4 *
//                                  (PointDomain[CellDomain[iter_].v1id].y + PointDomain[CellDomain[iter_].v2id].y +
//                                   PointDomain[CellDomain[iter_].v3id].y + PointDomain[CellDomain[iter_].v4id].y);
//
//            iter_++;
//    }
//
//}
//std::cout<<iter_<<std::endl;
std::cout<<"Collect boundary points"<<std::endl;

//Add bluff boundary id to struct
for(int i=0; i < Totalpoints; ++i){
    if((PointDomain[i].x - 3 * Ng * delx) >= 0 && (PointDomain[i].x - 4 * Ng * delx-1e-10) <= 0 && abs(PointDomain[i].y - Ng * dely)<1e-10){
        Square.Sid.push_back(i);
        std::cout<<"Sid"<<i<<std::endl;
    }
    if((PointDomain[i].x - 3 * Ng * delx) >= 0 && (PointDomain[i].x - 4 * Ng * delx-1e-10) <= 0 && abs(PointDomain[i].y - 2* Ng * dely)<1e-10){
        Square.Nid.push_back(i);
        std::cout<<"Sid"<<i<<std::endl;
    }
    if(abs(PointDomain[i].x - 3 * Ng * delx) < 1e-10 && (PointDomain[i].y -  Ng * dely+1e-10) >= 0 && PointDomain[i].y - 2* Ng * dely -1e-10 <= 0){
        Square.Wid.push_back(i);
    }
    if(abs(PointDomain[i].x - 4 * Ng * delx) < 1e-10 && (PointDomain[i].y -  Ng * dely+1e-10) >= 0 && PointDomain[i].y - 2* Ng * dely -1e-10 <= 0){
        Square.Eid.push_back(i);
    }
    if(abs(PointDomain[i].x - 0)<1e-10){
        Domain.Wid.push_back(i);
    }
    if((abs(PointDomain[i].x)-XL*Ng*delx)<1e-10){
        Domain.Eid.push_back(i);
    }
    if(abs(PointDomain[i].y-0)<1e-10){
        Domain.Sid.push_back(i);
    }
    if(abs(PointDomain[i].y-YL*Ng*dely)<1e-10){
        Domain.Nid.push_back(i);
    }
}





//Initialize BC and IC
for (int i = 0; i <Totalpoints ; ++i) {
   PointDomain[i].P = 1.013e+5;
   PointDomain[i].Pplus = 1.013e+5;
   PointDomain[i].U[0] = U_init;
   PointDomain[i].U[1] = 0;
   PointDomain[i].Uplus[0] = U_init;
   PointDomain[i].Uplus[1] = 0;
}
//Noslip condition
for (int i = 0; i <Totalpoints ; ++i) {
    if(abs(PointDomain[i].y-0) < 1e-10 || abs(PointDomain[i].y-YL*Ng*dely)<1e-10 ){
        PointDomain[i].U[0]=0;
        PointDomain[i].U[1]=0;
    }
}

for (auto point:Square.Sid){
    PointDomain[point].U[0]=0;
    PointDomain[point].U[1]=0;
}

for (auto point:Square.Eid){
    PointDomain[point].U[0]=0;
    PointDomain[point].U[1]=0;
}
for (auto point:Square.Wid){
    PointDomain[point].U[0]=0;
    PointDomain[point].U[1]=0;
}
for (auto point:Square.Nid){
    PointDomain[point].U[0]=0;
    PointDomain[point].U[1]=0;
}
//Pressure condition
//testzone
double **PMat = new double *[Totalpoints];
double **PPMat = new double *[size];
double *BMat = new double[Totalpoints];
for (int i = 0; i < Totalpoints; ++i) {
    PMat[i] = new double[Totalpoints];
    PPMat[i] = new double[Totalpoints];
};

for(int i = 0;i<size;++i){
     PMat[i][i] = 1;
}


std::cout<<"Building Mat"<<std::endl;
BuildLaplacePMat(PMat,BMat,PointDomain,Totalpoints,Square,Domain,delx,dely, rho, mu);
std::cout<<"Output PMat"<<std::endl;
std::cout<<PMat<<std::endl;



//Inital time loop
for(double step = 0 ;step<Totalsteps;++step){

    







}






}

int Coor2pid(double x,double y,std::vector<Points> PointDomain){//from Point coor to id ,
    int iter_ = 0;
    for(auto point : PointDomain){
        if(abs(point.x - x)<1e-10 && abs(point.y - y)<1e-10){
            return iter_;
        }
        iter_++;
    }

    std::cout << "Point not found"<<x<<","<<y<<std::endl;
    return -1;
};

int Coor2cid(double x,double y,std::vector<Cells> CellDomain){//from Point coor to id ,
    int iter_ = 0;
    for(auto cell : CellDomain){
        if(abs(cell.x - x)<1e-10 && abs(cell.y - y)<1e-10){
            return iter_;
        }
        iter_++;
    }

    std::cout << "Cell not found"<<std::endl;
    return -1;
};

void NablaP(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely) {
    double dpdx;
    double dpdy;
    dpdx = (CellDomain[cid + 1].P - CellDomain[cid - 1].P) / (2 * delx);
    int Nid = Coor2cid(CellDomain[cid].x, CellDomain[cid].y + dely, CellDomain);
    int Sid = Coor2cid(CellDomain[cid].x, CellDomain[cid].y - dely, CellDomain);
    dpdy = (CellDomain[Nid].y - CellDomain[Sid].y) / (2 * delx);
    double Nabla = dpdx + dpdy;
}
void NablaU(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely) {
    double ue;
    double uw;
    double vn;
    double vs;
    ue = (CellDomain[cid+1].U[0] + CellDomain[cid].U[0])/2;
    uw = (CellDomain[cid].U[0] + CellDomain[cid-1].U[0])/2;
    int Nid = Coor2cid(CellDomain[cid].x,CellDomain[cid].y+dely,CellDomain);
    int Sid = Coor2cid(CellDomain[cid].x,CellDomain[cid].y-dely,CellDomain);
    vn = (CellDomain[cid].U[1] + CellDomain[Nid].U[1])/2;
    vs = (CellDomain[cid].U[1] + CellDomain[Sid].U[1])/2;
    double NablaU = ue*dely-uw*dely+vn*dely-vs*delx;

}
void LaplaceU(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely) {
        double ue;
        double uw;
        double vn;
        double vs;
        ue = (CellDomain[cid+1].U[0] + CellDomain[cid].U[0])/2;
        uw = (CellDomain[cid].U[0] + CellDomain[cid-1].U[0])/2;
        int Nid = Coor2cid(CellDomain[cid].x,CellDomain[cid].y+dely,CellDomain);
        int Sid = Coor2cid(CellDomain[cid].x,CellDomain[cid].y-dely,CellDomain);
        vn = (CellDomain[cid].U[1] + CellDomain[Nid].U[1])/2;
        vs = (CellDomain[cid].U[1] + CellDomain[Sid].U[1])/2;
        double LaplaceU = (ue-uw)/delx + (vn-vs)/dely;
}

double* Predictor(std::vector<Points> PointDomain,int pid,double delx,double dely,double delt,double mu){
    double U_star[2];
    U_star[0] = SolveF(PointDomain,pid,delx,dely,mu)[0]*delt+PointDomain[pid].U[0];
    U_star[1] = SolveF(PointDomain,pid,delx,dely,mu)[1]*delt+PointDomain[pid].U[1];
    return U_star;
}

double* SolveF(std::vector<Points> PointDomain,int pid,double delx,double dely,double mu){
    static double F[2];
    F[0] = -PointDomain[pid].U[0]*divU(PointDomain,pid,delx,dely)+mu*nablaUx(PointDomain,pid,delx);
    F[1] = -PointDomain[pid].U[1]*divU(PointDomain,pid,delx,dely)+mu*nablaUy(PointDomain,pid,dely);
    return F;
}

double divU(std::vector<Points> PointDomain,int pid,double delx,double dely){
    double dudx;
    double dudy;
    double res;
    dudx = (PointDomain[pid-1].U[0] - PointDomain[pid+1].U[0])/(2*delx);
    int Nid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y+dely,PointDomain);
    int Sid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y-dely,PointDomain);
    dudy =(PointDomain[Nid].U[1]-PointDomain[Sid].U[1])/(2*dely);
    res = dudx + dudy;
    return res;
};
double nablaUx(std::vector<Points> PointDomain,int pid,double delphi){
    double dudphi;
    double dudy;
    double res;
    int Nid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y+delphi,PointDomain);
    int Sid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y-delphi,PointDomain);
    double nablaU;
    nablaU = (PointDomain[pid-1].U[0]-2*PointDomain[pid].U[0]+PointDomain[pid+1].U[0])/pow(delphi,2);

    return nablaU;
};
double nablaUy(std::vector<Points> PointDomain,int pid,double delphi){
    double dudphi;
    double dudy;
    double res;
    int Nid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y+delphi,PointDomain);
    int Sid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y-delphi,PointDomain);
    double nablaU;
    nablaU = (PointDomain[Nid].U[1]-2*PointDomain[pid].U[1]+PointDomain[Sid].U[0])/pow(delphi,2);
    return nablaU;
};

double* gradP(std::vector<Points> PointDomain,int pid,double delx,double dely){
    static double gradP[2];
    gradP[0] = (PointDomain[pid+1].P - PointDomain[pid-1].P)/(2*delx);
    int Nid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y+dely,PointDomain);
    int Sid = Coor2pid(PointDomain[pid].x,PointDomain[pid].y-dely,PointDomain);
    gradP[1] = (PointDomain[Nid].P-PointDomain[Sid].P)/(2*dely);
    return gradP;
}

//double* Corrector(){
//    double* U[2];
//    U[0] = U_star[0] - delt/rho * gradP();
//    U[1] = U_star[0] - delt/rho * gradP();
//}

void BuildLaplacePMat(double **PMat,double *BMat,std::vector<Points> PointDomain,int TotalPoints,Bluffbody Square,Bluffbody Domain,double delx,double dely,double rho,double mu) {

    double x, y, dpdx, dpdy;
    int Eid, Wid, Sid, Nid, SSid, NNid;
    //Create Mat
    for (int i = 0; i < TotalPoints; ++i) {
        for (int j = 0; j < TotalPoints; ++j) {

            if (InVector(i, Square.Wid)) {//Blunt West
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Nid = Coor2pid(x, y + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                dpdx = mu * (PointDomain[i].U[0] - 2 * PointDomain[i - 1].U[0] + PointDomain[i - 2].U[0]) /
                       pow(delx, 2) / delx;
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Wid = i - 1;
                Eid = i + 1;
                Nid = Coor2pid(x, y + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                PMat[i][i - 2] += -1 / 2 * pow(delx, 2);
                PMat[i][i] += 1 / 2 * pow(delx, 2);
                BMat[i] = rho * (-(pow(PointDomain[i].U[0], 2) - pow(PointDomain[i - 2].U[0], 2)) / (2 * delx) +
                                 mu * (PointDomain[i].U[0] - 2 * PointDomain[i - 1].U[0] + PointDomain[i - 2].U[0]) /
                                 pow(delx, 2)) - dpdx;
            } else if (InVector(i, Square.Eid)) {//Blunt East
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Nid = Coor2pid(x, y + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                dpdx = mu * (PointDomain[i].U[0] - 2 * PointDomain[i + 1].U[0] + PointDomain[i + 2].U[0]) /
                       pow(delx, 2) / delx;
                BMat[i] = rho * (-(pow(PointDomain[i + 2].U[0], 2) - pow(PointDomain[i].U[0], 2)) / (2 * delx) +
                                 mu * (PointDomain[i].U[0] - 2 * PointDomain[i + 1].U[0] + PointDomain[i + 2].U[0]) /
                                 pow(delx, 2)) - dpdx;
                PMat[i][i + 2] += -1 / 2 * pow(delx, 2);
                PMat[i][i] += 1 / 2 * pow(delx, 2);
            } else if (InVector(i, Square.Nid)) {//Blunt North
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Nid = Coor2pid(x, y + dely, PointDomain);
                NNid = Coor2pid(x, y + dely + dely, PointDomain);
                dpdy = mu * (PointDomain[i].U[1] - 2 * PointDomain[Nid].U[1] + PointDomain[NNid].U[1]) / pow(dely, 2) /
                       dely;
                BMat[i] = rho * (-(pow(PointDomain[NNid].U[1], 2) - pow(PointDomain[i].U[1], 2)) / (2 * dely) +
                                 mu * (PointDomain[i].U[1] - 2 * PointDomain[Nid].U[1] + PointDomain[NNid].U[1]) /
                                 pow(dely, 2)) - dpdy;
                PMat[i][NNid] += -1 / 2 * pow(dely, 2);
                PMat[i][i] += 1 / 2 * pow(dely, 2);
            } else if (InVector(i, Square.Sid)) {//Blunt South
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Sid = Coor2pid(x, y - dely, PointDomain);
                SSid = Coor2pid(x, y - dely - dely, PointDomain);
                dpdy = mu * (PointDomain[i].U[1] - 2 * PointDomain[Sid].U[1] + PointDomain[SSid].U[1]) / pow(dely, 2) /
                       dely;
                BMat[i] = rho * (-(pow(PointDomain[i].U[1], 2) - pow(PointDomain[SSid].U[1], 2)) / (2 * dely) +
                                 mu * (PointDomain[i].U[1] - 2 * PointDomain[Sid].U[1] + PointDomain[SSid].U[1]) /
                                 pow(dely, 2)) - dpdy;
                PMat[i][SSid] += -1 / 2 * pow(dely, 2);
                PMat[i][i] += 1 / 2 * pow(dely, 2);
            } else if (InVector(i, Domain.Nid)) {//Domain North
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                SSid = Coor2pid(x, y - dely - dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                dpdy = mu * (PointDomain[i].U[1] - 2 * PointDomain[Sid].U[1] + PointDomain[SSid].U[1]) / pow(dely, 2) /
                       dely;
                //dpdx = 0
                BMat[i] = rho * (-(pow(PointDomain[i].U[1], 2) - pow(PointDomain[SSid].U[1], 2)) / (2 * dely) +
                                 mu * (PointDomain[i].U[1] - 2 * PointDomain[Sid].U[1] + PointDomain[SSid].U[1]) /
                                 pow(dely, 2)) - dpdy;
                PMat[i][SSid] += -1 / 2 * pow(dely, 2);
                PMat[i][i] += 1 / 2 * pow(dely, 2);
            } else if (InVector(i, Domain.Sid)) {//Domain South
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Nid = Coor2pid(x, y + dely, PointDomain);
                NNid = Coor2pid(x, y + dely + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                dpdy = mu * (PointDomain[i].U[1] - 2 * PointDomain[Nid].U[1] + PointDomain[NNid].U[1]) / pow(dely, 2) /
                       dely;
                BMat[i] = rho * (-(pow(PointDomain[NNid].U[1], 2) - pow(PointDomain[i].U[1], 2)) / (2 * dely) +
                                 mu * (PointDomain[i].U[1] - 2 * PointDomain[Nid].U[1] + PointDomain[NNid].U[1]) /
                                 pow(dely, 2)) - dpdy;
                PMat[i][NNid] += -1 / 2 * pow(dely, 2);
                PMat[i][i] += 1 / 2 * pow(dely, 2);
            } else if (InVector(i, Domain.Wid)) {//Domain inlet
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Nid = Coor2pid(x, y + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                dpdx = (mu * (PointDomain[i].U[0] - 2 * PointDomain[i + 1].U[0] + PointDomain[i + 2].U[0]) /
                        pow(delx, 2) + (pow(PointDomain[i + 2].U[0], 2) - pow(PointDomain[i].U[0], 2)) / delx) / delx;
                BMat[i] = rho * (-(pow(PointDomain[i + 2].U[0], 2) - pow(PointDomain[i].U[0], 2)) / (2 * delx) +
                                 mu * (PointDomain[i].U[0] - 2 * PointDomain[i + 1].U[0] + PointDomain[i + 2].U[0]) /
                                 pow(delx, 2)) - dpdx;
                PMat[i][i + 2] += -1 / 2 * pow(delx, 2);
                PMat[i][i] += 1 / 2 * pow(delx, 2);
            } else if (InVector(i, Domain.Eid)) {//Domain outlet
                PMat[i][i] += 1;//equal to atm
                BMat[i] = 1.013e+5;//atm;
            } else {
                x = PointDomain[i].x;
                y = PointDomain[i].y;
                Wid = i - 1;
                Eid = i + 1;
                Nid = Coor2pid(x, y + dely, PointDomain);
                Sid = Coor2pid(x, y - dely, PointDomain);
                PMat[i][i] += -2 / pow(delx, 2) - 2 / pow(dely, 2);
                PMat[i][Wid] += 1 / pow(delx, 2);
                PMat[i][Eid] += 1 / pow(delx, 2);
                PMat[i][Nid] += 1 / pow(dely, 2);
                PMat[i][Sid] += 1 / pow(dely, 2);
                BMat[i] = rho * (-(pow(PointDomain[i + 1].U[0], 2) - pow(PointDomain[i - 1].U[0], 2)) / (2 * delx) -
                                 (pow(PointDomain[Nid].U[1], 2) - pow(PointDomain[Sid].U[1], 2)) / (2 * dely) +
                                 mu * (PointDomain[i - 1].U[0] - 2 * PointDomain[i].U[0] + PointDomain[i + 1].U[0]) /
                                 pow(delx, 2) / +mu *
                                 (PointDomain[Nid].U[1] - 2 * PointDomain[i].U[1] + PointDomain[Sid].U[1]) /
                                 pow(dely, 2));
            }
        }


    }

    //Solve Mat




}


bool InVector(int elem,std::vector<int> v){
    std::vector<int>::iterator it;
    it = std::find(v.begin(),v.end(),elem);
    if(it != v.end()){
        return true;
    }else{
        return false;
    }
}