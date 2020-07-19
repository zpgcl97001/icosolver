#include <iostream>
#include "math.h"
#include <vector>
#include <stdio.h>
#include "Deriv.h"
//#include "icoSolver.h";

struct Points{
    double x;
    double y;
    double U[2];
    double P;
    double Pplus;
    double Uplus[2];

};

struct Cells{
    double x;
    double y;
    int v1id;
    int v2id;
    int v3id;
    int v4id;
    double U[2];
    double Uplus[2];
    double P;
    double Pplus[2];
    double F;
    double Fplus;
};
struct Bluffbody{
    std::vector<int> Eid;
    std::vector<int> Wid;
    std::vector<int> Nid;
    std::vector<int> Sid;
};
int Coor2pid(double x,double y,std::vector<Points> PointDomain);
int Coor2cid(double x,double y,std::vector<Points> CellDomain);
void Nabla(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int pid,int type);
void NablaP(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely);
void NablaU(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely);
void LaplaceU(std::vector<Points> PointDomain,std::vector<Cells> CellDomain,int cid,double delx,double dely);
double* Predictor(std::vector<Points> PointDomain,int pid,double delx,double dely,double delt,double mu);

double* SolveF(std::vector<Points> PointDomain,int pid,double delx,double dely,double mu);
double divU(std::vector<Points> PointDomain,int pid,double delx,double dely);
double nablaUx(std::vector<Points> PointDomain,int pid,double delphi);
double nablaUy(std::vector<Points> PointDomain,int pid,double delphi);
double* gradP(std::vector<Points> PointDomain,int pid,double delx,double dely);
void BuildLaplacePMat(double **PMat,double* BMat,std::vector<Points> PointDomain,int TotalPoints,Bluffbody Square,Bluffbody Domain,double delx,double dely,double rho,double mu);
bool InVector(int elem,std::vector<int> v);
