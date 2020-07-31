#include <iostream>
#include "math.h"
#include <vector>
#include <stdio.h>
#include "Deriv.h"
#include <fstream>
#include <iomanip>
//#include "icoSolver.h";

//struct Points{
//    double x;
//    double y;
//    double U[2];
//    double P;
//    double Pplus;
//    double Uplus[2];
//
//};

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

class IcoSolver{


public:
    double Fx(std::vector<Points> u_Domain,std::vector<Points> v_Domain,int id,double mu,double delx,double dely,int Nx,std::string type){//uDomain
        double dudx,dvdy,ddudx;
        int vid,I,J;
        double uij,vij,vijm,vijp,uipj,uippj,uimj,uimmj,Nablau,Laplaceu,uijm,uijp;
        if(type =="P") {
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            vid = I+ J*(Nx-1);
            uij = u_Domain[id].Q;
            vij = (v_Domain[vid + Nx - 1].Q + v_Domain[vid + Nx - 2].Q + v_Domain[vid - 1].Q + v_Domain[vid].Q)/4;
            uipj = u_Domain[id+1].Q;
            uimj = u_Domain[id-1].Q;
            uijp = u_Domain[id+Nx].Q;
            uijm = u_Domain[id-Nx].Q;
            Nablau =  uij*(uipj-uimj)/(2*delx) + vij*(uijp-uijm)/(2*dely);
            Laplaceu = (uipj - 2*uij + uimj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely,2);
        }else if(type=="W"){//inlet
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            uij = u_Domain[id].Q;
            uijm = u_Domain[id-Nx].Q;
            uijp = u_Domain[id+Nx].Q;
            uipj = u_Domain[id+1].Q;
            uippj = u_Domain[id+2].Q;
            Nablau =  uij*(uipj-uij)/(delx) ;
            Laplaceu = (-2*uipj + uij + uippj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely,2);
        }else if(type=="E"){
            double dudx = 0;
            double dvdx = 0;
            uij = u_Domain[id].Q;
            vid = I+ J*(Nx-1);
            uijp = u_Domain[id + Nx].Q;
            uijm = u_Domain[id - Nx].Q;
            uimj = u_Domain[id-1].Q;
            uimmj =u_Domain[id-2].Q;
            Nablau =vij*(uijp-uijm)/(2*dely);
            Laplaceu = (-2*uimj + uij + uimmj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely,2);
        }else if(type=="N"){
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            vid = I+ J*(Nx-1);
            uij =  u_Domain[id].Q;
            if(I==0) {
                vij  =v_Domain[vid].Q/2;
                uijp = 0;
                uijm = (u_Domain[id].Q+u_Domain[id-Nx].Q)/2;
                uipj = u_Domain[id+1].Q;
                uippj = u_Domain[id + 2].Q;
                Nablau =  uij*(uipj-uij)/(delx) +vij*(uijp - uijm)/(dely/2);
                Laplaceu = (-2*uipj + uij + uippj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }else if(I==Nx-1){
                vij  =v_Domain[vid-1].Q/2;
                uijp = 0;
                uijm = (u_Domain[id].Q+u_Domain[id-Nx].Q)/2;
                uimj = u_Domain[id-1].Q;
                uimmj = u_Domain[id - 2].Q;
                Nablau =  uij*(uij-uimj)/(delx) +vij*(uijp - uijm)/(dely/2);
                Laplaceu = (-2*uimj + uij + uimmj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }else{
                vij = (0 + v_Domain[vid - 1].Q + v_Domain[vid].Q)/4;
                uijm = (u_Domain[id].Q+u_Domain[id-Nx].Q)/2;
                uijp = 0;
                uipj = (u_Domain[id+1].Q + u_Domain[id].Q)/2;
                uimj = (u_Domain[id-1].Q + u_Domain[id].Q)/2;
                Nablau =  uij*(uipj-uimj)/(delx) ;
                Laplaceu = (-2*uij + uimj + uipj)/pow(delx/2,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }

        }else if(type == "S"){
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            vid = I+ J*(Nx-1);
            uij =  u_Domain[id].Q;
            if(I==0) {
                vij  =v_Domain[vid+Nx-1].Q/2;
                uijm = 0;
                uijp = (u_Domain[id].Q+u_Domain[id+Nx].Q)/2;
                uipj = u_Domain[id+1].Q;
                uippj = u_Domain[id + 2].Q;
                Nablau =  uij*(uipj-uij)/(delx) +vij*(uijp - uijm)/(dely/2);
                Laplaceu = (-2*uipj + uij + uippj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }else if(I==Nx-1){
                vij  =v_Domain[vid-1].Q/2;
                uijm = 0;
                uijp = (u_Domain[id].Q+u_Domain[id+Nx].Q)/2;
                uimj = u_Domain[id-1].Q;
                uimmj = u_Domain[id - 2].Q;
                Nablau =  uij*(uij-uimj)/(delx) +vij*(uijp - uijm)/(dely/2);
                Laplaceu = (-2*uimj + uij + uimmj)/pow(delx,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }else{
                vij = (0 + v_Domain[vid - 1].Q + v_Domain[vid].Q)/4;
                uijp = (u_Domain[id].Q+u_Domain[id+Nx].Q)/2;
                uijm = 0;
                uipj = (u_Domain[id+1].Q + u_Domain[id].Q)/2;
                uimj = (u_Domain[id-1].Q + u_Domain[id].Q)/2;
                Nablau =  uij*(uipj-uij)/(delx) ;
                Laplaceu = (-2*uij + uimj + uipj)/pow(delx/2,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
            }
        }
        else if(type=="w") {
            double dudx = 0;
            double dvdx = 0;
            uij = 0 ;
            vid = I + J * (Nx - 1);
            uijp = 0;
            uijm = 0;
            uimj = u_Domain[id - 1].Q;
            uimmj = u_Domain[id - 2].Q;
            Nablau = 0;
            Laplaceu = (-2 * uimj + uij + uimmj) / pow(delx, 2) ;
        }else if(type=="e") {
            double dudx = 0;
            double dvdx = 0;
            uij = 0 ;
            vid = I + J * (Nx - 1);
            uijp = 0;
            uijm = 0;
            uipj = u_Domain[id + 1].Q;
            uippj = u_Domain[id + 2].Q;
            Nablau = 0;
            Laplaceu = (-2 * uipj + uij + uippj) / pow(delx, 2) ;
        }else if(type=="s") {
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            vid = I+ J*(Nx-1);
            uij =  u_Domain[id].Q;
            vij = (0 + v_Domain[vid - 1].Q + v_Domain[vid].Q)/4;
            uijm = (u_Domain[id].Q+u_Domain[id-Nx].Q)/2;
            uijp = 0;
            uipj = (u_Domain[id+1].Q + u_Domain[id].Q)/2;
            uimj = (u_Domain[id-1].Q + u_Domain[id].Q)/2;
            Nablau =  uij*(uipj-uimj)/(delx) ;
            Laplaceu = (-2*uij + uimj + uipj)/pow(delx/2,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
        }else if(type=="n") {
            I = u_Domain[id].I;
            J = u_Domain[id].J;
            vid = I+ J*(Nx-1);
            uij =  u_Domain[id].Q;
            vij = (0 + v_Domain[vid - 1].Q + v_Domain[vid].Q)/4;
            uijp = (u_Domain[id].Q+u_Domain[id+Nx].Q)/2;
            uijm = 0;
            uipj = (u_Domain[id+1].Q + u_Domain[id].Q)/2;
            uimj = (u_Domain[id-1].Q + u_Domain[id].Q)/2;
            Nablau =  uij*(uipj-uij)/(delx) ;
            Laplaceu = (-2*uij + uimj + uipj)/pow(delx/2,2) + (uijp-2*uij+uijm)/pow(dely/2,2);
        }else{
            std::cout<<"Claculating Fx:error type unknown"<<std::endl;
        }
        double Fx = -Nablau + mu *Laplaceu;
        return Fx;
    }
    double Fy(std::vector<Points> u_Domain,std::vector<Points> v_Domain,int id,double mu,double delx,double dely,int Nx,std::string type){//vDomain
        double vij,uij,vimj,vipj,vijp,vijm,Nablav,Laplacev,vijmm,vijpp,vippj,vimmj;
        if(type == "P") {
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J*Nx +I;
            vij = v_Domain[id].Q;
            uij =((u_Domain[uid+1].Q + u_Domain[uid - Nx+1].Q)+(u_Domain[uid].Q + u_Domain[uid - Nx].Q))/4;
            vimj =v_Domain[id-1].Q;
            vipj = v_Domain[id+1].Q;
            vijp = v_Domain[id+Nx-1].Q;
            vijm = v_Domain[id-Nx+1].Q;
            Nablav = uij * (vipj - vimj)/(2*delx) + vij * (vijp-vijm)/(2 *dely);
            Laplacev = (vipj -2*vij+vimj)/pow(delx,2) + (vijp - 2*vij+vijm)/pow(dely,2);
        }else if(type == "N"){
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J*Nx +I;
            vij= 0;
            vimj = 0;
            vipj = 0;
            vijm = v_Domain[id-Nx+1].Q;
            vijmm = v_Domain[id-2*(Nx-1)].Q;
            Nablav=0;
            Laplacev =  (vij - 2*vijm+vijmm)/pow(dely,2);
        }else if(type == "S") {
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = 0;
            vimj = 0;
            vipj = 0;
            vijp = v_Domain[id + Nx - 1].Q;
            vijpp = v_Domain[id + 2 * (Nx - 1)].Q;
            Nablav = 0;
            Laplacev = (vij - 2 * vijp + vijpp) / pow(dely, 2);
        }else if(type == "W") {
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = v_Domain[id].Q;
            uij = ((u_Domain[uid + 1].Q + u_Domain[uid - Nx + 1].Q) + (u_Domain[uid].Q + u_Domain[uid - Nx].Q)) / 4;
            vipj = v_Domain[id + 1].Q;
            vippj = v_Domain[id + 2].Q;
            vijp = v_Domain[id + Nx - 1].Q;
            vijm = v_Domain[id - Nx + 1].Q;
            Nablav = uij * (vipj - vij) / (delx) + vij * (vijp - vijm) / (2 * dely);
            Laplacev = (vippj - 2 * vipj + vij) / pow(delx, 2) + (vijp - 2 * vij + vijm) / pow(dely, 2);
        }else if(type == "E"){
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = v_Domain[id].Q;
            uij = ((u_Domain[uid + 1].Q + u_Domain[uid - Nx + 1].Q) + (u_Domain[uid].Q + u_Domain[uid - Nx].Q)) / 4;
            vimj = v_Domain[id - 1].Q;
            vimmj = v_Domain[id - 2].Q;
            vijp = v_Domain[id + Nx - 1].Q;
            vijm = v_Domain[id - Nx + 1].Q;
            Nablav = uij * (vij - vimj) / (delx) + vij * (vijp - vijm) / (2 * dely);
            Laplacev = (vimmj - 2 * vimj + vij) / pow(delx, 2) + (vijp - 2 * vij + vijm) / pow(dely, 2);

    }else if(type =="s"){
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J*Nx +I;
            vij= 0;
            vimj = 0;
            vipj = 0;
            vijm = v_Domain[id-Nx+1].Q;
            vijmm = v_Domain[id-2*(Nx-1)].Q;
            Nablav=0;
            Laplacev =  (vij - 2*vijm+vijmm)/pow(dely,2);
        }else if(type =="n") {
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = 0;
            vimj = 0;
            vipj = 0;
            vijp = v_Domain[id + Nx - 1].Q;
            vijpp = v_Domain[id + 2 * (Nx - 1)].Q;
            Nablav = 0;
            Laplacev = (vij - 2 * vijp + vijpp) / pow(dely, 2);
        }else if(type=="e"){
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = v_Domain[id].Q;
            uij = ((u_Domain[uid + 1].Q + u_Domain[uid - Nx + 1].Q) + (u_Domain[uid].Q + u_Domain[uid - Nx].Q)) / 4;
            vipj = v_Domain[id + 1].Q;
            vippj = v_Domain[id + 2].Q;
            vijp = v_Domain[id + Nx - 1].Q;
            vijm = v_Domain[id - Nx + 1].Q;
            Nablav = uij * (vipj - vij) / (delx) + vij * (vijp - vijm) / (2 * dely);
            Laplacev = (vippj - 2 * vipj + vij) / pow(delx, 2) + (vijp - 2 * vij + vijm) / pow(dely, 2);
        }else if(type=="w"){
            int I = v_Domain[id].I;
            int J = v_Domain[id].J;
            int uid = J * Nx + I;
            vij = v_Domain[id].Q;
            uij = ((u_Domain[uid + 1].Q + u_Domain[uid - Nx + 1].Q) + (u_Domain[uid].Q + u_Domain[uid - Nx].Q)) / 4;
            vimj = v_Domain[id - 1].Q;
            vimmj = v_Domain[id - 2].Q;
            vijp = v_Domain[id + Nx - 1].Q;
            vijm = v_Domain[id - Nx + 1].Q;
            Nablav = uij * (vij - vimj) / (delx) + vij * (vijp - vijm) / (2 * dely);
            Laplacev = (vimmj - 2 * vimj + vij) / pow(delx, 2) + (vijp - 2 * vij + vijm) / pow(dely, 2);
        }else{
            std::cout<<"Calculate Fy:Error:type not found"<<std::endl;
        }
        double Fy = -Nablav + mu *Laplacev;
//        std::cout<<"Fy   "<<id<<"   "<<Fy<<std::endl;
        return Fy;
    }
    void updateP(int P_points);


    void Post_staggered(std::string type,int time,std::vector<Points> Q_Domain ,int Nx,int Ny)
    {
        std::string filename = type + std::to_string(time) + ".tec";
        std::ofstream outfile;
        outfile.open(filename, std::ios::out);
        outfile << std::setprecision(14);
        outfile.setf(std::ios::scientific);
        outfile << "VARIABLES =\"x\" \n"
                << "\"y\" \n"
                << "\"Q\" \n";

        outfile<<"ZONE T = \"Rank" << 1<<"\" \n"
               <<"I= "<<Nx<<"\n"
               <<"J= "<<Ny<<"\n";
        outfile<<"DATAPACKING = POINT\n";

        for (auto i :Q_Domain) {
                outfile << i.x << " " << i.y << " " << i.Q<< "\n";
        };
        outfile.unsetf(std::ios::scientific);
        outfile << std::setprecision(6);
        outfile.close();
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



};