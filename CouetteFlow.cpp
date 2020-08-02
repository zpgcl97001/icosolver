#include "CouetteFlow.h"
#include "icoSolver.h"
#include <math.h>
#include "Matsolver.h"

//Create staggared Grid

//Initialize
int main() {

    double mu = 0.01;//Pa?¡ès
    double rho = 1; //kg/m3
    double Re = 200;
    double U_init;
    U_init = Re * mu / (rho);
//    U_init =0.01;
    double atm = 0;
    int Totalstep = 10000;
    double delt = 0.001;

    int Nx = 101;
    int Ny = 31;
    double delx = 0.1;
    double dely = 0.1;

    int Gx = Nx - 1;
    int Gy = Ny - 1;

    int P_points = (Nx - 1) * (Ny - 1);
    int u_points = Nx * (Ny - 1);
    int v_points = (Nx - 1) * (Ny);
    int Total_points = Nx * Ny;
    int NL = 11;

    std::vector<Points> P_Domain;
    std::vector<Points> u_Domain;
    std::vector<Points> v_Domain;
    std::vector<Points> Domain;
    Bluntbody Square;

    for (int i = 0; i < u_points; ++i) {
        u_Domain.push_back(Points());
        u_Domain[i].I = i % (Nx);
        u_Domain[i].J = i / (Nx);
        u_Domain[i].x = u_Domain[i].I * delx;
        u_Domain[i].y = u_Domain[i].J * dely + 0.5 * dely;
    }

    for (int i = 0; i < v_points; ++i) {
        v_Domain.push_back(Points());
        v_Domain[i].I = i % (Nx - 1);
        v_Domain[i].J = i / (Nx - 1);
        v_Domain[i].x = v_Domain[i].I * delx + 0.5 * delx;
        v_Domain[i].y = v_Domain[i].J * dely;
    }

    for (int i = 0; i < P_points; ++i) {
        P_Domain.push_back(Points());
        P_Domain[i].I = i % (Nx - 1);
        P_Domain[i].J = i / (Nx - 1);
        P_Domain[i].x = P_Domain[i].I * delx + 0.5 * delx;
        P_Domain[i].y = P_Domain[i].J * dely + 0.5 * dely;
        P_Domain[i].Q = atm;
//        std::cout<<i<<"    P"<<P_Domain[i].I<<","<<P_Domain[i].J<<std::endl;
    }

    for (int i = 0; i < Total_points; ++i) {
        Domain.push_back(Points());
        Domain[i].I = i % (Nx);
        Domain[i].J = i / (Nx);
        Domain[i].x = Domain[i].I * delx;
        Domain[i].y = Domain[i].J * dely;
    }

    //Store bluntbody id
    for (int i = 0; i <P_points ; ++i) {
        if(P_Domain[i].I==3*(NL-1)-1 && P_Domain[i].J>NL-2 && P_Domain[i].J< 2*NL-2){
            Square.PW_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<P_Domain[i].x<<","<<P_Domain[i].y<<std::endl;
        }else if(P_Domain[i].I==4*(NL-1) && P_Domain[i].J>NL-2 && P_Domain[i].J<2*NL-2){
            Square.PE_point.push_back(i);
//            std::cout<<"PE      "<<"I,J="<<P_Domain[i].x<<","<<P_Domain[i].y<<std::endl;
        }else if(P_Domain[i].J==NL-2 && P_Domain[i].I>3*(NL-1)-1 && P_Domain[i].I<4*(NL-1)) {
            Square.PS_point.push_back(i);
//            std::cout<<"PS      "<<"I,J="<<P_Domain[i].x<<","<<P_Domain[i].y<<std::endl;
        }else if(P_Domain[i].J==2*NL-2 && P_Domain[i].I>3*(NL-1)-1 && P_Domain[i].I<4*(NL-1)){
            Square.PN_point.push_back(i);
//            std::cout<<"PN      "<<"I,J="<<P_Domain[i].x<<","<<P_Domain[i].y<<std::endl;
        }else if(P_Domain[i].x>(3*(NL-1)-1)*delx && P_Domain[i].x< 4*(NL-1)*dely && P_Domain[i].y>(NL-1)*dely && P_Domain[i].y< 2*(NL-1)*dely){
            Square.PInner_point.push_back(i);
        }
    }

    for (int i = 0; i <u_points ; ++i) {
        if(abs(u_Domain[i].x - 3*(NL-1)*delx)<1e-5 && u_Domain[i].y>(NL-1)*dely && u_Domain[i].y< 2*(NL-1)*dely){
            Square.uW_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<u_Domain[i].x<<","<<u_Domain[i].y<<std::endl;
        }else if(abs(u_Domain[i].x - 4* (NL-1)*delx)<1e-5 && u_Domain[i].y>(NL-1)*dely && u_Domain[i].y< 2*(NL-1)*dely) {
            Square.uE_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<u_Domain[i].x<<","<<u_Domain[i].y<<std::endl;
        }else if(abs(u_Domain[i].y - 2* (NL-1)*dely-0.5*dely)<1e-5 && u_Domain[i].x>3*(NL-1)*delx && u_Domain[i].x< 4*(NL-1)*delx){
            Square.uN_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<u_Domain[i].x<<","<<u_Domain[i].y<<std::endl;
        }else if(abs(u_Domain[i].y -  (NL-2)*dely-0.5*dely)<1e-5 && u_Domain[i].x>3*(NL-1)*delx && u_Domain[i].x< 4*(NL-1)*delx){
            Square.uS_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<u_Domain[i].x<<","<<u_Domain[i].y<<std::endl;
        }else if(u_Domain[i].x>3*(NL-1)*delx && u_Domain[i].x< 4*(NL-1)*delx && u_Domain[i].y>(NL-1)*dely && u_Domain[i].y< 2*(NL-1)*dely){
            Square.uInner_point.push_back(i);
        }
    }

    for (int i = 0; i <v_points ; ++i) {
        if(abs(v_Domain[i].y - (NL-1)*dely)<1e-5 && v_Domain[i].x>3*(NL-1)*delx && v_Domain[i].x< 4*(NL-1)*delx){
            Square.vS_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<v_Domain[i].x<<","<<v_Domain[i].y<<std::endl;
        }else if(abs(v_Domain[i].y - 2*(NL-1)*dely)<1e-5 && v_Domain[i].x>3*(NL-1)*delx && v_Domain[i].x< 4*(NL-1)*delx){
            Square.vN_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<v_Domain[i].x<<","<<v_Domain[i].y<<std::endl;
        }else if(abs(v_Domain[i].x - 4* (NL-1)*delx-0.5*delx)<1e-5 && v_Domain[i].y>(NL-1)*delx && v_Domain[i].y< 2*(NL-1)*delx){
            Square.vE_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<v_Domain[i].x<<","<<v_Domain[i].y<<std::endl;
        }else if(abs(v_Domain[i].x - 3*(NL-1)*delx+0.5*delx)<1e-5 && v_Domain[i].y>(NL-1)*dely && v_Domain[i].y< 2*(NL-1)*dely){
            Square.vW_point.push_back(i);
//            std::cout<<"PW      "<<"I,J="<<v_Domain[i].x<<","<<v_Domain[i].y<<std::endl;
        }else if(v_Domain[i].x>3*(NL-1)*delx && v_Domain[i].x< 4*(NL-1)*delx && v_Domain[i].y>(NL-1)*dely && v_Domain[i].y< 2*(NL-1)*dely){
            Square.vInner_point.push_back(i);
        }
    }

//    for(auto i:Square.uE_point){
//        std::cout<<u_Domain[i].x<<"--"<<u_Domain[i].y<<std::endl;
//    }



    MatSol matsol;
//Initialize BC points


//initialize field

    for (auto point:P_Domain) {
        point.Q = atm;
    }


    for (int i=0;i<u_points;++i){
        if(u_Domain[i].I == 0) {
            u_Domain[i].Q = U_init;
        }else{
            u_Domain[i].Q = 0;
        }
    }

    for (auto point:v_Domain) {
        point.Q = 0;
    }

    std::vector<int> BC_S_u;
    std::vector<int> BC_N_u;


//Initialize F

    IcoSolver method;



    //Timeloop

    for(int step = 0;step<Totalstep;step++){
        std::cout<<"Step================================================================"<<step<<std::endl;
        std::cout<<"===================================================================="<<step<<std::endl;
        std::cout<<"===================================================================="<<step<<std::endl;
        //Revised No-slip condition
        for (auto p : v_Domain) {
            int I = p.I;
            int J = p.J;
            if(J==0 || J==Ny-1){
                p.Q = 0;
            }
        }

        for (auto i : Square.vS_point) {
            v_Domain[i].Q =  0 ;
        }
        for (auto i : Square.vN_point) {
            v_Domain[i].Q =  0 ;
        }
        for (auto i : Square.uW_point) {
            u_Domain[i].Q =  0 ;
        }
        for (auto i : Square.uE_point) {
            u_Domain[i].Q =  0 ;
        }
        //Calculate Fx,Fy
        for (int i = 0; i < u_points; ++i) {
            int I = u_Domain[i].I;
            int J = u_Domain[i].J;

            if(method.InVector(i,Square.uInner_point)){
                u_Domain[i].F = 0;
            }else if(method.InVector(i,Square.uW_point)){
                std::string type = "w";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.uE_point)){
                std::string type = "e";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.uS_point)){
                std::string type = "s";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.uN_point)){
                std::string type = "n";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }
            else if(J==0){
                std::string type = "S";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(J==Ny-2) {
                std::string type = "N";
                u_Domain[i].F = method.Fx(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else if(I==0){
                std::string type = "W";
                u_Domain[i].F = method.Fx(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else if(I==Nx-1){
                std::string type = "E";
                u_Domain[i].F = method.Fx(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else{
                std::string type ="P";
                u_Domain[i].F = method.Fx(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);

            }
            u_Domain[i].QPlus = u_Domain[i].Q + delt*u_Domain[i].F;
        }
        std::cout<<"Updated Fx"<<std::endl;

        for (int i = 0; i < v_points; ++i) {
            int I = v_Domain[i].I;
            int J = v_Domain[i].J;
            if(method.InVector(i,Square.vInner_point)){
                v_Domain[i].F = 0;
            }else if(method.InVector(i,Square.vW_point)){
                std::string type = "w";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.vE_point)){
                std::string type = "e";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.vS_point)){
                std::string type = "s";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(method.InVector(i,Square.vN_point)){
                std::string type = "n";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if (J==0){
                std::string type = "S";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }else if(J==Ny-1) {
                std::string type = "N";
                v_Domain[i].F = method.Fy(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else if(I==0){
                std::string type = "W";
                v_Domain[i].F = method.Fy(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else if(I==Nx-2){
                std::string type = "E";
                v_Domain[i].F = method.Fy(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
            }else{
                std::string type ="P";
                v_Domain[i].F = method.Fy(u_Domain,v_Domain,i,mu,delx,dely,Nx,type);
            }
            v_Domain[i].QPlus = v_Domain[i].Q + delt*v_Domain[i].F;
        }
        std::cout<<"Updated Fy"<<std::endl;


        //update P.rho * nablaF
        for (int i = 0; i < P_points; ++i) {
            int I = P_Domain[i].I;
            int J = P_Domain[i].J;
            int uid = J*Nx +I;
            int vid = J*(Nx-1)+I;
            P_Domain[i].F = rho/delt * ((-u_Domain[uid].QPlus + u_Domain[uid+1].QPlus)/delx + (-v_Domain[vid].QPlus + v_Domain[vid+Nx-1].QPlus)/dely);

        }

        std::cout<<"Update P.nablaF"<<std::endl;

        //UpdateP
        double TDMA_error = 1e10;
        while(TDMA_error>1e-10) {
            for (int i = 0; i < Ny - 1; ++i) {
                std::vector<double> TDMAa;
                std::vector<double> TDMAb;
                std::vector<double> TDMAc;
                std::vector<double> Matr;
                if (i == 0) {

                    for (int j = 0; j < Nx - 1; ++j) {
                        int Pid = j + (Nx - 1) * i;
                        if (j == 0) {
                            TDMAa.push_back(0);
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(2 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q -rho* dely * v_Domain[Pid].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid + Nx - 1].Q + Pghost) / pow(dely, 2));
                        } else if (j < Nx - 2) {
                            TDMAa.push_back(1 / pow(delx, 2));
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(1 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q -rho* dely * v_Domain[Pid].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid + Nx - 1].Q + Pghost) / pow(dely, 2));
                        } else {
                            TDMAa.push_back(0);
                            TDMAb.push_back(1);
                            TDMAc.push_back(0);
                            Matr.push_back(atm);
                        }
                    }
                }
                else if (i < Ny - 2) {
                    for (int j = 0; j < Nx - 1; ++j) {
                        int Pid = j + (Nx - 1) * i;
                        if (method.InVector(i*(Nx-1)+j,Square.PInner_point)){
                            TDMAa.push_back(0);
                            TDMAb.push_back(1);
                            TDMAc.push_back(0);
                            Matr.push_back(0);
                        }else if (method.InVector(i*(Nx-1)+j,Square.PS_point)){
                            TDMAa.push_back(1 / pow(delx, 2));
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(1 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q + dely * v_Domain[Pid + Nx - 1].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid - Nx + 1].Q + Pghost) / pow(dely, 2));
                        }else if (method.InVector(i*(Nx-1)+j,Square.PN_point)) {
                            TDMAa.push_back(1 / pow(delx, 2));
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(1 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q - dely * v_Domain[Pid].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid + Nx - 1].Q + Pghost) / pow(dely, 2));
                        }else if(method.InVector(i*(Nx-1)+j,Square.PW_point)) {
                            TDMAa.push_back(-1/delx);
                            TDMAb.push_back(1/delx);
                            TDMAc.push_back(0);
                            int uid = j+Nx*i;
                            Matr.push_back(u_Domain[uid].F);
                        }else if(method.InVector(i*(Nx-1)+j,Square.PE_point)) {
                            TDMAa.push_back(0);
                            TDMAb.push_back(-1/delx);
                            TDMAc.push_back(1/delx);
                            int uid = j+Nx*i;
                            Matr.push_back(u_Domain[uid+1].F);
                        }else if(j == 0) {
                            TDMAa.push_back(0);
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(2 / pow(delx, 2));
                            Matr.push_back(P_Domain[Pid].F -
                                           (P_Domain[Pid + Nx - 1].Q + P_Domain[Pid - Nx + 1].Q) / pow(dely, 2));
                        } else if (j < Nx - 2) {
                            TDMAa.push_back(1 / pow(delx, 2));
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(1 / pow(delx, 2));
                            Matr.push_back(
                                    P_Domain[Pid].F -
                                    (P_Domain[Pid + Nx - 1].Q + P_Domain[Pid - Nx + 1].Q) / pow(dely, 2));
                        } else {
                            TDMAa.push_back(0);
                            TDMAb.push_back(1);
                            TDMAc.push_back(0);
                            Matr.push_back(atm);
                        }
                    }

                }
                else {
                    for (int j = 0; j < Nx - 1; ++j) {
                        int Pid = j + (Nx - 1) * i;
                        if (j == 0) {
                            TDMAa.push_back(0);
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(2 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q + rho*dely * v_Domain[Pid + Nx - 1].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid - Nx + 1].Q + Pghost) / pow(dely, 2));
                        } else if (j < Nx - 2) {
                            TDMAa.push_back(1 / pow(delx, 2));
                            TDMAb.push_back(-2 / pow(delx, 2) - 2 / pow(dely, 2));
                            TDMAc.push_back(1 / pow(delx, 2));
                            double Pghost = P_Domain[Pid].Q +rho* dely * v_Domain[Pid + Nx - 1].F;
                            Matr.push_back(P_Domain[Pid].F - (P_Domain[Pid - Nx + 1].Q + Pghost) / pow(dely, 2));
                        } else {
                            TDMAa.push_back(0);
                            TDMAb.push_back(1);
                            TDMAc.push_back(0);
                            Matr.push_back(atm);
                        }
                    }
                }
                //Solve TDMA
                matsol.SolveTDMA(TDMAa, TDMAb, TDMAc, Matr, Nx);
                //Store Last step
                for (int j = 0; j < Nx - 1; ++j) {
                    int Pid = j + i * (Nx - 1);
                    P_Domain[Pid].P = P_Domain[Pid].Q;
                }
                //Update P
                for (int j = 0; j < Nx - 1; ++j) {
                    int Pid = j + i * (Nx - 1);
                    P_Domain[Pid].Q = Matr[j];
                }

            }

            //Calculate error
            TDMA_error = 0;
            for (auto point : P_Domain) {
                TDMA_error += pow(point.P - point.Q, 2);

            }
            TDMA_error = TDMA_error;
        }
        std::cout<<"Calculate P finish"<<std::endl;



        //Update u
        for (int i = 0; i < u_points; ++i) {
            int I = u_Domain[i].I;
            int J = u_Domain[i].J;
            if(I==0){
                int pid = J*(Nx-1)+I;//  P--U the Pid is left Uid
                u_Domain[i].Q = U_init;
                std::string type = "W";
                u_Domain[i].F = method.Fx(u_Domain, v_Domain, i, mu, delx, dely, Nx, type);
//                std::cout<<u_Domain[i].F<<"--"<<i<<std::endl;
                double dpdx =rho * u_Domain[i].F;
                u_Domain[i].P = dpdx;
                u_Domain[i].Q = u_Domain[i].QPlus - delt/rho * dpdx;
            }else if(I==Nx-1){
                int pid = J*(Nx-1)+I;//  P--U the Pid is left Uid
                double dpdx =u_Domain[i-1].P;
                u_Domain[i].Q = u_Domain[i].QPlus - delt/rho * dpdx;
                u_Domain[i].P = dpdx;
            }else{
                int pid = J*(Nx-1)+I;//  P--U the Pid is left Uid
                double dpdx = (P_Domain[pid].Q-P_Domain[pid-1].Q)/(delx);
                u_Domain[i].Q = u_Domain[i].QPlus - delt/rho * dpdx;
                u_Domain[i].P = dpdx;
            }


        }
        std::cout<<"updated u"<<std::endl;
        //update v
        for (int i = 0; i <v_points ; ++i) {
            int I = P_Domain[i].I;
            int J = P_Domain[i].J;
            if(J==0){
                v_Domain[i].Q = 0;
//
            }else if(J==Ny-1){
                v_Domain[i].Q = 0;
//
            }
            else{
                int pid = J*(Nx-1) + I;//P on the top y p///y
                double dpdy = (P_Domain[pid].Q - P_Domain[pid-Nx+1].Q)/(dely);
                v_Domain[i].P =dpdy;
                v_Domain[i].Q = v_Domain[i].QPlus - delt/rho * dpdy;
            }

        }
        std::cout<<"Updated v"<<std::endl;
        std::cout<<"===========================Pfield========================="<<std::endl;
        for(auto i:P_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.Q<<"\t";

        }
        std::cout<<std::endl;
        std::cout<<"===========================Nablavstarfield========================="<<std::endl;
        for(auto i:P_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.F<<"\t";

        }
        std::cout<<std::endl;


        //OutputP
        std::cout<<"==========================Ufield============================"<<std::endl;
        for(auto i:u_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.Q<<"\t";

        }

        std::cout<<std::endl;
        std::cout<<"==========================Vfield============================"<<std::endl;
        for(auto i:v_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.Q<<"\t";

        }

        std::cout<<std::endl;

        std::cout<<"==========================Ustarfield============================"<<std::endl;
        for(auto i:u_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.QPlus<<"\t";

        }

        std::cout<<std::endl;
        std::cout<<"==========================Vstarfield============================"<<std::endl;
        for(auto i:v_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.QPlus<<"\t";

        }

        std::cout<<std::endl;

        std::cout<<"==========================Fxfield============================"<<std::endl;
        for(auto i:u_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.F<<"\t";

        }

        std::cout<<std::endl;
        std::cout<<"==========================Fyfield============================"<<std::endl;
        for(auto i:v_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.F<<"\t";

        }

        std::cout<<std::endl;
        std::cout<<"==========================dpdyfield============================"<<std::endl;
        for(auto i:v_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.P<<"\t";

        }

        std::cout<<std::endl;
        std::cout<<"==========================dpdxfield============================"<<std::endl;
        for(auto i:u_Domain){
            int I=i.I;
            int J=i.J;
            if(I==0){
                std::cout<<std::endl;
            }
            std::cout<<i.P<<"\t";

        }

        std::cout<<std::endl;


        if(step % 10 == 0){
            //Combine Staggered Grid



            for (auto p : v_Domain) {
                int I = p.I;
                int J = p.J;
                if(J==0 || J==Ny-1){
                    p.Q = 0;
                }
            }

            for (auto i : Square.vS_point) {
                v_Domain[i].Q =  0 ;
            }
            for (auto i : Square.vN_point) {
                v_Domain[i].Q =  0 ;
            }
            for (auto i : Square.uW_point) {
                u_Domain[i].Q =  0 ;
            }
            for (auto i : Square.uE_point) {
                u_Domain[i].Q =  0 ;
            }


            for (int i=0; i<Total_points;++i) {
                int I = Domain[i].I;
                int J = Domain[i].J;
                int x = Domain[i].x;
                int y = Domain[i].y;
//
                if (I == 0 || J == 0 || I == Nx - 1 || J == Ny - 1) {
                    Domain[i].u = 0;
                    Domain[i].v = 0;
                } else {
                    int uid = J * Nx + I;
                    int vid = J * (Nx - 1) + I;
                    Domain[i].u = 0.5 * (u_Domain[uid].Q + u_Domain[uid - Nx].Q);
                    Domain[i].v = 0.5 * (v_Domain[vid].Q + v_Domain[vid - 1].Q);
//                std::cout<<Domain[i].u<<"--"<<Domain[i].v<<std::endl;
                }
            }



            method.Post_staggered("P",step,P_Domain,Nx-1,Ny-1);
            method.Post_staggered("U",step,u_Domain,Nx,Ny-1);
            method.Post_staggered("V",step,v_Domain,Nx-1,Ny);
            method.Post_Domain("Domain",step,Domain,Nx,Ny);
        }



    }

    //Postproess

    return 0;
}