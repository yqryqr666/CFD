#include <float.h>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <filesystem>
#include <string>
namespace fs = std::filesystem;

double Gamma = 1.4;
double CFL = 0.8;
double x0 = 0, x1 = 1, y0_new = 0, y1_new = 1;
double T_end = 0.25;
int nxCells = 300; int nyCells = 300;
//int nxPoints = nxCells + 1, nyPoints = nyCells + 1;
int ghost_cells = 2;
//double dx = (x1 - x0)/(nxPoints - 2), dy = (y1_new - y0_new)/(nyPoints - 2);
double dx = (x1 - x0)/nxCells, dy = (y1_new - y0_new)/nyCells;
double dt = 0;
//std::vector<double> v_rigid = {-1, 0};
std::vector<double> v_rigid = {0, 0};
double inf = std::numeric_limits<double>::infinity();


struct pvar {
    double rho, u, v, p;

    //  const 
    double& operator[](std::size_t i) {
        switch(i) {
            case 0: return rho;
            case 1: return u;
            case 2: return v;
            case 3: return p;
            default:
                throw std::out_of_range("Index out of range");
        }
    }

    // const 
    const double& operator[](std::size_t i) const {
        switch(i) {
            case 0: return rho;
            case 1: return u;
            case 2: return v;
            case 3: return p;
            default:
                throw std::out_of_range("Index out of range");
        }
    }
};
struct cvar {
    double rho, rhou, rhov, E;

    // 非 const 版本
    double& operator[](std::size_t i) {
        switch(i) {
            case 0: return rho;
            case 1: return rhou;
            case 2: return rhov;
            case 3: return E;
            default:
                throw std::out_of_range("Index out of range");
        }
    }

    // const 版本
    const double& operator[](std::size_t i) const {
        switch(i) {
            case 0: return rho;
            case 1: return rhou;
            case 2: return rhov;
            case 3: return E;
            default:
                throw std::out_of_range("Index out of range");
        }
    }
};


enum ShapeType {
    CIRCLE,
    SQUARE,
    TWO_CIRCLES,
    TWO_OVERLAPPING_CIRCLES,
    STATIC
    //MOVING
};



pvar con2pri(cvar cons, double Gamma){
    double rho = cons.rho; double rhou = cons.rhou; double rhov = cons.rhov; double E = cons.E;
    pvar pri;
    pri.rho = rho;
    pri.u = rhou / rho;
    pri.v = rhov / rho;
    pri.p = (Gamma - 1) * (E - 0.5 * (rhov * rhov + rhou * rhou)/ rho);
    return pri;
}
cvar pri2con(pvar prim, double Gamma){

    double rho = prim.rho; double u = prim.u; double v = prim.v; double p = prim.p;

    cvar con;
    con.rho = rho;
    con.rhou = rho * u;
    con.rhov = rho * v;
    con.E = p / (Gamma - 1) + 0.5 * rho * (v * v + u * u); //E = p/(Gamma-1) + 0.5 rho v^2
    return con;
}
std::vector<double> f(pvar prim){ //x-direction flux
    double rho = prim.rho; double u = prim.u;
    double v = prim.v; double p = prim.p;
    cvar cons = pri2con(prim,Gamma);
    double E = cons.E;
    std::vector<double> flux(4);
    flux[0] = rho * u;
    flux[1] = rho * u * u + p;
    flux[2] = rho * u * v;
    flux[3] = (E + p) * u;
    return flux;
}
std::vector<double> g(pvar prim){ // y-direction flux
    //std::vector<double> g = {rho * v, rho * u * v, rho * v * v + p, (E + p) * v};
    double rho = prim.rho; double u = prim.u;
    double v = prim.v; double p = prim.p;
    cvar cons = pri2con(prim,Gamma);
    double E = cons.E;
    std::vector<double> glux(4);
    glux[0] = rho * v;
    glux[1] = rho * u * v;
    glux[2] = rho * v * v + p;
    glux[3] = (E + p) * v;
    return glux;
}
double Minbee(double r){
    if(r <= 0) return 0;
    if(r <= 1) return r;
    return std::min(1.0, 2.0 / (1 + r));
}
void data_reconstruction(double uL, double u, double uR, double &uL_bar, double &uR_bar){
    double delta0, delta1, delta, r;
    delta0 = u - uL;//backward
    delta1 = uR- u;//forward
    delta = 0.5 * (delta0 + delta1);
    if (delta1 == 0){
        r = 9999;
    }else{
        r = delta0 / delta1;
    }
    double limiter;
    limiter = Minbee(r);
    uL_bar = u - 0.5 * limiter * delta;
    uR_bar = u + 0.5 * limiter * delta;
}
void half_time(cvar &uL_half, cvar &uR_half, cvar &uL, cvar &uR, double dt, double dx, double dy, double Gamma, char dir){
    pvar priL = con2pri(uL, Gamma);
    pvar priR = con2pri(uR, Gamma);
    if(dir == 'x'){
        std::vector<double> FL = f(priL);
        std::vector<double> FR = f(priR);

        uL_half.rho = uL.rho - 0.5 * (dt/dx) * (FR[0] - FL[0]);
        uL_half.rhou = uL.rhou - 0.5 * (dt/dx) * (FR[1] - FL[1]);
        uL_half.rhov = uL.rhov - 0.5 * (dt/dx) * (FR[2] - FL[2]);
        uL_half.E = uL.E - 0.5 * (dt/dx) * (FR[3] - FL[3]);

        uR_half.rho = uR.rho - 0.5 * (dt/dx) * (FR[0] - FL[0]);
        uR_half.rhou = uR.rhou- 0.5 * (dt/dx) * (FR[1] - FL[1]);
        uR_half.rhov = uR.rhov- 0.5 * (dt/dx) * (FR[2] - FL[2]);
        uR_half.E = uR.E - 0.5 * (dt/dx) * (FR[3] - FL[3]);

    }else if(dir == 'y'){
        std::vector<double> FL = g(priL);
        std::vector<double> FR = g(priR);

        uL_half.rho = uL.rho - 0.5 * (dt/dy) * (FR[0] - FL[0]);
        uL_half.rhou = uL.rhou - 0.5 * (dt/dy) * (FR[1] - FL[1]);
        uL_half.rhov = uL.rhov - 0.5 * (dt/dy) * (FR[2] - FL[2]);
        uL_half.E = uL.E - 0.5 * (dt/dy) * (FR[3] - FL[3]);

        uR_half.rho = uR.rho - 0.5 * (dt/dy) * (FR[0] - FL[0]);
        uR_half.rhou = uR.rhou- 0.5 * (dt/dy) * (FR[1] - FL[1]);
        uR_half.rhov = uR.rhov- 0.5 * (dt/dy) * (FR[2] - FL[2]);
        uR_half.E = uR.E - 0.5 * (dt/dy) * (FR[3] - FL[3]);
    }
}
std::vector<double> fri_flux(const cvar& uL, const cvar& uR, double dt, double dx, double dy, char dir){
    std::vector<double> flux_f(4);
    if(dir == 'x'){
        //left&right flux
        pvar uL_pri = con2pri(uL, Gamma);
        pvar uR_pri = con2pri(uR, Gamma);
        std::vector<double> FL = f(uL_pri);
        std::vector<double> FR = f(uR_pri);
        flux_f[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (dx / dt) * (uR.rho - uL.rho);
        flux_f[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (dx / dt) * (uR.rhou - uL.rhou);
        flux_f[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (dx / dt) * (uR.rhov - uL.rhov);
        flux_f[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (dx / dt) * (uR.E - uL.E);
    }else if(dir == 'y'){
        //left&right flux
        pvar uL_pri = con2pri(uL, Gamma);
        pvar uR_pri = con2pri(uR, Gamma);
        std::vector<double> FL = g(uL_pri);
        std::vector<double> FR = g(uR_pri);
        flux_f[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (dy / dt) * (uR.rho - uL.rho);
        flux_f[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (dy / dt) * (uR.rhou - uL.rhou);
        flux_f[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (dy / dt) * (uR.rhov - uL.rhov);
        flux_f[3] = 0.5 * (FL[3] + FR[3]) - 0.5 * (dy / dt) * (uR.E - uL.E);
   }
   return flux_f;
}
std::vector<double> RI_flux(const cvar& uL, const cvar& uR,  double dt, double dx, double dy, char dir){
    std::vector<double> flux_r;
    if(dir == 'x'){
        pvar uL_pri = con2pri(uL, Gamma);
        pvar uR_pri = con2pri(uR, Gamma);
        std::vector<double> FL = f(uL_pri);
        std::vector<double> FR = f(uR_pri);
        //intermediate state u*
        cvar u_star;
        u_star.rho = 0.5 * (uL.rho + uR.rho) - 0.5 * (dt/dx) * (FR[0] - FL[0]);
        u_star.rhou = 0.5 * (uL.rhou + uR.rhou) - 0.5 * (dt/dx) * (FR[1] - FL[1]);
        u_star.rhov = 0.5 * (uL.rhov + uR.rhov) - 0.5 * (dt/dx) * (FR[2] - FL[2]);
        u_star.E = 0.5 * (uL.E + uR.E) - 0.5 * (dt/dx) * (FR[3] - FL[3]);
        //u* flux
        pvar pri_star = con2pri(u_star, Gamma);
        flux_r = f(pri_star);
    }else if(dir == 'y'){
        pvar uL_pri = con2pri(uL, Gamma);
        pvar uR_pri = con2pri(uR, Gamma);
        //std::vector<double> flux(4);
        std::vector<double> FL = g(uL_pri);
        std::vector<double> FR = g(uR_pri);
        //intermediate state u*
        cvar u_star;
        u_star.rho = 0.5 * (uL.rho + uR.rho) - 0.5 * (dt/dy) * (FR[0] - FL[0]);
        u_star.rhou = 0.5 * (uL.rhou + uR.rhou) - 0.5 * (dt/dy) * (FR[1] - FL[1]);
        u_star.rhov = 0.5 * (uL.rhov + uR.rhov) - 0.5 * (dt/dy) * (FR[2] - FL[2]);
        u_star.E = 0.5 * (uL.E + uR.E) - 0.5 * (dt/dy) * (FR[3] - FL[3]);
        //u* flux
        pvar pri_star = con2pri(u_star, Gamma);
        flux_r = g(pri_star);
    }
    return flux_r;
}
std::vector<double> FORCE(const cvar& uL, const cvar& uR,  double dt, double dx, double dy, char dir){
    std::vector<double> flux(4);

    std::vector<double> flux_f = fri_flux(uL, uR, dt, dx, dy, dir);
    std::vector<double> flux_r = RI_flux(uL, uR, dt, dx, dy, dir);

    for(int i = 0; i < 4; i++){
        flux[i] = 0.5 * (flux_f[i] + flux_r[i]);
    }
    return flux;
}

std::vector<double> con2pri_hllc(std::vector<double> cons){
    std::vector<double> prim(3);
    prim[0] = cons[0];
    prim[1] = cons[1]/cons[0];
    prim[2] = (Gamma-1)*(cons[2]-((cons[1]*cons[1])/(2*cons[0])));
    return prim;
}
std::vector<double> pri2con_hllc(std::vector<double> prim){
    std::vector<double> cons(3);
    cons[0] = prim[0];
    cons[1] = prim[0]*prim[1];
    cons[2] = prim[2]/(Gamma-1)+(prim[0]*prim[1])*(prim[0]*prim[1])/(2*prim[0]);
    return cons;
}

double calcs(double p, double rho, double Gamma){
    return sqrt(Gamma * p / rho);
}
//L&R primitive pri
std::vector<double> HLLC(std::vector<double> L, std::vector<double> R){
    std::vector<double> L_cons = pri2con_hllc(L);
    std::vector<double> R_cons = pri2con_hllc(R);
    double cs_l = calcs(L[2], L[0], Gamma);
    double cs_R = calcs(R[2], R[0], Gamma);
    double S_R = std::max(fabs(L[1]) + cs_l, fabs(R[1]) + cs_R);
    double S_L = - S_R;
    double S_star = (R[2] - L[2] + L[0] * L[1] * (S_L - L[1]) - R[0] * R[1] * (S_R - R[1]))/(L[0] * (S_L - L[1]) - R[0] * (S_R - R[1]));

    double temp = L[0] * (S_L - L[1])/(S_L - S_star);
    std::vector<double> HLLC_L_cons(3);
    HLLC_L_cons[0] = temp;
    HLLC_L_cons[1] = temp*S_star;
    HLLC_L_cons[2] = temp * (L_cons[2]/L[0] + (S_star - L[1]) * (S_star + L[2]/(L[0] * (S_L - L[1]))));
    std::vector<double> HLLC_L = con2pri_hllc(HLLC_L_cons);
    return HLLC_L;

}
//L&R conservative, HLLC flux
std::vector<double> HLLC_flux(cvar L, cvar R, int dir){
    pvar L_pri = con2pri(L, Gamma);
    pvar R_pri = con2pri(R, Gamma);

    double cs_l = calcs(L_pri[3], L_pri[0], Gamma);
    double cs_R = calcs(R_pri[3], R_pri[0], Gamma);
    double S_R = std::max(fabs(L_pri[dir]) + cs_l, fabs(R_pri[dir]) + cs_R);
    double S_L = - S_R;
    double S_star = (R_pri[3] - L_pri[3] + L_pri[0] * L_pri[dir] * (S_L - L_pri[dir]) - R_pri[0] * R_pri[dir] * (S_R - R_pri[dir]))/(L_pri[0] * (S_L - L_pri[dir]) - R_pri[0] * (S_R - R_pri[dir]));
    if(S_L > 0){
        if(dir == 1){
            return f(L_pri);
        }else if(dir == 2){
            return g(L_pri);
        }
    }else if(S_star > 0){
        std::vector<double> f_L(4);
        std::vector<double> HLLC_f(4);
        if(dir == 1){
            f_L = f(L_pri);
            double temp = L_pri[0] * (S_L - L_pri[dir])/(S_L - S_star);
            std::vector<double> HLLC_L(4);
            HLLC_L[0] = temp;
            HLLC_L[1] = temp * S_star;
            HLLC_L[2] = temp * L_pri[2];
            HLLC_L[3] = temp * (L[3]/L_pri[0] + (S_star - L_pri[dir]) * (S_star + L_pri[3]/(L_pri[0] * (S_L - L_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_L[k] + S_L * (HLLC_L[k] - L[k]);
            }
        }else if(dir == 2){
            f_L = g(L_pri);
            double temp = L_pri[0] * (S_L - L_pri[dir])/(S_L - S_star);
            std::vector<double> HLLC_L(4);
            HLLC_L[0] = temp;
            HLLC_L[1] = temp * L_pri[1];
            HLLC_L[2] = temp * S_star;
            HLLC_L[3] = temp * (L[3]/L_pri[0] + (S_star - L_pri[dir]) * (S_star + L_pri[3]/(L_pri[0] * (S_L - L_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_L[k] + S_L * (HLLC_L[k] - L[k]);
            }
        }
        return HLLC_f;

    }else if(S_R > 0){
        std::vector<double> f_R(4);
        std::vector<double> HLLC_f(4);
        if(dir == 1){
            f_R = f(R_pri);
            double temp = R_pri[0] * (S_R - R_pri[dir])/(S_R - S_star);
            std::vector<double> HLLC_R(4);
            HLLC_R[0] = temp;
            HLLC_R[1] = temp * S_star;
            HLLC_R[2] = temp * R_pri[2];
            HLLC_R[3] = temp * (R[3]/R_pri[0] + (S_star - R_pri[dir]) * (S_star + R_pri[3]/(R_pri[0] * (S_R - R_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_R[k] + S_R * (HLLC_R[k] - R[k]);
            }
        }else if(dir == 2){
            f_R = g(R_pri);
            double temp = R_pri[0] * (S_R - R_pri[dir])/(S_R - S_star);
            std::vector<double> HLLC_R(4);
            HLLC_R[0] = temp;
            HLLC_R[1] = temp * R_pri[1];
            HLLC_R[2] = temp * S_star;
            HLLC_R[3] = temp * (R[3]/R_pri[0] + (S_star - R_pri[dir]) * (S_star + R_pri[3]/(R_pri[0] * (S_R - R_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_R[k] + S_R * (HLLC_R[k] - R[k]);
            }
        }
        return HLLC_f;
    }else{//S_R < 0
        if(dir == 1){
            return f(R_pri);
        }else if(dir == 2){
            return g(R_pri);
        }
    }   
}

//original slic
void evolve(std::vector<std::vector<std::vector<double>>>& flux, const std::vector<std::vector<cvar> >& grid, double dt, double dx, double dy, char dir){
    std::vector<std::vector<cvar>> halfL(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfR(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfB(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfT(grid.size(), std::vector<cvar>(grid[0].size()));
    if(dir == 'x'){
        for(int j = 1; j < nyCells + ghost_cells+1; j++){
            for(int i = 1; i < nxCells + ghost_cells+1; i++){

                cvar uL = grid[i - 1][j];
                cvar u_con = grid[i][j];//u
                cvar uR = grid[i + 1][j];

                //double rho, u, v, p, E;

                pvar priL = con2pri(uL, Gamma);
                pvar pri = con2pri(u_con, Gamma);
                pvar priR = con2pri(uR, Gamma);

                double rhoL_bar, rhoR_bar, uL_bar, uR_bar, vL_bar, vR_bar, pL_bar, pR_bar;
                data_reconstruction(priL.rho, pri.rho, priR.rho, rhoL_bar,rhoR_bar);
                data_reconstruction(priL.u, pri.u, priR.u, uL_bar,uR_bar);
                data_reconstruction(priL.v, pri.v, priR.v, vL_bar,vR_bar);
                data_reconstruction(priL.p, pri.p, priR.p, pL_bar,pR_bar);

                cvar uL_half, uR_half;
                cvar uL_con, uR_con;//uL_bar, uR_bar

                uL_con = pri2con({rhoL_bar, uL_bar, vL_bar, pL_bar}, Gamma);
                uR_con = pri2con({rhoR_bar, uR_bar, vR_bar, pR_bar}, Gamma);

                half_time(uL_half, uR_half, uL_con, uR_con, dt, dx, dy, Gamma, 'x');

                halfL[i][j] = uL_half;
                halfR[i][j] = uR_half;
            }
        }
        for(int j = 0; j < nyCells + 1; j++){
            for(int i = 0; i < nxCells + 1; i++){
                cvar uL = halfR[i+1][j+2];
                cvar uR = halfL[i+2][j+2];
                //flux[i][j] = FORCE(uL, uR, dt ,dx, dy, 'x');
                flux[i][j] = HLLC_flux(uL, uR, 1);
            }
        }
    }else if(dir =='y'){
        for (int j = 1; j < nyCells + ghost_cells+1; j++){
            for (int i = 1; i < nxCells + ghost_cells+1; i++){
                cvar uB = grid[i][j - 1];
                cvar u_con = grid[i][j];
                cvar uT = grid[i][j + 1];

                double rho, u, v, p, E;
                pvar priB = con2pri(uB, Gamma);
                pvar pri = con2pri(u_con, Gamma);
                pvar priT = con2pri(uT, Gamma);
                double rhoB_bar, rhoT_bar, uB_bar, uT_bar, vB_bar, vT_bar, pB_bar, pT_bar;

                data_reconstruction(priB.rho, pri.rho, priT.rho, rhoB_bar,rhoT_bar);
                data_reconstruction(priB.u, pri.u, priT.u, uB_bar,uT_bar);
                data_reconstruction(priB.v, pri.v, priT.v, vB_bar,vT_bar);
                data_reconstruction(priB.p, pri.p, priT.p, pB_bar,pT_bar);

                cvar uB_half, uT_half;
                cvar uB_con, uT_con;// uB_bar, uT_bar

                uB_con = pri2con({rhoB_bar, uB_bar, vB_bar, pB_bar}, Gamma);
                uT_con = pri2con({rhoT_bar, uT_bar, vT_bar, pT_bar}, Gamma);

                half_time(uB_half, uT_half, uB_con, uT_con, dt, dx, dy, Gamma, 'y');

                halfB[i][j] = uB_half;
                halfT[i][j] = uT_half;
            }
        }
        for (int j = 0; j < nyCells + 1; j++){
            for (int i = 0; i < nxCells + 1; i++){
                cvar uB = halfT[i+2][j+1];
                cvar uT = halfB[i+2][j+2];
                //flux[i][j] = FORCE(uB, uT, dt, dx, dy, 'y');
                flux[i][j] = HLLC_flux(uB,uT,2);

            }
        }

    }

}
double updatedt(const std::vector<std::vector<cvar> >& grid, double dx,double dy, double Gamma, double CFL){
    double max_speed = 0.0;
    for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
            pvar pri = con2pri(grid[i][j], Gamma);
            double v_ij = sqrt(pri.u * pri.u + pri.v * pri.v);
            double cs = calcs(pri.p, pri.rho, Gamma);
            // double speedx = fabs(pri.u) + cs;
            // double speedy = fabs(pri.v) + cs;
            double speed_diag = fabs(v_ij) + cs;

            if (speed_diag > max_speed){
                max_speed = speed_diag;
            }
        }
    }
    return CFL * fmin(dx,dy) / max_speed;
    //std::cout << "update: max_speed = " << max_speed << ", dt = " << CFL * dx / max_speed << std::endl;
}
void initializeGrid(std::vector<std::vector<cvar> >& grid){
    for(int i = 0; i < nxCells + 2*ghost_cells; i++){
        double x = (i - ghost_cells + 0.5) * dx;
        for(int j = 0; j < nyCells + 2*ghost_cells; j++){
            double y = (j - ghost_cells + 0.5) * dy;
            // if(x <= 0.2){
            //     grid[i][j] = pri2con({1,0,0,1},Gamma);
            // }
            // else{
            //     grid[i][j] = pri2con({0.125,0,0,0.1},Gamma);
            // }
            // if (x <= 0.2){
            //     cvar temp = pri2con ({1.3764, 0.394, 0.0, 1.5698}, Gamma);
            //     grid[i][j] = temp;
            // }else{
            //     cvar temp = pri2con({1.0, 0.0, 0.0, 1.0}, Gamma);
            //     grid[i][j] = temp;
            // }
            cvar temp = pri2con({1.0, 1.0, 0.0, 1.0}, Gamma);
            //cvar temp = pri2con({1.0, 0.0, 0.0, 1.0}, Gamma);
            grid[i][j] = temp;
            
            // if((x - 1) * (x - 1) + (y - 1) * (y - 1) < 0.4 * 0.4){
            //     grid[i][j] = pri2con({1.0, 0.0, 0.0, 1},Gamma);
            // }else{
            //     grid[i][j] = pri2con({0.125, 0.0, 0.0, 0.1},Gamma);
            // }
        }
    }
}
void initializeLevelSet(std::vector<std::vector<double>> &phi, ShapeType shape) {


    for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
        double x = x0 + (i - ghost_cells + 0.5) * dx;
        for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
            double y = y0_new + (j - ghost_cells + 0.5) * dy;


            switch (shape) {
            case CIRCLE:
            {
                // (0.6, 0.5)，r=0.2
                double cx = 0.6, cy = 0.5, R = 0.2;
                double dist = R - std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                phi[i][j] = dist;
                break;
            }
            case SQUARE:
            {
                // square (0.6, 0.5)，length 0.4 => x=[0.4,0.8]，y=[0.3,0.7]
                double half_len = 0.2; // half length
                double x_min = 0.6 - half_len;
                double x_max = 0.6 + half_len;
                double y_min = 0.5 - half_len;
                double y_max = 0.5 + half_len;
            
                double dist;
                bool inside = (x >= x_min && x <= x_max && y >= y_min && y <= y_max);//x-coord between 0.4 & 0.8, y-coord betwwn 0.3 & 0.7
                if (inside) {
                    // shortest dist from 
                    double dx_in = std::min(x - x_min, x_max - x);
                    double dy_in = std::min(y - y_min, y_max - y);

                    dist = std::min(dx_in, dy_in);
                    phi[i][j] = dist;  //positive inside
                } else {
                    // outside：
                    double dx_rect = 0.0;
                    if      (x < x_min) dx_rect = (x_min - x);
                    else if (x > x_max) dx_rect = (x - x_max);

                    double dy_rect = 0.0;
                    if      (y < y_min) dy_rect = (y_min - y);
                    else if (y > y_max) dy_rect = (y - y_max);

                    dist = std::sqrt(dx_rect * dx_rect + dy_rect * dy_rect);
                    phi[i][j] = - dist;
                }
                break;
            }
            case TWO_CIRCLES:
            {
                // 2 circles，radius0.2，(0.6, 0.25)  (0.6, 0.75)
                double r = 0.2;
                double cx1 = 0.6, cy1 = 0.25;
                double cx2 = 0.6, cy2 = 0.75;
                double dist1 = r - sqrt((x - cx1) * (x - cx1) + (y - cy1) * (y - cy1));
                double dist2 = r - sqrt((x - cx2) * (x - cx2) + (y - cy2) * (y - cy2));
                // union => =maximum dist od 2 distance
                phi[i][j] = std::max(dist1, dist2);
                break;
            }
            case TWO_OVERLAPPING_CIRCLES:
            {
                // 2circles，radius0.2，(0.6, 0.35)  (0.6, 0.65)
                double r = 0.2;
                double cx1 = 0.6, cy1 = 0.35;
                double cx2 = 0.6, cy2 = 0.65;
                double dist1 = r - std::sqrt((x - cx1) * (x - cx1) + (y - cy1) * (y - cy1));
                double dist2 = r - std::sqrt((x - cx2) * (x - cx2) + (y - cy2) * (y - cy2));
                phi[i][j] = std::max(dist1, dist2);
                break;
            }
            case STATIC:
            {
                //phi[i][j] = 0.2 - sqrt(((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)));
                double cx = 0.5, cy = 0.5, R = 0.2;
                double dist = R - std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
                phi[i][j] = dist;
                break;
            }
            

            
            default:

                phi[i][j] = 999.0;
                break;

        }
            
      }
    }

//trans_phi(phi);
}
void compute_normal_vector(std::vector<std::vector<std::vector<double>>> &normal_vec, std::vector<std::vector<double>> phi){
    for (int i = ghost_cells; i < nxCells + 2 * ghost_cells - 1; i++) {
        double x = x0 + (i - ghost_cells + 0.5) * dx;
        for (int j = ghost_cells; j < nyCells + 2 * ghost_cells - 1; j++) {
            double y = y0_new + (j - ghost_cells + 0.5) * dy;
                    double phi_x = (phi[i+1][j] - phi[i-1][j]) / (2.0 * dx) + 1e-16;
                    double phi_y = (phi[i][j+1] - phi[i][j-1]) / (2.0 * dy) + 1e-16;
                    double norm = sqrt(phi_x*phi_x + phi_y*phi_y) + 1e-16;
                    normal_vec[i][j][0] = phi_x  / norm;
                    normal_vec[i][j][1] = phi_y  / norm;
                   
            }
        }
    }

void trans(std::vector<std::vector<cvar> >& grid){
    //BC
    for (int j = ghost_cells; j < ghost_cells + nyCells; j++) {
        for (int i = 0; i < ghost_cells; i++) {
            grid[i][j] = grid[ghost_cells][j]; // 左边界
            grid[nxCells + ghost_cells + i][j] = grid[nxCells + ghost_cells - 1][j]; // 右边界
        }
    }

    for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
        for (int j = 0; j < ghost_cells; j++) {
            grid[i][j] = grid[i][ghost_cells]; // 底边界
            grid[i][nyCells + ghost_cells + j] = grid[i][nyCells + ghost_cells - 1]; // 顶边界
        }
    }
}

pvar normal2xy(std::vector<double> L_star, std::vector<double> n, std::vector<double> Vt){
    pvar ans;
    ans.rho = L_star[0];
    ans.u = Vt[0] + n[0] * L_star[1];
    ans.v = Vt[1] + n[1] * L_star[1];
    ans.p = L_star[2];
    return ans;

}
double getVn(std::vector<double> n, pvar pri){
    return n[0] * pri[1] + n[1] * pri[2];
}
std::vector<double> getVt(std::vector<double> n , pvar pri){
    double Vn = getVn(n, pri);
    std::vector<double> Vt(2);
    Vt[0] = pri[1] - Vn * n[0];
    Vt[1] = pri[2] - Vn * n[1];
    return Vt;
}
std::vector<double> get_R_state(std::vector<double> L_state_pri_normal, std::vector<double> n, std::vector<double> v_rigid){
    return {L_state_pri_normal[0], - L_state_pri_normal[1] + 2 * (v_rigid[0] * n[0] + v_rigid[1] * n[1]) ,L_state_pri_normal[2]};
}


double getQijk(double dx, double dy, double Qx, double Qy, std::vector<double> n){
    double up = dy * fabs(n[0]) * Qx + dx * fabs(n[1]) * Qy;
    double bottom = dy * fabs(n[0]) + dx * fabs(n[1]);
    return up / bottom ;
}

void constant_extrap(std::vector<std::vector<cvar>> &grid, std::vector<std::vector<std::vector<double>>> n, std::vector<std::vector<int>> interface, std::vector<std::vector<double>> phi ){
    for(int i = 2; i < ghost_cells + nxCells; i++){
        for(int j = 2; j < ghost_cells + nyCells; j++){
            if(phi[i][j] >= 0 && interface[i][j] == 0){
                grid[i][j] = {inf, inf, inf, inf};
            }
        }
    }
    
    for(int i = 2; i < ghost_cells + nxCells; i++){
        for(int j = 2; j < nyCells + ghost_cells; j++){
            if(phi[i][j] >= 0 && interface[i][j] == 0){
                cvar Qx;
                cvar Qy;
                if(n[i][j][0] > 0){
                    Qx = grid[i-1][j];
                }else{
                    Qx = grid[i+1][j];
                }
                if(n[i][j][1] > 0){
                    Qy = grid[i][j-1];
                }else{
                    Qy = grid[i][j+1];
                }
                for(int k = 0; k < 4; k++){
                    double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                    if(fabs(Qijk) < fabs(grid[i][j][k])){
                        grid[i][j][k] = Qijk;
                    }
                    
                }
            }

        }
    }
    
    for(int i = ghost_cells + nxCells - 1; i > 1; i--){
        for(int j = ghost_cells + nyCells - 1; j > 1; j--){
            if(phi[i][j] >= 0 && interface[i][j] == 0){
                cvar Qx;
                cvar Qy;
                if(n[i][j][0] > 0){
                    Qx = grid[i-1][j];
                }else{
                    Qx = grid[i+1][j];
                }
                if(n[i][j][1] > 0){
                    Qy = grid[i][j-1];
                }else{
                    Qy = grid[i][j+1];
                }
                for(int k = 0; k < 4; k++){
                    double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                    if(fabs(Qijk) < fabs(grid[i][j][k])){
                        grid[i][j][k] = Qijk;
                    }
                    
                }
            }
        }
    }
    
    for(int i = ghost_cells + nxCells - 1; i > 1; i--){
        for(int j = 2; j < ghost_cells + nyCells; j++){
            if(phi[i][j] >= 0 && interface[i][j] == 0){
                cvar Qx;
                cvar Qy;
                if(n[i][j][0] > 0){
                    Qx = grid[i-1][j];
                }else{
                    Qx = grid[i+1][j];
                }
                if(n[i][j][1] > 0){
                    Qy = grid[i][j-1];
                }else{
                    Qy = grid[i][j+1];
                }
                for(int k = 0; k < 4; k++){
                    double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                    if(fabs(Qijk) < fabs(grid[i][j][k])){
                        grid[i][j][k] = Qijk;
                    }
                    
                }
            }
        }
    }
    
    for(int i = 2; i < ghost_cells + nxCells; i++){
        for(int j = ghost_cells + nyCells - 1; j > 1; j--){
            if(phi[i][j] >= 0 && interface[i][j] == 0){
                cvar Qx;
                cvar Qy;
                if(n[i][j][0] > 0){
                    Qx = grid[i-1][j];
                }else{
                    Qx = grid[i+1][j];
                }
                if(n[i][j][1] > 0){
                    Qy = grid[i][j-1];
                }else{
                    Qy = grid[i][j+1];
                }
                for(int k = 0; k < 4; k++){
                    double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                    if(fabs(Qijk) < fabs(grid[i][j][k])){
                        grid[i][j][k] = Qijk;
                    }
                    
                }
            }
        }
    }

}

int main(){
    std::vector<std::vector<cvar>> grid(nxCells + 2 * ghost_cells, std::vector<cvar>(nyCells + 2 * ghost_cells));
    initializeGrid(grid);
    std::vector<std::vector<cvar>> gridplus(nxCells + 2 * ghost_cells, std::vector<cvar>(nyCells + 2 * ghost_cells));
    std::vector<std::vector<double>> phi(nxCells + 2 * ghost_cells, std::vector<double> (nyCells + 2 * ghost_cells));
    initializeLevelSet(phi, STATIC);
    std::vector<std::vector<std::vector<double>>> normal_vec(nxCells + 2 * ghost_cells, std::vector<std::vector<double>> (nyCells + 2 * ghost_cells, std::vector<double>(2)));
    compute_normal_vector(normal_vec, phi);
    std::vector<std::vector<std::vector<double>>> flux_x(nxCells + 1, std::vector<std::vector<double>>(nyCells + 1, std::vector<double>(4)));
    std::vector<std::vector<std::vector<double>>> flux_y(nxCells + 1, std::vector<std::vector<double>>(nyCells + 1, std::vector<double>(4)));
    std::vector<std::vector<int>> interface(nxCells + 2 * ghost_cells, std::vector<int>(nyCells + 2 * ghost_cells));




    int iter = 0; double t = 0.0;
    do {
        

        for(int i = ghost_cells; i < nxCells +  ghost_cells; i++){
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                double right_nei = phi[i][j] * phi[i + 1][j];
                double left_nei = phi[i][j] * phi[i - 1][j];
                double top_nei = phi[i][j] * phi[i][j + 1];
                double bot_nei = phi[i][j] * phi[i][j - 1];
                if (right_nei <= 0.0||left_nei <= 0.0||top_nei <= 0.0||bot_nei <= 0.0){
                    interface[i][j] = 1;
                }else{
                    interface[i][j] = 0;
                }
            }
        }

        //gfm
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++){
            double x = x0 + (i - ghost_cells + 0.5) * dx;
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                double y = y0_new + (j - ghost_cells + 0.5) * dy;
                if(phi[i][j] >= 0 && interface[i][j] == 1){

                    std::vector<double> p = {x - phi[i][j] * normal_vec[i][j][0], y - phi[i][j] * normal_vec[i][j][1]};
                    std::vector<double> p2 = {p[0] - 1.5 * dx * normal_vec[i][j][0], p[1] - 1.5 * dy * normal_vec[i][j][1]};
                    std::vector<double> p2_ij = {((p2[0] - x0)/dx) + ghost_cells - 0.5, ((p2[1] - y0_new)/dy) + ghost_cells - 0.5};
                    int temp1 = std::floor(p2_ij[0]);
                    int temp2 = std::floor(p2_ij[1]);
                    std::vector<int> Q11 = {temp1, temp2};
                    temp1 = std::floor(p2_ij[0]);
                    temp2 = std::ceil(p2_ij[1]);
                    std::vector<int> Q12 = {temp1, temp2};
                    temp1 = std::ceil(p2_ij[0]);
                    temp2 = std::floor(p2_ij[1]);
                    std::vector<int> Q21 = {temp1, temp2};
                    temp1 = std::ceil(p2_ij[0]);
                    temp2 = std::ceil(p2_ij[1]);
                    std::vector<int> Q22 = {temp1, temp2};
                    double x1 = x0 + (Q11[0] - ghost_cells + 0.5) * dx;
                    double x2 = x0 + (Q21[0] - ghost_cells + 0.5) * dx;
                    double y1 = y0_new + (Q11[1] - ghost_cells + 0.5) * dy;
                    double y2 = y0_new + (Q12[1] - ghost_cells + 0.5) * dy;
                    cvar L_state;
                    
                    if(p2_ij[0] == std::floor(p2_ij[0]) && p2_ij[1] == std::floor(p2_ij[1])){//if int, floor = ceil
                        L_state.rho  = grid[Q11[0]][Q11[1]][0];
                        L_state.rhou = grid[Q11[0]][Q11[1]][1];
                        L_state.rhov = grid[Q11[0]][Q11[1]][2];
                        L_state.E    = grid[Q11[0]][Q11[1]][3];

                    }else if(p2_ij[0] == std::floor(p2_ij[0])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (y2 - p2[1])/(y2 - y1) * grid [Q11[0]][Q11[1]][k] + (p2[1] - y1)/(y2 - y1) * grid[Q12[0]][Q12[1]][k];
                            L_state[k] = R3;
                        }
                    }else if(p2_ij[1] == std::floor(p2_ij[1])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (x2 - p2[0])/(x2 - x1) * grid[Q11[0]][Q11[1]][k] + (p2[0] - x1)/(x2 - x1) * grid[Q21[0]][Q21[1]][k];
                            L_state[k] = R3;
                        }
                    }else{
                        for(int k = 0; k < 4; k++){
                            double R1 = (x2 - p2[0])/(x2 - x1) * grid[Q11[0]][Q11[1]][k] + (p2[0] - x1)/(x2 - x1) * grid[Q21[0]][Q21[1]][k];
                            double R2 = (x2 - p2[0])/(x2 - x1) * grid[Q12[0]][Q12[1]][k] + (p2[0] - x1)/(x2 - x1) * grid[Q22[0]][Q22[1]][k];
                            double R3 = (y2 - p2[1])/(y2 - y1) * R1 + (p2[1] - y1)/(y2 - y1) * R2;
                            L_state[k] = R3;
                        }
                    }
                    L_state = grid[Q11[0]][Q11[1]];

                    pvar L_state_pri = con2pri(L_state, Gamma);
                    std::vector<double> L_state_pri_normal = {L_state_pri[0], getVn(normal_vec[i][j], L_state_pri), L_state_pri[3]};
                    std::vector<double> Vt = getVt(normal_vec[i][j], L_state_pri);
                    std::vector<double> R_state_pri_normal = get_R_state(L_state_pri_normal, normal_vec[i][j], v_rigid);

                    std::vector<double> L_star_pri_normal = HLLC(L_state_pri_normal, R_state_pri_normal);
                    pvar L_star_pri = normal2xy(L_star_pri_normal, normal_vec[i][j], Vt);
                    cvar L_star_con = pri2con(L_star_pri, Gamma);

                    grid[i][j] = L_star_con;
                }

            }
        }

        constant_extrap(grid, normal_vec, interface, phi);


        iter += 1;
        dt = updatedt(grid,dx,dy,Gamma, CFL);
        t += dt;
        if(t > T_end){
            t = T_end;
        }
        std::cout<< std::round(t/T_end * 1000)/1000 * 100 << "%" <<" dt: "<<dt<< std::endl;



        evolve(flux_x, grid, dt, dx, dy, 'x');
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
                gridplus[i][j].rho  = grid[i][j].rho-(dt / dx) * (flux_x[i-1][j-2][0] - flux_x[i-2][j-2][0]);
                gridplus[i][j].rhou = grid[i][j].rhou-(dt / dx) * (flux_x[i-1][j-2][1] - flux_x[i-2][j-2][1]);
                gridplus[i][j].rhov = grid[i][j].rhov-(dt / dx) * (flux_x[i-1][j-2][2] - flux_x[i-2][j-2][2]);
                gridplus[i][j].E    = grid[i][j].E-(dt / dx) * (flux_x[i-1][j-2][3] - flux_x[i-2][j-2][3]);
            }
        }
        grid = gridplus;
        trans(grid);



        evolve(flux_y, grid, dt, dx, dy, 'y');
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
                gridplus[i][j].rho  = grid[i][j].rho - (dt / dy) * (flux_y[i-2][j-1][0] - flux_y[i-2][j-2][0]);
                gridplus[i][j].rhou = grid[i][j].rhou - (dt / dy) * (flux_y[i-2][j-1][1] - flux_y[i-2][j-2][1]);
                gridplus[i][j].rhov = grid[i][j].rhov - (dt / dy) * (flux_y[i-2][j-1][2] - flux_y[i-2][j-2][2]);
                gridplus[i][j].E    = grid[i][j].E - (dt / dy) * (flux_y[i-2][j-1][3] - flux_y[i-2][j-2][3]);
            }
        }
        grid = gridplus;
        trans(grid);














        std::string path = "data/";
        std::string file = path + "data" + std::to_string(iter) + ".dat";
        std::ofstream writefile(file);

        for(int j = 2; j < nyCells + ghost_cells; j++){
            double y = y0_new + (j - 1.5) * dy;
            for(int i = 2; i < nxCells + ghost_cells; i++){
                double x  = x0 + (i - 1.5) * dx;
                pvar prim = con2pri(grid[i][j], Gamma);
                double speed_mag = sqrt(prim.u * prim.u + prim.v * prim.v);
                writefile << x << " " << y << " " << prim.rho << " " << speed_mag << " " << prim.p << " " \
                       << phi[i][j] <<  std::endl;
            }
        }




    } while(t < T_end);



    // std::ofstream writefile2("final.dat");

    // for(int j = 2; j < nyCells + ghost_cells; j++){
    //     double y = y0_new + (j - 1.5) * dy;
    //     for(int i = 2; i < nxCells + ghost_cells; i++){
    //         double x  = x0 + (i - 1.5) * dx;
    //         pvar prim = con2pri(grid[i][j], Gamma);
    //         double speed_mag = sqrt(prim.u * prim.u + prim.v * prim.v);
    //         writefile2 << x << " " << y << " " << prim.rho << " " << speed_mag << " " << prim.p << " " \
    //                << phi[i][j] <<  std::endl;
    //     }
    // }


    std::ofstream writefile2("final.dat");

    for(int j = 2; j < nyCells + ghost_cells - 1; j++){
        double y = y0_new + (j - 1.5) * dy;

        for(int i = 2; i < nxCells + ghost_cells - 1; i++){
            double x = x0 + (i - 1.5) * dx;

            
            pvar prim       = con2pri(grid[i][j], Gamma);
            double rho      = prim.rho;
            double u        = prim.u;
            double v        = prim.v;
            double p        = prim.p;
            double speed_mag = std::sqrt(u*u + v*v);
            double phi_val  = phi[i][j];

            
            pvar prim_left  = con2pri(grid[i-1][j],   Gamma);
            pvar prim_right = con2pri(grid[i+1][j],   Gamma);
            pvar prim_down  = con2pri(grid[i][j-1],   Gamma);
            pvar prim_up    = con2pri(grid[i][j+1],   Gamma);

            double rho_left  = prim_left.rho;
            double rho_right = prim_right.rho;
            double rho_down  = prim_down.rho;
            double rho_up    = prim_up.rho;

            
            double drhodx = (rho_right - rho_left) / (2.0 * dx);
            double drhody = (rho_up    - rho_down) / (2.0 * dy);

            double grad_rho = std::sqrt(drhodx * drhodx + drhody * drhody);

            // mock-schlieren: 
            double mock_schlieren = std::exp(-20.0 * grad_rho / (1000.0 * rho));

            // ：x, y, rho, speed_mag, p, phi, mock_schlieren
            if(phi[i][j] < 0){
                writefile2 << x << " " << y << " "
                        << rho << " "
                        << speed_mag << " "
                        << p << " "
                        << phi_val << " "
                        << mock_schlieren << std::endl;
            }else{
                writefile2 << x << " " << y << " " << 0 << " " << 0<< " " << 0<< " " << 0 << " "<< 0<< std::endl;
            }
        }
    }

    writefile2.close();
    











    return 0;
}