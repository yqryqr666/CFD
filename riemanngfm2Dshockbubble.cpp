#include <float.h>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <string>
#include <omp.h>


enum cons{rho=0, rhou=1, rhov=2, Ene=3};
enum pri{u=1, v=2, pre=3};

constexpr double Gamma1 = 1.4;
constexpr double Gamma2 = 1.67;
constexpr double CFL = 0.8;
constexpr double x0 = 0, x1 = 0.225, y0_new = -0.0445, y1_new = 0.0445;
//constexpr double x0 = 0, x1 = 0.225, y0_new = -0.445, y1_new = 0.445;
constexpr double T_end = 0.000293525; //t1=0.000293525, t2=0.00058705, t3=0.000880575, t4=0.0011741
constexpr int nxCells = 1000; constexpr int nyCells = 394; //500*197 放大2倍
//int nxPoints = nxCells + 1, nyPoints = nyCells + 1;
constexpr int ghost_cells = 2;
//double dx = (x1 - x0)/(nxPoints - 2), dy = (y1_new - y0_new)/(nyPoints - 2);
constexpr double dx = (x1 - x0)/nxCells, dy = (y1_new - y0_new)/nyCells;
double dt = 0;
constexpr double weight = 0;


constexpr std::array<double,4> w_shock = {1.77551,110.6273,0,159059.99}; //rou vx vy p
constexpr std::array<double,4> w_air = {1.29,0,0,101325}; //rou vx vy p
constexpr std::array<double,4> w_helium = {0.214,0,0,101325}; //rou vx vy p



std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> grid1 = {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> grid2 = {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> gridplus1 = {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> gridplus2 = {};

std::array<std::array<double,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> phi = {};
std::array<std::array<double,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> phiplus = {};

std::array<std::array<std::array<double, 2>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> normal_vec = {};
std::array<std::array<int,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> interface = {};

std::array<std::array<std::array<double, 4>,nyCells + 1>,nxCells + 1> flux1_x = {};
std::array<std::array<std::array<double, 4>,nyCells + 1>,nxCells + 1> flux1_y = {};
std::array<std::array<std::array<double, 4>,nyCells + 1>,nxCells + 1> flux2_x = {};
std::array<std::array<std::array<double, 4>,nyCells + 1>,nxCells + 1> flux2_y = {};

std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u1_L_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u1_R_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u1_B_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u1_T_bar_half= {};

std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u2_L_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u2_R_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u2_B_bar_half= {};
std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> u2_T_bar_half= {};

//std::array<std::array<std::array<double, 2>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> normal_vec = {};




std::array<double, 4> con2pri(const std::array<double, 4> &cons, const double &Gamma){
    std::array<double, 4> pri = {};
    pri[rho] = cons[rho];
    pri[u] = cons[rhou] / cons[rho];
    pri[v] = cons[rhov] / cons[rho];
    pri[pre] = (Gamma - 1) * (cons[Ene] - 0.5 * (cons[rhov] * cons[rhov] + cons[rhou] * cons[rhou])/ cons[rho]);
    return pri;
}
std::array<double, 4> pri2con(const std::array<double, 4> &prim, const double &Gamma){

    std::array<double, 4> con = {};
    con[rho] = prim[rho];
    con[rhou] = prim[rho] * prim[u];
    con[rhov] = prim[rho] * prim[v];
    con[Ene] = prim[pre] / (Gamma - 1) + 0.5 * prim[rho] * (prim[v] * prim[v] + prim[u] * prim[u]); //E = p/(Gamma-1) + 0.5 rho v^2
    return con;

}


std::array<double, 3> con2pri_hllc(std::array<double, 3> cons, double Gamma){
    std::array<double, 3> prim;
    prim[0] = cons[0];
    prim[1] = cons[1]/cons[0];
    prim[2] = (Gamma-1)*(cons[2]-((cons[1]*cons[1])/(2*cons[0])));
    return prim;
}
std::array<double, 3> pri2con_hllc(std::array<double, 3> prim, double Gamma){
    std::array<double, 3> cons;
    cons[0] = prim[0];
    cons[1] = prim[0]*prim[1];
    cons[2] = prim[2]/(Gamma-1)+(prim[0]*prim[1])*(prim[0]*prim[1])/(2*prim[0]);
    return cons;
}
////////////////////////////////////////////////////////////////////////
double calcs(const double &p, const double &rho, const double &Gamma){
    return sqrt(Gamma * p / rho);
}
////////////////////////////////////////////////////////////////////////
double updatedt(const std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> &grid1,
                const std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> &grid2){
    double max_speed = 0.0;
    for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
            std::array<double,4> pri1 = con2pri(grid1[i][j], Gamma1);
            std::array<double,4> pri2 = con2pri(grid2[i][j], Gamma2);
            // double v_ij1 = sqrt(pri1.u * pri1.u + pri1.v * pri1.v);
            // double v_ij2 = sqrt(pri2.u * pri2.u + pri2.v * pri2.v);
            double cs1 = calcs(pri1[3], pri1[0], Gamma1);
            double cs2 = calcs(pri2[3], pri2[0], Gamma2);
            // double speedx = fabs(pri.u) + cs;
            // double speedy = fabs(pri.v) + cs;
            // double speed_diag1 = fabs(v_ij1) + cs1;
            // double speed_diag2 = fabs(v_ij2) + cs2;
            // double speed_diag = std::max(speed_diag1, speed_diag2);
            // if (speed_diag > max_speed){
            //     max_speed = speed_diag;
            // }
            double speed_x1 = fabs(pri1[1]) + cs1;
            double speed_y1 = fabs(pri1[2]) + cs1;
            double speed_x2 = fabs(pri2[1]) + cs2;
            double speed_y2 = fabs(pri2[2]) + cs2;
            double speed_x = std::max(speed_x1, speed_x2);
            double speed_y = std::max(speed_y1, speed_y2);
            double current_speed = std::max(speed_x, speed_y);
            if(current_speed > max_speed){
                max_speed = current_speed;
            }
        }
    }
    return CFL * std::min(dx,dy) / max_speed;
    //std::cout << "update: max_speed = " << max_speed << ", dt = " << CFL * dx / max_speed << std::endl;
}

////////////////////////////////////////////////////////////////////////
double Minbee(const double &r){

    if(r <= 0) return 0;
    if(r <= 1) return r;
    return std::min(1.0, 2.0 / (1 + r));
}
////////////////////////////////////////////////////////////////////////
std::array<double, 4> g(const std::array<double, 4>  &cons, const double &Gamma){ // y-direction flux
    //std::vector<double> g = {rho * v, rho * u * v, rho * v * v + p, (E + p) * v};

    std::array<double, 4>  prim = con2pri(cons,Gamma);

    std::array<double, 4>  glux = {};
    glux[0] = prim[rho]* prim[v];
    glux[1] = prim[rho]* prim[u] * prim[v];
    glux[2] = prim[rho] * prim[v] * prim[v] + prim[pre];
    glux[3] = (cons[Ene] + prim[pre]) * prim[v];
    return glux;
}

std::array<double, 4> f(const std::array<double, 4>  &cons, const double &Gamma){ // x-direction flux
    //std::vector<double> f = {rho * u, rho * u * u + p, rho * u * v, (E + p) * u};
    std::array<double, 4>  prim = con2pri(cons,Gamma);
    std::array<double, 4>  flux = {};
    flux[0] = prim[rho]* prim[u];
    flux[1] = prim[rho]* prim[u] * prim[u] + prim[pre];
    flux[2] = prim[rho]* prim[u] * prim[v];
    flux[3] = (cons[Ene] + prim[pre]) * prim[u];
    return flux;
}
////////////////////////////////////////////////////////////////////////
std::array<double,4> HLLC_flux(const std::array<double,4> &L, const std::array<double,4> &R, const double &Gamma, int dir){
    std::array<double,4> L_pri = con2pri(L, Gamma);
    std::array<double,4> R_pri = con2pri(R, Gamma);

    double cs_l = calcs(L_pri[pre], L_pri[rho], Gamma);
    double cs_R = calcs(R_pri[pre], R_pri[rho], Gamma);
    double S_R = std::max(fabs(L_pri[dir]) + cs_l, fabs(R_pri[dir]) + cs_R);
    double S_L = - S_R;
    double S_star = (R_pri[pre] - L_pri[pre] + L_pri[rho] * L_pri[dir] * (S_L - L_pri[dir]) - R_pri[rho] * R_pri[dir] * (S_R - R_pri[dir]))/(L_pri[rho] * (S_L - L_pri[dir]) - R_pri[rho] * (S_R - R_pri[dir]));
    if(S_L > 0){
        if(dir == 1){
            return f(L, Gamma);
        }else if(dir == 2){
            return g(L, Gamma);
        }
    }else if(S_star > 0){
        std::array<double,4> f_L = {};
        std::array<double,4> HLLC_f = {};
        if(dir == 1){
            f_L = f(L, Gamma);
            double temp = L_pri[0] * (S_L - L_pri[dir])/(S_L - S_star);
            std::array<double,4> HLLC_L = {};
            HLLC_L[0] = temp;
            HLLC_L[1] = temp * S_star;
            HLLC_L[2] = temp * L_pri[2];
            HLLC_L[3] = temp * (L[3]/L_pri[0] + (S_star - L_pri[dir]) * (S_star + L_pri[3]/(L_pri[0] * (S_L - L_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_L[k] + S_L * (HLLC_L[k] - L[k]);
            }
        }else if(dir == 2){
            f_L = g(L, Gamma);
            double temp = L_pri[0] * (S_L - L_pri[dir])/(S_L - S_star);
            std::array<double,4> HLLC_L = {};
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
        std::array<double,4> f_R = {};
        std::array<double,4> HLLC_f = {};
        if(dir == 1){
            f_R = f(R, Gamma);
            double temp = R_pri[0] * (S_R - R_pri[dir])/(S_R - S_star);
            std::array<double,4> HLLC_R = {};
            HLLC_R[0] = temp;
            HLLC_R[1] = temp * S_star;
            HLLC_R[2] = temp * R_pri[2];
            HLLC_R[3] = temp * (R[3]/R_pri[0] + (S_star - R_pri[dir]) * (S_star + R_pri[3]/(R_pri[0] * (S_R - R_pri[dir]))));
            for (int k = 0; k < 4; k++){
                HLLC_f[k] = f_R[k] + S_R * (HLLC_R[k] - R[k]);
            }
        }else if(dir == 2){
            f_R = g(R, Gamma);
            double temp = R_pri[0] * (S_R - R_pri[dir])/(S_R - S_star);
            std::array<double,4> HLLC_R = {};
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
            return f(R, Gamma);
        }else if(dir == 2){
            return g(R, Gamma);
        }
    }
}
////////////////////////////////////////////////////////////////////////
void evolve(const double &dt, char dir){
    if(dir == 'x'){
        for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
            for(int j = ghost_cells - 1; j < nyCells + ghost_cells + 1; j++){
                std::array<double,4> u1_L_bar = {};
                std::array<double,4> u1_R_bar = {};
                std::array<double,4> u2_L_bar = {};
                std::array<double,4> u2_R_bar = {};
                for(int k = 0; k < 4; k++){
                    //m1
                    double deltaiminus1 = grid1[i][j][k] - grid1[i-1][j][k];
                    double deltaiplus1 = grid1[i+1][j][k] - grid1[i][j][k];
                    double deltai1 = 0.5 * ((1 + weight) * deltaiminus1 + (1 - weight) * deltaiplus1);
                    double r1;
                    if(deltaiminus1 == 0 && deltaiplus1 == 0){
                        r1 = 0;
                    }else if (deltaiminus1 == 0){
                        r1 = 0;
                    }else if (deltaiplus1 == 0){
                        r1 = 9999;
                    }else{
                        r1 = deltaiminus1 / deltaiplus1;
                    }
                    double limiter1 = Minbee(r1);
                    u1_L_bar[k] = grid1[i][j][k] - 0.5 * limiter1 * deltai1;
                    u1_R_bar[k] = grid1[i][j][k] + 0.5 * limiter1 * deltai1;

                    //m2
                    double deltaiminus2 = grid2[i][j][k] - grid2[i-1][j][k];
                    double deltaiplus2 = grid2[i+1][j][k] - grid2[i][j][k];
                    double deltai2 = 0.5 * ((1 + weight) * deltaiminus2 + (1 - weight) * deltaiplus2);
                    double r2;
                    if(deltaiminus2 == 0 && deltaiplus2 == 0){
                        r2 = 0;
                    }else if (deltaiminus2 == 0){
                        r2 = 0;
                    }else if (deltaiplus2 == 0){
                        r2 = 9999;
                    }else{
                        r2 = deltaiminus2 / deltaiplus2;
                    }
                    double limiter2 = Minbee(r2);
                    u2_L_bar[k] = grid2[i][j][k] - 0.5 * limiter2 * deltai2;
                    u2_R_bar[k] = grid2[i][j][k] + 0.5 * limiter2 * deltai2;
                }
                std::array<double, 4> flux1_L = f(u1_L_bar, Gamma1);
                std::array<double, 4> flux1_R = f(u1_R_bar, Gamma1);
                std::array<double, 4> flux2_L = f(u2_L_bar, Gamma2);
                std::array<double, 4> flux2_R = f(u2_R_bar, Gamma2);
                for(int k = 0; k < 4; k++){
                    u1_L_bar_half[i][j][k] = u1_L_bar[k] - 0.5 * dt / dx * (flux1_R[k] - flux1_L[k]);
                    u1_R_bar_half[i][j][k] = u1_R_bar[k] - 0.5 * dt / dx * (flux1_R[k] - flux1_L[k]);
                    u2_L_bar_half[i][j][k] = u2_L_bar[k] - 0.5 * dt / dx * (flux2_R[k] - flux2_L[k]);
                    u2_R_bar_half[i][j][k] = u2_R_bar[k] - 0.5 * dt / dx * (flux2_R[k] - flux2_L[k]);
                }
            }
        }
        for(int i = 0; i < nxCells + 1; i++){
            for(int j = 0; j < nyCells + 1; j++){
                flux1_x[i][j] = HLLC_flux(u1_R_bar_half[i + 1][j + 2], u1_L_bar_half[i + 2][j + 2], Gamma1, 1);
                flux2_x[i][j] = HLLC_flux(u2_R_bar_half[i + 1][j + 2], u2_L_bar_half[i + 2][j + 2], Gamma2, 1);
            }
        }
        for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
            for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
                for(int k = 0; k < 4; k++){
                    gridplus1[i][j][k] = grid1[i][j][k] - dt / dx * (flux1_x[i - 1][j - 2][k] - flux1_x[i - 2][j - 2][k]);
                    gridplus2[i][j][k] = grid2[i][j][k] - dt / dx * (flux2_x[i - 1][j - 2][k] - flux2_x[i - 2][j - 2][k]);
                }
            }
        }
        grid1 = gridplus1;
        grid2 = gridplus2;

    }else{
        for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
            for(int j = ghost_cells - 1; j < nyCells + ghost_cells + 1; j++){
                std::array<double,4> u1_B_bar = {};
                std::array<double,4> u1_T_bar = {};
                std::array<double,4> u2_B_bar = {};
                std::array<double,4> u2_T_bar = {};
                for(int k = 0; k < 4; k++){
                    //m1
                    double deltaiminus1 = grid1[i][j][k] - grid1[i][j - 1][k];
                    double deltaiplus1 = grid1[i][j + 1][k] - grid1[i][j][k];
                    double deltai1 = 0.5 * ((1 + weight) * deltaiminus1 + (1 - weight) * deltaiplus1);
                    double r1;
                    if(deltaiminus1 == 0 && deltaiplus1 == 0){
                        r1 = 0;
                    }else if (deltaiminus1 == 0){
                        r1 = 0;
                    }else if (deltaiplus1 == 0){
                        r1 = 9999;
                    }else{
                        r1 = deltaiminus1 / deltaiplus1;
                    }
                    double limiter1 = Minbee(r1);
                    u1_B_bar[k] = grid1[i][j][k] - 0.5 * limiter1 * deltai1;
                    u1_T_bar[k] = grid1[i][j][k] + 0.5 * limiter1 * deltai1;

                    //m2
                    double deltaiminus2 = grid2[i][j][k] - grid2[i][j - 1][k];
                    double deltaiplus2 = grid2[i][j + 1][k] - grid2[i][j][k];
                    double deltai2 = 0.5 * ((1 + weight) * deltaiminus2 + (1 - weight) * deltaiplus2);
                    double r2;
                    if(deltaiminus2 == 0 && deltaiplus2 == 0){
                        r2 = 0;
                    }else if (deltaiminus2 == 0){
                        r2 = 0;
                    }else if (deltaiplus2 == 0){
                        r2 = 9999;
                    }else{
                        r2 = deltaiminus2 / deltaiplus2;
                    }
                    double limiter2 = Minbee(r2);
                    u2_B_bar[k] = grid2[i][j][k] - 0.5 * limiter2 * deltai2;
                    u2_T_bar[k] = grid2[i][j][k] + 0.5 * limiter2 * deltai2;



                }
                std::array<double, 4> flux1_B = g(u1_B_bar, Gamma1);
                std::array<double, 4> flux1_T = g(u1_T_bar, Gamma1);
                std::array<double, 4> flux2_B = g(u2_B_bar, Gamma2);
                std::array<double, 4> flux2_T = g(u2_T_bar, Gamma2);
                for(int k = 0; k < 4; k++){
                    u1_B_bar_half[i][j][k] = u1_B_bar[k] - 0.5 * dt / dy * (flux1_T[k] - flux1_B[k]);
                    u1_T_bar_half[i][j][k] = u1_T_bar[k] - 0.5 * dt / dy * (flux1_T[k] - flux1_B[k]);
                    u2_B_bar_half[i][j][k] = u2_B_bar[k] - 0.5 * dt / dy * (flux2_T[k] - flux2_B[k]);
                    u2_T_bar_half[i][j][k] = u2_T_bar[k] - 0.5 * dt / dy * (flux2_T[k] - flux2_B[k]);
                }
            }
        }
        for(int i = 0; i < nxCells + 1; i++){
            for(int j = 0; j < nyCells + 1; j++){
                flux1_y[i][j] = HLLC_flux(u1_T_bar_half[i + 2][j + 1], u1_B_bar_half[i + 2][j + 2], Gamma1, 2);
                flux2_y[i][j] = HLLC_flux(u2_T_bar_half[i + 2][j + 1], u2_B_bar_half[i + 2][j + 2], Gamma2, 2);
            }
        }
        for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
            for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
                for(int k = 0; k < 4; k++){
                    gridplus1[i][j][k] = grid1[i][j][k] - dt / dy * (flux1_y[i - 2][j - 1][k] - flux1_y[i - 2][j - 2][k]);
                    gridplus2[i][j][k] = grid2[i][j][k] - dt / dy * (flux2_y[i - 2][j - 1][k] - flux2_y[i - 2][j - 2][k]);
                }
            }
        }
        grid1 = gridplus1;
        grid2 = gridplus2;
    }
}
////////////////////////////////////////////////////////////////////////
void trans(std::array<std::array<std::array<double, 4>,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells>& grid){
    //BC
    for (int j = ghost_cells; j < ghost_cells + nyCells; j++) {
        for (int i = 0; i < ghost_cells; i++) {
            grid[i][j] = grid[ghost_cells][j]; // left
            grid[nxCells + ghost_cells + i][j] = grid[nxCells + ghost_cells - 1][j]; // right
        }
    }


    for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
        for (int j = 0; j < ghost_cells; j++) {
            grid[i][j] = grid[i][ghost_cells]; // bottom
            grid[i][nyCells + ghost_cells + j] = grid[i][nyCells + ghost_cells - 1]; // top
        }
    }
}

// phi: size [nxCells+2*ghost_cells][nyCells+2*ghost_cells]
void transPhi(
  std::array<std::array<double, nyCells + 2*ghost_cells>,nxCells + 2*ghost_cells>& phi){
    for(int j = ghost_cells; j < ghost_cells + nyCells; ++j){
        for(int i = 0; i < ghost_cells; ++i){
            phi[i][j] = phi[ghost_cells][j];
            phi[nxCells + ghost_cells + i][j] = phi[nxCells + ghost_cells - 1][j];
        }
    }
    for(int i = 0; i < nxCells + 2*ghost_cells; ++i){
        for(int j = 0; j < ghost_cells; ++j){
            phi[i][j] = phi[i][ghost_cells];
            phi[i][nyCells + ghost_cells + j] = phi[i][nyCells + ghost_cells - 1];
        }
    }
}

////////////////////////////////////////////////////////////////////////
void initializephiGrid(){
    for(int i = 0; i < nxCells+ghost_cells*2; i++){
        double x = x0 + (i-ghost_cells+0.5)*dx;
        for(int j = 0; j < nyCells+ghost_cells*2; j++){
            double y = y0_new + (j-ghost_cells+0.5)*dy;


            //helium shock bubble test
            double distox = fabs(x-0.035);
            double distoy = fabs(y-0);
            double distocentre = sqrt(distox*distox+distoy*distoy);
            phi[i][j] = 0.025-distocentre;

            // mach test
            //Initial u1
            if(x<=0.005) {
                grid1[i][j] = pri2con(w_shock,Gamma1);
                grid2[i][j] = pri2con(w_shock,Gamma1);
            }else{
                if(distocentre<=0.025){
                    grid1[i][j] = pri2con(w_helium,Gamma1);
                    grid2[i][j] = pri2con(w_helium,Gamma2);
                }
                else{
                    grid1[i][j] = pri2con(w_air,Gamma1);
                    grid2[i][j] = pri2con(w_air,Gamma2);
                }
            }
        }
    }
}

//  void initializephiGrid(){
//         for(int i = 0; i < nxCells+ghost_cells*2; i++){
//             double x = x0 + (i-ghost_cells+0.5)*dx;
//             for(int j = 0; j < nyCells+ghost_cells*2; j++){
//                 double y = y0_new + (j-ghost_cells+0.5)*dy;
//
//                 //double r = std::sqrt((x-1.0)*(x-1.0) + (y-1.0)*(y-1.0));
//
//
//                 phi[i][j] = x - 0.4;
//
//                 if(x<=0.4) {
//                     grid1[i][j] = pri2con({1.0, 0.0, 0.0, 1.0},Gamma1);
//                     grid2[i][j] = pri2con({1.0, 0.0, 0.0, 1.0},Gamma2);
//                 }else{
//                     grid1[i][j] = pri2con({0.1, 0.0 ,0.0 , 0.125},Gamma1);
//                     grid2[i][j] = pri2con({0.1, 0.0 ,0.0 , 0.125}, Gamma2);
//
//                 }
//             }
//         }
// }
////////////////////////////////////////////////////////////////////////

void compute_normal_vector(){
    //紫色
    #pragma omp parallel for
    for (int i = ghost_cells + 1; i < nxCells + ghost_cells - 1; i++) {
        //double x = x0 + (i - ghost_cells + 0.5) * dx;
        for (int j = ghost_cells + 1; j < nyCells + ghost_cells - 1; j++) {
            //double y = y0_new + (j - ghost_cells + 0.5) * dy;
            double phi_x = (phi[i+1][j] - phi[i-1][j]) / (2.0 * dx);
            double phi_y = (phi[i][j+1] - phi[i][j-1]) / (2.0 * dy);
            double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
            normal_vec[i][j][0] = phi_x  / norm;
            normal_vec[i][j][1] = phi_y  / norm;
        }
    }
    //红色下面
    #pragma omp parallel for
    for(int i = ghost_cells + 1; i < nxCells + ghost_cells - 1; i++){
        double phi_x = (phi[i+1][ghost_cells] - phi[i-1][ghost_cells]) / (2.0 * dx);
        double phi_y = (phi[i][ghost_cells + 1] - phi[i][ghost_cells]) / (dy);
        double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
        normal_vec[i][ghost_cells][0] = phi_x  / norm;
        normal_vec[i][ghost_cells][1] = phi_y  / norm;
    }
    //红色上面
    #pragma omp parallel for
    for(int i = ghost_cells + 1; i < nxCells + ghost_cells - 1; i++){
        double phi_x = (phi[i+1][nyCells + ghost_cells - 1] - phi[i-1][nyCells + ghost_cells - 1]) / (2.0 * dx);
        double phi_y = (phi[i][nyCells + ghost_cells - 1] - phi[i][nyCells + ghost_cells - 2]) / (dy);
        double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
        normal_vec[i][nyCells + ghost_cells - 1][0] = phi_x  / norm;
        normal_vec[i][nyCells + ghost_cells - 1][1] = phi_y  / norm;
    }
    //红色左边
    #pragma omp parallel for
    for(int j = ghost_cells + 1; j < nyCells + ghost_cells - 1; j++){
        double phi_x = (phi[ghost_cells+1][j] - phi[ghost_cells][j]) / (dx);
        double phi_y = (phi[ghost_cells][j + 1] - phi[ghost_cells][j - 1]) / (2 * dy);
        double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
        normal_vec[ghost_cells][j][0] = phi_x  / norm;
        normal_vec[ghost_cells][j][1] = phi_y  / norm;
    }
    //红色右边
    #pragma omp parallel for
    for(int j = ghost_cells + 1; j < nyCells + ghost_cells - 1; j++){
        double phi_x = (phi[nxCells + ghost_cells - 1][j] - phi[nxCells + ghost_cells - 1 - 1][j]) / (dx);
        double phi_y = (phi[nxCells + ghost_cells - 1][j + 1] - phi[nxCells + ghost_cells - 1][j - 1]) / (2 * dy);
        double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
        normal_vec[nxCells + ghost_cells - 1][j][0] = phi_x  / norm;
        normal_vec[nxCells + ghost_cells - 1][j][1] = phi_y  / norm;
    }
    //左下角
    double phi_x = (phi[ghost_cells+1][ghost_cells] - phi[ghost_cells][ghost_cells]) / (dx);
    double phi_y = (phi[ghost_cells][ghost_cells+1] - phi[ghost_cells][ghost_cells]) / (dy);
    double norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
    normal_vec[ghost_cells][ghost_cells][0] = phi_x  / norm;
    normal_vec[ghost_cells][ghost_cells][1] = phi_y  / norm;
    //左上角
    phi_x = (phi[ghost_cells + 1][nyCells + ghost_cells - 1] - phi[ghost_cells][nyCells + ghost_cells - 1]) / (dx);
    phi_y = (phi[ghost_cells][nyCells + ghost_cells - 1] - phi[ghost_cells][nyCells + ghost_cells - 1 - 1]) / (dy);
    norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
    normal_vec[ghost_cells][nyCells + ghost_cells - 1][0] = phi_x  / norm;
    normal_vec[ghost_cells][nyCells + ghost_cells - 1][1] = phi_y  / norm;
    //右下角
    phi_x = (phi[nxCells + ghost_cells - 1][ghost_cells] - phi[nxCells + ghost_cells - 1 - 1][ghost_cells]) / (dx);
    phi_y = (phi[nxCells + ghost_cells - 1][ghost_cells+1] - phi[nxCells + ghost_cells - 1][ghost_cells]) / (dy);
    norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
    normal_vec[nxCells + ghost_cells - 1][ghost_cells][0] = phi_x  / norm;
    normal_vec[nxCells + ghost_cells - 1][ghost_cells][1] = phi_y  / norm;
    //右上角
    phi_x = (phi[nxCells + ghost_cells - 1][nyCells + ghost_cells - 1] - phi[nxCells + ghost_cells - 1 - 1][nyCells + ghost_cells - 1]) / (dx);
    phi_y = (phi[nxCells + ghost_cells - 1][nyCells + ghost_cells - 1] - phi[nxCells + ghost_cells - 1][nyCells +  ghost_cells - 1 - 1]) / (dy);
    norm = std::sqrt(phi_x*phi_x + phi_y*phi_y);
    normal_vec[nxCells + ghost_cells - 1][nyCells + ghost_cells - 1][0] = phi_x  / norm;
    normal_vec[nxCells + ghost_cells - 1][nyCells + ghost_cells - 1][1] = phi_y  / norm;
}

double getVn(std::array<double, 2> n, std::array<double, 4> prim){
    return n[0] * prim[1] + n[1] * prim[2];
}

std::array<double, 2> getVt(std::array<double, 2> n , std::array<double, 4> prim){
    double Vn = getVn(n, prim);
    //std::vector<double> Vt(2);
    std::array<double, 2> Vt;
    Vt[0] = prim[1] - Vn * n[0];
    Vt[1] = prim[2] - Vn * n[1];
    return Vt;
}

std::array<double, 3> HLLC_L_star(std::array<double, 3> L, std::array<double, 3> R, double Gamma1, double Gamma2){
    std::array<double, 3> L_cons = pri2con_hllc(L, Gamma1);
    std::array<double, 3> R_cons = pri2con_hllc(R, Gamma2);
    double cs_l = calcs(L[2], L[0], Gamma1);
    double cs_R = calcs(R[2], R[0], Gamma2);
    double S_R = std::max(fabs(L[1]) + cs_l, fabs(R[1]) + cs_R);
    double S_L = - S_R;
    double S_star = (R[2] - L[2] + L[0] * L[1] * (S_L - L[1]) - R[0] * R[1] * (S_R - R[1]))/(L[0] * (S_L - L[1]) - R[0] * (S_R - R[1]));

    double temp = L[0] * (S_L - L[1])/(S_L - S_star);
    std::array<double, 3> HLLC_L_cons;
    HLLC_L_cons[0] = temp;
    HLLC_L_cons[1] = temp*S_star;
    HLLC_L_cons[2] = temp * (L_cons[2]/L[0] + (S_star - L[1]) * (S_star + L[2]/(L[0] * (S_L - L[1]))));
    std::array<double, 3> HLLC_L = con2pri_hllc(HLLC_L_cons, Gamma1);
    return HLLC_L;
}
std::array<double, 3> HLLC_R_star(std::array<double, 3> L, std::array<double, 3> R, double Gamma1, double Gamma2){
    std::array<double, 3> L_cons = pri2con_hllc(L, Gamma1);
    std::array<double, 3> R_cons = pri2con_hllc(R, Gamma2);
    double cs_l = calcs(L[2], L[0], Gamma1);
    double cs_R = calcs(R[2], R[0], Gamma2);
    double S_R = std::max(fabs(L[1]) + cs_l, fabs(R[1]) + cs_R);
    double S_L = - S_R;
    double S_star = (R[2] - L[2] + L[0] * L[1] * (S_L - L[1]) - R[0] * R[1] * (S_R - R[1]))/(L[0] * (S_L - L[1]) - R[0] * (S_R - R[1]));
    double temp = R[0] * (S_R - R[1])/(S_R - S_star);
    std::array<double, 3> HLLC_R_cons;
    HLLC_R_cons[0] = temp;
    HLLC_R_cons[1] = temp*S_star;
    HLLC_R_cons[2] = temp * (R_cons[2]/R[0] + (S_star - R[1]) * (S_star + R[2]/(R[0] * (S_R - R[1]))));
    std::array<double, 3> HLLC_R = con2pri_hllc(HLLC_R_cons, Gamma2);
    return HLLC_R;
}

std::array<double, 4> normal2xy(std::array<double, 3> L_star, std::array<double, 2> n, std::array<double, 2> Vt){
    std::array<double, 4> ans;
    ans[rho] = L_star[0];
    ans[u] = Vt[0] + n[0] * L_star[1];
    ans[v] = Vt[1] + n[1] * L_star[1];
    ans[pre] = L_star[2];
    return ans;

}


double getQijk(double dx, double dy, double Qx, double Qy, std::array<double, 2> n){
    double up = dy * fabs(n[0]) * Qx + dx * fabs(n[1]) * Qy;
    double bottom = dy * fabs(n[0]) + dx * fabs(n[1]);
    if(n[0] == 0){
        return Qy;
    }else if(n[1] == 0){
        return Qx;
    }else{
        return up / bottom ;
    }

}

void constant_extrap(std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &grid1,
                     std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &grid2,
                     const std::array<std::array<std::array<double, 2>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &n,
                     const std::array<std::array<int,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> &interface,
                     const std::array<std::array<double,nyCells + 2 * ghost_cells>,nxCells + 2 * ghost_cells> &phi,
                     int niter = 10){
    for(int i = 2; i < ghost_cells + nxCells; i++){
        for( int j = 2; j < ghost_cells + nyCells; j++){
            if(interface[i][j] == 0){
                if(phi[i][j] >= 0){//grid1的ghost cells
                    grid1[i][j] = {1e9, 1e9, 1e9, 1e9};
                }else{//grid2的ghost cells
                    grid2[i][j] = {1e9, 1e9, 1e9, 1e9};
                }
            }
        }
    }

    for(int iter = 0; iter < niter; iter++){
        //l down corner
        for(int i = 2; i < ghost_cells + nxCells; i++){
            for(int j = 2; j < nyCells + ghost_cells; j++){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        std::array<double, 4> Qx = {};
                        std::array<double, 4> Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid1[i-1][j];
                        }else{
                            Qx = grid1[i+1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid1[i][j-1];
                        }else{
                            Qy = grid1[i][j+1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid1[i][j][k])){
                                grid1[i][j][k] = Qijk;
                            }
                        }
                    }else{
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid2[i+1][j];
                        }else{
                            Qx = grid2[i-1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid2[i][j+1];
                        }else{
                            Qy = grid2[i][j-1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid2[i][j][k])){
                                grid2[i][j][k] = Qijk;
                            }
                        }

                    }
                }

            }
        }
        //r up corner
        for(int i = ghost_cells + nxCells - 1; i > 1; i--){
            for(int j = ghost_cells + nyCells - 1; j > 1; j--){
            if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid1[i-1][j];
                        }else{
                            Qx = grid1[i+1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid1[i][j-1];
                        }else{
                            Qy = grid1[i][j+1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid1[i][j][k])){
                                grid1[i][j][k] = Qijk;
                            }
                        }
                    }else{
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid2[i+1][j];
                        }else{
                            Qx = grid2[i-1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid2[i][j+1];
                        }else{
                            Qy = grid2[i][j-1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid2[i][j][k])){
                                grid2[i][j][k] = Qijk;
                            }
                        }

                    }
                }

            }
        }
        //r down corner
        for(int i = ghost_cells + nxCells - 1; i > 1; i--){
            for(int j = 2; j < ghost_cells + nyCells; j++){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid1[i-1][j];
                        }else{
                            Qx = grid1[i+1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid1[i][j-1];
                        }else{
                            Qy = grid1[i][j+1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid1[i][j][k])){
                                grid1[i][j][k] = Qijk;
                            }
                        }
                    }else{
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid2[i+1][j];
                        }else{
                            Qx = grid2[i-1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid2[i][j+1];
                        }else{
                            Qy = grid2[i][j-1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid2[i][j][k])){
                                grid2[i][j][k] = Qijk;
                            }
                        }

                    }
                }

            }
        }
        //l up corner
        for(int i = 2; i < ghost_cells + nxCells; i++){
            for(int j = ghost_cells + nyCells - 1; j > 1; j--){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid1[i-1][j];
                        }else{
                            Qx = grid1[i+1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid1[i][j-1];
                        }else{
                            Qy = grid1[i][j+1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid1[i][j][k])){
                                grid1[i][j][k] = Qijk;
                            }
                        }
                    }else{
                        std::array<double, 4>  Qx = {};
                        std::array<double, 4>  Qy = {};
                        if(n[i][j][0] > 0){
                            Qx = grid2[i+1][j];
                        }else{
                            Qx = grid2[i-1][j];
                        }
                        if(n[i][j][1] > 0){
                            Qy = grid2[i][j+1];
                        }else{
                            Qy = grid2[i][j-1];
                        }
                        for(int k = 0; k < 4; k++){
                            double Qijk = getQijk(dx, dy, Qx[k], Qy[k], n[i][j]);
                            if(fabs(Qijk) < fabs(grid2[i][j][k])){
                                grid2[i][j][k] = Qijk;
                            }
                        }

                }
            }


        }
    }
    }
}


void advectPhi(std::array<std::array<double, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &phi,
            std::array<std::array<double, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &phiplus ,
            const std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &grid1,
            const std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> &grid2,
            double dt, double dx, double dy, double Gamma1, double Gamma2) {

    // 更新内部点
    #pragma omp parallel for
    for (int i = ghost_cells; i < ghost_cells + nxCells; i++) {
        for (int j = ghost_cells; j < ghost_cells + nyCells; j++) {
            double u_val, v_val;
            if(phi[i][j] <= 0){
                std::array<double, 4> prim1 = con2pri(grid1[i][j], Gamma1);
                u_val = prim1[u];
                v_val = prim1[v];
            }else{
                std::array<double, 4> prim2 = con2pri(grid2[i][j], Gamma2);
                u_val = prim2[u];
                v_val = prim2[v];
            }

            // 计算 φ_x
            double phi_x;
            if (u_val > 0)
                phi_x = (phi[i][j] - phi[i-1][j]) / dx;
            else
                phi_x = (phi[i+1][j] - phi[i][j]) / dx;

            // 计算 φ_y
            double phi_y;
            if (v_val > 0)
                phi_y = (phi[i][j] - phi[i][j-1]) / dy;
            else
                phi_y = (phi[i][j+1] - phi[i][j]) / dy;

            // 更新 φ 前向
            phiplus[i][j] = phi[i][j] - dt * (u_val * phi_x + v_val * phi_y);
        }
    }
    phi = phiplus;
}




std::array<double, 2> getphiijk(double dx, double dy, double phix, double phiy){
    double a = dy * dy + dx * dx;
    double b = -2 * (phix * dy * dy + phiy * dx * dx);
    double c = phix * phix * dy * dy + phiy * phiy * dx * dx - dx * dx * dy * dy;
    double delta = b * b - 4 * a * c;

    if(delta < 0){
        if(fabs(phix) > fabs(phiy)){
            a = 1;
            b = -2 * phiy;
            c = phiy * phiy - dy * dy;
            delta = b * b - 4 * a * c;

        }else{
            a = 1;
            b = -2 * phix;
            c = phix * phix - dx * dx;
            delta = b * b - 4 * a * c;

        }
    }
    if(delta < 0){
        return {0, 0};
    }

    double phiijk;
    if(phix > 0){
        phiijk = (-b + sqrt(delta)) / (2 * a);
    }else{
        phiijk = (-b - sqrt(delta)) / (2 * a);
    }
    return {1, phiijk};
}

void reinitialisation( int niter = 1){
    for(int i = 2; i < ghost_cells + nxCells; i++){
        for( int j = 2; j < ghost_cells + nyCells; j++){
            if(interface[i][j] == 0){
                if(phi[i][j] >= 0){//grid1的ghost cells
                    phi[i][j] = 999;
                }else{//grid2的ghost cells
                    phi[i][j] = -999;
                }
            }
        }
    }

    for(int iter = 0; iter < niter; iter++){
        //l down corner
        for(int i = 2; i < ghost_cells + nxCells; i++){
            for(int j = 2; j < nyCells + ghost_cells; j++){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        double phix;
                        double phiy;

                        phix = std::min(phi[i-1][j], phi[i+1][j]);
                        phiy = std::min(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }

                    }else{
                        double phix;
                        double phiy;

                        phix = std::max(phi[i-1][j], phi[i+1][j]);
                        phiy = std::max(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }

                    }
                }
            }

        }

        //r up corner
        for(int i = ghost_cells + nxCells - 1; i > 1; i--){
            for(int j = ghost_cells + nyCells - 1; j > 1; j--){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        double phix;
                        double phiy;

                        phix = std::min(phi[i-1][j], phi[i+1][j]);
                        phiy = std::min(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }
                    }else{
                        double phix;
                        double phiy;

                        phix = std::max(phi[i-1][j], phi[i+1][j]);
                        phiy = std::max(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }

                    }
                }
            }
        }
        //r down corner
        for(int i = ghost_cells + nxCells - 1; i > 1; i--){
            for(int j = 2; j < ghost_cells + nyCells; j++){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        double phix;
                        double phiy;

                        phix = std::min(phi[i-1][j], phi[i+1][j]);
                        phiy = std::min(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }
                    }else{
                        double phix;
                        double phiy;

                        phix = std::max(phi[i-1][j], phi[i+1][j]);
                        phiy = std::max(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }

                    }
                }
            }
        }
        //l up corner
        for(int i = 2; i < ghost_cells + nxCells; i++){
            for(int j = ghost_cells + nyCells - 1; j > 1; j--){
                if(interface[i][j] == 0){
                    if(phi[i][j] >= 0){
                        double phix;
                        double phiy;

                        phix = std::min(phi[i-1][j], phi[i+1][j]);
                        phiy = std::min(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }
                    }else{
                        double phix;
                        double phiy;

                        phix = std::max(phi[i-1][j], phi[i+1][j]);
                        phiy = std::max(phi[i][j-1], phi[i][j+1]);

                        std::array<double, 2> phiijk = getphiijk(dx, dy, phix, phiy);
                        if(phiijk[0] == 1){//不是ill defined
                            if(fabs(phiijk[1]) < fabs(phi[i][j])){
                                phi[i][j] = phiijk[1];
                            }
                        }

                    }

                }
            }
        }
    }
}





void store_data1(){
    std::ofstream grid1_file("grid1.dat");
        for(int i = 0; i < nxCells + 2*ghost_cells; i++) {
            double x = x0 + (i - ghost_cells + 0.5) * dx;
            for(int j = 0; j < nyCells + 2*ghost_cells; j++) {
                double y = y0_new + (j - ghost_cells + 0.5) * dy;
                grid1_file << i << " " << j << " " << x << " " << y << " "
                        << grid1[i][j][rho] << " "
                        << phi[i][j] << " " << interface[i][j] << std::endl;

            }
        }

    std::ofstream grid2_file("grid2.dat");// Print grid2 values
    // std::cout << "\nGrid2 values:" << std::endl;
    for(int i = 0; i < nxCells +  2*ghost_cells; i++) {
        double x = x0 + (i - ghost_cells + 0.5) * dx;
        for(int j = 0; j < nyCells + 2*ghost_cells; j++) {
            double y = y0_new + (j - ghost_cells + 0.5) * dy;
            grid2_file << i << " " << j << " " << x << " " << y << " "
                    << grid2[i][j][rho] << " "
                    << phi[i][j] << " " << interface[i][j] << std::endl;

        }
    }
}

void store_data2(){
    std::ofstream grid1_file("grid.dat");
    grid1_file<<std::fixed<<std::setprecision(16);
    for(int i = 0; i < nxCells + 2*ghost_cells; i++) {
        double x = x0 + (i - ghost_cells + 0.5) * dx;
        for(int j = 0; j < nyCells + 2*ghost_cells; j++) {
            double y = y0_new + (j - ghost_cells + 0.5) * dy;
            if (phi[i][j] <= 0) {
                std::array<double, 4> prim1 = con2pri(grid1[i][j], Gamma1);
                double internal_ene = prim1[pre] / ((Gamma1 - 1) * prim1[rho]);
                grid1_file << x << " " << y << " " << prim1[rho] << " " << prim1[u] << " "
                << prim1[v] << " " << prim1[pre]<< " " << internal_ene << " "<< phi[i][j] << " " << std::endl;
            }else {
                std::array<double, 4> prim2 = con2pri(grid2[i][j], Gamma1);
                double internal_ene = prim2[pre] / ((Gamma1 - 1) * prim2[rho]);
                grid1_file << x << " " << y << " " << prim2[rho] << " " << prim2[u] << " "
                << prim2[v] << " " << prim2[pre] << " " << internal_ene << " "<< phi[i][j] << " " << std::endl;
            }
        }
    }
}




int main(){
    omp_set_num_threads(16);
    //initialize interface
    initializephiGrid();
    // store_data2();
    // exit(0);
    //find interface
    #pragma omp parallel for
    for(int i = ghost_cells; i < nxCells +  ghost_cells; i++){
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                double right_nei = phi[i][j] * phi[i + 1][j];
                double left_nei = phi[i][j] * phi[i - 1][j];
                double top_nei = phi[i][j] * phi[i][j + 1];
                double bot_nei = phi[i][j] * phi[i][j - 1];
                if (right_nei < 0.0||left_nei < 0.0||top_nei < 0.0||bot_nei < 0.0){
                    interface[i][j] = 1;
                }else{
                    interface[i][j] = 0;
                }
            }
    }//已看出0和1

    // for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
    //     for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
    //         std::cout << interface[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // exit(0);
    //find normal vector
    compute_normal_vector();

    //  都对



    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> orig1;
    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> orig2;
    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> tempA1;
    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> tempA2;
    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> tempB1;
    static std::array<std::array<std::array<double, 4>, nyCells + 2 * ghost_cells>, nxCells + 2 * ghost_cells> tempB2;

    int iter = 0; double t = 0.0;
    double t0 = omp_get_wtime();
    do {




        //riemann gfm
        #pragma omp parallel for
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++){
            double x = x0 + (i - ghost_cells + 0.5) * dx;
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                double y = y0_new + (j - ghost_cells + 0.5) * dy;

                if(interface[i][j] == 1){

                    //left 双线性插值的目的是在已有的网格数据上估计出某个非整数网格位置处的流体状态
                    std::array<double,2> p = {x - phi[i][j] * normal_vec[i][j][0], y - phi[i][j] * normal_vec[i][j][1]};
                    std::array<double,2> p1 = {p[0] + 1.5 * dx * normal_vec[i][j][0], p[1] + 1.5 * dy * normal_vec[i][j][1]};
                    std::array<double,2> p2 = {p[0] - 1.5 * dx * normal_vec[i][j][0], p[1] - 1.5 * dy * normal_vec[i][j][1]};
                    std::array<double,2> p1_ij = {((p1[0] - x0)/dx) + ghost_cells - 0.5, ((p1[1] - y0_new)/dy) + ghost_cells - 0.5};
                    std::array<double,2> p2_ij = {((p2[0] - x0)/dx) + ghost_cells - 0.5, ((p2[1] - y0_new)/dy) + ghost_cells - 0.5};



                    ///////////////////////////////////插p2
                    int temp1_p2 = std::floor(p2_ij[0]);
                    int temp2_p2 = std::floor(p2_ij[1]);
                    std::vector<int> Q11_p2 = {temp1_p2, temp2_p2};


                    temp1_p2 = std::floor(p2_ij[0]);
                    temp2_p2 = std::ceil(p2_ij[1]);
                    std::vector<int> Q12_p2 = {temp1_p2, temp2_p2};
                    // std::vector<int> Q12_p2 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p2)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p2))
                    // };

                    temp1_p2 = std::ceil(p2_ij[0]);
                    temp2_p2 = std::floor(p2_ij[1]);
                    std::vector<int> Q21_p2 = {temp1_p2, temp2_p2};
                    // std::vector<int> Q21_p2 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p2)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p2))
                    // };

                    temp1_p2 = std::ceil(p2_ij[0]);
                    temp2_p2 = std::ceil(p2_ij[1]);
                    std::vector<int> Q22_p2 = {temp1_p2, temp2_p2};
                    // std::vector<int> Q22_p2 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p2)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p2))
                    // };

                    double x1_p2 = x0 + (Q11_p2[0] - ghost_cells + 0.5) * dx;
                    double x2_p2 = x0 + (Q21_p2[0] - ghost_cells + 0.5) * dx;
                    double y1_p2 = y0_new + (Q11_p2[1] - ghost_cells + 0.5) * dy;
                    double y2_p2 = y0_new + (Q12_p2[1] - ghost_cells + 0.5) * dy;
                    std::array<double, 4> L_state;
                    //当他刚好在一个点上
                    if(p2_ij[0] == std::floor(p2_ij[0]) && p2_ij[1] == std::floor(p2_ij[1])){//if int, floor = ceil
                        L_state[rho]  = grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][0];
                        L_state[rhou] = grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][1];
                        L_state[rhov] = grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][2];
                        L_state[Ene]    = grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][3];

                    }else if(p2_ij[0] == std::floor(p2_ij[0])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (y2_p2 - p2[1])/(y2_p2 - y1_p2) * grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][k] + (p2[1] - y1_p2)/(y2_p2 - y1_p2) * grid1[std::floor(p2_ij[0])][std::ceil(p2_ij[1])][k];
                            L_state[k] = R3;
                        }
                    }else if(p2_ij[1] == std::floor(p2_ij[1])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (x2_p2 - p2[0])/(x2_p2 - x1_p2) * grid1[std::floor(p2_ij[0])][std::floor(p2_ij[1])][k] + (p2[0] - x1_p2)/(x2_p2 - x1_p2) * grid1[std::ceil(p2_ij[0])][std::floor(p2_ij[1])][k];
                            L_state[k] = R3;
                        }
                    }else{
                        for(int k = 0; k < 4; k++){
                            double R1 = (x2_p2 - p2[0])/(x2_p2 - x1_p2) * grid1[Q11_p2[0]][Q11_p2[1]][k] + (p2[0] - x1_p2)/(x2_p2 - x1_p2) * grid1[Q21_p2[0]][Q21_p2[1]][k];
                            double R2 = (x2_p2 - p2[0])/(x2_p2 - x1_p2) * grid1[Q12_p2[0]][Q12_p2[1]][k] + (p2[0] - x1_p2)/(x2_p2 - x1_p2) * grid1[Q22_p2[0]][Q22_p2[1]][k];
                            double R3 = (y2_p2 - p2[1])/(y2_p2 - y1_p2) * R1 + (p2[1] - y1_p2)/(y2_p2 - y1_p2) * R2;
                            L_state[k] = R3;
                        }
                    }








                    ///////////////////////////////////叉p1
                    int temp1_p1 = std::floor(p1_ij[0]);
                    int temp2_p1 = std::floor(p1_ij[1]);
                    std::vector<int> Q11_p1 = {temp1_p1, temp2_p1};
                    // std::vector<int> Q11_p1 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p1)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p1))
                    // };
                    temp1_p1 = std::floor(p1_ij[0]);
                    temp2_p1 = std::ceil(p1_ij[1]);
                    std::vector<int> Q12_p1 = {temp1_p1, temp2_p1};
                    // std::vector<int> Q12_p1 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p1)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p1))
                    // };
                    temp1_p1 = std::ceil(p1_ij[0]);
                    temp2_p1 = std::floor(p1_ij[1]);
                    std::vector<int> Q21_p1 = {temp1_p1, temp2_p1};
                    // std::vector<int> Q21_p1 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p1)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p1))
                    // };
                    temp1_p1 = std::ceil(p1_ij[0]);
                    temp2_p1 = std::ceil(p1_ij[1]);
                    std::vector<int> Q22_p1 = {temp1_p1, temp2_p1};
                    // std::vector<int> Q22_p1 = {
                    //     std::max(0, std::min(nxCells + 2 * ghost_cells - 1, temp1_p1)),
                    //     std::max(0, std::min(nyCells + 2 * ghost_cells - 1, temp2_p1))
                    // };

                    double x1_p1 = x0 + (Q11_p1[0] - ghost_cells + 0.5) * dx;
                    double x2_p1 = x0 + (Q21_p1[0] - ghost_cells + 0.5) * dx;
                    double y1_p1 = y0_new + (Q11_p1[1] - ghost_cells + 0.5) * dy;
                    double y2_p1 = y0_new + (Q12_p1[1] - ghost_cells + 0.5) * dy;
                    std::array<double, 4> R_state;
                    //当他刚好在一个点上
                    if(p1_ij[0] == std::floor(p1_ij[0]) && p1_ij[1] == std::floor(p1_ij[1])){//if int, floor = ceil
                        R_state[rho]  = grid2[std::floor(p1_ij[0])][std::floor(p1_ij[1])][0];
                        R_state[rhou] = grid2[std::floor(p1_ij[0])][std::floor(p1_ij[1])][1];
                        R_state[rhov] = grid2[std::floor(p1_ij[0])][std::floor(p1_ij[1])][2];
                        R_state[Ene]    = grid2[std::floor(p1_ij[0])][std::floor(p1_ij[1])][3];

                    }else if(p1_ij[0] == std::floor(p1_ij[0])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (y2_p1 - p1[1])/(y2_p1 - y1_p1) * grid2 [std::floor(p1_ij[0])][std::floor(p1_ij[1])][k] + (p1[1] - y1_p1)/(y2_p1 - y1_p1) * grid2[std::floor(p1_ij[0])][std::ceil(p1_ij[1])][k];
                            R_state[k] = R3;
                        }
                    }else if(p1_ij[1] == std::floor(p1_ij[1])){
                        for(int k = 0; k < 4; k++){
                            double R3 = (x2_p1 - p1[0])/(x2_p1 - x1_p1) * grid2[std::floor(p1_ij[0])][std::floor(p1_ij[1])][k] + (p1[0] - x1_p1)/(x2_p1 - x1_p1) * grid2[std::ceil(p1_ij[0])][std::floor(p1_ij[1])][k];
                            R_state[k] = R3;
                        }
                    }else{
                        for(int k = 0; k < 4; k++){
                            double R1 = (x2_p1 - p1[0])/(x2_p1 - x1_p1) * grid2[Q11_p1[0]][Q11_p1[1]][k] + (p1[0] - x1_p1)/(x2_p1 - x1_p1) * grid2[Q21_p1[0]][Q21_p1[1]][k];
                            double R2 = (x2_p1 - p1[0])/(x2_p1 - x1_p1) * grid2[Q12_p1[0]][Q12_p1[1]][k] + (p1[0] - x1_p1)/(x2_p1 - x1_p1) * grid2[Q22_p1[0]][Q22_p1[1]][k];
                            double R3 = (y2_p1 - p1[1])/(y2_p1 - y1_p1) * R1 + (p1[1] - y1_p1)/(y2_p1 - y1_p1) * R2;
                            R_state[k] = R3;
                        }
                    }

                    if(phi[i][j] >= 0){
                        std::array<double, 4> L_state_pri = con2pri(L_state, Gamma1);
                        std::array<double, 4> R_state_pri = con2pri(R_state, Gamma2);
                        std::array<double, 3> L_state_pri_normal = {L_state_pri[0], getVn(normal_vec[i][j], L_state_pri), L_state_pri[3]};
                        std::array<double, 3> R_state_pri_normal = {R_state_pri[0], getVn(normal_vec[i][j], R_state_pri), R_state_pri[3]};
                        std::array<double, 2> Vt = getVt(normal_vec[i][j], L_state_pri);
                        //std::vector<double> R_state_pri_normal = get_R_state(L_state_pri_normal, normal_vec[i][j], v_rigid);

                        std::array<double, 3> L_star_pri_normal = HLLC_L_star(L_state_pri_normal, R_state_pri_normal, Gamma1, Gamma2);
                        std::array<double, 4> L_star_pri = normal2xy(L_star_pri_normal, normal_vec[i][j], Vt);
                        std::array<double, 4> L_star_con = pri2con(L_star_pri, Gamma1);

                        grid1[i][j] = L_star_con;
                    }else{
                        std::array<double, 4> R_state_pri = con2pri(R_state, Gamma2);
                        std::array<double, 4> L_state_pri = con2pri(L_state, Gamma1);
                        std::array<double, 3> L_state_pri_normal = {L_state_pri[0], getVn(normal_vec[i][j], L_state_pri), L_state_pri[3]};
                        std::array<double, 3> R_state_pri_normal = {R_state_pri[0], getVn(normal_vec[i][j], R_state_pri), R_state_pri[3]};
                        //std::vector<double> R_state_pri_normal = {R_state_pri[0], getVn(normal_vec[i][j], R_state_pri), R_state_pri[3]};
                        std::array<double, 2> Vt = getVt(normal_vec[i][j], R_state_pri);
                        //state_pri_normal = get_R_state(R_state_pri_normal, normal_vec[i][j], v_rigid);

                        std::array<double, 3> R_star_pri_normal = HLLC_R_star(L_state_pri_normal, R_state_pri_normal, Gamma1, Gamma2);
                        std::array<double, 4> R_star_pri = normal2xy(R_star_pri_normal, normal_vec[i][j], Vt);
                        std::array<double, 4> R_star_con = pri2con(R_star_pri, Gamma2);

                        grid2[i][j] = R_star_con;
                    }
                }
            }
        }





        constant_extrap(grid1, grid2, normal_vec, interface, phi, 5);
        trans(grid1);
        trans(grid2);






        // Print grid1 values
        iter += 1;
        dt = updatedt(grid1, grid2);
        t += dt;
        if(t > T_end){
            t = T_end;
        }
        std::cout<< std::round(t/T_end * 10000)/10000 * 100 << "%" << std::endl;




        advectPhi(phi, phiplus, grid1, grid2, dt, dx, dy, Gamma1, Gamma2);
        //transPhi(phi);



        // store_data();


        #pragma omp parallel for
        for(int i = ghost_cells; i < nxCells +  ghost_cells; i++){
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                double right_nei = phi[i][j] * phi[i + 1][j];
                double left_nei = phi[i][j] * phi[i - 1][j];
                double top_nei = phi[i][j] * phi[i][j + 1];
                double bot_nei = phi[i][j] * phi[i][j - 1];
                if (right_nei < 0.0||left_nei < 0.0||top_nei < 0.0||bot_nei < 0.0){
                    interface[i][j] = 1;
                }else{
                    interface[i][j] = 0;
                }
            }
        }



        reinitialisation(1);

        compute_normal_vector();






        // dt = 0.5*dt;
        // evolve(dt, 'x');
        // trans(grid1);
        // trans(grid2);
        // dt = 2*dt;
        // evolve(dt, 'y');
        // trans(grid1);
        // trans(grid2);
        // dt = 0.5*dt;
        // evolve(dt, 'x');
        // trans(grid1);
        // trans(grid2);
        // dt = 2*dt;




        trans(grid1); trans(grid2);
        orig1 = grid1;
        orig2 = grid2;
        evolve(dt, 'x');
        trans(grid1); trans(grid2);
        evolve(dt, 'y');
        tempA1 = grid1;
        tempA2 = grid2;

        grid1 = orig1;
        grid2 = orig2;

        evolve(dt, 'y');
        trans(grid1); trans(grid2);
        evolve(dt, 'x');
        tempB1 = grid1;
        tempB2 = grid2;
        for (int j = ghost_cells; j < nyCells+ghost_cells; ++j) {
            for (int i = ghost_cells; i < nxCells+ghost_cells; ++i) {
                for (int k = 0; k < 4; ++k) {
                    grid1[i][j][k] = 0.5*(tempA1[i][j][k] + tempB1[i][j][k]);
                    grid2[i][j][k] = 0.5*(tempA2[i][j][k] + tempB2[i][j][k]);
                }
            }
        }








    }while(t < T_end);
    double t1 = omp_get_wtime();

    double elapsed = t1 - t0;
    std::printf("Elapsed wall‐clock time: %.3f s\n", elapsed);



    store_data2();


    return 0;
}


