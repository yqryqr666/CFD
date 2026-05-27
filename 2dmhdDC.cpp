#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <cmath>

struct pvar{
    double rho, v_x, v_y, v_z, p, B_x, B_y, B_z, psi;
};

struct cvar{
    double rho, rhov_x, rhov_y, rhov_z, U, B_x, B_y, B_z, psi;
    double values[9];
};
//1d to 2d
double Gamma = 5.0 / 3.0; 
double CFL = 0.9;
double x0 = 0, x1 = 1, y0_new = 0, y1_new = 1; 
double T_end = 0.5;
int nxCells = 256, nyCells = 256; 
int nxPoints = nxCells + 1, nyPoints = nyCells + 1; 
int ghost_cells = 3;
double dx = (x1 - x0)/(nxPoints - 2), dy = (y1_new - y0_new)/(nyPoints - 2); 

void initializeGrid(std::vector<pvar> & grid);
pvar con2pri(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y, double B_z, double Gamma);
cvar pri2con(double rho, double v_x, double v_y, double v_z, double p , double B_x, double B_y, double B_z, double Gamma);
std::vector<double> f(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y,  double B_z, double psi, double p, double Gamma, double c_h);
std::vector<double> g(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y,  double B_z, double psi, double p, double Gamma, double c_h);
double Minbee(double r);
//double vanLeer(double r);
std::vector<double> fri_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double c_h, char dir);
std::vector<double> RI_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double Gamma, double c_h, char dir);
std::vector<double> FORCE(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double Gamma, char dir);
void data_reconstruction(double uL, double u, double uR, double &uL_bar, double &uR_bar);
void half_time(cvar &uL_half, cvar &uR_half, cvar &uL, cvar& uR, double dt, double dx, double dy, double Gamma, double c_h, char dir);
void SLIC(std::vector<std::vector<std::vector<double> > >& flux, const std::vector<std::vector<cvar>>& grid, double dt, double dx, double dy, double Gamma, double c_h, char dir);
double compute_Cf(double rho, double p, double B_x, double B_y, double B_z, double Gamma);
double update(const std::vector<std::vector<cvar> >& grid, double dx,double dy, double Gamma, double CFL);
void trans(std::vector<std::vector<cvar> >& grid);
void periodic(std::vector<std::vector<cvar> >& grid);
double ch(std::vector<std::vector<cvar> >& grid);
double compute_cp(double c_h);
void compute_source_term(std::vector<std::vector<cvar>>& source,  std::vector<std::vector<cvar>>& grid);
void update_source(std::vector<std::vector<cvar>> &source, std::vector<std::vector<cvar>> &u, std::vector<std::vector<cvar>> &u_con, double c_h, double c_p, double dt);


void initializeGrid(std::vector<std::vector<pvar> >& grid){
    for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
            double x = x0 + (i - ghost_cells + 0.5) * dx;
            double y = y0_new + (j - ghost_cells + 0.5) * dy;
//          if (y < 0.5){
//              grid[i][j] = (pvar){1.0, 0.0, 0.0, 1};
//          }else{
//              grid[i][j] = (pvar){0.125, 0.0, 0.0, 0.1};
//          }
            
            // if((x - 1) * (x - 1) + (y - 1) * (y - 1) < 0.4 * 0.4){
            //     grid[i][j] = (pvar){Gamma * Gamma, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
            // }else{
            //     grid[i][j] = (pvar){0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};
            // }
            grid[i][j] = (pvar){Gamma * Gamma, -sin(2 * M_PI * y), sin(2 * M_PI * x), 0,Gamma, -sin(2 * M_PI * y), sin(4 * M_PI * x), 0.0};
        }
    }
}

pvar con2pri(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y, double B_z, double psi, double Gamma){
    pvar pri;
    pri.rho = rho;
    pri.v_x = rhov_x / rho;
    pri.v_y = rhov_y / rho;
    pri.v_z = rhov_z / rho;
    double v_square = (pri.v_x * pri.v_x + pri.v_y * pri.v_y + pri.v_z * pri.v_z);
    double B_square = (B_x * B_x + B_y * B_y + B_z * B_z);
    pri.p = (Gamma - 1) * (U - 0.5 * rho * v_square - 0.5 * B_square);
    pri.B_x = B_x;
    pri.B_y = B_y;
    pri.B_z = B_z;
    pri.psi = psi;
    return pri;
}

cvar pri2con(double rho, double v_x, double v_y, double v_z, double p , double B_x, double B_y, double B_z, double psi, double Gamma){
    cvar con;
    con.rho = rho;
    con.rhov_x = rho * v_x;
    con.rhov_y = rho * v_y;
    con.rhov_z = rho * v_z;
    double v_square = (v_x * v_x + v_y * v_y + v_z * v_z);
    double B_square = (B_x * B_x + B_y * B_y + B_z * B_z);
    con.U = p / (Gamma - 1) + 0.5 * rho * v_square + 0.5 * B_square;
    con.B_x = B_x;
    con.B_y = B_y;
    con.B_z = B_z;
    con.psi = psi;
    return con;
}

std::vector<double> f(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y,  double B_z, double psi, double p, double Gamma, double c_h){
    std::vector<double> flux(9);
    double v_x = rhov_x / rho;
    double v_y = rhov_y / rho;
    double v_z = rhov_z / rho;

    double B_square = B_x * B_x + B_y * B_y + B_z * B_z;
    double v_dot_B = v_x * B_x + v_y * B_y + v_z * B_z;
    flux[0] = rhov_x;
    flux[1] = rho * v_x * v_x + p + 0.5 * B_square - B_x * B_x;
    flux[2] = rho * v_x * v_y - B_x * B_y;
    flux[3] = rho * v_x * v_z - B_x * B_z;
    flux[4] = (U + p + 0.5 * B_square) * v_x - v_dot_B * B_x;
    flux[5] = psi;
    flux[6] = B_y * v_x - B_x * v_y;
    flux[7] = B_z * v_x - B_x * v_z;
    flux[8] = c_h * c_h * B_x;
    return flux;
}

std::vector<double> g(double rho, double rhov_x, double rhov_y, double rhov_z, double U, double B_x, double B_y,  double B_z, double psi, double p, double Gamma, double c_h){
    std::vector<double> glux(9);
    double v_x = rhov_x / rho;
    double v_y = rhov_y / rho;
    double v_z = rhov_z / rho;

    double B_square = B_x * B_x + B_y * B_y + B_z * B_z;
    double v_dot_B = v_x * B_x + v_y * B_y + v_z * B_z;
    glux[0] = rho * v_y;
    glux[1] = rho * v_x * v_y - B_x * B_y;
    glux[2] = rho * v_y * v_y + p + 0.5 * B_square - B_y * B_y;
    glux[3] = rho * v_y * v_z - B_y * B_z;
    glux[4] = (U + p + 0.5 * B_square) * v_y - B_y * v_dot_B;
    glux[5] = B_x * v_y - B_y * v_x;
    glux[6] = psi;
    glux[7] = B_z * v_y - B_y * v_z;
    glux[8] = c_h * c_h * B_y; 
    return glux;
}

double Minbee(double r){
    if(r <= 0) return 0;
    if(r <= 1) return r;
    return std::min(1.0, 2.0 / (1 + r));
 }
// double vanLeer(double r) {
//     return (r + std::fabs(r)) / (1.0 + std::fabs(r));
// }



std::vector<double> fri_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double c_h, char dir){
    std::vector<double> flux_f(9);
    
    std::vector<double> FL(9), FR(9);
    double delta = (dir == 'x') ? (dx / dt) : (dy / dt);
    
    if (dir == 'x') {
        FL = f(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, pL, Gamma, c_h);
        FR = f(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, pR, Gamma, c_h);
    } else if (dir == 'y') {
        FL = g(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, pL, Gamma, c_h);
        FR = g(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, pR, Gamma, c_h);
    }
    
    std::vector<double> uL_values = {uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi};
    std::vector<double> uR_values = {uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi};
    
    for (size_t i = 0; i < 9; ++i) {
        flux_f[i] = 0.5 * (FL[i] + FR[i]) - 0.5 * delta * (uR_values[i] - uL_values[i]);
    }
    
    return flux_f;
}


std::vector<double> RI_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double Gamma, double c_h, char dir){
    std::vector<double> flux_r;
    
    std::vector<double> FL(9), FR(9);
    double delta = (dir == 'x') ? (dt / dx) : (dt / dy);
    
    if (dir == 'x') {
        FL = f(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, pL, Gamma, c_h);
        FR = f(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, pR, Gamma, c_h);
    } else if (dir == 'y') {
        FL = g(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, pL, Gamma, c_h);
        FR = g(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, pR, Gamma, c_h);
    }
    
    cvar u_star;
    std::vector<double> uL_values = {uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi};
    std::vector<double> uR_values = {uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi};
    std::vector<double*> u_star_values = {&u_star.rho, &u_star.rhov_x, &u_star.rhov_y, &u_star.rhov_z, &u_star.U, &u_star.B_x, &u_star.B_y, &u_star.B_z, &u_star.psi};
    
    for (size_t i = 0; i < 9; ++i) {
        *u_star_values[i] = 0.5 * (uL_values[i] + uR_values[i]) - 0.5 * delta * (FR[i] - FL[i]);
    }
    
    pvar pri_star = con2pri(u_star.rho, u_star.rhov_x, u_star.rhov_y, u_star.rhov_z, u_star.U, u_star.B_x, u_star.B_y, u_star.B_z, u_star.psi, Gamma);
    
    flux_r = (dir == 'x') ? f(u_star.rho, u_star.rhov_x, u_star.rhov_y, u_star.rhov_z, u_star.U, u_star.B_x, u_star.B_y, u_star.B_z, u_star.psi, pri_star.p, Gamma, c_h)
                            : g(u_star.rho, u_star.rhov_x, u_star.rhov_y, u_star.rhov_z, u_star.U, u_star.B_x, u_star.B_y, u_star.B_z, u_star.psi, pri_star.p, Gamma, c_h);
    
    return flux_r;
}

std::vector<double> FORCE(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double dy, double Gamma, double c_h, char dir){
    std::vector<double> flux(9);
    std::vector<double> flux_f = fri_flux(uL, uR, pL, pR, dt, dx, dy, c_h, dir);
    std::vector<double> flux_r = RI_flux(uL, uR, pL, pR, dt, dx, dy, Gamma, c_h, dir);
    for(int i = 0; i < 9; i++){
        flux[i] = 0.5 * (flux_f[i] + flux_r[i]);
    }
    return flux;
}

void data_reconstruction(double uL, double u, double uR, double &uL_bar, double &uR_bar){
    double delta0, delta1, delta, r;
    delta0 = u - uL;//backward
    delta1 = uR- u;//forward
    delta = 0.5 * (delta0 + delta1);
    if (delta1 == 0){
        r = 0;
    }else{
        r = delta0 / delta1;
    }
    double phi;
    phi = Minbee(r);
    uL_bar = u - 0.5 * phi * delta;
    uR_bar = u + 0.5 * phi * delta;
}

void half_time(cvar &uL_half, cvar &uR_half, cvar &uL, cvar& uR, double dt, double dx, double dy, double Gamma, double c_h, char dir){
    pvar priL = con2pri(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, Gamma);
    pvar priR = con2pri(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, Gamma);
    std::vector<double> FL(9), FR(9);
    double delta = (dir == 'x') ? (dt / dx) : (dt / dy);
    
    if (dir == 'x') {
        FL = f(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, priL.p, Gamma, c_h);
        FR = f(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, priR.p, Gamma, c_h);
    } else if (dir == 'y') {
        FL = g(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, priL.p, Gamma, c_h);
        FR = g(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, priR.p, Gamma, c_h);
    }
    
    std::vector<double> uL_values = {uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi};
    std::vector<double> uR_values = {uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi};
    std::vector<double*> uL_half_values = {&uL_half.rho, &uL_half.rhov_x, &uL_half.rhov_y, &uL_half.rhov_z, &uL_half.U, &uL_half.B_x, &uL_half.B_y, &uL_half.B_z, &uL_half.psi};
    std::vector<double*> uR_half_values = {&uR_half.rho, &uR_half.rhov_x, &uR_half.rhov_y, &uR_half.rhov_z, &uR_half.U, &uR_half.B_x, &uR_half.B_y, &uR_half.B_z, &uR_half.psi};
    
    for (size_t i = 0; i < 9; ++i) {
        *uL_half_values[i] = uL_values[i] - 0.5 * delta * (FR[i] - FL[i]);
        *uR_half_values[i] = uR_values[i] - 0.5 * delta * (FR[i] - FL[i]);
    }
}


void SLIC(std::vector<std::vector<std::vector<double>>>& flux, const std::vector<std::vector<cvar> >& grid, double dt, double dx, double dy, double Gamma, double c_h, char dir){
    std::vector<std::vector<cvar>> halfL(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfR(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfB(grid.size(), std::vector<cvar>(grid[0].size()));
    std::vector<std::vector<cvar>> halfT(grid.size(), std::vector<cvar>(grid[0].size()));
    if(dir == 'x'){
        for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
            for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
            
                cvar uL = grid[i - 1][j];
                cvar u_con = grid[i][j];//u
                cvar uR = grid[i + 1][j];
                
                pvar priL = con2pri(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, Gamma);
                pvar pri = con2pri(u_con.rho, u_con.rhov_x, u_con.rhov_y, u_con.rhov_z, u_con.U, u_con.B_x, u_con.B_y, u_con.B_z, u_con.psi, Gamma);
                pvar priR = con2pri(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, Gamma);
                double rhoL_bar, rhoR_bar, v_xL_bar, v_xR_bar, v_yR_bar, v_yL_bar, v_zL_bar, v_zR_bar, pL_bar, pR_bar, B_xL_bar, B_xR_bar, B_yL_bar, B_yR_bar, B_zL_bar, B_zR_bar, psiL_bar, psiR_bar;
                data_reconstruction(priL.rho, pri.rho, priR.rho, rhoL_bar, rhoR_bar);
                data_reconstruction(priL.v_x, pri.v_x, priR.v_x, v_xL_bar, v_xR_bar);
                data_reconstruction(priL.v_y, pri.v_y, priR.v_y, v_yL_bar, v_yR_bar);
                data_reconstruction(priL.v_z, pri.v_z, priR.v_z, v_zL_bar, v_zR_bar);
                data_reconstruction(priL.p, pri.p, priR.p, pL_bar, pR_bar);
                data_reconstruction(priL.B_x, pri.B_x, priR.B_x, B_xL_bar, B_xR_bar);
                data_reconstruction(priL.B_y, pri.B_y, priR.B_y, B_yL_bar, B_yR_bar);
                data_reconstruction(priL.B_z, pri.B_z, priR.B_z, B_zL_bar, B_zR_bar);
                data_reconstruction(priL.psi, pri.psi, priR.psi, psiL_bar, psiR_bar);

                
                cvar uL_half, uR_half;
                cvar uL_con, uR_con;//uL_bar, uR_bar
                
                uL_con = pri2con(rhoL_bar, v_xL_bar, v_yL_bar, v_zL_bar, pL_bar, B_xL_bar, B_yL_bar, B_zL_bar, psiL_bar, Gamma);
                uR_con = pri2con(rhoR_bar, v_xR_bar, v_yR_bar, v_zR_bar, pR_bar, B_xR_bar, B_yR_bar, B_zR_bar, psiR_bar, Gamma);
                
                half_time(uL_half, uR_half, uL_con, uR_con, dt, dx, dy, Gamma, c_h, 'x');
            
                halfL[i][j] = uL_half;
                halfR[i][j] = uR_half;
            }
        }
        for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
            for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
                cvar uL = halfR[i][j];
                cvar uR = halfL[i + 1][j];
                pvar priL = con2pri(uL.rho, uL.rhov_x, uL.rhov_y, uL.rhov_z, uL.U, uL.B_x, uL.B_y, uL.B_z, uL.psi, Gamma);
                pvar priR = con2pri(uR.rho, uR.rhov_x, uR.rhov_y, uR.rhov_z, uR.U, uR.B_x, uR.B_y, uR.B_z, uR.psi, Gamma);
                flux[i][j] = FORCE(uL, uR, priL.p, priR.p, dt, dx, dy, Gamma, c_h, 'x');
            }
        }
    }else if(dir =='y'){
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++){
            for (int j = ghost_cells - 1; j < nyCells + ghost_cells + 1; j++){
                cvar uB = grid[i][j - 1];
                cvar u_con = grid[i][j];
                cvar uT = grid[i][j + 1];
                
                pvar priB = con2pri(uB.rho, uB.rhov_x, uB.rhov_y, uB.rhov_z, uB.U, uB.B_x, uB.B_y, uB.B_z, uB.psi, Gamma);
                pvar pri = con2pri(u_con.rho, u_con.rhov_x, u_con.rhov_y, u_con.rhov_z, u_con.U, u_con.B_x, u_con.B_y, u_con.B_z, u_con.psi, Gamma);
                pvar priT = con2pri(uT.rho, uT.rhov_x, uT.rhov_y, uT.rhov_z, uT.U, uT.B_x, uT.B_y, uT.B_z, uT.psi, Gamma);
                // double rhoB_bar, rhoT_bar, uB_bar, uT_bar, vB_bar, vT_bar, pB_bar, pT_bar;
                double rhoB_bar, rhoT_bar, v_xB_bar, v_xT_bar, v_yB_bar, v_yT_bar, v_zB_bar, v_zT_bar, pB_bar, pT_bar, B_xB_bar, B_xT_bar, B_yB_bar, B_yT_bar, B_zB_bar, B_zT_bar, psiB_bar, psiT_bar;

            
                data_reconstruction(priB.rho, pri.rho, priT.rho, rhoB_bar, rhoT_bar);
                data_reconstruction(priB.v_x, pri.v_x, priT.v_x, v_xB_bar, v_xT_bar);
                data_reconstruction(priB.v_y, pri.v_y, priT.v_y, v_yB_bar, v_yT_bar);
                data_reconstruction(priB.v_z, pri.v_z, priT.v_z, v_zB_bar, v_zT_bar);
                data_reconstruction(priB.p, pri.p, priT.p, pB_bar,pT_bar);
                data_reconstruction(priB.B_x, pri.B_x, priT.B_x, B_xB_bar, B_xT_bar);
                data_reconstruction(priB.B_y, pri.B_y, priT.B_y, B_yB_bar, B_yT_bar);
                data_reconstruction(priB.B_z, pri.B_z, priT.B_z, B_zB_bar, B_zT_bar);
                data_reconstruction(priB.psi, pri.psi, priT.psi, psiB_bar, psiT_bar);
            
                cvar uB_half, uT_half;
                cvar uB_con, uT_con;// uB_bar, uT_bar
            
                uB_con = pri2con(rhoB_bar, v_xB_bar, v_yB_bar, v_zB_bar, pB_bar, B_xB_bar, B_yB_bar, B_zB_bar, psiB_bar, Gamma);
                uT_con = pri2con(rhoT_bar, v_xT_bar, v_yT_bar, v_zT_bar, pT_bar, B_xT_bar, B_yT_bar, B_zT_bar, psiT_bar, Gamma);
                
                half_time(uB_half, uT_half, uB_con, uT_con, dt, dx, dy, Gamma, c_h, 'y');
            
                halfB[i][j] = uB_half;
                halfT[i][j] = uT_half;
            }
        }
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++){
            for (int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
                cvar uB = halfT[i][j];
                cvar uT = halfB[i][j + 1];

                pvar priB = con2pri(uB.rho, uB.rhov_x, uB.rhov_y, uB.rhov_z, uB.U, uB.B_x, uB.B_y, uB.B_z, uB.psi, Gamma);
                pvar priT = con2pri(uT.rho, uT.rhov_x, uT.rhov_y, uT.rhov_z, uT.U, uT.B_x, uT.B_y, uT.B_z, uT.psi, Gamma);
                flux[i][j] = FORCE(uB, uT, priB.p, priT.p, dt, dx, dy, Gamma, c_h, 'y');
            }
        }
    }
}

double compute_Cf(double rho, double p, double B_x, double B_y, double B_z, double Gamma) {
    
    double Cs = sqrt(Gamma * p / rho);  // Cs
    //double Ca = fabs(B_x) / sqrt(rho); // alfven
    double Ca = sqrt((B_x * B_x + B_y * B_y + B_z * B_z) / rho);
    // double term1 = sqrt((Cs * Cs + Ca * Ca) * (Cs * Cs + Ca * Ca) - 4 * (Cs * Cs * (B_x * B_x) / rho));
    // //double term1 = sqrt((Cs * Cs + Ca * Ca) * (Cs * Cs + Ca * Ca) - 4 * (Cs * Cs * B_x * B_x) / rho);
    // if(term1 = 0) term1 = 0;
    double disc_x = (Cs * Cs + Ca * Ca) * (Cs * Cs + Ca * Ca) - 4 * (Cs * Cs * (B_x * B_x) / rho);
    if (disc_x < 0) disc_x = 0;  // avoid sqrt < 0
    double term1 = sqrt(disc_x);

    double Cf = sqrt(0.5 * (Cs * Cs + Ca * Ca) + term1);
    return Cf;
}

double update(const std::vector<std::vector<cvar> >& grid, double dx,double dy, double Gamma, double CFL){
    double max_speed = 0.0;
    for(int i = 0; i < grid.size(); i++){
        for(int j = 0; j < grid.size(); j++){
            pvar pri = con2pri(grid[i][j].rho, grid[i][j].rhov_x, grid[i][j].rhov_y, grid[i][j].rhov_z, grid[i][j].U, grid[i][j].B_x, grid[i][j].B_y, grid[i][j].B_z, grid[i][j].psi, Gamma);
            double c_f = compute_Cf(pri.rho, pri.p, pri.B_x, pri.B_y, pri.B_z, Gamma);
            double Cs = sqrt(Gamma * pri.p / pri.rho);
            double speed_x = fabs(pri.v_x) + c_f;
            double speed_y = fabs(pri.v_y) + c_f;
            double speed = std::max(speed_x, speed_y);

            if (speed > max_speed){
                max_speed = speed;
            }
        }
    }
    return CFL * fmin(dx,dy) / max_speed;
    //std::cout << "update: max_speed = " << max_speed << ", dt = " << CFL * dx / max_speed << std::endl;
}

double ch(std::vector<std::vector<cvar> >& grid){
    double c_h = 1e-8;
    std::vector<std::vector<pvar>> pvar_grid(nxCells + 2 * ghost_cells, std::vector<pvar>(nyCells + 2 * ghost_cells));
    for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
        for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
            pvar_grid[i][j] = con2pri(grid[i][j].rho, grid[i][j].rhov_x, grid[i][j].rhov_y, grid[i][j].rhov_z, 
                                    grid[i][j].U, grid[i][j].B_x, grid[i][j].B_y, grid[i][j].B_z, grid[i][j].psi, Gamma);
            double rho = pvar_grid[i][j].rho;
            double p = pvar_grid[i][j].p;
            double B_x = pvar_grid[i][j].B_x;
            double B_y = pvar_grid[i][j].B_y;
            double B_z = pvar_grid[i][j].B_z;
            double B2 = B_x * B_x + B_y * B_y + B_z * B_z;
            double cf_x = sqrt(0.5 * (((Gamma * p + B2)/(rho)) + sqrt( ((Gamma * p + B2)/rho) * ((Gamma * p + B2)/rho) - 4 * ((Gamma * p * B_x * B_x)/(rho * rho)) )));
			double cf_y = sqrt(0.5 * (((Gamma * p + B2)/(rho)) + sqrt( ((Gamma * p + B2)/rho) * ((Gamma * p + B2)/rho) - 4 * ((Gamma * p * B_y * B_y)/(rho * rho)) )));
			double cf_z = sqrt(0.5 * (((Gamma * p + B2)/(rho)) + sqrt( ((Gamma * p + B2)/rho) * ((Gamma * p + B2)/rho) - 4 * ((Gamma * p * B_z * B_z)/(rho * rho)) )));
            if(c_h < std::max({fabs(pvar_grid[i][j].v_x) + cf_x, std::max(fabs(pvar_grid[i][j].v_y) + cf_y, fabs(pvar_grid[i][j].v_z) + cf_z)})){
				c_h = std::max({fabs(pvar_grid[i][j].v_x) + cf_x, std::max(fabs(pvar_grid[i][j].v_y) + cf_y, fabs(pvar_grid[i][j].v_z) + cf_z)});
			}
        }
    }
    return c_h;
}

double compute_cp(double c_h) {
    return sqrt(0.18 * c_h);
}

void compute_source_term(std::vector<std::vector<cvar>>& source,  std::vector<std::vector<cvar>>& grid) {
    //double factor = - (ch * ch) / (cp * cp);  // 计算 Source Term 系数
    double c_h = ch(grid);
    double c_p = compute_cp(c_h);
    double factor = - (c_h * c_h) / (c_p * c_p);
    for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
        for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
            source[i][j].rho = 0.0;
            source[i][j].rhov_x = 0.0;
            source[i][j].rhov_y = 0.0;
            source[i][j].rhov_z = 0.0;
            source[i][j].U = 0.0;
            source[i][j].B_x = 0.0;
            source[i][j].B_y = 0.0;
            source[i][j].B_z = 0.0;
            source[i][j].psi = factor * grid[i][j].psi;  // S_psi = - (ch^2 / cp^2) * psi
        }
    }
}

void update_source(std::vector<std::vector<cvar>> &source, std::vector<std::vector<cvar>> &u, std::vector<std::vector<cvar>> &u_con, double c_h, double c_p, double dt) {
    std::vector<std::vector<std::vector<double>>> K1(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K2(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K3(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K4(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    
    // 第一阶段：计算 K1 并更新中间状态 u = u_con + 0.5*K1
    compute_source_term(source, u_con);
    for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int k = 0; k < 9; k++) {
                K1[i][j][k] = dt * source[i][j].values[k];
                u[i][j].values[k] = u_con[i][j].values[k] + 0.5 * K1[i][j][k];
            }
        }
    }
    
    // 第二阶段：以更新后的 u 计算源项，得到 K2，然后 u = u_con + 0.5*K2
    compute_source_term(source, u);
    for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int k = 0; k < 9; k++) {
                K2[i][j][k] = dt * source[i][j].values[k];
                u[i][j].values[k] = u_con[i][j].values[k] + 0.5 * K2[i][j][k];
            }
        }
    }
    
    // 第三阶段：以更新后的 u 计算源项，得到 K3，然后 u = u_con + K3
    compute_source_term(source, u);
    for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int k = 0; k < 9; k++) {
                K3[i][j][k] = dt * source[i][j].values[k];
                u[i][j].values[k] = u_con[i][j].values[k] + 0.5 * K3[i][j][k];
            }
        }
    }
    
    // 第四阶段：以更新后的 u 计算源项，得到 K4
    compute_source_term(source, u);
    for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            for (int k = 0; k < 9; k++) {
                K4[i][j][k] = dt * source[i][j].values[k];
                // u_con^(n+1) = u_con + 1/6*(K1 + 2*K2 + 2*K3 + K4)
                u_con[i][j].values[k] = u_con[i][j].values[k] + (1.0/6.0) * (K1[i][j][k] + 2 * K2[i][j][k] + 2 * K3[i][j][k] + K4[i][j][k]);
            }
        }
    }
}




void trans(std::vector<std::vector<cvar> >& grid){
    // BC
    for (int j = ghost_cells; j < ghost_cells + nyCells; j++) {
        for (int i = 0; i < ghost_cells; i++) {
            grid[i][j] = grid[ghost_cells][j]; // 左边界
            grid[nxCells + ghost_cells + i][j] = grid[nxCells + ghost_cells - 1][j]; // 右边界
        }
    }
    for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        for (int j = 0; j < ghost_cells; j++) {
            grid[i][j] = grid[i][ghost_cells]; // 底边界
            grid[i][nyCells + ghost_cells + j] = grid[i][nyCells + ghost_cells - 1]; // 顶边界
        }
    }
}

void periodic(std::vector<std::vector<cvar> >& grid){
    for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
        for(int i = 0; i < ghost_cells; i++){
            grid[i][j] = grid[nxCells + i][j]; // left
            grid[nxCells + ghost_cells + i][j] = grid[ghost_cells + i][j];// right
        }
    }
    for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
        for(int j = 0; j < ghost_cells; j++){
            grid[i][j] = grid[i][nyCells + j]; //top
            grid[i][nyCells + ghost_cells + j] = grid[i][ghost_cells + j]; // bottom
        }
    }
}

int main(){
    std::vector<std::vector<pvar> > primitive_grid(nxCells + 2 * ghost_cells, std::vector<pvar>(nyCells + 2 * ghost_cells));
    std::vector<std::vector<cvar>> source(nxCells + 2 * ghost_cells, std::vector<cvar>(nyCells + 2 * ghost_cells));
    std::vector<std::vector<cvar>> u(nxCells + 2 * ghost_cells, std::vector<cvar>(nyCells + 2 * ghost_cells));

    std::vector<std::vector<std::vector<double>>> K1(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K2(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K3(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    std::vector<std::vector<std::vector<double>>> K4(nxCells + 2*ghost_cells, std::vector<std::vector<double>>(nyCells + 2*ghost_cells, std::vector<double>(9, 0.0)));
    

    initializeGrid(primitive_grid);
    

    std::vector<std::vector<cvar> > grid(nxCells + 2 * ghost_cells, std::vector<cvar>(nyCells + 2 * ghost_cells));
    for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
        for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
            grid[i][j] = pri2con(primitive_grid[i][j].rho, primitive_grid[i][j].v_x, primitive_grid[i][j].v_y, primitive_grid[i][j].v_z,  primitive_grid[i][j].p, primitive_grid[i][j].B_x, primitive_grid[i][j].B_y, primitive_grid[i][j].B_z, primitive_grid[i][j].psi, Gamma);
        }
    }
    periodic(grid);


    double t = 0.0;


    
    for (;;){
        // for (int j = ghost_cells; j < ghost_cells + nyCells; j++) {
        //     for (int i = 0; i < ghost_cells; i++) {
        //         grid[i][j] = grid[ghost_cells][j]; // 左边界
        //         grid[nxCells + ghost_cells + i][j] = grid[nxCells + ghost_cells - 1][j]; // 右边界
        //     }
        // }
        // for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        //     for (int j = 0; j < ghost_cells; j++) {
        //         grid[i][j] = grid[i][ghost_cells]; // 底边界
        //         grid[i][nyCells + ghost_cells + j] = grid[i][nyCells + ghost_cells - 1]; // 顶边界
        //     }
        // }
        periodic(grid);
        if (t > T_end){
            break;
        }

        
        //double dt = update(grid, dx, dy, Gamma, CFL);
         double dt = update(grid, dx, dy, Gamma, CFL);
        // if (dt > (T_end - t)) {
        //     dt = T_end - t;  
        // }
        double c_h = ch(grid);
        double c_p = compute_cp(c_h);
        

        //std::cout <<"update done" << std::endl;

        
        std::vector<std::vector<std::vector<double>>> flux_x(nxCells + 2 * ghost_cells - 1, std::vector<std::vector<double>>(nyCells + 2 * ghost_cells, std::vector<double>(9)));

        
        SLIC(flux_x, grid, dt, dx, dy, Gamma, c_h, 'x');
        //std::cout << "slic x done" << std::endl;

        
        //std::vector<std::vector<cvar> > grid_x = grid;

        
        // for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        //     for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
        //         grid[i][j].rho  -= (dt / dx) * (flux_x[i][j][0] - flux_x[i - 1][j][0]);
        //         grid[i][j].rhou -= (dt / dx) * (flux_x[i][j][1] - flux_x[i - 1][j][1]);
        //         grid[i][j].rhov -= (dt / dx) * (flux_x[i][j][2] - flux_x[i - 1][j][2]);
        //         grid[i][j].E    -= (dt / dx) * (flux_x[i][j][3] - flux_x[i - 1][j][3]);
        //     }
        // }
        //std::cout << "update x done" << std::endl;
//         for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
//             for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
//                 for (size_t k = 0; k < 8; ++k) {
//                     grid[i][j].values[k] -= (dt / dx) * (flux_x[i][j][k] - flux_x[i - 1][j][k]);
//         }
//     }
// }
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++){
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++){
                
                std::vector<double*> grid_values = {&grid[i][j].rho, &grid[i][j].rhov_x, &grid[i][j].rhov_y, 
                                                    &grid[i][j].rhov_z, &grid[i][j].U, &grid[i][j].B_x, 
                                                    &grid[i][j].B_y, &grid[i][j].B_z, &grid[i][j].psi};
                for (size_t k = 0; k < 9; ++k) {
                    *grid_values[k] -= (dt / dx) * (flux_x[i][j][k] - flux_x[i-1][j][k]);
                }
            }
        }


        
        // for (int j = ghost_cells; j < ghost_cells + nyCells; j++){
        //     for (int i = 0; i < ghost_cells; i++) {
        //         grid[i][j] = grid[ghost_cells][j]; // 左边界
        //         grid[nxCells + ghost_cells + i][j] = grid[nxCells + ghost_cells - 1][j]; // 右边界
        //     }
        // }
        // for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
        //     for (int j = 0; j < ghost_cells; j++) {
        //         grid[i][j] = grid[i][ghost_cells]; // 底边界
        //         grid[i][nyCells + ghost_cells + j] = grid[i][nyCells + ghost_cells - 1]; // 顶边界
        //     }
        // }        //std::cout << "bc done" << std::endl;
        periodic(grid);

        
        std::vector<std::vector<std::vector<double>>> flux_y(nxCells + 2 * ghost_cells, std::vector<std::vector<double>>(nyCells + 2 * ghost_cells - 1, std::vector<double>(9)));

        
        SLIC(flux_y, grid, dt, dx, dy, Gamma, c_h,'y');
    
        
        for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
            for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
                std::vector<double*> grid_values = {&grid[i][j].rho, &grid[i][j].rhov_x, &grid[i][j].rhov_y, &grid[i][j].rhov_z, &grid[i][j].U, &grid[i][j].B_x, &grid[i][j].B_y, &grid[i][j].B_z, &grid[i][j].psi};
                for (size_t k = 0; k < 9; ++k) {
                    *grid_values[k] -= (dt / dy) * (flux_y[i][j][k] - flux_y[i][j - 1][k]);
                }
            }
        }

        periodic(grid);

        update_source(source, u, grid, c_h, c_p, dt);
        periodic(grid);

        
        
        t += dt;

        
        std::cout << "Time: " << t << " / " << T_end << ", dt: " << dt << std::endl;
    }

    
    std::ofstream outFile("DCmhd2D_results.dat");
    for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
        for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
            pvar pri = con2pri(grid[i][j].rho, grid[i][j].rhov_x, grid[i][j].rhov_y, grid[i][j].rhov_z, grid[i][j].U, grid[i][j].B_x, grid[i][j].B_y, grid[i][j].B_z, grid[i][j].psi, Gamma);
            //outFile << i * dx << " " << j * dy << " " << pri.rho << " " << " " << pri.p << " " << std::endl;
            outFile << i * dx << " " << j * dy << " " << pri.rho  << " " << std::endl;
        }
    }
    outFile.close();

    std::cout << "Simulation completed. Results saved in 'DCmhd2D_results.dat'." << std::endl;
    return 0;
}