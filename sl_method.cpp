#include "sl_method.h"
#include <iostream>
#include <cmath>
#include <omp.h>

SL_method::SL_method()
{

}

void SL_method::set_velX(CF_2 &velocity_x_t0){
    velocity_x = &velocity_x_t0;
}

void SL_method::set_velY(CF_2 &velocity_y_t0){
    velocity_y = &velocity_y_t0;
}

void SL_method::set_xd(std::vector<double> & xd_){
    xd=xd_;
}

void SL_method::set_yd(std::vector<double> & yd_){
    yd=yd_;
}

void SL_method::set_grid(Grid2D & grid_t0){
    grid = grid_t0;
}

void SL_method::set_solution(std::vector<double> & solution_t0){
    solution = solution_t0;

}

void SL_method::set_True(std::vector<double> & solution_true){
    True = solution_true;

}


void SL_method::Runge_Kutta2(double dt){
     int Size = grid.get_N()*grid.get_M();
     std::vector<double> Solution(Size);

#pragma omp parallel for
     for ( int n = 0; n<Size;n++ ){
         double x = grid.x_from_n(n);
         double y = grid.y_from_n(n);
         double x1 = x-(dt/2)*(*velocity_x)(x,y);
         double y1 = y-(dt/2)*(*velocity_y)(x,y);
         xd[n] = x - dt*(*velocity_x)(x1,y1);
         yd[n] = y - dt*(*velocity_y)(x1,y1);
     }
}

void SL_method::bilinear_interpolationENO(){

std::vector<double> TempSol(grid.get_N()*grid.get_M());
#pragma omp parallel for
    for ( int n = 0; n<grid.get_N()*grid.get_M();n++ ){

        //Get gridpoints to the right and left of x
        int i_min = (int) floor( (xd[n] - grid.get_xmin())/grid.get_dx() );
            i_min = std::max(0,i_min);
        int i_max = (int) ceil ((xd[n] - grid.get_xmin())/grid.get_dx());
            i_max = std::min(i_max, grid.get_N()-1);

        //Check if x is a gridpoint
        if (i_min == i_max){
                if (i_min==0){
                    i_max = i_min+1;
                }
                else {
                    i_min = i_max-1;
                }
        }

        //Get gridpoints above and below y
        int j_min = (int) floor( (yd[n] - grid.get_ymin())/grid.get_dy() );
            j_min = std::max(0,j_min);
        int j_max = (int) ceil ( (yd[n] - grid.get_ymin())/grid.get_dy() );
          j_max = std::min(j_max, grid.get_M()-1);

        //Check if y is a gridpoint
        if (j_min == j_max){
                if (j_min==0){
                    j_max = j_min+1;
                }
                else {
                    j_min = j_max-1;
                }
        }

        //Generate the location of the corners
        int Corner_00 = grid.n_from_ij(i_min, j_min);
        int Corner_01 = grid.n_from_ij(i_min, j_max);
        int Corner_10 = grid.n_from_ij(i_max, j_min);
        int Corner_11 = grid.n_from_ij(i_max, j_max);

        double x_min = grid.x_from_n(Corner_00);
        double y_min = grid.y_from_n(Corner_00);

        double x_max = grid.x_from_n(Corner_11);
        double y_max = grid.y_from_n(Corner_11);

        //Get dx and dy
        double dx = grid.get_dx();
        double dy = grid.get_dy();



        //Interpolate!
        TempSol[n] = (1./(dx*dy))*(solution[Corner_00]*(x_max-xd[n])*(y_max-yd[n]) +
                                    solution[Corner_01]*(x_max-xd[n])*(yd[n]-y_min) +
                                    solution[Corner_10]*(xd[n]-x_min)*(y_max-yd[n]) +
                                    solution[Corner_11]*(xd[n]-x_min)*(yd[n]-y_min));

    }
    solution=TempSol;
}

void SL_method::save_vtk(std::string file_name,std::string file_name2){
//    Grid2D grid(N,M, xmin, xmax, ymin, ymax);
    //save as .vtk file
    grid.initialize_VTK_file(file_name);
    grid.print_VTK_Format(True,"True",file_name);
    grid.initialize_VTK_file(file_name2);
    grid.print_VTK_Format(solution,"EstSol",file_name2);
}
