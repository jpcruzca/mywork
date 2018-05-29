#include <iostream>
#include "rk4-sim.h"
#include <cmath>
#include <cstdlib>
#include <vector>

void initial_condition(std::vector<double> & pos, std::vector<double> & vel)
{
  //posiciones
  pos = {0.994, 0.0, 0.0};
  
  //velocidades
  vel = {0.0, -2.0015851063790825224053786224, 0.0};
}

void print(const std::vector<double> & v)
{
  for (int i = 0; i < v.size(); ++i)
    {
      std::cout << v[i] << " ";
    }
}

double rvar(const std::vector<double> & pos, const std::vector<double> & vel)
{
  return std::sqrt((pos[0]+u)*(pos[0]+u)+(pos[1]*pos[1])+(pos[2]*pos[2]));
}

double svar(const std::vector<double> & pos, const std::vector<double> & vel)
{
  return std::sqrt((pos[0]+u-1)*(pos[0]+u-1)+(pos[1]*pos[1])+(pos[2]*pos[2]));
}

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id)
{
  double r = rvar(pos,vel);
  double s = svar(pos,vel);
  
  if(id==1){ //returna funcion para calcular la velocidad en x
    return (pos[0] + 2*vel[1] - ((1-u)*(pos[0]+u)/(r*r*r)) - (u*(pos[0]-1+u)/(s*s*s)));
  }
  if(id==2){ //returna funcion para calcular la velocidad en y
    return (pos[1] - 2*vel[0] - ((1-u)*pos[1]/(r*r*r)) - (u*pos[1]/(s*s*s)));
  }
  if(id==3){ //returna funcion para calcular la velocidad en z
    return (((u-1)*pos[2]/(r*r*r)) - (u*pos[2]/(s*s*s)));
  }
  if(id==4){ //retorna la velocidad en x
    return vel[0];
  }
  if(id==5){ //retorna la velocidad en y
    return vel[1];
  }
  if(id==6){ //retorna la velocidad en z
    return vel[2];
  }
  else{ //advierte de un error
    std::cerr << "Error!!!!" << id << std::endl;
    exit(1);
  }
}

void rk4(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend)
{
  std::vector<double> k1, k2, k3, k4, vaux, paux;
  k1.resize(6);
  k2.resize(6);
  k3.resize(6);
  k4.resize(6);
  vaux.resize(vel.size());
  paux.resize(vel.size());
  const double dt=T/6000.0;
  const int N =(int)(tend-tini)/dt;

  //calculo
  for (int nt = 0; nt < 1*N; ++nt) {
    double t = tini + dt*nt;
    
    // k1
    for(int ii = 0; ii < vel.size(); ++ii) {
      k1[ii+3] = dt*compute(pos, vel,t, ii+1);
      k1[ii] = dt*compute(pos, vel,t, ii+4);
    }
    // k2 aux
    for(int ii = 0; ii < vel.size(); ++ii) {
      vaux[ii] = vel[ii] + k1[ii+3]/2;
      paux[ii] = pos[ii] + k1[ii]/2;
    }
    //k2
    for(int ii = 0; ii < vel.size(); ++ii) {
     k2[ii+3] = dt*compute(paux,vaux,t+dt/2, ii+1);
     k2[ii] = dt*compute(paux,vaux,t+dt/2, ii+4);
    }
    // k3 aux
    for(int ii = 0; ii < vel.size(); ++ii) {
      vaux[ii] = vel[ii] + k2[ii+3]/2;
      paux[ii] = pos[ii] + k2[ii]/2;
    }
    //k3
    for(int ii = 0; ii < vel.size(); ++ii) {
      k3[ii+3] = dt*compute(paux,vaux,t+dt/2, ii+1);
      k3[ii] = dt*compute(paux,vaux,t+dt/2, ii+4);
    }
    // k4 aux
    for(int ii = 0; ii < vel.size(); ++ii) {
      vaux[ii] = vel[ii] + k3[ii+3];
      paux[ii] = pos[ii] + k3[ii];
    }
    //k4
    for(int ii = 0; ii < vel.size(); ++ii) {
      k4[ii+3] = dt*compute(paux, vaux, t + dt, ii+1);
      k4[ii] = dt*compute(paux, vaux, t + dt, ii+4);
    }
    // new vel
    for(int ii = 0; ii < vel.size(); ++ii) {
      vel[ii] = vel[ii] + (k1[ii+3] + 2*k2[ii+3] + 2*k3[ii+3] + k4[ii+3])/6.0;
      pos[ii] = pos[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }
    
    std::cout << t  << " ";
    print(pos);
    print(vel);
    //plot_gnuplot(pos);
    // plot_gnuplot(vel);
    std::cout << "\n";
  }
  
} 
void init_gnuplot(void)
{
  //std::cout << "set contour " << std::endl;
  //std::cout << "set terminal gif animate " << std::endl;
  std::cout << "set out 'anim.gif' " << std::endl;
}

void plot_gnuplot(const std::vector<double>& m)
{
  std::cout << "plot '-' u 2:3 " << std::endl;
  print(m);
  std::cout << "e" << std::endl;  
}



