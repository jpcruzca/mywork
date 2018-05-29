#include <iostream>
#include "dp54-sim.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

void dp54(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend)
{
  std::vector<double> k1, k2, k3, k4, k5, k6, k7, vaux, paux;
  double pzy = 0.0;
  double vzy = 0.0;
  k1.resize(6);
  k2.resize(6);
  k3.resize(6);
  k4.resize(6);
  k5.resize(6);
  k6.resize(6);
  k7.resize(6);
  vaux.resize(vel.size());
  paux.resize(vel.size());
  double dt=0.001;
  
  //calculo
  for (int tt=0;tt<137;tt++)
    {
      double t = tini + dt*tt;
      
      // k1
      for(int ii = 0; ii < vel.size(); ++ii) {
	k1[ii+3] = dt*compute(pos, vel,t, ii+1);//deberia ser diferente.
	k1[ii] = dt*compute(pos, vel,t, ii+4);
      }
      // k2 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + k1[ii+3]/5.0;
	paux[ii] = pos[ii] + k1[ii]/5.0;
      }
      //k2
      for(int ii = 0; ii < vel.size(); ++ii) {
	k2[ii+3] = dt*compute(paux,vaux,t+dt/5, ii+1);
	k2[ii] = dt*compute(paux,vaux,t+dt/5, ii+4);
      }
      // k3 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 3.0*k1[ii+3]/40.0 + 9.0*k2[ii+3]/40.0;
	paux[ii] = pos[ii] + 3.0*k1[ii]/40.0 + 9*k2[ii]/40.0;
      }
      //k3
      for(int ii = 0; ii < vel.size(); ++ii) {
	k3[ii+3] = dt*compute(paux,vaux,3*(t+dt)/10, ii+1);
	k3[ii] = dt*compute(paux,vaux,3*(t+dt)/10, ii+4);
      }
      // k4 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 44*k1[ii+3]/45 - 56*k2[ii+3]/15 + 32*k3[ii+3]/9;
	paux[ii] = pos[ii] + 44*k1[ii]/45 - 56*k2[ii]/15 + 32*k3[ii]/9;
      }
      //k4
      for(int ii = 0; ii < vel.size(); ++ii) {
	k4[ii+3] = dt*compute(paux,vaux,4*(t+dt)/5, ii+1);
	k4[ii] = dt*compute(paux,vaux,4*(t+dt)/5, ii+4);
      }
      // k5 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 38744.0*k1[ii+3]/13122 - 25360.0*k2[ii+3]/2187 + 64448.0*k3[ii+3]/6561 - 212.0*k4[ii+3]/729;
	paux[ii] = pos[ii] + 38744*k1[ii]/13122 - 25360*k2[ii]/2187 + 64448*k3[ii]/6561 - 212*k4[ii]/729;
      }
      //k5
      for(int ii = 0; ii < vel.size(); ++ii) {
	k5[ii+3] = dt*compute(paux,vaux,8*(t+dt)/9, ii+1);
	k5[ii] = dt*compute(paux,vaux,8*(t+dt)/9, ii+4);
      }      
      // k6 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 9017*k1[ii+3]/3168 - 355*k2[ii+3]/33 + 46732*k3[ii+3]/5247 + 49*k4[ii+3]/176 - 5103*k5[ii+3]/18656;
	paux[ii] = pos[ii] + 9017*k1[ii]/3168 - 355*k2[ii]/33 + 46732*k3[ii]/5247 + 49*k4[ii]/176 - 5103*k5[ii]/18656;
      }
      //k6
      for(int ii = 0; ii < vel.size(); ++ii) {
	k6[ii+3] = dt*compute(paux, vaux, t + dt, ii+1);
	k6[ii] = dt*compute(paux, vaux, t + dt, ii+4);
      }
      // k7 aux
      for(int ii = 0; ii < vel.size(); ++ii) {
	vaux[ii] = vel[ii] + 35*k1[ii+3]/384 + 500*k3[ii+3]/1113 + 125*k4[ii+3]/192 - 2187*k5[ii+3]/6784 + 11*k6[ii+3]/84;
	paux[ii] = pos[ii] + 35*k1[ii]/384 + 500*k3[ii]/1113 + 125*k4[ii]/192 - 2187*k5[ii]/6784 + 11*k6[ii]/84;
      }
      //k7
      for(int ii = 0; ii < vel.size(); ++ii) {
	k7[ii+3] = dt*compute(paux, vaux, t + dt, ii+1);
	k7[ii] = dt*compute(paux, vaux, t + dt, ii+4);
      }
      // new vel,pos
      for(int ii = 0; ii < vel.size(); ++ii) {
	pos[ii] = pos[ii] + 35*k1[ii]/384 + 500*k3[ii]/1113 + 125*k4[ii]/192 - 2187*k5[ii]/6784 + 11*k6[ii]/84;
	vel[ii] = vel[ii] + 35*k1[ii+3]/384 + 500*k3[ii+3]/1113 + 125*k4[ii+3]/192 - 2187*k5[ii+3]/6784 + 11*k6[ii+3]/84;
      }
      
            //errores
      for (int ii = 0; ii < vel.size(); ++ii)
	{
	  vaux[ii] = 71*k1[ii+3]/57600 - 71*k3[ii+3]/16695 + 71*k4[ii+3]/1920 - 17253*k5[ii+3]/339200 + 22*k6[ii+3]/525 - k7[ii+3]/40;
	  paux[ii] = 71*k1[ii]/57600 - 71*k3[ii]/16695 + 71*k4[ii]/1920 - 17253*k5[ii]/339200 + 22*k6[ii]/525 - k7[ii]/40;
	}

      pzy = norm(paux);
      vzy = norm(vaux);
      
      //new dt
      dt=dtnew(pzy, vzy, dt); 

      //print
      std::cout <<dt<<" "<< t  << " ";
      print(pos);
      print(vel);
      std::cout << "\n";
      
    }
}  
  

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id)
{
  double r = rvar(pos,vel);
  double s = svar(pos,vel);
  
  if(id==1){ //retorna funcion para calcular la velocidad en x
    return (pos[0] + 2*vel[1] - ((1-u)*(pos[0]+u)/(r*r*r)) - (u*(pos[0]-1+u)/(s*s*s)));
  }
  if(id==2){ //retorna funcion para calcular la velocidad en y
    return (pos[1] - 2*vel[0] - ((1-u)*pos[1]/(r*r*r)) - (u*pos[1]/(s*s*s)));
  }
  if(id==3){ //retorna funcion para calcular la velocidad en z
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
    std::cerr << "Error en la funcion compute pelmazo!!!!" << id << std::endl;
    exit(1);
  }
}

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

double saux(const double x, const double dt)
{
  double aux = eps*dt/(2*std::fabs(x));
  return std::pow(aux, 1.0/5.0); 
}

double dtnew(const double p, const double v, const double dt)
{
  double sprom = (saux(p,dt) + saux(v,dt))/2.0;
  double hp=0.001;
  double hg=0.138;
  if(dt*sprom <=hp)
    return hp;
  else if (dt*sprom >=hg)
    return hg;
  else
    return dt*sprom;
  //*/
  //return sprom*dt;
}

double norm(const std::vector<double> & x)
{
  double res = 0;
  for (int i = 0; i < x.size(); ++i)
  {
    res += x[i]*x[i];
  }
  return std::sqrt(res);
}
