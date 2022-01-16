#include <iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>

using namespace std;

double drand48(){
return rand()/(RAND_MAX+1.0);
}

int main()
{
    //simulation parameters
int therms = 10000000;
int sweeps = 10000000;
double delta = 0.00001;
int tau;

    // physics parameters
int T = 64; // number of time slices
double omega = 1.0; // frequency omega
double m = 1.0; // mass m

double x[T+1];
double x_new[T+1];
double dS=0;
double corr[T+1];

//Initialization
for(int i = 0;i<=T;i++){
    corr[i] = 0;
    x[i] = min(i,T-i);
    x_new[i] = x[i];
}
double u;
//thermalization
int counter = 0;
while(counter<=therms){
    //selecting a random t to change x
    tau = int((T)*(drand48()));
    x_new[tau] = x_new[tau] + 2*delta*(drand48()-0.5);
    x_new[T] = x_new[0];
    //computing change in action
    for(int i = 0;i<T;i++){

    dS = dS - (m/2.0)*((omega*omega/4.0)*(x[i+1]+x[i])*(x[i+1]+x[i]) + (x[i+1]-x[i])*(x[i+1]-x[i]));
    dS = dS + (m/2.0)*((omega*omega/4.0)*(x_new[i+1]+x_new[i])*(x_new[i+1]+x_new[i]) + (x_new[i+1]-x_new[i])*(x_new[i+1]-x_new[i]));
    }
   u = drand48();
   if(u<=exp(-dS)){
    x[tau] = x_new[tau];
    x[T] = x_new[T];
   }
   x_new[tau] = x[tau];
   x_new[T] = x[T];

   counter++;
}
//fstream acceptance;
//acceptance.open("Acceptance_values.txt",ios::out);
counter = 0;
int accept = 0;
//MCMC
while(counter<=sweeps){
    //selecting a random t to change x
    tau = int((T)*(drand48()));
    x_new[tau] = x_new[tau] + 2*delta*(drand48()-0.5);
     x_new[T] = x_new[0];
    //computing change in action
    for(int i = 0;i<T;i++){

    dS = dS - (m/2.0)*((omega*omega/4.0)*(x[i+1]+x[i])*(x[i+1]+x[i]) + (x[i+1]-x[i])*(x[i+1]-x[i]));
    dS = dS + (m/2.0)*((omega*omega/4.0)*(x_new[i+1]+x_new[i])*(x_new[i+1]+x_new[i]) + (x_new[i+1]-x_new[i])*(x_new[i+1]-x_new[i]));
    }
   u = drand48();
   if(u<=exp(-dS)){
       //cout<<accept<<'\n';
        accept++;
        //acceptance<<counter<<"\t"<<100.0*accept/counter<<"\n";
    x[tau] = x_new[tau];
    x[T] = x_new[T];
   }
   x_new[tau] = x[tau];
   x_new[T] = x[T];

   //computing correlator
   for(int t =0;t<=T;t++){
    corr[t] = corr[t] + (1.0/sweeps)*exp(-omega*abs(x[t]-x[0]))/(2.0*m*omega);
   }

   counter++;
}
//acceptance.close();
cout<<"acceptance rate is:"<<100.0*accept/sweeps<<"%\n";
fstream f_data;
f_data.open("Correlator_64.txt",ios::out);
for(int t = 0;t<=T;t++){
f_data<<t<<"\t"<<corr[t]<<"\n";
}
f_data.close();
    return 0;
}
