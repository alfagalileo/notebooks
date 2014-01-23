//Ising 2D por Wolff
#include<cmath>
#include<iostream>
#include "Random.h"
using namespace std;

const int L=18;
const int L2=L*L;
const int k=1, J=1;

const int teq=(int)(10*pow(L/8.0,1.2));
const int tcorr=(int)(1*pow(L/8.0,1.2));
const double TC=2/log(1+sqrt(2.0));//TC=2.269;

//GRUPO1
//~ const double Tmin=2.2;
//~ const double Tmax=2.4;
//~ const double Tstep=(Tmax-Tmin)/20;
//~ const int Nmuestras=1000000;

//GRUPO2
const double Tmin=TC-1.3/pow(L/8.0,1.0);
const double Tmax=TC+3.0/pow(L/8.0,1.0);
const double Tstep=(Tmax-Tmin)/500;
const int Nmuestras=10000;

class SpinSystem{
private:
  int s[L][L]; int E,M; 
  int Cluster[L2][2]; int write;
public:
  void InicieArriba(void);
  void Muestre(void);
  void MuestreCluster(void);
  bool NoEstaEnElCLuster(int i,int j);
  int WolffCluster(Crandom & ran2,double p);
  double GetE(void){return (double)E;};
  double GetM(void){return fabs((double)M);};
};
void SpinSystem::InicieArriba(void){
  for(int i=0;i<L;i++)
    for(int j=0;j<L;j++)
      s[i][j]=1;
  M=L2; E=-2*L2;
  write=0;
}
void SpinSystem::Muestre(void){
  for(int i=0;i<L;i++){
    for(int j=0;j<L;j++)
      cout<<s[i][j]<<" ";
    cout<<endl;
  }
  cout<<"M="<<M<<" ; E="<<E<<endl<<endl;
}
bool SpinSystem::NoEstaEnElCLuster(int i,int j){
  int search; bool NoEsta=true;
  for(search=0;search<write;search++)
    if(Cluster[search][0]==i && Cluster[search][1]==j){
      NoEsta=false; break;
    }
  return NoEsta;
}

void SpinSystem::MuestreCluster(void){
  int i,j,search; bool NoEsta;
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      if(NoEstaEnElCLuster(i,j))
	cout<<"0 ";
      else
	cout<<"1 ";
    }cout<<endl;
  }
  cout<<"Tamaño del Cluster="<<write<<endl<<endl;
}
int SpinSystem::WolffCluster(Crandom & ran2,double p){
  int n,i,j,south,north,east,west,dE,sref,read;
  //Escojo el primer espín, y lo añado a la lista de leer;
  write=read=0;
  n=(int) (L2*ran2.r()); i=n/L; j=n%L;
  sref=s[i][j]; Cluster[write][0]=i; Cluster[write][1]=j; write++;
  while(read<write){
    //voltearlo
    i=Cluster[read][0];  j=Cluster[read][1];
    south=(i+1)%L; north=(i+L-1)%L; east=(j+1)%L; west=(j+L-1)%L; 
    dE=2*s[i][j]*(s[south][j]+s[north][j]+s[i][east]+s[i][west]);
    s[i][j]*=-1; E+=dE; M+=2*s[i][j];
    read++;
    //introducir nuevos en la lista 
    //North
    if(s[north][j]==sref) //Si está en la misma dirección,
      if(ran2.r()<p) //y los dados nos favorecen para agregarlo,
	if(NoEstaEnElCLuster(north,j)){//y no está en el cluster
	  Cluster[write][0]=north; Cluster[write][1]=j; write++; //agregarlo.
	}    
    //South
    if(s[south][j]==sref) //Si está en la misma dirección,
      if(ran2.r()<p) //y los dados nos favorecen para agregarlo,
	if(NoEstaEnElCLuster(south,j)){//y no está en el cluster
	  Cluster[write][0]=south; Cluster[write][1]=j; write++; //agregarlo.
	}
    //East
    if(s[i][east]==sref) //Si está en la misma dirección,
      if(ran2.r()<p) //y los dados nos favorecen para agregarlo,
	if(NoEstaEnElCLuster(i,east)){//y no está en el cluster
	  Cluster[write][0]=i; Cluster[write][1]=east; write++; //agregarlo.
	}
    //West
    if(s[i][west]==sref) //Si está en la misma dirección,
      if(ran2.r()<p) //y los dados nos favorecen para agregarlo,
	if(NoEstaEnElCLuster(i,west)){//y no está en el cluster
	  Cluster[write][0]=i; Cluster[write][1]=west; write++; //agregarlo al culster	
	}
  }
  return write;
}


int main(){
  SpinSystem Ising;
  Crandom ran2(0);
  int mcs,t,ccorr,c;
  double Mprom,M2prom,M4prom,Eprom,E2prom,E;
  double Cv,Xs,Ub;
  double T,p;

  /*
  //PARA PROBAR EL FUNCIONAMIENTO DEL ALGORITMO DE WOLFF
  T=TC; p=1-exp(-(2*J)/(k*T)); t=0;
  //Inicio el sistema
  Ising.InicieArriba();
  //Mostrar
  Ising.Muestre();
  Ising.MuestreCluster();
  //Hacer un paso de cluster
  t+=Ising.WolffCluster(ran2,p);
  t+=Ising.WolffCluster(ran2,p);
  t+=Ising.WolffCluster(ran2,p);
  //Mostrar
  Ising.Muestre();
  Ising.MuestreCluster();
  */
  /*
  //PARA ESTUDIAR EL TIEMPO DE EQUILIBRIO
  double M[teq];
  T=TC;  p=1-exp(-(2*J)/(k*T));
  //Inicio el sistema
  Ising.InicieArriba();
  //Iniciar los contadores;
  for(t=0;t<teq;t++)
    M[t]=0;
  //Cargarlos
  for(c=0;c<Nmuestras;c++){
    for(t=0;t<teq;t++){
      M[t]+=Ising.GetM();
      for(mcs=0;mcs<L2;)
	mcs+=Ising.WolffCluster(ran2,p);
    }
  }
  //Promediarlos
  for(t=0;t<teq;t++)
    M[t]/=Nmuestras;
  //Mostrarlos
  for(t=0;t<teq;t++)
    cout<<t<<" "<<M[t]<<endl;
  */
  
  //PARA CALCULAR LOS VALORES MEDIOS EN FUNCION DE LA TEMPERATURA
  double M;
  for(T=Tmin;T<Tmax;T+=Tstep){
    p=1-exp(-(2*J)/(k*T));
    //Inicio el sistema
    Ising.InicieArriba();
    //Equilibro
    for(t=0;t<teq;t++){
      for(mcs=0;mcs<L2;) //Un MCSS
	mcs+=Ising.WolffCluster(ran2,p);
    }
    //Inicio Acumuladores en cero
    Mprom=M2prom=M4prom=Eprom=E2prom=0;
    //Tomo muestras
    for(t=0;t<Nmuestras;t++){
      E=Ising.GetE();  M=Ising.GetM();
      Mprom+=M; M2prom+=M*M; M4prom+=M*M*M*M; Eprom+=E; E2prom+=E*E;
      for(ccorr=0;ccorr<tcorr;ccorr++)
	for(mcs=0;mcs<L2;) //Un MCSS
	mcs+=Ising.WolffCluster(ran2,p);
    }
    //Normalizo los acumuladores
    Mprom/=Nmuestras; M2prom/=Nmuestras; M4prom/=Nmuestras;
    Eprom/=Nmuestras; E2prom/=Nmuestras;
    //Calculo lo que me interesa
    Cv=1/(k*T*T)*(E2prom-Eprom*Eprom);
    Xs=1/(k*T)*(M2prom-Mprom*Mprom);
    Ub=1.0-1.0/3.0*(M4prom/(M2prom*M2prom));
    //Imprimo
    cout<<T<<" "<<Eprom<<" "<<Mprom<<" "<<Cv<<" "<<Xs<<" "<<Ub<<" "<<M2prom<<endl;
  }
  
  return 0;
}
