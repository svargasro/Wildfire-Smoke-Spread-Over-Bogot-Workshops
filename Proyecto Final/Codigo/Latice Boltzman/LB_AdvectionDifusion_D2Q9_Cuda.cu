//LB for Waves D2Q5 on CUDA
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 128
#define Ly 128
#define N 32 //Threads per Block
const int M=(Lx*Ly+N-1)/N; //Blocks per Grid
#define Q 5
const int ArraySize=Lx*Ly*Q;

const float W0=1.0/3;

const float C=0.5; // C<0.707 cells/click
const float C2=C*C;
const float AUX0=1-3*C2*(1-W0);

const float tau=0.5;
const float Utau=1.0/tau;
const float UmUtau=1-Utau;

//------------ PROGRAMMING ON THE DEVICE ----------------
//---------------Constants (Symbols)----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];
__constant__ float d_C[3];   // d_C[0]=C,  d_C[1]=C2,  d_C[2]=AUX, 
__constant__ float d_tau[3]; // d_tau[0]=tau,  d_tau[1]=Utau,  d_tau[2]=UmUtau, 
//----------Functions called by the device itself
//Data index
__device__ int d_n(int ix,int iy,int i){
  return ix*Ly*Q+iy*Q+i;  
}
//Macroscopic Fields
__device__ float d_rho(int ix,int iy,float *d_f){
  float sum=0; int i,n0;
  for(i=0;i<Q;i++){
    n0=d_n(ix,iy,i); sum+=d_f[n0];
  }
  return sum;
}
__device__ float d_Jx(int ix,int iy,float *d_f){
  float sum=0; int i,n0;
  for(i=0;i<Q;i++){
    n0=d_n(ix,iy,i); sum+=d_Vx[i]*d_f[n0];
  }
  return sum;
}
__device__ float d_Jy(int ix,int iy,float *d_f){
  float sum=0; int i,n0;
  for(i=0;i<Q;i++){
    n0=d_n(ix,iy,i); sum+=d_Vy[i]*d_f[n0];
  }
  return sum;
}
//Equilibrium Functions
__device__ float d_feq(float rho0,float Jx0,float Jy0,int i){
  return 3*d_w[i]*(d_C[1]*rho0+d_Vx[i]*Jx0+d_Vy[i]*Jy0);
}
__device__ float d_f0eq(float rho0,float Jx0,float Jy0){
  return rho0*d_C[2];
}
//---------------------KERNELS----------------------------
__global__ void d_Collision(float *d_f,float *d_fnew,float *d_test){
  //Define internal registers
  int icell,ix,iy,i,n0;  float rho0,Jx0,Jy0;
  //Find which thread an which cell should I work
  icell=blockIdx.x*blockDim.x+threadIdx.x;
  ix=icell/Ly; iy=icell%Ly;
  //Compute the macroscopic fields
  rho0=d_rho(ix,iy,d_f); //rho
  Jx0=d_Jx(ix,iy,d_f);   //Jx0
  Jy0=d_Jy(ix,iy,d_f);   //Jy0
  if (ix==Lx/4 && iy== Lx/3) d_test[0]=Jy0; //OJO
  //Collide and compute fnew
  n0=d_n(ix,iy,0); d_fnew[n0]=d_tau[2]*d_f[n0]+d_tau[1]*d_f0eq(rho0,Jx0,Jy0);
  for(i=1;i<Q;i++){ //on each direction
    n0=d_n(ix,iy,i); d_fnew[n0]=d_tau[2]*d_f[n0]+d_tau[1]*d_feq(rho0,Jx0,Jy0,i);
  }
}
__global__ void d_ImposeFields(float *d_f,float *d_fnew,float RhoSource){
  //Define internal registers
  int ix,iy,i,n0;  float rho0,Jx0,Jy0;
  //There is only one thread and for one cell in the center
  ix=Lx/2; iy=Ly/2;
  //Compute the macroscopic fields
  rho0=RhoSource; //rho
  Jx0=d_Jx(ix,iy,d_f);   //Jx0
  Jy0=d_Jy(ix,iy,d_f);   //Jy0
  //Collide and compute fnew
  n0=d_n(ix,iy,0); d_fnew[n0]=d_f0eq(rho0,Jx0,Jy0);
  for(i=1;i<Q;i++){ //on each direction
    n0=d_n(ix,iy,i); d_fnew[n0]=d_feq(rho0,Jx0,Jy0,i);
  }
}
__global__ void d_Advection(float *d_f,float *d_fnew){
  //Define internal registers
  int icell,ix,iy,i,ixnext,iynext,n0,n0next;
  //Find which thread an which cell should I work
  icell=blockIdx.x*blockDim.x+threadIdx.x;
  ix=icell/Ly; iy=icell%Ly;
  //Move the contents to the neighboring cells
  for(i=0;i<Q;i++){ //on each direction
    ixnext=(ix+d_Vx[i]+Lx)%Lx; iynext=(iy+d_Vy[i]+Ly)%Ly;//periodic boundaries
    n0=d_n(ix,iy,i); n0next=d_n(ixnext,iynext,i);
    d_f[n0next]=d_fnew[n0]; 
  }
}
//------------ PROGRAMMING ON THE HOST ----------------
//-------------LatticeBoltzmann class------------
class LatticeBoltzmann{
private:
  float h_C[3];   // h_C[0]=C,  h_C[1]=C2,  h_C[2]=AUX, 
  float h_tau[3]; // h_tau[0]=tau,  h_tau[1]=Utau,  h_tau[2]=UmUtau, 
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]
  float *h_f,*h_fnew;  float *d_f,*d_fnew;// f[ix][iy][i]
  float *h_Test,*d_Test; //Just for tests
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  float h_rho(int ix,int iy);
  float h_feq(float rho0,float Jx0,float Jy0,int i);
  void Start(float rho0,float Jx0,float Jy0);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(const char * NameFile);
  void ShowTest(void);
};  
LatticeBoltzmann::LatticeBoltzmann(void){
  //CONSTANTS(d_Symbols)
  //---Charge constantes on the Host-----------------
  //running constants
  h_C[0]=C;  h_C[1]=C2;  h_C[2]=AUX0;
  h_tau[0]=tau;  h_tau[1]=Utau;  h_tau[2]=UmUtau;
  //weights
  h_w[0]=W0; h_w[1]=h_w[2]=h_w[3]=h_w[4]=(1.0-W0)/4;
  //velocity vectors
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;  
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //------Send to the Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_C,h_C,3*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_tau,h_tau,3*sizeof(float),0,cudaMemcpyHostToDevice);
  //DISTRIBUTION FUNCTIONS
  //Build the dynamic matrices on the host
  h_f=new float [ArraySize];  h_fnew=new float [ArraySize];
  //Build the dynamic matrices on the device
  cudaMalloc((void**) &d_f,ArraySize*sizeof(float));
  cudaMalloc((void**) &d_fnew,ArraySize*sizeof(float));
  //Test variables
   h_Test=new float [1]; cudaMalloc((void**) &d_Test,sizeof(float));
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] h_f;  delete[] h_fnew;
  cudaFree(d_f);  cudaFree(d_fnew);
  //Test variables
  delete[] h_Test; cudaFree(d_Test);
}
float LatticeBoltzmann::h_rho(int ix,int iy){
  //Note: Please import data from device before running
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i); sum+=h_fnew[n0];
  }
  return sum;
}
float LatticeBoltzmann::h_feq(float rho0,float Jx0,float Jy0,int i){
  if(i>0)
    return 3*h_w[i]*(C2*rho0+h_Vx[i]*Jx0+h_Vy[i]*Jy0);
  else
    return rho0*AUX0;
}
void LatticeBoltzmann::Start(float rho0,float Jx0,float Jy0){
  int ix,iy,i,n0;
  //Charge on the Host
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i); h_f[n0]=h_feq(rho0,Jx0,Jy0,i);
      }
  //Send to the Device
  cudaMemcpy(d_f,h_f,ArraySize*sizeof(float),cudaMemcpyHostToDevice);
}  
void LatticeBoltzmann::Collision(void){
  //Do everything on the Device
  dim3 ThreadsPerBlock(N,1,1);
  dim3 BlocksPerGrid(M,1,1);
  d_Collision<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew,d_test); //OJO, quitar test
  cudaMemcpy(h_Test,d_Test,sizeof(float),cudaMemcpyDeviceToHost); //OJO
  cout<<"Test="<<h_Test[0]<<endl; //OJO
}
void LatticeBoltzmann::ImposeFields(int t){
  float lambda=10, omega=2*M_PI/lambda*C;
  float RhoSource=10*sin(omega*t);
  dim3 ThreadsPerBlock(1,1,1); //A single thread (in this case)
  dim3 BlocksPerGrid(1,1,1);
  d_ImposeFields<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew,RhoSource);
}
void LatticeBoltzmann::Advection(void){
  //Do everything on the Device
  dim3 ThreadsPerBlock(N,1,1);
  dim3 BlocksPerGrid(M,1,1);
  d_Advection<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew);
}
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix,iy;
  //Bring back the data from Device to Host
  cudaMemcpy(h_fnew,d_fnew,ArraySize*sizeof(float),cudaMemcpyDeviceToHost);
  //Print for gnuplot splot
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=h_rho(ix,iy);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}
void LatticeBoltzmann::ShowTest(void){
  //Bring back test data from Device to Host
  cudaMemcpy(h_Test,d_Test,sizeof(float),cudaMemcpyDeviceToHost);
  cout<<"Test="<<h_Test[0]<<endl;
}

//--------------- GLOBAL FUNCTIONS ------------

int main(void){
  LatticeBoltzmann Waves;
  int t,tmax=100;
  float rho0=0,Jx0=0,Jy0=0;

  //Start
  Waves.Start(rho0,Jx0,Jy0);
  //Run
  for(t=0;t<tmax;t++){
    Waves.Collision();
    Waves.ImposeFields(t);
    Waves.Advection();
  }
  //Print Results
  Waves.Print("Waves_CUDA.dat");
 
  return 0;
}  
