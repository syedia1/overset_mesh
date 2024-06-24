#include <iostream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

class DataStructure {
public: 
    int Nx, Ny;
    double dx, dy, Lx, Ly;
    double dt;
    int nVar;
    double ***Var;
    double ***VarOld;
    double ***Ff;
    double *Init;
    double *residual;
    double ulid;
    int Dim;
    double volp;
    double nu;
    double rho;
}; 

void CopyNewtoOld(DataStructure *rect)
{
    for(int k=0; k<rect->nVar; k++){
      for(int i = 0; i<rect->Nx+2; i++)
       {
           for(int j = 0; j<rect->Ny+2; j++)
                {
                     rect->VarOld[k][i][j] = rect->Var[k][i][j];
                }
       }
    }
}

void allocate(DataStructure *rect)
{
   rect->Nx = 100;
   rect->Ny = 100;
   rect->Lx = 1.0;
   rect->Ly = 1.0;
   rect->dx = (rect->Lx/rect->Nx);
   rect->dy = (rect->Ly/rect->Ny);
   rect->dt =0.001;//0.5*(rect->dx*rect->dx * rect->dy*rect->dy)/(rect->alpha * (rect->dx*rect->dx + rect->dy*rect->dy));

   rect->nVar = 3;
   rect->Init = new double[rect->nVar];
   rect->ulid = 1.0;
   rect->Var = new double**[rect->nVar];
   rect->VarOld = new double**[rect->nVar];
   rect->Dim = 2;
   rect->volp = rect->dx*rect->dy;
   rect->nu = 0.01;
   rect->rho = 1.0;
   for(int i = 0; i<rect->nVar; i++)
   {
       rect->Init[i] = 0.0;
       rect->Var[i] =  new double*[rect->Nx+2];
       rect->VarOld[i] =  new double*[rect->Nx+2];
       for(int j = 0; j<(rect->Nx+2); j++)
       {
            rect->Var[i][j]=new double[rect->Ny+2];
            rect->VarOld[i][j]=new double[rect->Ny+2];
       }
   }
   rect->residual = new double[rect->nVar];
   rect->Ff = new double**[2*rect->Dim]; //No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S
   for(int k = 0; k<(2*rect->Dim); k++)
   {
       rect->Ff[k] = new double*[rect->Nx+2];
           for(int j = 0; j<(rect->Nx+2); j++)
             {
                rect->Ff[k][j]=new double[rect->Ny+2];
             }
   }
   
   
}
void ApplyBC(DataStructure *rect, int k)
{
    switch (k)
    {
        case 0: //U
    
                    for(int j = 1; j<rect->Ny+1; j++)
                    {
                        rect->Var[k][0][j]= 2*(0.0) - rect->Var[k][1][j];  //Left
                        rect->Var[k][rect->Nx+1][j]= 2*(0.0) - rect->Var[k][rect->Nx][j];  //Right
                    }
                    for(int i = 1; i<rect->Nx+1; i++)
                    {
                        rect->Var[k][i][rect->Ny+1]= 2*rect->ulid - rect->Var[k][i][rect->Ny]; //Top
                        rect->Var[k][i][0]= 2*(0.0) - rect->Var[k][i][1]; //Bottom
                    }
                    break;

        case 1: //V
     
                    for(int j = 1; j<rect->Ny+1; j++)
                    {
                        rect->Var[k][0][j]= 2*(0.0) - rect->Var[k][1][j];  //Left
                        rect->Var[k][rect->Nx+1][j]= 2*(0.0) - rect->Var[k][rect->Nx][j];  //Right
                    }
                    for(int i = 1; i<rect->Nx+1; i++)
                    {
                        rect->Var[k][i][rect->Ny+1]= 2*(0.0) - rect->Var[k][i][rect->Ny]; //Top
                        rect->Var[k][i][0]= 2*(0.0) - rect->Var[k][i][1]; //Bottom
                    }
                    break;
             
        case 2: //P (All Neumann Condition)
     
                    for(int j = 1; j<rect->Ny+1; j++)
                    {
                        rect->Var[k][0][j]= rect->Var[k][1][j];   //Left
                        rect->Var[k][rect->Nx+1][j]= rect->Var[k][rect->Nx][j];  //Right
                    }
                    for(int i = 1; i<rect->Nx+1; i++)
                    {
                        rect->Var[k][i][rect->Ny+1]= rect->Var[k][i][rect->Ny]; //Top
                        rect->Var[k][i][0]= rect->Var[k][i][1]; //Bottom
                    }
                    break;
    }
    
}

void LinearInterpolation(DataStructure *rect)
{
   for(int i = 1; i<rect->Nx+1; i++)
            {
                for(int j = 1; j<rect->Ny+1; j++)
                        {
                            
                                 rect->Ff[0][i][j] = (rect->Var[0][i][j] + rect->Var[0][i+1][j])*rect->dy*0.5;  //East Face
                                 rect->Ff[1][i][j] = (rect->Var[1][i][j] + rect->Var[1][i][j+1])*rect->dx*0.5;  //North Face
                                 rect->Ff[2][i][j] = -(rect->Var[0][i][j] + rect->Var[0][i-1][j])*rect->dy*0.5;  //West Face
                                 rect->Ff[3][i][j] = -(rect->Var[1][i][j] + rect->Var[1][i][j-1])*rect->dx*0.5;  //South Face                  
                        }
            } 
}

void initialize(DataStructure *rect)
{
   //initializing interior values
   for(int k=0; k<rect->nVar; k++){

            for(int i = 1; i<rect->Nx+1; i++)
            {
                for(int j = 1; j<rect->Ny+1; j++)
                        {
                            rect->Var[k][i][j] = rect->Init[k];
                        }
            }
    ApplyBC(rect, k);
   }
   
   CopyNewtoOld(rect);
   LinearInterpolation(rect);
}


void GetOutput(DataStructure *rect, std::string input)
{
    FILE *fp;
    fp = fopen(input.c_str(),"w");

     for(int k=0; k<rect->nVar; k++){

         fprintf(fp,"\n ########## Data for k = %d ############ \n", k);
            for(int i = 0; i<rect->Nx+2; i++)
            {
                for(int j = 0; j<rect->Ny+2; j++)
                        {
                            fprintf(fp, "%lf \t", rect->Var[k][i][j]); 
                        }
                        fprintf(fp,"\n");
            }

    }
    fclose(fp);
} 
void ExplicitSolve(DataStructure *rect)
{ 
   rect->residual = 0;
   for(int i = 1; i<rect->Nx+1; i++)
   {
       for(int j = 1; j<rect->Ny+1; j++)
            {
            }
   }
}
void SimpleUpwind(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k)
{
    double ue, uw, un, us;   
    double sum =0;
    if(rect->Ff[0][i][j]>=0)
    {
        ue= rect->Var[k][i][j];
        sum += rect->Ff[0][i][j];
    }
    else
        ue=rect->Var[k][i+1][j];
    
    if(rect->Ff[2][i][j]>=0){
        uw= rect->Var[k][i][j];
        sum += rect->Ff[2][i][j];
    }
    else
        uw=rect->Var[k][i-1][j];
    
    if(rect->Ff[1][i][j]>=0){
        un= rect->Var[k][i][j];
        sum += rect->Ff[1][i][j];
    }
    else
        un= rect->Var[k][i][j+1];
    
    if(rect->Ff[3][i][j]>=0){
        us= rect->Var[k][i][j];
        sum += rect->Ff[3][i][j];
    }
    else
        us= rect->Var[k][i][j-1];

    *Fc = (ue*rect->Ff[0][i][j] + uw*rect->Ff[2][i][j]  + un*rect->Ff[1][i][j]  +us*rect->Ff[3][i][j] );
    *ap_c = sum*rect->volp;
}


void Quick(DataStructure *rect, double *Fc, double *ap_c, int i, int j, int k)
{ 
    double ue, uw, un, us;   
    double sum =0;
    //East
    if(rect->Ff[0][i][j]>=0)
    {
        ue= 0.75*rect->Var[k][i][j] + 0.375*rect->Var[k][i+1][j] - 0.125*rect->Var[k][i-1][j];
        sum += 0.75*rect->Ff[0][i][j];
    }
    else
    {
        ue=0.75*rect->Var[k][i+1][j] + 0.375*rect->Var[k][i][j] - 0.125*rect->Var[k][i+2][j];
        sum += 0.375*rect->Ff[0][i][j];
    }
    //West
    if(rect->Ff[2][i][j]>=0){
        uw= 0.75*rect->Var[k][i][j] + 0.375*rect->Var[k][i-1][j] - 0.125*rect->Var[k][i+1][j];
        sum += 0.75*rect->Ff[2][i][j];
    }
    else
    {
        uw=0.75*rect->Var[k][i-1][j] + 0.375*rect->Var[k][i][j] - 0.125*rect->Var[k][i-2][j];
        sum += 0.375*rect->Ff[2][i][j];
    }
    //North
    if(rect->Ff[1][i][j]>=0){
        un= 0.75*rect->Var[k][i][j] + 0.375*rect->Var[k][i][j+1] - 0.125*rect->Var[k][i][j-1];
        sum += 0.75*rect->Ff[1][i][j];
    }
    else
    {
        un= 0.75*rect->Var[k][i][j+1] + 0.375*rect->Var[k][i][j] - 0.125*rect->Var[k][i][j+2];
        sum += 0.375*rect->Ff[1][i][j];
    }
    if(rect->Ff[3][i][j]>=0){
        us= 0.75*rect->Var[k][i][j] + 0.375*rect->Var[k][i][j-1] - 0.125*rect->Var[k][i][j+1];
        sum += 0.75*rect->Ff[3][i][j];
    }
    else
    {
        us= 0.75*rect->Var[k][i][j-1] + 0.375*rect->Var[k][i][j] - 0.125*rect->Var[k][i][j-2];
        sum += 0.375*rect->Ff[3][i][j];
    }
    *Fc = (ue*rect->Ff[0][i][j] + uw*rect->Ff[2][i][j]  + un*rect->Ff[1][i][j]  +us*rect->Ff[3][i][j] );
    *ap_c = sum*rect->volp;
}


void DiffusiveFlux(DataStructure *rect, double *Fd, double *ap_d, int i, int j, int k)
{
     *Fd=rect->volp*((rect->Var[k][i+1][j]-2.0*rect->Var[k][i][j]+rect->Var[k][i-1][j])/(rect->dx*rect->dx) + (rect->Var[k][i][j+1]-2.0*rect->Var[k][i][j]+rect->Var[k][i][j-1])/(rect->dy*rect->dy));
     *ap_d = -rect->volp*(2.0/(rect->dx * rect->dx) + 2.0/(rect->dy * rect->dy));
}

void UpdateFlux(DataStructure *rect)
{
   for(int i = 1; i<rect->Nx+1; i++)
            {
                for(int j = 1; j<rect->Ny+1; j++)
                        {
                                 rect->Ff[0][i][j] += -rect->dt/rect->rho*(rect->Var[2][i+1][j] - rect->Var[2][i][j])*rect->dy/rect->dx;  //East Face
                                 rect->Ff[1][i][j] += -rect->dt/rect->rho*(rect->Var[2][i][j+1] - rect->Var[2][i][j])*rect->dx/rect->dy;  //North Face
                                 rect->Ff[2][i][j] += -rect->dt/rect->rho*(rect->Var[2][i-1][j] - rect->Var[2][i][j])*rect->dy/rect->dx;  //West Face
                                 rect->Ff[3][i][j] += -rect->dt/rect->rho*(rect->Var[2][i][j-1] - rect->Var[2][i][j])*rect->dx/rect->dy;  //South Face                  
                        }
            } 
}

void PressureFlux(DataStructure *rect,double *Fp, int i, int j, int k)
{
    if(k==0)
     *Fp=-rect->dy/(rect->rho*2.0) * (rect->Var[2][i+1][j] - rect->Var[2][i-1][j]);
    else if(k==1)
     *Fp=-rect->dx/(rect->rho*2.0) * (rect->Var[2][i][j+1] - rect->Var[2][i][j-1]);   
    //rect->volp*((rect->Var[k][i+1][j]-2.0*rect->Var[k][i][j]+rect->Var[k][i-1][j])/(rect->dx*rect->dx) + (rect->Var[k][i][j+1]-2.0*rect->Var[k][i][j]+rect->Var[k][i][j-1])/(rect->dy*rect->dy));
     
}

void ImplicitSolve(DataStructure *rect)
{
    for(int k=0; k<rect->nVar; k++){
            rect->residual[k] = 0.0; 
    }
   double Fc, ap_c, Fd, ap_d, ap, R, rms, Fp;  //No of faces per cell = 4  i=0,1,2,3 -> E,N,W,S

  //Solving for U and V
    for(int k=0; k<2; k++){
  
      do
      {
          rms =0.0; 
           for(int i = 1; i<rect->Nx+1; i++)
          {
              for(int j = 1; j<rect->Ny+1; j++)
                   {
                      //  SimpleUpwind(rect, &Fc, &ap_c, i,j, k);  
                        Quick(rect, &Fc, &ap_c, i,j, k);  
                        DiffusiveFlux(rect, &Fd, &ap_d, i, j, k); 
                        R = -(rect->volp/rect->dt * (rect->Var[k][i][j] - rect->VarOld[k][i][j]) + Fc + (-rect->nu)*Fd) + Fp;
                        ap = rect->volp/rect->dt + ap_c + (-rect->nu)*ap_d;
                        rect->Var[k][i][j] = rect->Var[k][i][j] + R/ap;
                        rms = rms + R*R;
                   }
          }
         ApplyBC(rect, k);
         rms = sqrt(rms/(rect->Nx*rect->Ny));
      }while(rms > 1e-6);
    }
   
   LinearInterpolation(rect);
   
   //Solving for P
   int k=2;
      do
      {
           rms =0.0; 
           for(int i = 1; i<rect->Nx+1; i++)
          {
              for(int j = 1; j<rect->Ny+1; j++)
                   {
                        DiffusiveFlux(rect, &Fd, &ap_d, i, j, k); 
                        double LHS = Fd;
                        double RHS = rect->rho/rect->dt *(rect->Ff[0][i][j] + rect->Ff[1][i][j] + rect->Ff[2][i][j] +  rect->Ff[3][i][j]);
                        R = RHS-LHS;
                        ap = ap_d;
                        rect->Var[k][i][j] = rect->Var[k][i][j] + R/ap;
                        rms = rms + R*R;
                   }
          }
         ApplyBC(rect, k);
         rms = sqrt(rms/(rect->Nx*rect->Ny));
      }while(rms > 1e-6);
         
    // Correcting Velocity
           for(int i = 1; i<rect->Nx+1; i++)
          {
              for(int j = 1; j<rect->Ny+1; j++)
                   {
                     k=0;
                     rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt/rect->rho *(rect->Var[2][i+1][j]-rect->Var[2][i-1][j])/(2*rect->dx); 
                     k=1;
                     rect->Var[k][i][j] = rect->Var[k][i][j] - rect->dt/rect->rho *(rect->Var[2][i][j+1]-rect->Var[2][i][j-1])/(2*rect->dy); 
                     
                     rect->residual[0] += (rect->Var[0][i][j] - rect->VarOld[0][i][j])*(rect->Var[0][i][j] - rect->VarOld[0][i][j]);
                     rect->residual[1] += (rect->Var[1][i][j] - rect->VarOld[1][i][j])*(rect->Var[1][i][j] - rect->VarOld[1][i][j]);
                     rect->residual[2] += (rect->Var[2][i][j] - rect->VarOld[2][i][j])*(rect->Var[2][i][j] - rect->VarOld[2][i][j]);
                   }
          }
    ApplyBC(rect, 0);
    ApplyBC(rect, 1);
    UpdateFlux(rect);     
}

bool ConvergenceCheck(DataStructure *rect)
{
   double rms[rect->nVar];
   for(int k=0; k<rect->nVar; k++)
   {
   rms[k] = sqrt(rect->residual[k]/(rect->Nx*rect->Ny));
        rms[k] = rms[k]/rect->dt;
        cout << "\t" << rms[k];
   }
   cout<<endl;
   if(rms[0] > 1e-12 || rms[1] > 1e-6 || rms[2] > 1e-6)
   {
       CopyNewtoOld(rect);
       return false;
   }
   else
   {
       return true;
   }
   
}
void Solve(DataStructure *rect)
{
    int count =0;
    do
    {
        /* code */
        count++;
        ImplicitSolve(rect);
       // ExplicitSolve(rect);
      //  ApplyBC(rect,);
        cout<<count;
    }while (!ConvergenceCheck(rect));
}

int main()
{
    DataStructure rect;
    allocate(&rect);
    initialize(&rect);
    Solve(&rect);
    GetOutput(&rect, "output_Quick.dat");
    return 0;
}