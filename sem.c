/*

A 1D SEM code for elastic wave propagation in a homogeneous medium. By Afsaneh Mohammadzaheri, 2014
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

static const int NGNOD=2 ;
static const double vpvs = 1.78;

void drive_lagrange(int , double *, double *, double **);
void ZELEGL (int , double *, double *);
void VALEPO(int , double , double *, double *, double *);
void DELEGL(int ,double * , double *, double *, double *);
void shapefunc(int , double *, double **, double **);


typedef struct{double *mass; 
               double mue,rho, vs;
               double *Jacobi, *Jacobian;
               int *cnvt; }spelement;

int main()
{
 FILE *out , *in;
  spelement *element;

  double **dldxi, **shape, **dshape;  
  double  *wgll, *xi, *U, *Unew, *Uold,srct;//[400] ;
  double  len_medium, Rho, vp, sigma, tmp, Invjacobi, dt,courant;
  double time, xa, dxmin, dx, vmax,dlen;
  int    it,tstp, tpeak;
  double *Mass, *fe, *xmsh, *x, *tmpx1, *Gforc; 
  int N_ord, num_elm, num_glo, N, ne;
  int i, ie, ic, j, ng, src_n,src_ne,rec_ne, rec_n,rec_ix;
  char c, srctyp;


  in=fopen("parameters","r");
  out=fopen("seismogram","w");

  fscanf(in,"%lf,%lf,%lf,%lf,%lf", &len_medium,&Rho,&vp,&courant,&time);
  c=fgetc(in);

  fscanf(in, "%d,%d,%d,%d,%d,%d\n",&N_ord,&num_elm,&src_ne,&src_n, &rec_ne,&rec_n);
fclose(in);
if(src_n>N_ord) src_n = N_ord;

  num_glo= N_ord*num_elm + 1	;			  //total number of nodes

 // printf("%d\n",N_ord);

   element = (spelement *) calloc(num_elm,sizeof(spelement));
  wgll = (double *)calloc((N_ord+1), sizeof(double));
  xi = (double *)calloc((N_ord+1), sizeof(double));
  fe = (double *)calloc((N_ord+1),sizeof(double));
  Mass=(double *)calloc(num_glo,sizeof(double));
  U=(double *)calloc(num_glo,sizeof(double));
  Uold=(double *)calloc(num_glo,sizeof(double));
  Unew=(double *)calloc(num_glo,sizeof(double));
  xmsh=(double *) malloc(num_glo*sizeof(double));
  x=(double *)calloc(num_glo,sizeof(double));
  tmpx1 = (double *)calloc((N_ord+1),sizeof(double));
  Gforc = (double *)calloc(num_glo,sizeof(double));
for(i=0;i<num_glo;i++){
   U[i]=0;
   Unew[i]=0;
   Uold[i]=0;
}

  // regular grid  
  i=0;
  xmsh[i] = 0;
printf("length of the medium= %lf\n",len_medium);
  for(i=1;i<num_elm+1;i++){ 
 //   i=i+1;
     dlen= len_medium/num_elm;
     xmsh[i]= xmsh[i-1] + dlen;
     if(xmsh[i]> len_medium)
   xmsh[i]= len_medium;

  }  // end of while

/****************************/
  
/* elasticity coefficient */
  
 for (ie=0;ie<num_elm; ie++){
     for(j=0;j<N_ord;j++){ 
       element[ie].vs = vp/ vpvs;
       element[ie].rho = Rho;
       element[ie].mue = Rho * (vp/ vpvs)*(vp/ vpvs);
     }
  }

  /* Lagrange derivation */

  N = N_ord+1;
  dldxi= (double **)malloc(N * sizeof(double*));
  for (i = 0 ; i < N; i++){
      dldxi[i] = (double *)malloc(N * sizeof(double));
  }
  //printf("%d\n",N_ord);
  drive_lagrange(N_ord, wgll, xi, dldxi);

/***************************************************/
/*           1D shape function             */

  shape= (double **)malloc(NGNOD * sizeof(double*));
  dshape= (double **)malloc(NGNOD * sizeof(double*));
  for (i = 0 ; i < NGNOD; i++){
      shape[i] = (double *)malloc (N * sizeof(double));
      dshape[i] = (double *)malloc (N * sizeof(double));
  }
  
  shapefunc(N, xi, shape, dshape);
  /****************************************/
  /*                Jacobi                */


  for(i = 0; i < num_elm; i++){
    element[i].Jacobi =  (double *) calloc(N,sizeof(double));
    element[i].Jacobian =  (double *) calloc(N,sizeof(double));
    element[i].mass =  (double *) calloc(N,sizeof(double)); 
    element[i].cnvt =  (int *) calloc(N,sizeof(int));
  }
  for(ie=0;ie<num_elm;ie++){
      for(i=0;i<N;i++){
          for(j=0;j<2;j++){
             xa=xmsh[ie+j];
             element[ie].Jacobi[i]=element[ie].Jacobi[i]+dshape[j][i]*xa;
         }
         element[ie].Jacobian[i]=element[ie].Jacobi[i];
     } 

  }
//exit(1);
/****************************************/

/* creating grid in transformed coordinates*/

N=N_ord-1;
for(ie=0; ie < num_elm;ie++){
   for(i=0;i< N+1;i++){
      for(j=0;j<NGNOD;j++){
         x[ie*(N+1)+i]=x[ie*(N+1)+i]+shape[j][i]*xmsh[ie+j];
      }
    }
}

for(j=0;j<NGNOD;j++){
  x[num_glo-1] = x[num_glo-1] + shape[j][N+1]*xmsh[num_elm-1+j];
}
dxmin=100000;
for(i=0;i<(num_glo-1);i++){
   dx=x[i+1]-x[i];
if (dx<= dxmin) dxmin=dx;
}


/*           ******************************                 */
             /* defining the time vector */
vmax= vp/vpvs;
  dt = courant * dxmin/vmax ;  //Stability Criterion
/* time step */

tstp= time/dt;
 //printf("%lf\t%d\t%lf\n",vmax,tstp,dt);

/*      defining source :Delta peak or exp */
  //if(srctyp == 'del'){

    tpeak = 1 ;
   // if(tpeak< tstp){
       srct=1/dt;
   //printf("%lf\n",srct);
//}
   /* else if{
       printf("tpeak must be smaller than %lf",time);
    }*/
//}
/*else if(srctyp == 'exp'){
t0=(1/frq);
for(i=0;i<tstp;i++){
   srct[i]=-exp(-((i-t0)*dt)**2);
}

} */ //end of else if 

/* source position : not necessary*/
//src_ix= (src_ne-1)*N_ord+src_n;
//srcloc=x[src_ix];

/* elemental mass matrix*/

N=N_ord+1;
for(ie=0;ie<num_elm;ie++){
   for(i=0;i<N;i++){
      element[ie].mass[i]= element[ie].Jacobian[i]*element[ie].rho*wgll[i] ;
    }
    //printf("\n");
}
//exit(1);
/*           ******************************                 */
/*                1D connectivity matrix                    */
 
for(ie=0;ie<num_elm;ie++){
    for(j=0;j<N;j++){
       element[ie].cnvt[j] = ie*(N-1)+j ;
    }
}

/*          ********************************              */
/*           assembling the elemental mass matrix */

  for(ie=0;ie<num_elm;ie++){
     for(i=0;i<N;i++){
       ic = element[ie].cnvt[i];
       Mass[ic]= Mass[ic]+ element[ie].mass[i];
       
    }
}

/*          ********************************              */
/*         calculating the forces at every point          */ 
it=0;
while(it < tstp){

for(i=0;i<num_glo;i++){
   Gforc[i]=0;
}
ic=0;
N=N_ord+1;
for(ie=0;ie<num_elm;ie++){

    for(i=0;i<N;i++){
      tmp = 0;
      fe[i]=0;
      tmpx1[i]=0;
        for(j=0;j<N;j++){
           ic = element[ie].cnvt[j];
           tmp= tmp + U[ic] * dldxi[j][i];     //du/dxi
        }
      Invjacobi= 1/(element[ie].Jacobian[i]);
      sigma = element[ie].mue* tmp *Invjacobi ;
      tmpx1[i] = element[ie].Jacobi[i]* sigma * Invjacobi ;
      }     // for i
for(i=0;i<N;i++){
  tmp = 0;
  for(j=0;j<N; j++){

     tmp = tmp + tmpx1[j] * dldxi[i][j] * wgll[j];
    
  }//
  
  fe[i] = - tmp;
  if(ie == src_ne && i == src_n && it == tpeak){
      fe[i] = fe[i] + 2 *(vp/vpvs) * srct;//[it]; 
  }
}
  for(i=0;i<N;i++){
  ic = element[ie].cnvt[i];
  Gforc[ic] = Gforc[ic]+fe[i];
 }
}   //ie
             /******************************************/
/* rigid and priodic boundary condition*/
             /*****************************************/
/*    solution to wave equation       */

for(i=0;i<num_glo;i++){
  Unew[i]= dt*dt * (1/Mass[i]) * Gforc[i] + 2*U[i] - Uold[i];
 if(Unew[i] != 0.000000 && Unew[i]!= -0.000000){
 }
}

//
rec_ix = (rec_ne)*N + rec_n;
 
for(i=0;i<num_glo;i++){
  Uold[i]=U[i];
  U[i]=0;
  U[i]=Unew[i];
  if(i == rec_ix)
    fprintf(out,"%lf\n",Unew[i]);
   Unew[i]=0;
}

it=it+1;
} // time loop

fclose(out);
 return 0;
}
/////////////////////////////////
void shapefunc(int NGLL, double *xigll, double **shape, double **dshape)
{

// generate the 1D shape functions and their derivatives (2 nodes)
  double xi, l1xi,l2xi, l1pxi, l2pxi ;
  int i;

  for(i=0; i<NGLL; i++)
  {

      xi=xigll[i];
      l1xi= 0.5*(1.0 - xi);
      l2xi= 0.5*(1.0 + xi);

      l1pxi= -0.5 ;
      l2pxi= 0.5 ;

      shape[0][i] = l1xi ;
      shape[1][i] = l2xi ;

      dshape[0][i] = l1pxi ;
      dshape[1][i] = l2pxi ;

}  // end of for

}  // end of shape1D function

void drive_lagrange(int N, double *wgllx, double *xi, double **dl)
{

  int i,j;
  double TEMP[N+1];
  double ET[N+1], VN[N+1], LAG[N+1];
  
  ZELEGL(N,ET,VN);
  for(i=0;i< N+1; i++){
     xi[i]= ET[i];
     wgllx[i] = 2/((N+1)*N* VN[i]*VN[i]) ;
  }

  for(i=0; i<N+1; i++){
     LAG[i]=1;
    DELEGL(N,xi,VN,LAG,TEMP);
    LAG[i]=0;
     for(j=0;j<N+1;j++){
      dl[i][j]=TEMP[j];
     TEMP[j]=0;
    }
  }
   


  } // end of drive_lagranger 

//////////////////////
void ZELEGL (int N, double *ET, double *VN)
{
      int N2, i, it ;
      double SN, X, Y, DY, D2Y, PI, C, ETX;
   
      N2 = 0 ;
      if (N == 0) exit(1);
                                                                        
         N2 = (N - 1) / 2 ;                                                  
         SN = 1.0 * (2 * N - 4 * N2 - 3) ;                                       
         ET[0] = -1.0;                                                  
         ET[N] = 1;                                                   
         VN[0] = SN ;                                                    
         VN[N] = 1 ;                                                  
      if (N == 1) exit(1);                                              
         ET[N2+1] = 0;  
         X = 0;                                                       
      VALEPO(N, X, &Y, &DY, &D2Y);                                         
         VN[N2+1] = Y  ; 
      if(N == 2)  exit(1);                                               
                                                                        
         PI = 3.14159265358979323846 ;                                  
         C  = PI/N ;                                            
      for(i=1; i<= N2;i++){                                              
         ETX = cos(C*i);                                        
      for(it=1;it<9; it++){                                                       
         VALEPO(N,ETX,&Y,&DY,&D2Y);                                       
         ETX = ETX-DY/D2Y ;                                              
      }                                                          
         ET[i] = -ETX;                                                   
         ET[N-i] = ETX ;                                                 
         VN[i] = Y*SN;                                                   
         VN[N-i] = Y;  
//printf("%lf\t%lf\n",ET[i],ET[N-i]);                                                   
     }                                                         
}    // end of ZELEGL     
                                                  
     
///////////////////////////////////////////////////////////////
void VALEPO(int N, double X, double* Y, double* DY, double * D2Y)
{
 int I ;
 double YP,DYP,D2YP,C1,C2,C4,YM,DYM,D2YM;

         *Y   = 1;                                                     
         *DY  = 0;                                                     
         *D2Y = 0;                                                     
      if (N == 0) exit(1);                                              
                                                                        
         *Y   = X ;                                                        
         *DY  = 1 ;                                                     
         *D2Y = 0 ;                                                     
      if(N == 1) exit(1);                                               
                                                                  
         YP   = 1;                                                    
         DYP  = 0;                                                    
         D2YP = 0;                                                    
      for(I=2;I<= N;I++){                                                        
         C1 = I ;                                                 
         C2 = 2*C1-1 ;                                              
         C4 = C1-1 ;                                                   
         YM = *Y ;                                                         
         *Y  = (C2*X*YM-C4*YP)/C1 ;                                        
         YP = YM ;                                                       
         DYM  = *DY ;                                                     
         *DY   = (C2*X*DYM-C4*DYP+C2*YP)/C1 ;                               
         DYP  = DYM ;                                                     
         D2YM = *D2Y ;                                                    
         *D2Y  = (C2*X*D2YM-C4*D2YP+2*C2*DYP)/C1 ;                      
         D2YP = D2YM ;                                                   
      }     
                                                     
 }  // end of VALEPO     

////////////////////////////////////////////////

void DELEGL(int N,double * ET, double *VN, double *QN, double *DQN)
{
  int I,J;
  double SU,VI,EI,VJ,EJ,DN,C ;

   DQN[0] = 0.0 ;   
/* for(I=0;I<N+1;I++){
              printf("%lf\n",QN[I]); 
}  */                                             
      if (N == 0) exit(1);                                             
                                                                        
      for(I=0; I<N+1; I++) {                                           
          SU = 0.0 ;                                                     
          VI = VN[I];                                                    
          EI = ET[I];                                                    
        for(J=0; J<N+1; J++){                                                       
          if (I == J) continue ;
          VJ = VN[J] ;                                                    
          EJ = ET[J] ;                                                    
          SU = SU+QN[J]/(VJ*(EI-EJ)) ;                                   
        }                                                          
        DQN[I] = VI*SU ;                                               
     }     // end for I                                                         
                                                                        
          DN = 1.0*N ;                                               
          C  = .25*DN*(DN+1.0) ; 
          DQN[0] = DQN[0]-C*QN[0] ;                                       
          DQN[N] = DQN[N]+C*QN[N];  
                                                                                                                                 
 }   // end of DELEGL

