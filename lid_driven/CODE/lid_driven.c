#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main(){
    FILE *file1 = fopen("streamlines_Re400.dat", "w");
    FILE *file2 = fopen("velocity_Re400.dat", "w");
    FILE *file3 = fopen("vorticity_Re400.dat", "w");
    FILE *file4 = fopen("u_velocity_Re400.dat", "w");
    FILE *file5 = fopen("v_velocity_Re400.dat", "w");

   int i, j, m, n,p;
   double dx,dy,x,y;

   m=128; //Number of points along x-direction
   n=128; //Number of points along y-direction

   p=(m-2)*(n-2); //Total number of interior points 

   dx= 1.0/(m-1); 
   dy= 1.0/(n-1);

   double beta=(dx/dy);
   double Re =400.0; //Reynold Number

   double  psi[m][n],omega[m][n],u[m][n],v[m][n],V[m][n];
   double  psi_pre[m][n], omega_pre[m][n];
   double  error_psi=0.0;
   double  error_omega=0.0;
   int iteration=0;

 //Initialization and Boundary Conditions   
    for ( j=0; j<n; j++)
    {
        for ( i=0; i<m; i++)
        {
            if (j==0)
            {
            u[i][j]=0.0;
            v[i][j]=0.0;
            psi[i][j]=0.0;
            }
            else if (j==(n-1))
            {
            u[i][j]=1.0;
            v[i][j]=0.0;
            psi[i][j]=0.0;                
            }
             else if (i==0)
            {
            u[i][j]=0.0;
            v[i][j]=0.0;
            psi[i][j]=0.0;                
            }
             else if (i==(m-1))
            {
            u[i][j]=0.0;
            v[i][j]=0.0;
            psi[i][j]=0.0;                
            }
            else
            {
            u[i][j]=0.0;
            v[i][j]=0.0;
            psi[i][j]=0.0;
            }           
        }       
    }
    for ( j=0; j<n; j++)
    {
        for (i=0; i<m; i++)
        {
            if (j==0)
            {
                omega[i][j]=(2.0/pow(dy,2))*(psi[i][j]-psi[i][j+1]);
            }
            else if (j==(n-1))
            {
                omega[i][j]=((-2.0)/pow(dy,2))*(psi[i][j-1]-psi[i][j]+dy*u[i][j]);
            }
            else if (i==0)
            {
                omega[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i+1][j]);
            }
            else if (i==(m-1))
            {
                omega[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i-1][j]);
            }
            else
            {
                omega[i][j]=0.0;
            }           
        }
        
    }
//Gauss-Siedel method
    do
    {
        for ( j= 0; j<n; j++)
        {  
           for ( i=0; i<m; i++)
           {
            psi_pre[i][j]=psi[i][j];
            omega_pre[i][j]=omega[i][j];
           }        
        }       
       for ( j=1; j<(n-1); j++)
        {
        for ( i=1; i<(m-1); i++)
        {
            //Formula for Stream function
            psi[i][j]=(1.0/(2.0*(1+pow(beta,2))))*(psi[i+1][j]+psi[i-1][j] + (psi[i][j+1]+psi[i][j-1]+pow(dx,2)*omega[i][j])*pow(beta,2));
                                                                                             
        }
        
       }

       for ( j=1; j<(n-1); j++)
       {
          for ( i=1; i<(m-1); i++)
          {
            //Formula for vorticity calculaions
            double T_1=(1.0-(psi[i][j+1]-psi[i][j-1])*((beta*Re)/4))*omega[i+1][j];
            double T_2=(1.0+(psi[i][j+1]-psi[i][j-1])*((beta*Re)/4))*omega[i-1][j];
            double T_3=(1.0+(psi[i+1][j]-psi[i-1][j])*((Re)/(4*beta)))*pow(beta,2)*omega[i][j+1];
            double T_4=(1.0-(psi[i+1][j]-psi[i-1][j])*(Re/(4*beta)))*pow(beta,2)*omega[i][j-1];

            omega[i][j]=(0.5/(1+pow(beta,2)))*(T_1 + T_2 + T_3 + T_4);
          }
          
       } 
 //update vorticy at boundaries      
      for ( j=0; j<n; j++)
        {
        for (i=0; i<m; i++)
         {
            if (j==0)
            {
                omega[i][j]=(2.0/pow(dy,2))*(psi[i][j]-psi[i][j+1]);
            }
            else if (j==(n-1))
            {
                omega[i][j]=((-2.0)/pow(dy,2))*(psi[i][j-1]-psi[i][j]+dy*u[i][j]);
            }
            else if (i==0)
            {
                omega[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i+1][j]);
            }
            else if (i==(m-1))
            {
                omega[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i-1][j]);
            } 
          }
        }
//Calculation of errors
       error_psi=0.0;
       error_omega=0.0;
       
       for ( j=1; j<(n-1); j++)
       {
          for ( i=1; i<(m-1); i++)
          {

            error_psi=error_psi+pow((psi[i][j]-psi_pre[i][j]),2.0);
            error_omega=error_omega+pow((omega[i][j]-omega_pre[i][j]),2.0);

           
          }
          
       }
       error_psi=sqrt(error_psi/p);
       error_omega=sqrt(error_omega/p);
//print of iteration and errors
       printf("iteration=%d\t",iteration);
       printf("error_psi=%.10lf\terror_omega=%.10lf\n",error_psi,error_omega);
       iteration++;
    }     
    while (error_psi>1.0e-6 || error_omega>1.0e-6);
    
    for ( j=1; j<(n-1); j++)
    {
        for ( i=1; i<(m-1); i++)
        {
            //Calculation for velocity
            u[i][j]=(0.5/dy)*(psi[i][j+1]-psi[i][j-1]);
            v[i][j]=(-0.5/dx)*(psi[i+1][j]-psi[i-1][j]); 
            V[i][j]=sqrt(pow(u[i][j],2)+pow(v[i][j],2));
        }
        
    }
//print in file
    fprintf(file1,"TITLE =STREAMLINES \n");
    fprintf(file1,"VARIABLES= X,Y,PSI\n");  
    fprintf(file1,"ZONE T=BLOCK I=128 J=128,F=POINT\n");

    fprintf(file2,"TITLE =VELOCITY \n");
    fprintf(file2,"VARIABLES=X,Y,U,V \n"); 
    fprintf(file2,"ZONE T=BLOCK I=128 J=128,F=POINT\n");
     

    fprintf(file3,"TITLE = VORTICITY\n");
    fprintf(file3,"VARIABLES= X,Y,VORTICITY\n"); 
    fprintf(file3,"ZONE T=BLOCK I=128 J=128,F=POINT\n"); 

    fprintf(file4,"TITLE = U-VELOCITY\n");
    fprintf(file4,"VARIABLES=Y,U \n");  
    fprintf(file4,"ZONE T=BLOCK I=128 J=128,F=POINT\n");

    fprintf(file5,"TITLE = V-VELOCITY\n");
    fprintf(file5,"VARIABLES=X,V \n");  
    fprintf(file5,"ZONE T=BLOCK I=128 J=128,F=POINT\n");
    
       for ( j=0; j<n; j++)
        { 
             y=j*dy;
        for ( i=0; i<m; i++)
        {
            x=i*dx;
//Results
             // printf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi[i][j],omega[i][j]);
             fprintf(file1,"%lf\t%lf\t%lf\n",x,y,psi[i][j]);
             fprintf(file2,"%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j]);
             fprintf(file3,"%lf\t%lf\t%lf\n",x,y,omega[i][j]);
             fprintf(file4,"%lf\t%lf\n",y,u[63][j]);
             fprintf(file5,"%lf\t%lf\n",x,v[i][63]);
                          
        }        
    }

    fclose(file1);
    fclose(file2);
    fclose(file3);
    fclose(file4);
    fclose(file5);
    return 0;    
}