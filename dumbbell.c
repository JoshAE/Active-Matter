* Active 2D particles with a box search method*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define N 5000 // Number of particles

#define wsizex 1000
#define wsizey 1000
#define bx 99
#define by 99
#define del_t 0.01 //time step

#define k 50
#define r 1


void force(void);

double x[N],xold[N],phi;
double y[N];   	
double vx[N];
double sizex,sizey;
char fname[100];

double vy[N], yold[N];
double Fx[N];
double Fy[N], vxout,vyout;

int n,box,nn[4][bx*by],bsrch,iy,ix,bn;


double head[bx*by],list[N],nx,actbox[N],ny; //Setup for linked list


int main(void)
{	int id,iout=0; //graphics
	int i,j,a,b,c; //integers for loops

	phi=0.5;
    sizex=sqrt(N*M_PI*r*r/phi);
    sizey=sizex;



	//Partuicle speeds, forces and positions

	double dx,dy,d;
	int bn;
	double bin[100];
	double Fdrive, Pe;
	double Vo=0,rin;
	double scale, DT, DR,D;


srand48(1234);
FILE *fp;





scale=wsizex/sizex;

nx=(int)((sizex-1)/2);
ny=(int)(2*nx);
rin=(sizex-2)/(2*nx);
n=0;




//Main time loop
Pe=50;

    DR=Vo/(Pe*r);
    D=1/Pe;

    for(ix=0;ix<bx;ix++)
    {
        for(iy=0;iy<by;iy++)
        {
            i=ix+bx*iy;

            nn[0][i]=((ix-1+bx)%bx)+bx*((iy+1)%by);
            nn[1][i]=ix+bx*((iy+1)%by);
            nn[2][i]=((ix+1)%bx)+bx*((iy+1)%by);
            nn[3][i]=((ix+1)%bx)+bx*iy;

        }
    }


//initialising conditions
    for (a=0;a<N;a++)
    {
        vx[a]=0;
        vy[a]=0;

    }

//initialising particle positions

    for (b=0;b<ny;b++)
    {

        for (a=0;a<nx;a++)
        {
            if(n < N)
            {


                x[n]=rin+2*rin*a;
                y[n]=rin+2*rin*b;
                n++;
            }
        }

    }


    for (i = 1; i <= 1000000; i++) {

        if (i > 1000) {
            Vo = 1;
        }

        force();

        for (a = 0; a < N; a += 2) {
            b = a + 1;
            dx = x[b] - x[a];
            if (dx > sizex / 2) dx = dx - sizex;
            else if (dx < -sizex / 2) dx = dx + sizex;

            dy = y[b] - y[a];

            if (dy > sizey / 2) dy = dy - sizey;
            else if (dy < -sizey / 2) dy = dy + sizey;

            d = sqrt((dx * dx) + (dy * dy));

            vx[a] = Vo * (dx / d) + Fx[a];

            vy[a] = Vo * (dy / d) + Fy[a];


            x[a] += vx[a] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);
            y[a] += vy[a] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);

            if (x[a] < 0) {
                x[a] = x[a] + sizex;
            }
            if (x[a] > sizex) {
                x[a] = x[a] - sizex;
            }
            if (y[a] < 0) {
                y[a] = y[a] + sizey;
            }
            if (y[a] > sizey) {
                y[a] = y[a] - sizey;
            }

            vx[b] = Vo * (dx / d) + Fx[b];
            vy[b] = Vo * (dy / d) + Fy[b];
            x[b] += vx[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);
            y[b] += vy[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);

            if (x[b] < 0) {
                x[b] = x[b] + sizex;
            }
            if (x[b] > sizex) {
                x[b] = x[b] - sizex;
            }
            if (y[b] < 0) {
                y[b] = y[b] + sizey;
            }
            if (y[b] > sizey) {
                y[b] = y[b] - sizey;
            }
        }


        if (i == 999999) {

            sprintf(fname, "DATA_N=2500/data_%lf.txt", Pe);
            fp=fopen(fname,"w");
            for (c=0;c<N;c+=2)
            {
                fprintf(fp,"%010.5f %010.5f %010.5f %010.5f\n",x[c],y[c],x[c+1],y[c+1]);
            }
        }

    }

}









void force(void)
{
double dx,dy,d;	
double phi;
double Fdrive;
int a,b,c;

phi=0.5;
sizex=sqrt(N*M_PI*r*r/phi);
sizey=sizex;


for (c = 0; c < N; c++) {
            Fx[c] = 0;
            Fy[c] = 0;
        }


        for (a = 0; a < by * bx; a++) {
            head[a] = N;
        }


        for (a = 0; a < N; a++) {
            box = ((int) (x[a] * bx / (sizex))) + bx * ((int) (y[a] * by / (sizey)));
            //printf("%d %d %lf %lf\n",box,a,x[a],y[a]);
            list[a] = head[box];
            head[box] = a;
            actbox[a] = box;
        }

        for (b = 0; b < N; b++) {
            a = list[b];

            while (a < N) {

                dx = x[b] - x[a];
                if (dx > sizex / 2) dx = dx - sizex;
                else if (dx < -sizex / 2) dx = dx + sizex;

                dy = y[b] - y[a];

                if (dy > sizey / 2) dy = dy - sizey;
                else if (dy < -sizey / 2) dy = dy + sizey;

                d = sqrt((dx * dx) + (dy * dy));

                if (d < 2 * r) {
                    Fx[b] += k * (2 * r - d) * (dx / d);
                    Fx[a] += -k * (2 * r - d) * (dx / d);
                    Fy[b] += k * (2 * r - d) * (dy / d);
                    Fy[a] += -k * (2 * r - d) * (dy / d);

                }

                a = list[a];
            }


            box = actbox[b];
            for (bsrch = 0; bsrch < 4; bsrch++) {
                a = head[nn[bsrch][box]];
                while (a < N) {
                    dx = x[b] - x[a];
                    if (dx > sizex / 2) dx = dx - sizex;
                    else if (dx < -sizex / 2) dx = dx + sizex;

                    dy = y[b] - y[a];

                    if (dy > sizey / 2) dy = dy - sizey;
                    else if (dy < -sizey / 2) dy = dy + sizey;

                    d = sqrt((dx * dx) + (dy * dy));

                    if (d < 2 * r) {
                        Fx[b] += k * (2 * r - d) * (dx / d);
                        Fx[a] += -k * (2 * r - d) * (dx / d);
                        Fy[b] += k * (2 * r - d) * (dy / d);
                        Fy[a] += -k * (2 * r - d) * (dy / d);

                    }
                    a = list[a];
                }
            }

        }

        for (a = 0; a < N; a += 2) {
            b = a + 1;
            dx = x[b] - x[a];
            if (dx > sizex / 2) dx = dx - sizex;
            else if (dx < -sizex / 2) dx = dx + sizex;

            dy = y[b] - y[a];

            if (dy > sizey / 2) dy = dy - sizey;
            else if (dy < -sizey / 2) dy = dy + sizey;

            d = sqrt((dx * dx) + (dy * dy));

            if (d > 2 * r) {
                Fx[b] += k * (2 * r - d) * (dx / d);
                Fx[a] += -k * (2 * r - d) * (dx / d);
                Fy[b] += k * (2 * r - d) * (dy / d);
                Fy[a] += -k * (2 * r - d) * (dy / d);
            }
        }
}

