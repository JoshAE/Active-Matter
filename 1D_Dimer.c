/* Active Psudeo 1D particles with a box search method*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2.h>
#include <g2_X11.h>
#include <unistd.h>

#define N 1000 // Number of particles
#define sizex 200
#define sizey 40
#define wsizex 1000
#define wsizey 200
#define bx 49
#define by 9
#define del_t 0.01 //time step
#define k 50
#define r 1


void force(void);

double x[N],xold[N];
double y[N];
double vx[N];
double vy[N], yold[N];
double Fx[N];
double Fy[N], vxout,vyout;
int n,box,nn[4][bx*by],bsrch,iy,ix,bn;


double head[bx*by],list[N],nx,actbox[N],ny; //Setup for linked list


int main(void)
{	int id,iout=0; //graphics
    int i,j,a,b,c; //integers for loops





    //Partuicle speeds, forces and positions

    double dx,dy,d;
    int bn;
    double phi,bin[100];
    double Fdrive, Pe;
    double Vo=0,rin;
    double scale, DT, DR,D;

    srand48(1234);
    FILE *fp;


    Pe=200;

    phi=N*M_PI*r*r/(sizex*sizey);
    printf("%lf\n",phi);
    DR=Vo/(Pe*r);
    D=1/Pe;


    scale=wsizex/sizex;

    nx=(int)(5*N/6);
    ny=(int)(N/6);
    rin=sizex/(2*nx);
    n=0;


    id=g2_open_X11(wsizex,wsizey);
    g2_set_auto_flush(id,0);


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



//Main time loop
    for (i=1;i<=1000000; i++)
    {

        if(i>1000) {
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

        if (i > iout) {
            for (a = 0; a < N; a += 2) {
                b = a + 1;
                g2_pen(id, 0);
                g2_filled_circle(id, scale * xold[a], scale * yold[a], scale * r);
                g2_pen(id, 3);
                g2_filled_circle(id, scale * x[a], scale * y[a], scale * r);
                xold[a] = x[a];
                yold[a] = y[a];

                g2_pen(id, 0);
                g2_filled_circle(id, scale * xold[b], scale * yold[b], scale * r);
                g2_pen(id, 19);
                g2_filled_circle(id, scale * x[b], scale * y[b], scale * r);
                xold[b] = x[b];
                yold[b] = y[b];
            }

            g2_flush(id);
            iout = iout + 10;
        }

    }
    for (a=0;a<100;a++)
    {
        bin[a]=0;
    }

    for (a=0;a<N;a++)
    {
        bn=(int)(x[a]-0.5);
        bin[bn]+=1;

    }

    fp=fopen("Dimer_Velocity","w");
    vxout=0;
    vyout=0;

    for (a=0;a<N;a++)
    {
        vxout+=vx[a];
        vyout+=vy[a];

    }
    fprintf(fp,"%lf %lf\n",vxout,vyout);

}





void force(void)
{
    double dx,dy,d;
    double phi;
    double Fdrive, Pe;
    int a,b,c;



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





