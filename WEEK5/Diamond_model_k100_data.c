#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#define N 16

#define wsizex 1000
#define wsizey 1000
#define bx 99
#define by 99
#define del_t 0.001 //time step
#define k 100
#define r 1
#define r2 0.25 //r2 must be less than r1

void force(void);

double x[N],xold[N],phi,ext;
double radii[N];

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
{	int id,iout=0,iout2,count=0; //graphics
    int i,j,a,b,c; //integers for loops

    phi=0.5;
    sizex=sqrt((((N/2)*M_PI*r*r)+((N/2)*M_PI*r2*r2))/phi);
    sizey=sizex;

    ext=sqrt(((r2+r)*(r2+r))-(r*r));

    for (a=0;a<N;a+=4){
        radii[a]=r;
        radii[a+1]=r;
        radii[a+2]=r2;
        radii[a+3]=r2;
    }


    iout2=0.2*(10000/del_t);
    //Particle speeds, forces and positions

    double dx,dy,d;
    int bn;
    double bin[100];
    double Fdrive, Pe;
    double Vo=0,rinx,riny;
    double scale, DT, DR,D;


    srand48(1234);
    FILE *fp;




    scale=wsizex/sizex;

    nx=(int)(sizex-1);
    ny=(int)(sizey-1);
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

    for (b=0;b<ny;b+=3)
    {

        for (a=0;a<nx;a+=4)
        {
            if(n < N)
            {

                x[n]=a+1;
                x[n+1]=a+3;
                x[n+2]=a+2;
                x[n+3]=a+2;

                y[n]=b+1.5;
                y[n+1]=b+1.5;
                y[n+2]=b+1.5+ext;
                y[n+3]=b+1.5-ext;
                n+=4;
            }
        }

    }





    for (i = 0; i < 10000/del_t; i++) {


        if (i > 1000000 * del_t) {
            Vo = 1;
        }

        force();

        for (a = 0; a < N; a += 4) {
            b = a + 1;

            //main particles
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
            } else if (x[a] > sizex) {
                x[a] = x[a] - sizex;
            }
            if (y[a] < 0) {
                y[a] = y[a] + sizey;
            } else if (y[a] > sizey) {
                y[a] = y[a] - sizey;
            }

            vx[b] = Vo * (dx / d) + Fx[b];
            vy[b] = Vo * (dy / d) + Fy[b];
            x[b] += vx[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);
            y[b] += vy[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);

            if (x[b] < 0) {
                x[b] = x[b] + sizex;
            } else if (x[b] > sizex) {
                x[b] = x[b] - sizex;
            }
            if (y[b] < 0) {
                y[b] = y[b] + sizey;
            } else if (y[b] > sizey) {
                y[b] = y[b] - sizey;
            }





            //Small particles


            vx[a + 2] = Vo * (dx / d) + Fx[a + 2];

            vy[a + 2] = Vo * (dy / d) + Fy[a + 2];


            x[a + 2] += vx[a + 2] * del_t;
            y[a + 2] += vy[a + 2] * del_t;

            if (x[a + 2] < 0) {
                x[a + 2] = x[a + 2] + sizex;
            } else if (x[a + 2] > sizex) {
                x[a + 2] = x[a + 2] - sizex;
            }
            if (y[a + 2] < 0) {
                y[a + 2] = y[a + 2] + sizey;
            } else if (y[a + 2] > sizey) {
                y[a + 2] = y[a + 2] - sizey;
            }

            vx[b + 2] = Vo * (dx / d) + Fx[b + 2];
            vy[b + 2] = Vo * (dy / d) + Fy[b + 2];
            x[b + 2] += vx[b + 2] * del_t;
            y[b + 2] += vy[b + 2] * del_t;

            if (x[b + 2] < 0) {
                x[b + 2] = x[b + 2] + sizex;
            } else if (x[b + 2] > sizex) {
                x[b + 2] = x[b + 2] - sizex;
            }
            if (y[b + 2] < 0) {
                y[b + 2] = y[b + 2] + sizey;
            } else if (y[b + 2] > sizey) {
                y[b + 2] = y[b + 2] - sizey;
            }


        }

    }

        sprintf(fname, "yeet.txt");
        fp = fopen(fname, "w");
        for (a=0;a<N;a++){

            fprintf(fp, "%010.5f %010.5f \n", x[a], y[a]);
        }


}






void force(void) {
    double dx, dy, d;
    double sep;
    int a, b, c;


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
            sep = radii[a] + radii[b];

            dx = x[b] - x[a];
            if (dx > sizex / 2) dx = dx - sizex;
            else if (dx < -sizex / 2) dx = dx + sizex;

            dy = y[b] - y[a];

            if (dy > sizey / 2) dy = dy - sizey;
            else if (dy < -sizey / 2) dy = dy + sizey;

            d = sqrt((dx * dx) + (dy * dy));

            if (d < sep) {
                Fx[b] += k * (sep - d) * (dx / d);
                Fx[a] += -k * (sep - d) * (dx / d);
                Fy[b] += k * (sep - d) * (dy / d);
                Fy[a] += -k * (sep - d) * (dy / d);

            }

            a = list[a];
        }


        box = actbox[b];
        for (bsrch = 0; bsrch < 4; bsrch++) {
            a = head[nn[bsrch][box]];
            while (a < N) {

                sep = radii[a] + radii[b];

                dx = x[b] - x[a];
                if (dx > sizex / 2) dx = dx - sizex;
                else if (dx < -sizex / 2) dx = dx + sizex;

                dy = y[b] - y[a];

                if (dy > sizey / 2) dy = dy - sizey;
                else if (dy < -sizey / 2) dy = dy + sizey;

                d = sqrt((dx * dx) + (dy * dy));

                if (d < sep) {
                    Fx[b] += k * (sep - d) * (dx / d);
                    Fx[a] += -k * (sep - d) * (dx / d);
                    Fy[b] += k * (sep - d) * (dy / d);
                    Fy[a] += -k * (sep - d) * (dy / d);

                }

                a = list[a];
            }
        }

    }


    for (a = 0; a < N; a += 4) {
        b = a + 1;

        dx = x[b] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > 2 * r) {
            Fx[b] += -k * (2 * r - d) * (dx / d);
            Fx[a] += k * (2 * r - d) * (dx / d);
            Fy[b] += -k * (2 * r - d) * (dy / d);
            Fy[a] += k * (2 * r - d) * (dy / d);

        }

        dx = x[b] - x[a+1];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b] - y[a+1];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx[b] += -k * (r+r2- d) * (dx / d);
            Fx[a+1] += k * (r+r2 - d) * (dx / d);
            Fy[b] += -k * (r+r2 - d) * (dy / d);
            Fy[a+1] += k * (r+r2 - d) * (dy / d);

        }

        dx = x[b+1] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b+1] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx[b+1] += -k * (r+r2 - d) * (dx / d);
            Fx[a] += k * (r+r2- d) * (dx / d);
            Fy[b+1] += -k * (r+r2 - d) * (dy / d);
            Fy[a] += k * (r+r2 - d) * (dy / d);

        }

        dx = x[a+1] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[a+1] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx[a+1] += -k * (r+r2 - d) * (dx / d);
            Fx[a] += k * (r+r2 - d) * (dx / d);
            Fy[a+1] += -k * (r+r2 - d) * (dy / d);
            Fy[a] += k * (r+r2 - d) * (dy / d);

        }

        dx = x[b+1] - x[b];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b+1] - y[b];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r + r2) {
            Fx[b+1] += -k * (r+r2 - d) * (dx / d);
            Fx[b] += k * (r+r2 - d) * (dx / d);
            Fy[b+1] += -k * (r+r2 - d) * (dy / d);
            Fy[b] += k * (r+r2 - d) * (dy / d);

        }

    }
}