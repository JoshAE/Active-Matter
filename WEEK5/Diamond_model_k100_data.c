/* Active 2D particles with a box search method*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2.h>
#include <g2_X11.h>
#include <unistd.h>
#include <g2_gd.h>
#define N 5000

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


double x2[N],y2[N],Fx2[N],Fy2[N],x2old[N],y2old[N],vx2[N],vy2[N],xoutput[N*100],youtput[N*100];
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
    sizex=sqrt(((N*M_PI*r*r)+(N*M_PI*r2*r2))/phi);
    sizey=sizex;

    ext=sqrt((2*r*r2)-(r2*r2));




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

    nx=(int)(88);
    ny=(int)(80);
    rinx=10/88;
    riny=18/80;
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
        vx2[a]=0;
        vy2[a]=0;

    }

//initialising particle positions

    for (b=0;b<ny;b++)
    {

        for (a=0;a<nx;a+=2)
        {
            if(n < N)
            {


                x[n]=rinx+(rinx*a/2)+(a*r*2);
                y[n]=riny+(riny*b)+(2*r*b);



                x[n+1]=x[n]+2*r;
                y[n+1]=y[n];

                n+=2;
            }
        }

    }
    for (a=0;a<N;a+=2){

        x2[a]=x[a]+r;
        y2[a]=y[a]+ext;
        x2[a+1]=x[a]+r;
        y2[a+1]=y[a]-ext;
    }



    for (i = 0; i < 10000/del_t; i++) {


        if (i > 1000000 * del_t) {
            Vo = 1;
        }

        force();

        for (a = 0; a < N; a += 2) {
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
            }
            else if (x[a] > sizex) {
                x[a] = x[a] - sizex;
            }
            if (y[a] < 0) {
                y[a] = y[a] + sizey;
            }
            else if (y[a] > sizey) {
                y[a] = y[a] - sizey;
            }

            vx[b] = Vo * (dx / d) + Fx[b];
            vy[b] = Vo * (dy / d) + Fy[b];
            x[b] += vx[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);
            y[b] += vy[b] * del_t + sqrt(2 * del_t * D * 3) * (2 * drand48() - 1);

            if (x[b] < 0) {
                x[b] = x[b] + sizex;
            }
            else if (x[b] > sizex) {
                x[b] = x[b] - sizex;
            }
            if (y[b] < 0) {
                y[b] = y[b] + sizey;
            }
            else if (y[b] > sizey) {
                y[b] = y[b] - sizey;
            }





            //Small particles


            vx2[a] = Vo * (dx / d) + Fx2[a];

            vy2[a] = Vo * (dy / d) + Fy2[a];


            x2[a] += vx2[a] * del_t;
            y2[a] += vy2[a] * del_t;

            if (x2[a] < 0) {
                x2[a] = x2[a] + sizex;
            }
            else if (x2[a] > sizex) {
                x2[a] = x2[a] - sizex;
            }
            if (y2[a] < 0) {
                y2[a] = y2[a] + sizey;
            }
            else if (y2[a] > sizey) {
                y2[a] = y2[a] - sizey;
            }

            vx2[b] = Vo * (dx / d) + Fx2[b];
            vy2[b] = Vo * (dy / d) + Fy2[b];
            x2[b] += vx2[b] * del_t;
            y2[b] += vy2[b] * del_t;

            if (x2[b] < 0) {
                x2[b] = x2[b] + sizex;
            }
            else if (x2[b] > sizex) {
                x2[b] = x2[b] - sizex;
            }
            if (y2[b] < 0) {
                y2[b] = y2[b] + sizey;
            }
            else if (y2[b] > sizey) {
                y2[b] = y2[b] - sizey;
            }


        }




    }

    sprintf(fname, "Diamondshape_packdata/Diamondpositions_Pe=50_r2=%f.txt", r2);
    fp = fopen(fname, "w");
    for (a=0;a<N;a++){

        fprintf(fp, "%010.5f %010.5f %010.5f %010.5f \n", x[a], y[a],x2[a],y2[a]);
    }
}






void force(void)
{
    double dx,dy,d;
    double phi;
    double Fdrive;
    int a,b,c;



    for (c = 0; c < N; c++) {
        Fx[c] = 0;
        Fy[c] = 0;
        Fx2[c] = 0;
        Fy2[c] = 0;
    }


    //Force for repulsion and attraction dumbbells

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


            dx = x2[b] - x[a];
            if (dx > sizex / 2) dx = dx - sizex;
            else if (dx < -sizex / 2) dx = dx + sizex;

            dy = y2[b] - y[a];

            if (dy > sizey / 2) dy = dy - sizey;
            else if (dy < -sizey / 2) dy = dy + sizey;

            d = sqrt((dx * dx) + (dy * dy));

            if (d < r+r2) {
                Fx2[b] += k * (r+r2 - d) * (dx / d);
                Fx[a] += -k * (r+r2 - d) * (dx / d);
                Fy2[b] += k * (r+r2 - d) * (dy / d);
                Fy[a] += -k * (r+r2 - d) * (dy / d);

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

                dx = x2[b] - x[a];
                if (dx > sizex / 2) dx = dx - sizex;
                else if (dx < -sizex / 2) dx = dx + sizex;

                dy = y2[b] - y[a];

                if (dy > sizey / 2) dy = dy - sizey;
                else if (dy < -sizey / 2) dy = dy + sizey;

                d = sqrt((dx * dx) + (dy * dy));

                if (d < r+r2) {
                    Fx2[b] += k * (r+r2 - d) * (dx / d);
                    Fx[a] += -k * (r+r2 - d) * (dx / d);
                    Fy2[b] += k * (r+r2 - d) * (dy / d);
                    Fy[a] += -k * (r+r2 - d) * (dy / d);

                }


                a = list[a];
            }
        }

    }

    //attractive big same dimer
    for (a = 0; a < N; a += 2) {
        b = a + 1;

        dx = x[b] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > 2 * (r)) {
            Fx[b] += k * (2 * (r) - d) * (dx / d);
            Fx[a] += -k * (2 * (r) - d) * (dx / d);
            Fy[b] += k * (2 * (r) - d) * (dy / d);
            Fy[a] += -k * (2 * (r) - d) * (dy / d);
        }






        //big small same dimer

        dx = x2[b] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y2[b] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx2[b] += k * (r+r2 - d) * (dx / d);
            Fx[a] += -k * (r+r2 - d) * (dx / d);
            Fy2[b] += k * (r+r2 - d) * (dy / d);
            Fy[a] += -k * (r+r2 - d) * (dy / d);

        }






        dx = x[b] - x2[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y[b] - y2[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx[b] += k * (r+r2 - d) * (dx / d);
            Fx2[a] += -k * (r+r2 - d) * (dx / d);
            Fy[b] += k * (r+r2 - d) * (dy / d);
            Fy2[a] += -k * (r+r2 - d) * (dy / d);

        }




        dx = x2[a] - x[a];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y2[a] - y[a];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx2[a] += k * (r+r2 - d) * (dx / d);
            Fx[a] += -k * (r+r2 - d) * (dx / d);
            Fy2[a] += k * (r+r2 - d) * (dy / d);
            Fy[a] += -k * (r+r2 - d) * (dy / d);

        }



        dx = x2[b] - x[b];
        if (dx > sizex / 2) dx = dx - sizex;
        else if (dx < -sizex / 2) dx = dx + sizex;

        dy = y2[b] - y[b];

        if (dy > sizey / 2) dy = dy - sizey;
        else if (dy < -sizey / 2) dy = dy + sizey;

        d = sqrt((dx * dx) + (dy * dy));

        if (d > r+r2) {
            Fx2[b] += k * (r+r2 - d) * (dx / d);
            Fx[b] += -k * (r+r2 - d) * (dx / d);
            Fy2[b] += k * (r+r2 - d) * (dy / d);
            Fy[b] += -k * (r+r2 - d) * (dy / d);

        }


    }








//box search for small repulsive
    for (a = 0; a < by * bx; a++) {
        head[a] = N;
    }


    for (a = 0; a < N; a++) {
        box = ((int) (x2[a] * bx / (sizex))) + bx * ((int) (y2[a] * by / (sizey)));
//printf("%d %d %lf %lf\n",box,a,x[a],y[a]);
        list[a] = head[box];
        head[box] = a;
        actbox[a] = box;
    }

    for (b = 0; b < N; b++) {
        a = list[b];

        while (a < N) {

            dx = x2[b] - x2[a];
            if (dx > sizex / 2) dx = dx - sizex;
            else if (dx < -sizex / 2) dx = dx + sizex;

            dy = y2[b] - y2[a];

            if (dy > sizey / 2) dy = dy - sizey;
            else if (dy < -sizey / 2) dy = dy + sizey;

            d = sqrt((dx * dx) + (dy * dy));

            if (d < 2 * r2) {
                Fx2[b] += k * (2 * r2 - d) * (dx / d);
                Fx2[a] += -k * (2 * r2 - d) * (dx / d);
                Fy2[b] += k * (2 * r2 - d) * (dy / d);
                Fy2[a] += -k * (2 * r2 - d) * (dy / d);

            }

            a = list[a];
        }


        box = actbox[b];
        for (bsrch = 0; bsrch < 4; bsrch++) {
            a = head[nn[bsrch][box]];
            while (a < N) {
                dx = x2[b] - x2[a];
                if (dx > sizex / 2) dx = dx - sizex;
                else if (dx < -sizex / 2) dx = dx + sizex;

                dy = y2[b] - y2[a];

                if (dy > sizey / 2) dy = dy - sizey;
                else if (dy < -sizey / 2) dy = dy + sizey;

                d = sqrt((dx * dx) + (dy * dy));

                if (d < 2 * r2) {
                    Fx2[b] += k * (2 * r2 - d) * (dx / d);
                    Fx2[a] += -k * (2 * r2 - d) * (dx / d);
                    Fy2[b] += k * (2 * r2 - d) * (dy / d);
                    Fy2[a] += -k * (2 * r2 - d) * (dy / d);

                }
                a = list[a];
            }
        }

    }




}
