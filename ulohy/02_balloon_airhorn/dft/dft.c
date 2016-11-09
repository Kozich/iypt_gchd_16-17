#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

typedef struct c_{
    double r;
    double i;
}   c_double;

typedef struct p_{
    double x;
    struct p_ *n;
}   point;

typedef struct pf_{
    unsigned long s;
    point *fp;
}   pf;

c_double c_exp(c_double c){
    c_double out;
    out.i=exp(c.r)*sin(c.i);
    out.r=exp(c.r)*cos(c.i);
    return out;
}

c_double c_add(c_double a, c_double b){
    a.i+=b.i;
    a.r+=b.r;
    return a;
}

c_double c_t_c(c_double a, c_double b){
    a.i=a.r*b.i+a.i*b.r;
    a.r=a.r*b.r-a.i*b.i;
    return a;
}

c_double c_t_r(c_double a, double b){
    a.i*=b;
    a.r*=b;
    return a;
}

int c_gt_c(c_double a, c_double b){
    if(a.r*a.r+a.i*a.i<b.r*b.r+b.i*b.i){
        return 0;
    }
    else{
        return 1;
    }
}

pf *input(FILE *s){
    char *line=malloc(40);
    if(!fgets(line, 40, s)){
        printf("Empty file");
        return (pf *) NULL;
    }
    pf *f=(pf *) malloc(sizeof(pf));
    f->s=0;
    f->fp=(point *) NULL;
    printf("a\n");
    fflush(stdout);
    do{
        point *np=(point *) malloc(sizeof(point));
        sscanf(line,"%lf",&np->x);
        np->n=f->fp;
        f->fp=np;
        f->s++;
    } while(fgets(line, 40, s));
    printf("Size of field is %lu\n", f->s);
    fflush(stdout);
    return f;
}

c_double fft_coe(pf *f, double o){
    c_double sum;
    sum.i=0;
    sum.r=0;
    c_double e;
    double ec=-2*PI*o/(f->s);
    e.r=0;
    e.i=ec;
    point *p=f->fp;
    for(unsigned long i=0; i<f->s; i++){
        sum=c_add(sum,c_t_r(c_exp(c_t_r(e, i)),p->x));
        p=p->n;
    }
    return sum;
}

int main(int argc, char *argv[]){
    fprintf(stdout,"Initializing Data\n");
    FILE *t;
    if(t=fopen(argv[1],"r+")){
        printf("%s\n",argv[1]);
    }
    pf *f=input(t);
    c_double out;
    out.r=0;
    out.i=0;
    c_double a;
    a.r=1;
    a.i=1;
    FILE *of=fopen("of","w+");
    fprintf(stdout,"Starting DFT\n");
    time_t ot=time(NULL);
    for(unsigned long u=1;u<f->s/2;u++){
        c_double k = fft_coe(f,u);
        fprintf(of,"%lu %lf\n", u*8000/(f->s), k.r*k.r + k.i*k.i);
    }
    time_t rt=time(NULL);
    fprintf(stdout,"Finished\n");
    fprintf(stdout,"Runtime was %i s\n",(int) difftime(rt,ot));
    return 0;
}

