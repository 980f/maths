/*  ALHALHALHALHALHALHALHALHALHALHALHALHALHALHALHALHALHALHALHALH
    Integral preserving Quadratic compression:
        useful for xrf data.

    QUADCOMP input_file  output_file

*/
#pragma hdrfile "quadcomp.sym"
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <stdio.h>
#include <ctype.h>
#include <io.h>
#include <fcntl.h>
#pragma hdrstop

const char *ownedby="(C) 1992 by Andy Heilveil";
int noisy=0;
int outform=2; //0 binary,1 raws,2 xrayspex (8194)
int opmode=1;  //0 linear, 1 quadratic
// declare type that may change from long to int or vice versa:
typedef unsigned int datum;     // data in and out
typedef int coord;    // calculation coordinate system


/*make a parameter block that can be shared between compression
    and expansion so that one day we can dynamically slide from one
    to the other.
*/

typedef struct {
	coord perout;       /*y gain per step of X*/
	coord perinp;       /*x gain per step of Y, rate at which it changes*/
	coord twiddle;		/*for dithering calculations*/
	datum *din,*dout;	/*point at arrays for in and out,
						   caller worries about overlap effects */
	int incnt,outcnt,outnominal; /*number of points in,out,and quess at out*/
} qcpak;

//#include <lowmath.c>
//extracted just the one needed:
/* arcane math routines:
    mpdiv multiplies a by b then divides by c. it differs from anything
        reasonable in C by maintaining temporary extended precision.
*/
unsigned mpdiv(unsigned a,unsigned b,unsigned c){
    asm{
        mov ax,a;       //16 bit in ax
        mul word ptr b; //*16= 32 bit in dx:ax
        mov cx,c
        cmp dx,cx     ; // pre check for overflows, /0
        jae errout    ;//overflow if dx>=cx.,if dx<cx then cx !=0
        div cx; // div 16=>16 bit in ax
    }
    return _AX;         // which is where return values go
errout:
    return 0xFFFF;
}


char yeses[]="Yy\r ";
void breakpoint (int why){ //debug tool
int k;
     printf("\r\n  BREAKPOINT %d reached, Continue?  ",why);
     if(!(k=getch())) k=256+getch();
     if(strchr(yeses,k)) return; else exit(-1);
}

/* for expansion routines
    macro BITE is what we do each time we look at a new input:
        fetch next cell  ;  count it ; precalc full size output
*/
#define BITE unit=*small++; step--; chunk=mpdiv(unit,icnt,ocnt)

int lexpand(datum *small,datum *big,int icnt,int ocnt){//linear expander
coord twiddle;
datum unit,chunk,chomp;
int step;
int ochk=0; //will count outputs
    step=icnt;
    twiddle=-ocnt;
    BITE;
    do {
        twiddle+=icnt;/* gen next test value, placed here means we have
            optimized to output every cycle */

        if( twiddle<0){ //a piece of the input
            *big++=chunk; //output a chunk
            unit-=chunk;  //this much left
        }else{
            if(twiddle>0){ //split across boundary
                *big=unit ;//dregs of last cell
                //read new cell
                BITE;
                *big++=(chomp=mpdiv(unit,twiddle,ocnt));//this much of new
                unit-=chomp;
            }else{ //perfectly divided
                *big++=unit; //last chunk, uses 'unit' to handle round off
                BITE;
            }
            twiddle-=ocnt;
        }
        ochk++; //actual count of outputs
    } while(step);
    while(ochk<ocnt){
        *big++=chunk;
        ochk++;
    }
    return ochk;
}

/*  for expansion q.perinp is a large and growing larger amount
    while q.perout is fixed. q.perout is the denominator of the
    expansion constant, q.a is the numerator. The local expansion
    factor at any input cell boundary is q.perinp/q.perout.
*/

qcpak qexpand(qcpak q){//quadratic expander
datum unit,chunk,chomp;
int step;
int ochk=0; //will count outputs
    step=q.perout;
    q.twiddle-=q.perinp; //(negative) fraction remaining part of input cell
    BITE;
    do {
        q.twiddle+=q.perout;/* gen next test value, placed here means we have
            optimized to output every cycle */
        if( q.twiddle<0){ //a piece of the input
            *big++=chunk; //output a chunk
            unit-=chunk;  //this much left
        }else{
            if(q.twiddle>0){ //split across boundary
                *big=unit ;//dregs of last cell
                //read new cell
                BITE;
                *big++=(chomp=mpdiv(unit,q.twiddle,q.perinp));//this much of new
                unit-=chomp;
            }else{ //perfectly divided
                *big++=unit; //last chunk, uses 'unit' to handle round off
                BITE;
            }
            q.perinp += q.a ;//quadraticafier
            q.twiddle-=q.perinp; //takes more fixed size
        }
        ochk++; //actual count of outputs
    } while(step);
    while(ochk<q.perinp){
        *big++=chunk;
        ochk++;
    }
    return ochk;
}


qcpack lsquish(qcpack q){//optimized linear squeeze
    return q;
}

/*  for expansion q.perout is a large and growing larger amount
    while q.perinp is fixed. q.perinp is the denominator of the
    expansion constant, q.a is the numerator. The local compression
    factor at any output cell boundary is q.perinp/q.perout.
*/


qcpak qsquish(qcpak q){//squeezer
int step;
datum xword,nxword,yword;
	xword=0; q.outcnt=0; step=q.incnt;
	if(noisy){
        printf("\r\n Starting with dx=% 5d, dy=% 5d, steps=%d",        
                                q.perinp,       q.perout,   step);
        printf("\r\n step  twiddle  dx  xdatum   outcnt   ");
    }
//bias loop to pick up first data point:
    q.twiddle-=q.perinp;//(negative) fractional part of output remaining to be built
	while(step){
        if (q.twiddle<0){ //do output
            xword+=(yword=*q.din++); //this is what preserves integrals
            q.twiddle+= q.perinp;
            step--; //consume all inputs
        }else{ //next input
            if(dirf){ //split
                nxword= mpdiv(yword,q.twiddle,q.perinp); //belongs to next output cell
                *q.dout++=xword-nxword;
                if(noisy) printf(" outputting %5u",q.dout[-1]);
                q.outcnt++;
                xword=nxword; //fractional start on next output

            } //else stretching..
            q.perout+=q.a; //without this one statement we are linear.
            q.twiddle-= q.perout;
		}
		if(noisy){
			printf("\r\n %5d  % 5d  % 5d  %5u  %5d  %d",
                    step,q.twiddle,q.perout,xword,q.outcnt,dirf);
			if(q.perout<=0) breakpoint(2);
		}

	}

    *q.dout++=xword-nxword; q.outcnt++;
    printf("\r\n Final output %5u (%u/%u ths of normal)",
           q.dout[-1] ,q.perout+q.twiddle,q.perinp);
	return q;
}
/*
	Set up for compression run.
*/
int quadcomp(FILE *bigfile,FILE *smallfile,double nin,double nout,double astep){
qcpak q;
    q.incnt=nin; /*we will take care to not run off the end of
                the actual data irregardless of other considerations*/
    q.outnominal=nout;


if (0.==astep){ //linear
	q.perout=nin;         //yes, the output gets the input count and vice versa
	q.perinp=nout;
	q.a= 0;
}else{
		//need to find best integer pair instead of truncating num and denom
		q.perinp= q.outnominal-1;   //to make a be an integer for integral cfactor
		q.a= (2.*nin/nout)-2;
		q.perout=q.perinp; //start 1:1 /*** arbitrary (but typical) decision !!***/
	}
}

	printf("\r\n Magic number is %d/%d for %d input points compressed %lg:%lg ",
							  q.a, q.perinp,   q.incnt,              nin,nout);

	if(!(q.din=(datum *)calloc(q.incnt+10,sizeof(datum)))) return 2;
	if(!(q.dout=(datum *)calloc(q.outnominal+10,sizeof(datum)))) return 3; //extra big for safety

	q.incnt=fread(q.din,sizeof(datum),q.incnt,bigfile);

	datum *din=q.din;datum *dout=q.dout; //preserve pointers
    q.twiddle=0; //initialize fractional part of coord so that we may
                 //slickly jump from expansion to compression
                 //without a blip.
    q=qsquish(q); //low level go ahead and do it.

	//accept returned size
	switch(outform){
	case 0:
		if (q.outcnt!=fwrite(dout,sizeof(datum),q.outcnt,smallfile)) return 4;
		break;
	case 1:
		if (q.outcnt!=fwrite(dout,sizeof(datum),q.outcnt,smallfile)) return 4;
		break;
	case 2: //re-expand
        if (4096!=lexpand(dout,din,q.outcnt,4096)) return 5;
		if (4096!=fwrite(din,2,4096,smallfile)) return 4;
		break;
	}
    return 0; //success
}

char *ermsg[]={"Success",           //0
    "File problems with input",     //1
    "Couldn't make room for input", //2
    "Couldn't find room for output buffer", //3
    "File problems with output",     //4
    "Didn't expand right for XRAYSPEX", //5
    ""
};
int ermsglast=sizeof(ermsg)/sizeof(int);
void erexit(int ecode){
	cputs("\r\n Exit Status:");
    printf( ecode<ermsglast?ermsg[ecode]:"Unusual Error #%d  ",ecode);
    exit(ecode);  //for batch file 'errlevel'
}

void main(int argc,char *argv[]){
FILE *fin,*fout;
char *p,c;

double numinp,numout; //number of points in
double cfactor=1.; //inp/out (linked to above), total compression desired.
double astep=-1.; //flag that this needs to be calculated
char *iname= "";
char *oname= "";
char aname[66];
    clrscr();
    printf("\r You are running version %s \r\n    of QUADCOMP quadratic data scaler by ALH.",argv[0]);
    printf("\r\n   Run with ? as operand to get full command line syntax.");
    while(--argc){  //scan all arguments
        switch (*(p=argv[argc])++){ //branch on first character of each argument
            case '?': // help
                cputs("\r\n quadcomp options");
                cputs("\r\n /factor  shrink by (real) factor  *factor  expand by factor,if missing: -4");
                cputs("\r\n  filename (no special character) defines input file, if missing: QUADCOMP.IN");
                cputs("\r\n =filename defines output file, if missing: input.OUT");
                cputs("\r\n -mode 0=linear else (real) quadratic rate, if missing:quadratic calculated");
                cputs("\r\n @format output file 0=simple binary, 1=raws, 2=xrayspex");
                cputs("\r\n !noise level  control screen printouts, smaller is quieter.");
                break;
            case '@':
                outform=atoi(p);
                break;
            case '=': //
                oname=p;
                break;
            case '*': //scale up factor
                cfactor*=atof(p);
                break;
            case '/': //scale down factor
                cfactor/=atof(p);
                break;
            case '-': //fudged in curve
                astep=atof(p);
                break;
            case '!': // noise
                noisy=atoi(p);
                break;
            default:
                if(isalnum(*--p)){
                    iname=p;
                }else{
                    printf("<%s> is not a valid command, perhaps it needs some spaces? " ,p);
                    argv[argc++]="?"; //slip in ? command
                }
            }
    }//end all args scan


    if(cfactor==1.) cfactor=.25;
    if(!*iname){
        iname="quadcomp.in";
    }
    printf("\r\n Using %s for input data file.",iname);
    fin=  fopen(iname,"rb");
    if (fin==NULL) erexit(1);

    if(!*oname){
        oname= aname;
        p= strrchr(iname,'.');
        if(p) {int cut=p-iname;
            strncpy(oname,iname,cut);
            oname[cut]='\0';
        }else{
            strcpy(oname,iname);
        }
        strcat(oname,".out");
    }
    printf("\r\n Using %s for output data file.",oname);
    fout= fopen(oname,"wb");
    if (fout==NULL) erexit(4);
    long along=  filelength(fileno(fin));
    if( 0L>along) erexit(1);
    while('\n'!=(c=fgetc(fin))){
        --along; //strip ascii header
        if (outform==1){   //copy header
            fputc(c,fout);
        }
	}
	--along;
    if(outform==1) fputc('\n',fout);
    if(outform==2) { fputc('\r',fout);fputc('\n',fout);}
//    if(EOF==fflush(fout)) erexit(4); //was getting wierd behaviour
    numinp= along/sizeof(datum); //number of input points
    numout= numinp*cfactor;
	printf("\r\n Instructions are to create %lg outputs from %lg inputs..",
                                             numout,          numinp);
	erexit(quadcomp(fin,fout,numinp,numout,astep));

}

/*
    Derivation:

    First one needs to understand Integral preserving compression
using a ratio of integers as the scale factor. This is of course
a little bit more complicated then a compression that is n:1. For
compression we are adding together multiple inputs to generate
each output. To allow for arbitrary rescaling we must make sure
that when we split an two input cell into parts of two output
cells that the sum of the split parts equals the whole. In this
program this is done by calculating one part to the nearest
integer, then getting the other part by subtract the first from
the whole. The previous is a picky point but if one were to
iterate compression, expansion,compression,expansion any loss of
counts would accumulate.

    TO avoid explicit fractional representations and the loss of
exactness that the finite precision results of computer division
cause we maintain our fractions as ratios of integers.This may be
viewed as deciding on a denominator and scaling all numbers
according to it. Eventually we will have to really do divides to
apportion data in to multipl data outs BUT with proper attention
we can follow our rounding to preserve integrals.


*/
