#include <cstdio>
#include <math.h>

#define PI 3.14159265358979323846
#define TWO_PI (2*PI)
#define rad2deg (180.0 / PI)
#define deg2rad (PI / 180.0)


class cmplx
{
    public:
    double re;
    double im;

    cmplx(double real = 0, double imaginary = 0);


    cmplx operator+ (const cmplx r) const;
    cmplx operator+= (const cmplx r);

    double getlength();
    double getangle();

};

cmplx::cmplx(double real, double imaginary)
{
    re = real;
    im = imaginary;
}


cmplx cmplx::operator+ (const cmplx r) const
{
    return cmplx(re+r.re,im+r.im);
}

cmplx cmplx::operator+= (const cmplx r)
{
    re+=r.re;
    im+=r.im;
    return *this;
}

double cmplx::getlength()
{
    return(sqrt( (re*re) + (im*im) ));
}

double cmplx::getangle()
{
    return atan2(im,re);
}

void fouriertransform(double* input, cmplx* output, int length)
{
    for(int k = 0;k<length;k++)
    {
        cmplx kval;
        for(int n = 0;n<length;n++)
        {
            double twopiknoverlength = (TWO_PI * k * n)/length;
            kval+=cmplx(input[n]*cos(twopiknoverlength), -(input[n]*sin(twopiknoverlength)));
        }
        kval.re/=length;
        kval.im/=length;
        output[k] = kval;
    }
}


void printdoubles(double* thedoubles, int length)
{
    for(int i = 0;i<length;i++)
    {
        printf("%f, ",thedoubles[i]);
    }
    printf("\r\n");
}

void printall(double* theinput, int length)
{
    printf("input: \r\n");
    printdoubles(theinput,length);

    cmplx output[length];
    fouriertransform(theinput,output,length);
    double amplitudes[length];
    double phases[length];
    for(int i = 0;i<length;i++)
    {
        amplitudes[i] = output[i].getlength();
        phases[i] = output[i].getangle()*rad2deg;
    }
    printf("amplitudes: \r\n");
    printdoubles(amplitudes,length);
    printf("phases: \r\n");
    printdoubles(phases,length);

}


int main (int argc, const char* argv[])
{
    double testvals[100];// = {1,1,1,0,0,0,1,1,1,0,0,0};
    for(int i = 0;i<100;i++)
    {
        testvals[i] = sin((float)(i/100.0)*TWO_PI);
    }

    printall(testvals,100);
    return 0;
}