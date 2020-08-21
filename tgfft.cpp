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



    cmplx operator+ (const cmplx &r) const;
    cmplx operator+= (const cmplx &r);
    cmplx operator* (const cmplx &r) const;
    cmplx operator*= (const cmplx &r);

    void copyfrom(const cmplx &tocopy);
    cmplx operator= (const cmplx toassign);
    cmplx(const cmplx &tocopy);

    double getlength();
    double getangle();

};


cmplx::cmplx(double real, double imaginary)
{
    re = real;
    im = imaginary;
}

void cmplx::copyfrom(const cmplx &tocopy)
{
    re = tocopy.re;
    im = tocopy.im;
}

cmplx cmplx::operator= (const cmplx toassign)
{
    copyfrom(toassign);
    return *this;
}

cmplx::cmplx(const cmplx &tocopy)
{
    copyfrom(tocopy);
}


cmplx cmplx::operator+ (const cmplx &r) const
{
    return cmplx(re+r.re,im+r.im);
}

cmplx cmplx::operator+= (const cmplx &r)
{
    re+=r.re;
    im+=r.im;
    return *this;
}

cmplx cmplx::operator* (const cmplx &r) const
{
    return cmplx( (re*r.re)-(im*r.im), (re*r.im)+(im*r.re) );
}

cmplx cmplx::operator*= (const cmplx &r)
{
    copyfrom(*this*r);
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




void fouriertransform(cmplx* input, cmplx* output, int length)
{
    double norm = (1.0/length);
    for(int k = 0;k<length;k++)
    {
        cmplx kval;
        for(int n = 0;n<length;n++)
        {
            double twopiknoverlength = (TWO_PI * k * n)/length;
            kval += input[n] * cmplx(cos(twopiknoverlength), -(sin(twopiknoverlength)));
            //kval+=cmplx(input[n]*cos(twopiknoverlength), -(input[n]*sin(twopiknoverlength)));
        }
        output[k] = kval*norm;
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



cmplx* getcomplex(double* theinput, int length)
{
    cmplx* output = new cmplx[length];
    for(int i = 0;i<length;i++)
    {
        output[i].re = theinput[i];
    }
    return output;
}


void printall(double* theinput, int length)
{
    printf("input: \r\n");
    printdoubles(theinput,length);

    cmplx* complexinput = getcomplex(theinput,length);

    cmplx output[length];
    fouriertransform(complexinput,output,length);
    double amplitudes[length];
    double phases[length];
    for(int i = 0;i<length;i++)
    {
        amplitudes[i] = output[i].getlength();
        if(amplitudes[i] > 0.0000001)
        {
            phases[i] = output[i].getangle()*rad2deg;
        }
        else
        {
            phases[i] = 0;
        }        
    }
    printf("amplitudes: \r\n");
    printdoubles(amplitudes,length);
    printf("phases: \r\n");
    printdoubles(phases,length);

    cmplx backthroughoutput[length];
    //mer??
    fouriertransform(output,backthroughoutput,length);
    //mer?
    for(int i = 0;i<length;i++)
    {
        printf("%f, ",backthroughoutput[i].re);//, backthroughoutput[i].im);
    }
    printf("\r\n");

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