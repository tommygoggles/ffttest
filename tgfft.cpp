#include <cstdio>
#include <math.h>

#include <fstream>//ifstream
#include <cstring>//strncmp

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






void addcosfloat(double* signal, int length, double frequency, double level, double phase)
{
    for(int i = 0;i<length;i++)
    {
        signal[i] += cos(phase+((float)(i)/length)*TWO_PI*frequency )*level;
    }
}

void addcos(double* signal, int length, int frequency, double level, double phase)
{
    for(int i = 0;i<length;i++)
    {
        signal[i] += cos(phase+((float)(i)/length)*TWO_PI*(double)frequency )*level;
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


void combine(double* signal, double* amplitudes, double* phases, int length)
{
    for(int i = 0;i<length;i++)
    {
        addcos(signal,length,i,amplitudes[i],phases[i]);
    }
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
        if(amplitudes[i] > 0.00000000001)
        {
            phases[i] = output[i].getangle();//*rad2deg;
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

    /*
    cmplx backthroughoutput[length];
    //mer??
    for(int i = 0;i<length;i++)
    {
        output[i].im = -output[i].im;
    }
    fouriertransform(output,backthroughoutput,length);
    //mer?
    for(int i = 0;i<length;i++)
    {
        //im flip, too.
        backthroughoutput[i].re*=length;
    }
    for(int i = 0;i<length;i++)
    {
        printf("%f, ",backthroughoutput[i].re);//, backthroughoutput[i].im);
    }
    printf("\r\n");*/

    double combined[length];
    for(int i = 0;i<length;i++)
    {
        combined[i] = 0;
    }
    combine(combined,amplitudes,phases,length);
    printf("recombined: \r\n");
    printdoubles(combined,length);

    printf("error?: \r\n");
    for(int i = 0;i<length;i++)
    {
        printf("%f, ",theinput[i]-combined[i]);
    }
}


double getfrequency(int length, double sourcefrequency, int cyclesperwhole)
{
    double seconds = (double)length/sourcefrequency;
    double frequency = (double)cyclesperwhole/seconds;
    int other = length-cyclesperwhole;
    double otherfrequency = (double)other/seconds;
    printf("seconds: %f, frequency: %f, other: %f\r\n",seconds,frequency,otherfrequency);
    return frequency;
}



class wavdata
{
    public:
    int filesize;
    char* data;

    short channels;
    unsigned int samplerate;
    short bitspersample;
    short bytespersample;

    unsigned int startofaudio;
    unsigned int bytesofaudio;
    unsigned int bufferlength;

    wavdata();
    ~wavdata();

    int savedata(const char* filename);
    int loaddata(const char* filename);
    int load(const char* filename);
};

wavdata::wavdata() : filesize(0), data(0)
{

}

wavdata::~wavdata()
{
    if(data)
    {
        delete(data);
    }
}


//16-bit samples are stored as
//2's-complement signed integers, ranging from -32768 to 32767


int wavdata::load(const char* filename)
{
    if(loaddata(filename) != 0){return 1;}
    if(filesize <= 44){return 2;}
    if(strncmp("RIFF",(const char*)&(data[0]),4)){return 3;}
    if( *((int*)&(data[4])) != (filesize-8)  ){return 4;}
    if(strncmp("WAVE",(const char*)&(data[8]),4)){return 5;}
    if(strncmp("fmt ",(const char*)&(data[12]),4)){return 6;}
    if( *((int*)&(data[16])) != 16 ){return 7;}//16 bytes follow for PCM. others not supported yet.
    if( *((short*)&(data[20])) != 1 ){return 8;}

    channels = *((short*)&(data[22]));
    samplerate = *((unsigned int*)&(data[24]));
    bitspersample = *((short*)&(data[34]));
    bytespersample = bitspersample/8;
    /*28        4   ByteRate         == SampleRate * NumChannels * BitsPerSample/8
        32        2   BlockAlign       == NumChannels * BitsPerSample/8
                        The number of bytes for one sample including
                        all channels. I wonder what happens when
                        this number isn't an integer?*/

    //floataudio = false;

    printf("%s",(const char*)&(data[36]));
    if(strncmp("data",(const char*)&(data[36]),4)){return 9;}
    bytesofaudio = *((unsigned int*)&(data[40]));
    bufferlength = ((bytesofaudio*8)/channels)/bitspersample;//samplesofaudio

    //assert( (bytesofaudio+44) == filelength);// actually 4 x 44??
    startofaudio = 44;

    return 0;
}

int wavdata::loaddata(const char* filename)
{
    std::ifstream stream(filename, std::ios_base::in | std::ios_base::binary);
    if(!stream.is_open())
    {
        return 1;
    }
    stream.seekg(0, std::ios::end);
    filesize = stream.tellg();
    stream.seekg(0, std::ios::beg);
    data = new char[filesize];
    stream.read(data,filesize);
    return 0;
}

int wavdata::savedata(const char* filename)
{
    std::fstream outputfile;
    outputfile.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    if(!outputfile.is_open())
    {
        return 1;
    }
    outputfile.seekg(0, std::ios::beg);
    outputfile.write(data, filesize);
    outputfile.close();
    return 0;
}


int main (int argc, const char* argv[])
{
    wavdata mong;
    int l = mong.load("iine.wav");
    printf("%i\r\n",l);

    double testvals[100];// = {1,0,0,0,0,0,0,0,0,0,0,1};
    for(int i = 0;i<100;i++)
    {
        testvals[i] = 0.0;// ((double)rand()/RAND_MAX*2)-1.0;
    }

    /*for(int i = 0;i<20;i++)
    {
        testvals[i] = ((double)rand()/RAND_MAX*2)-1.0;
    }*/

    addcos(testvals,100,1,1.0,0);
    addcos(testvals,100,27,1.0,0);

    printall(testvals,100);

    /*printf("\r\n48000kHz, 200 samples\r\n");
    for(int i = 1;i<=100;i++)
    {
        getfrequency(200,48000,i);
    }*/
    //these are summed together, at the points of the sampling rate. what does that look like?
    //getfrequency(2000,48000,1000);
    return 0;
}