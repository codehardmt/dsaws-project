/**
    @file
    simplemsp~: a simple audio object for Max
    original by: jeremy bernstein, jeremy@bootsquad.com
    @ingroup examples
*/

#include "ext.h"            // standard Max include, always required (except in Jitter)
#include "ext_obex.h"        // required for "new" style objects
#include "z_dsp.h"            // required for MSP objects
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MIDC_OFFSET (261.62556530059868 / 256.0)



// struct to represent the object's state
typedef struct _dsaws {
    t_pxobject        ob;            // the object itself (t_pxobject in MSP instead of t_object)
    float si[1024]; //storing every si in each voice
    float phase[1024]; // storing every phase in each voice
    short w_connected[1];
} t_dsaws;


// method prototypes
double calculateWL(double sr, double freq);
double octcps(double cps);
double cpsoct(double oct);
void *dsaws_new(t_symbol *s, long argc, t_atom *argv);
void dsaws_free(t_dsaws *x);
void dsaws_assist(t_dsaws *x, void *b, long m, long a, char *s);
void process_saw(t_dsaws* sptr, int index);
void dsawsz_float(t_dsaws* x, double freq);
void dsaws_dsp64(t_dsaws *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void dsaws_perform64(t_dsaws *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);


// global class pointer variable
static t_class *dsaws_class = NULL;


//***********************************************************************************************

// calculate how many samples in one wavelength given the frequency and sampling rate
// for example: 440 hz in 44100 sampling rate (in one sec), each wavelength is 44100/440, each wavelength is 100.227 samples long, and afterward we will calculate the sample increment within -1 to 1(saw wave).
double calculateWL(double sr, double freq){
    float wavelength = sr / freq; // Calculate the number of samples in one wavelength.
    return wavelength;
}

double octcps(double cps)
{
    return log(cps / MIDC_OFFSET) / M_LN2;
}

double cpsoct(double oct)
{
    return pow(2.0, oct) * MIDC_OFFSET;
}

void process_saw(t_dsaws* x, int index){ // arguement is a pointer to struct, we want to change the "real" value stored in the memory address, int index is how many voices
    
    x->phase[index] = x->phase[index] + x->si[index]; //logically *s.phase (syntax thing) pointer to the real value, calculating each voice's different phase + different si
    if (x->phase[index] >= 1.0){
        x->phase[index] = -1.0;
    }
}

void dsaws_detune(t_dsaws* x, double base, double detune){ // calculate sampling increment for add-up saw waves
    
    
    base = octcps(base); // taking base frequency converted to e.g. 8.00
    double base_detune = detune; // set the base_detune for increment later
    
    //calculating each freq of spacing-out detune when the saw becomes more, and further calculating the sampling increment
    for(int i = 0; i < 1024; i++){
        
        if(i == 0){ // first initial value
            x->si[i] = 2.0 / calculateWL(sys_getsr(), cpsoct(base));
        }
        
        else if((i + 1) % 2 == 0)// if it's odd num
        {
            x->si[i] = 2.0 / calculateWL(sys_getsr(), cpsoct(base - detune));
        }
        
        else // even num
        {
            x->si[i] = 2.0 / calculateWL(sys_getsr(), cpsoct(base + detune));
            
        }
        detune = detune + base_detune; // instead of adding the base num, add the base_detune num to make it equally wider.
    }
    
}


void ext_main(void *r)
{
    // object initialization, note the use of dsp_free for the freemethod, which is required
    // unless you need to free allocated memory, in which case you should call dsp_free from
    // your custom free function.

    t_class *c = class_new("dsaws~", (method)dsaws_new, (method)dsp_free, (long)sizeof(t_dsaws), 0L, A_GIMME, 0);

    class_addmethod(c, (method)dsaws_dsp64,     "dsp64",    A_CANT, 0);
    class_addmethod(c, (method)dsaws_assist,    "assist",    A_CANT, 0);
    class_addmethod(c, (method)dsawsz_float,    "float",    A_FLOAT, 0);
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    dsaws_class = c;
}


void *dsaws_new(t_symbol *s, long argc, t_atom *argv) // creating object
{
    t_dsaws *x = (t_dsaws *)object_alloc(dsaws_class);
    
    if (x) {
        dsp_setup((t_pxobject *)x, 1);    // MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        
        
        
        outlet_new(x, "signal");
        outlet_new(x, "signal"); // signal outlet (note "signal" rather than NULL)
        
        
        x->w_connected[0] = 0;
        
        int i;
        for(i = 0; i < 1024; i++){
            x->phase[i] = -1;
        } // set the phase to -1 everytime
        float frequency = atom_getfloat(argv);
        if(argc >= 1){ // argument take in object as freq
            
            if(frequency > 0)
            {
                x->si[i] = 2.0 / calculateWL(sys_getsr(), frequency); // -1 to 1 divide samplerate and frequecy = sample increment
                dsaws_detune(x, frequency, 0.001);
            }
            else
            {
                error("please enter frequency > 0.\n");
                return;
            }
        }
        else
        {
            float defaultFreq = 440.0;
            x->si[i] = 2.0 / calculateWL(sys_getsr(), defaultFreq);
            dsaws_detune(x, frequency, 0.001);
            
        }
    }
    return (x);
}



// NOT CALLED!, we use dsp_free for a generic free function
void dsaws_free(t_dsaws *x) // called when object is deleted
{
    ;
}


void dsaws_assist(t_dsaws *x, void *b, long m, long a, char *s) // when we move mouse, it pops up information
{
    if (m == ASSIST_INLET) { //inlet
        sprintf(s, "I am inlet %ld", a);
    }
    else {    // outlet
        sprintf(s, "I am outlet %ld", a);
    }
}






// float function handle changing the floats that going in freq inlet into frequency
void dsawsz_float(t_dsaws* x, double freq){
    if(freq > 0){
        x->si = 2.0 / calculateWL(sys_getsr(), freq);
    }
    else{
        error("please enter frequency > 0.\n");
        
    }
    
}

/*int dsawsz_integer(t_dsaws* x, int num){
    if(num > 0 | num < 1024){
        return num;
    }
    else{
        error("please enter a valid number > 0.\n");
    }
    
    
}
*/









// registers a function for the signal chain in Max
void dsaws_dsp64(t_dsaws *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{

    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling function
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method

    x->w_connected[0] = count[0];
    
    //0 disconnected, 1 connected (inlet) decided by the users
   

  


    object_method(dsp64, gensym("dsp_add64"), x, dsaws_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void dsaws_perform64(t_dsaws* x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) //double **ins is a pointer to an array of pointers, where each pointer represents a signal inlet.
{
    t_double *inL = ins[0];        // we get audio for each inlet of the object from the **ins argument
    
    t_double *outL = outs[0];    // we get audio for each outlet of the object from the **outs argument
    t_double *outR = outs[1];
    
    
    
    
    long n = sampleframes;
    
    
    
    
    for (int time=0; time < sampleframes; time++){ // how many samples user/maxmsp are going to take at a time
        
        if (x->w_connected[0]){
            t_double inputL = *inL; // freq
            if(inputL > 0)
            {
                dsawsz_float(x, inputL);
            }
            else{
                error("please enter frequency > 0.\n");
            }
        }
        
        inL++;
        float allsaws = 0;
        float sum = 0; // initiallized the sum in phase[index]'s data type
        for (int index=0; index < 4; index++) //之後要改回1024 index: how many voices
        {
            if (x->w_connected[0]){
                t_double inputL = *inL; // freq
                process_saw(x, index);
                allsaws = x->phase[index] + dsaws_detune(x, inputL, 0.001);
                    
                }
                sum = sum + allsaws; // adding all the saw waves
                
            }
            *outL++= sum / 4; // amplitude divided by the array phase[1024] /改回1024
            *outR++= sum / 4;
        }
    }
}


   
