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
    short w_connected[3]; //([]means how many inlets)
    float freqp; //storing users' detune parameter
    float voicenump; //storing users' voices parameter
    float detunep; // storing users' detune parameter
} t_dsaws;


// method prototypes
double calculateWL(double sr, double freq);
double octcps(double cps);
double cpsoct(double oct);
void *dsaws_new(t_symbol *s, long argc, t_atom *argv);
void dsaws_free(t_dsaws *x);
void dsaws_assist(t_dsaws *x, void *b, long m, long a, char *s);
void process_saw(t_dsaws* sptr, int index);
void dsaws_changefreq(t_dsaws* x, double freq);
void dsaws_changedetune(t_dsaws* x, double detune);
void dsawsz_float(t_dsaws* x, double input);
void dsaws_int(t_dsaws* x, long n);
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


void dsaws_detune(t_dsaws* x, double base, double detune){ // calculate sampling increment for add-up saw waves, base freq & add or minus number (0, 1, -1, 2, -2...)
  
    base = octcps(base); // taking base frequency converted to e.g. 8.00
    double base_detune = detune; // set the base_detune for increment later
    
    //calculating each freq of spacing-out detune when the saw becomes more, and further calculating the different sampling increment
    for(int i = 0; i < x->voicenump; i++){ //keep calling this func that it should not calculate 1024 everytime, just calculate what users' entering

        
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

        detune = detune + base_detune; // instead of adding the base num, add the base_detune num to make it equally wider. base_detune will not change, it's the original value that user set, but the detune will add up every time.

    }
    
}


void ext_main(void *r)
{
    // object initialization, note the use of dsp_free for the freemethod, which is required
    // unless you need to free allocated memory, in which case you should call dsp_free from
    // your custom free function.

    t_class *c = class_new("dsaws~", (method)dsaws_new, (method)dsp_free, (long)sizeof(t_dsaws), 0L, A_GIMME, 0);


    class_addmethod(c, (method)dsaws_dsp64,     "dsp64",    A_CANT, 0); //set up signal chains in max
    class_addmethod(c, (method)dsaws_assist,    "assist",    A_CANT, 0);
    class_addmethod(c, (method)dsaws_changefreq,    "float",    A_FLOAT, 0);
    class_addmethod(c, (method)dsaws_changedetune,    "float",    A_FLOAT, 0);
    class_addmethod(c, (method)dsawsz_float,    "float",    A_FLOAT, 0);
    class_addmethod(c, (method)dsaws_int, "int", A_LONG, 0);
    // for the non-signal inlet - connection between users's messages to the written function (binding)
    
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    dsaws_class = c;
}




void *dsaws_new(t_symbol *s, long argc, t_atom *argv) // creating object // parameter -> takes argument // argv is an array of arguments


{
    t_dsaws *x = (t_dsaws *)object_alloc(dsaws_class);
    
    if (x) {

        dsp_setup((t_pxobject *)x, 3);    // MSP inlets: arg is # of inlets and is REQUIRED!
        // use 0 if you don't need inlets
        //intin(x, 2); // for the non-signal inlet, the num is for 3rd inlet (0:1st, 1:2nd, 2: 3rd. index of the inlets

        
        outlet_new(x, "signal");
        outlet_new(x, "signal"); // signal outlet (note "signal" rather than NULL)
        

        //initializtion to not connect!
        x->w_connected[0] = 0;
        x->w_connected[1] = 0;
        x->w_connected[2] = 0;
        
        int i;
        for(i = 0; i < 1024; i++){ // initialized every voices to -1, 1024 is the highest num the func can calculate
            x->phase[i] = -1;
        }
        
        
        //SET UP THE ARGUMENT OF THE OBJECT
            x->freqp = 440.0;
            x->detunep = 0.0;
        if(argc >= 1){// 1st argument in object: freq
            
            x->freqp = atom_getfloat(argv);//All atom_getfloat does is turn the pointer that gets passed to it into a float)
            if(x->freqp > 0)
            {
                // x->si[i] = 2.0 / calculateWL(sys_getsr(), frequency); // -1 to 1 divide samplerate and frequecy = sample increment
                dsaws_detune(x, x->freqp, x->detunep);

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

           // x->si[i] = 2.0 / calculateWL(sys_getsr(), defaultFreq);
            dsaws_detune(x, defaultFreq, x->detunep);
        
        }
    }
        
       
        
        if(argc >= 2){ // 2nd argument in object: num of voices
            
            x->voicenump = atom_getfloat(argv+1); // (argv+1) pointer to the 2nd argument of the array
            
        }
        else{
            x->voicenump = 1; //if user didn't type anything, default will be 1(to make sound, one saw wave)
            }
    
    
        x->detunep = 0.0;
         
        if(argc >= 3){ // 3rd argument in object: detune (in oct.pc)
            x->detunep = atom_getfloat(argv+2);
            
            if(x->detunep > 0){
                dsaws_detune(x, x->freqp, x->detunep);
            }
            else
            {
                error("please enter detune parameter > 0.\n");
                return;
            }
        }
        else
        {
            float defaultdetunep = 0.0;
            dsaws_detune(x, x->freqp, defaultdetunep);
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







// function that handle changing the floats that going in freq inlet into frequency
void dsaws_changefreq(t_dsaws* x, double freq){
    if(freq > 0){
        
        dsaws_detune(x, freq, x->detunep); // calculating every voice's si
        x->freqp = freq;
        // x->val = arg; arg is replacing the value in x->val. (i used to write backward...) read from right to left
        
    }
    else{
        error("please enter float frequency > 0.\n");
        
    }
    
}

// function that handle changing the floats that going in detune inlet into detune parameter
void dsaws_changedetune(t_dsaws* x, double detune){
    if(detune > 0)
    {
        dsaws_detune(x, x->freqp, detune);
        x->detunep = detune;
    
    }
    else
    {
        error("please enter a valid number > 0.\n");
    }
    
}
 

void dsawsz_float(t_dsaws* x, double input) // making the message that users send into different inlet works. (only message)
{
    long inlet_number = proxy_getinlet((t_object *)x);
    post("this is inlet %d.\n", inlet_number);
    if (inlet_number == 0){ //first inlet
        dsaws_changefreq(x, input);
    
    }
    else if (inlet_number == 1){ //second inlet
        
        x->voicenump = round(input); //round() 四捨五入 .0 （rounding）
    }
    else if (inlet_number == 2){ //third inlet
        
        dsaws_changedetune(x, input);
    }
    else
        object_error((t_object *)x, "please enter valid num to freq, voices, and detune num.");
}


void dsaws_int(t_dsaws* x, long n){ //allowing the users sending both integeter and float, this func is registered in ext_main method.
    dsawsz_float(x, n);
}



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
    x->w_connected[1] = count[1];
    x->w_connected[2] = count[2];
    
    //0 disconnected, 1 connected (inlet) decided by the users
   

  


    object_method(dsp64, gensym("dsp_add64"), x, dsaws_perform64, 0, NULL);
}


// this is the 64-bit perform method audio vectors
void dsaws_perform64(t_dsaws* x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) //double **ins is a pointer to an array of pointers, where each pointer represents a signal inlet. // paramater->takes signal
{
    t_double *inL = ins[0];        // we get audio for each inlet of the object from the **ins argument
    t_double *inM = ins[1];
    t_double *inR = ins[2];
    
    t_double *outL = outs[0];    // we get audio for each outlet of the object from the **outs argument
    t_double *outR = outs[1];
    
    long n = sampleframes;
    
    t_double inputL = *inL; // freq
    t_double inputM = *inM; // num of voices
    t_double inputR = *inR; // detune parameter

    
    
    for (int time=0; time < sampleframes; time++){ // how many samples user/maxmsp are going to take at a time

        // CALL THE FUNCTION CALCULATES EACH VOICE'S SI//
        if (x->w_connected[0]){ // if the 1st inlet is connected by a signal
            if(inputL > 0)
            {
                dsaws_changefreq(x, inputL); // it calls dsaws_detune func that calculating every voice's si
          
            }
            else
            {
                error("please enter frequency > 0.\n");
            }
        }
        //SET THE VOICE OF NUM INLET IF CONNECTED//
        if (x->w_connected[1]){
            if(inputM > 0)
            {
                
                x->voicenump = round(inputM);
          
            }
            else
            {
                error("please enter a valid number > 0.\n");
            }
        
        }
        //SET THE DETUNE INLET IF CONNECTED//
        if (x->w_connected[2]){
            if(inputR > 0){
                
                dsaws_changedetune(x, inputR);// taking input M as the detune parameter
                
            }
            else
            {
                error("please enter a valid number > 0, using oct.pc to think.\n");
            }
        }
        
        
        inL++;
        inM++;
        inR++;
       
        //ADD UP ALL THE SAW WAVES//
        float sum = 0; // initiallized the sum in phase[index]'s data type
        for (int index=0; index < x->voicenump; index++) //index: how many voices
        {
            
                process_saw(x, index);
                sum = sum + x->phase[index]; // adding all the saw waves -- when the for loop loops multiple times(depends on the voice num), it add up every saw wave
                
        }
            *outL++= sum / x->voicenump; // amplitude divided by the num of the voice to reduce amp.
            *outR++= sum / x->voicenump;
        }
    }




   
