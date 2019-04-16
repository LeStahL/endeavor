#version 130
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0.,1e-3,x); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq(float a) { return clip(50.*_sin(a)); }
float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float _saw(float a) { return (2.*fract(a) - 1.); }
float freqC1(float note){ return 32.7 * pow(2., note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return round(sin(.5*PI*float(n))); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }
float fhelp(float x) { return 1. + .333*x; } // 1. + .33333*x + .1*x*x + .02381*x*x*x + .00463*x*x*x*x;

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

const float BPM = 35.;
const float BPS = BPM/60.;
const float SPB = 60./BPM;

const float Fsample = 48000.; // PRODUCTION: CHANGE THIS BACK TO 44100.
const float Tsample = 1./Fsample;

const float stereo_delay = 2e-4; //enhance the stereo feel - this is experimental since I included the stereo functionality

const float filterthreshold = 1e-3;

//TEXCODE

float doubleslope(float t, float a, float d, float s)
{
    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);
}

float drop_phase(float time, float t1, float f0, float f1)
{
    float t = min(time, t1);
    float phi = f0*t + .5*(f1-f0)/t1*t*t;

    if(time > t1)
    {
        phi += f1 * (time - t1);
    }
    return phi;
}



// One-dimensional value noise from https://www.shadertoy.com/view/wdj3D1 (NR4)

float lpnoise(float t, float fq) // kudos to Dmitry Andreev - and'2014!
{
    t *= fq;

    float tt = fract(t);
    float tn = t - tt;
    tt = smoothstep(0.0, 1.0, tt);

    // does pseudorandom(...) somehow equal hash22 noise?
    float n0 = pseudorandom(floor(tn + 0.0) / fq);
    float n1 = pseudorandom(floor(tn + 1.0) / fq);

    return mix(n0, n1, tt);
}

float reverb_phase(float t, float amt)
{
    float r = lpnoise(t, 100.0) + 0.2*lpnoise(t, 550.0) + 0.1*lpnoise(t, 1050.0)*exp(-5.*t);
    return amt * r;
}

float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)
{
    float att = x/A;
    float dec = 1. - (1.-S)*(x-H-A)/D;
    float rel = (x <= L-R) ? 1. : (L-x)/R;
    return (x<A ? att : x<A+H ? 1 : x<A+H+D ? dec : x<=L-R ? S : x<=L ? (L-x)/R : 0.);
}


float s_atan(float a) { return 2./PI * atan(a); }
float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } 

float supershape(float s, float amt, float A, float B, float C, float D, float E)
{
    float w;
    float m = sign(s);
    s = abs(s);

    if(s<A) w = B * smoothstep(0.,A,s);
    else if(s<C) w = C + (B-C) * smoothstep(C,A,s);
    else if(s<=D) w = s;
    else if(s<=1.)
    {
        float _s = (s-D)/(1.-D);
        w = D + (E-D) * (1.5*_s*(1.-.33*_s*_s));
    }
    else return 1.;
    
    return m*mix(s,w,amt);
}

float GAC(float t, float offset, float a, float b, float c, float d, float e, float f, float g)
{
    t = t - offset;
    return t<0. ? 0. : a + b*t + c*t*t + d*_sin(e*t) + f*exp(-g*t);
}

float comp_SAW(int N, float inv_N) {return inv_N * minus1hochN(N);}
float comp_TRI(int N, float inv_N) {return N % 2 == 0 ? 0. : inv_N * inv_N * minus1hochNminus1halbe(N);}
float comp_SQU(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (minus1hochN(N)*_sin(PW*float(N)+.25) - 1.);}

float MADD(float t, float f, float phase, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
{
    float ret = 0.;
    float INR = keyF==1 ? 1./CO : f/CO;
    float IRESQ = keyF==1 ? 1./RES_Q : 1./(RES_Q*f);
    
    float p = f*t + phase;
    for(int N=1; N<=NMAX; N+=NINC)
    {
        float float_N = float(N);
        float inv_N = 1./float_N;
        float comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N)    +  (-MIX)  * comp_SAW(N,inv_N)
                       : MIX < 1. ?   MIX    * comp_TRI(N,inv_N)    + (1.-MIX) * comp_SQU(N,inv_N,PW)
                                  : (MIX-1.) * comp_HAE(N,inv_N,PW) + (2.-MIX) * comp_SQU(N,inv_N,PW);

        float filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));
        
        if(abs(filter_N*comp_mix) < 1e-6) break; //or is it wise to break already?
        
        ret += comp_mix * filter_N * (_sin(float_N * p) + _sin(float_N * p * (1.+DET)));
    }
    return s_atan(ret);
}

float QFM_FB(float PH, float FB) // my guessing of feedback coefficients, FB>0 'saw', FB<0 'sq'
{
    if(FB > 0.) return abs(FB) * .8*_sin(PH + .35*_sin(PH));
    else return abs(FB) * _sin(PH + .5*PI);
}

float QFM(float t, float f, float phase, float LV1, float LV2, float LV3, float LV4, float FR1, float FR2, float FR3, float FR4, float FB1, float FB2, float FB3, float FB4, float ALGO)
{
    int iALGO = int(ALGO);
    float PH1 = FR1 * f * t + phase;
    float PH2 = FR2 * f * t + phase;
    float PH3 = FR3 * f * t + phase;
    float PH4 = FR4 * f * t + phase;
    
    float LINK41 = 0., LINK42 = 0., LINK43 = 0., LINK32 = 0., LINK31 = 0., LINK21 = 0.; 
    if(iALGO == 1)       {LINK43 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 2)  {LINK42 = 1.; LINK32 = 1.; LINK21 = 1.;}    
    else if(iALGO == 3)  {LINK41 = 1.; LINK32 = 1.; LINK21 = 1.;}
    else if(iALGO == 4)  {LINK42 = 1.; LINK43 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 5)  {LINK41 = 1.; LINK31 = 1.; LINK21 = 1.;}
    else if(iALGO == 6)  {LINK43 = 1.; LINK32 = 1.;}
    else if(iALGO == 7)  {LINK43 = 1.; LINK32 = 1.; LINK31 = 1.;}
    else if(iALGO == 8)  {LINK21 = 1.; LINK43 = 1.;}
    else if(iALGO == 9)  {LINK43 = 1.; LINK42 = 1.; LINK41 = 1.;}
    else if(iALGO == 10) {LINK43 = 1.; LINK42 = 1.;}
    else if(iALGO == 11) {LINK43 = 1.;}

    float OP4 = LV4 * _sin(PH4 + QFM_FB(PH4, FB4));
    float OP3 = LV3 * _sin(PH3 + QFM_FB(PH3, FB3) + LINK43*OP4);
    float OP2 = LV2 * _sin(PH2 + QFM_FB(PH2, FB2) + LINK42*OP4 + LINK32*OP3);
    float OP1 = LV1 * _sin(PH1 + QFM_FB(PH1, FB1) + LINK41*OP4 + LINK31*OP3 + LINK32*OP2);
    
    float sum = OP1;
    if(LINK21 > 0.) sum += OP2;
    if(LINK31 + LINK32 > 0.) sum += OP3;
    if(LINK41 + LINK42 + LINK43 > 0.) sum += OP4;
    
    return s_atan(sum);
}

float reverbFsaw3_IIR(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4)
{
    int imax = int(log(filterthreshold)/log(IIRgain));
    float delay[4] = float[4](IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    
    float sum = 0.;
    
    // 4 IIR comb filters
    for(int d=0; d<4; d++)
    {
        float fac = 1.;
        
        for(int i=0; i<imax; i++)
        {
            float _TIME = time - float(i)*delay[d] * (.8 + .4*pseudorandom(sum));
            sum += fac*(theta(_TIME*SPB)*exp(-8.*_TIME*SPB)*((.5+(.5*_psq(8.*_TIME*SPB)))*(2.*fract(f*_TIME+0.)-1.)));
            fac *= -IIRgain;
        }
    }
    return .25*sum;
}

float reverbFsaw3_AP1(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1)
{
    // first allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbFsaw3_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    float fac = 1. - APgain * APgain;
    
    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));
    
    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel1 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbFsaw3_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}

float reverbFsaw3(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1, float APdel2)
{   // // based on this Schroeder Reverb from Paul Wittschen: http://www.paulwittschen.com/files/schroeder_paper.pdf
    // todo: add some noise...
    // second allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbFsaw3_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
    float fac = 1. - APgain * APgain;

    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));

    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel2 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbFsaw3_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}
float bandpassBPsaw1(float time, float f, float tL, float vel, float fcenter, float bw, float M)
{
    float y = 0.;
        
    float facM = 2.*PI/M;
    float facL = 2.*PI*Tsample * (fcenter - bw);
    float facH = 2.*PI*Tsample * (fcenter + bw);
    
    if(facL < 0.) facL = 0.;
    if(facH > PI) facH = PI;
    
    float _TIME, mm, w, h;
    
    M--;
    for(float m=1.; m<=M; m++)
    {
        mm = m - .5*M;
        w = .42 - .5 * cos(mm*facM) - .08 * cos(2.*mm*facM);
        h = 1./(PI*mm) * (sin(mm*facH) - sin(mm*facL));
        
        _TIME = time - m*Tsample;
        y += w*h*(2.*fract(f*_TIME+0.)-1.);
    }
    
    return s_atan(M*M*y); // I DO NOT CARE ANYMORE
}
float avglpBDbody3f(float time, float f, float tL, float vel, float N)
{    
    int iN = int(N);

    float _TIME = time;
    float avg = 0.;
    
    for(int i = 0; i < iN; i++)
    {
          _TIME = time - float(i)*Tsample;
          avg += s_atan(smoothstep(0.,.01,_TIME)*smoothstep(.3+.1,.1,_TIME)*MADD(_TIME,(150.+(60.-150.)*smoothstep(-.2, 0.,-_TIME)),5.,10,1,.8,1.,1.,1.,.1,0.,0.,1) + 1.5*.5*step(_TIME,.05)*_sin(_TIME*1100.*5.*_saw(_TIME*800.*5.)) + 1.5*(1.-exp(-1000.*_TIME))*exp(-40.*_TIME)*_sin((400.-200.*_TIME)*_TIME*_sin(1.*(150.+(60.-150.)*smoothstep(-.2, 0.,-_TIME))*_TIME))) / N;
    }
    return avg;
}
float reverbsnrrev_IIR(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4)
{
    int imax = int(log(filterthreshold)/log(IIRgain));
    float delay[4] = float[4](IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    
    float sum = 0.;
    
    // 4 IIR comb filters
    for(int d=0; d<4; d++)
    {
        float fac = 1.;
        
        for(int i=0; i<imax; i++)
        {
            float _TIME = time - float(i)*delay[d] * (.8 + .4*pseudorandom(sum));
            sum += fac*clip((1.+1.6)*(_tri(_TIME*(350.+(6000.-800.)*smoothstep(-.01,0.,-_TIME)+(800.-350.)*smoothstep(-.01-.01,-.01,-_TIME)))*smoothstep(-.1,-.01-.01,-_TIME) + .7*fract(sin(_TIME*90.)*4.5e4)*doubleslope(_TIME,.05,.3,.3),-1., 1.)*doubleslope(_TIME,0.,.25,.3));
            fac *= -IIRgain;
        }
    }
    return .25*sum;
}

float reverbsnrrev_AP1(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1)
{
    // first allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbsnrrev_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    float fac = 1. - APgain * APgain;
    
    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));
    
    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel1 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbsnrrev_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}

float reverbsnrrev(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1, float APdel2)
{   // // based on this Schroeder Reverb from Paul Wittschen: http://www.paulwittschen.com/files/schroeder_paper.pdf
    // todo: add some noise...
    // second allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbsnrrev_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
    float fac = 1. - APgain * APgain;

    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));

    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel2 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbsnrrev_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}


float protokick(float t, float f_start, float f_end, float fdecay, float hold, float decay, float drive, float detune, float rev_amount, float rev_hold, float rev_decay, float rev_drive)
{
    float phi = drop_phase(t, fdecay, f_start, f_end);
    float rev_phi = phi + reverb_phase(t, rev_amount);
    return clamp(drive*.5*(_sin(phi)+_sin((1.-detune)*phi)),-1.,1.) * exp(-max(t-hold, 0.)/decay)
         + clamp(rev_drive*.5*(_sin(rev_phi)+_sin((1.-detune)*rev_phi)),-1.,1.) * exp(-max(t-rev_hold, 0.)/rev_decay);
} 



uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;
uniform int iTexSize;
uniform sampler2D iSequence;
uniform float iSequenceWidth;

// Read short value from texture at index off
float rshort(float off)
{
    float hilo = mod(off, 2.);
    off *= .5;
    vec2 ind = (vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth))+.05)/iSequenceWidth;
    vec4 block = texture(iSequence, ind);
    vec2 data = mix(block.rg, block.ba, hilo);
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloat(int off)
{
    float d = rshort(float(off));
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

#define NTRK 8
#define NMOD 19
#define NPTN 7
#define NNOT 463

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float trk_slide(int index)  {return     rfloat(index+1+4*NTRK);}
float mod_on(int index)     {return     rfloat(index+1+5*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+5*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+5*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+5*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+5*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_pan(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+3*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+4*NNOT);}
float note_slide(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+5*NNOT);}

vec2 mainSynth(float time)
{
    float max_mod_off = 53.;
    int drum_index = 35;
    time += 10.2857;
    
    float sL = 0.;
    float sR = 0.;
    float dL = 0.;
    float dR = 0.;

    // mod for looping
    float BT = mod(BPS * time, max_mod_off);
    time = SPB * BT;
    
    float time2 = time - stereo_delay;
    float sidechain = 1.;

    float amaysynL, amaysynR, amaydrumL, amaydrumR, B, Bon, Boff, Bprog, Bproc, L, tL, _t, _t2, vel, rel, f, amtL, amtR;
    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, syn, drum;

    for(int trk = 0; trk < NTRK; trk++)
    {
        tsep0 = trk_sep(trk);
        tsep1 = trk_sep(trk + 1);

        syn = trk_syn(trk);
        rel = trk_rel(trk);
 
        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1)); _modU++);             
        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + rel); _modL++);

        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            B = BT - mod_on(_mod);

            ptn   = mod_ptn(_mod);
            psep0 = ptn_sep(ptn);
            psep1 = ptn_sep(ptn + 1);
                         
            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1)); _noteU++);
            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + rel); _noteL++);
            //here: could introduce "monosynth" mode that sets _noteL = _noteU

            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                amaysynL  = 0.;
                amaysynR  = 0.;
                amaydrumL = 0.;
                amaydrumR = 0.;

                Bon   = note_on(_note);
                Boff  = note_off(_note) + rel;
                L     = Boff - Bon;
                tL    = L * SPB;
                Bprog = B - Bon;
                Bproc = Bprog / L;
                _t    = Bprog * SPB;
                _t2   = _t - stereo_delay;
                vel   = note_vel(_note);
                amtL  = clamp(1. - note_pan(_note), 0., 1.);
                amtR  = clamp(1. + note_pan(_note), 0., 1.); 

                if(syn == drum_index)
                {
                    drum = int(note_pitch(_note));
                    if(drum == 0) { sidechain = min(sidechain, 1. - clamp(1e4 * Bprog,0.,1.) + pow(Bprog/(L-rel),8.)); }
                    else if(drum == 4){
                        amaydrumL = .2*fract(sin(_t*100.*.9)*50000.*.9)*doubleslope(_t,.03,.1,.1);
                        amaydrumR = .2*fract(sin(_t2*100.*.9)*50000.*.9)*doubleslope(_t2,.03,.1,.1);
                    }
                    else if(drum == 5){
                        amaydrumL = .4*(.6+(.25*_psq(4.*B)))*fract(sin(_t*100.*.3)*50000.*2.)*doubleslope(_t,0.,.05,0.);
                        amaydrumR = .4*(.6+(.25*_psq(4.*B)))*fract(sin(_t2*100.*.3)*50000.*2.)*doubleslope(_t2,0.,.05,0.);
                    }
                    
                    dL += amtL * s_atan(trk_norm(trk) * amaydrumL * theta(Bprog));
                    dR += amtR * s_atan(trk_norm(trk) * amaydrumR * theta(Bprog));
                }
                else
                {
                    f = freqC1(note_pitch(_note) + mod_transp(_mod));

                    if(abs(note_slide(_note)) > 1e-3) // THIS IS SLIDEY BIZ
                    {
                        float Bslide = trk_slide(trk);
                        float fac = note_slide(_note) * log(2.)/12.;
                        if (Bprog <= Bslide)
                        {
                            float help = 1. - Bprog/Bslide;
                            f *= Bslide * (fhelp(fac) - help * fhelp(fac*help*help)) / Bprog;
                        }
                        else
                        {
                            f *= 1. + (Bslide * (fhelp(fac)-1.)) / Bprog;
                        }
                    }

                    float env = theta(Bprog) * (1. - smoothstep(Boff-rel, Boff, B));
                    if(syn == 0){amaysynL = _sin(f*_t); amaysynR = _sin(f*_t2);}
                    else if(syn == 4){
                        amaysynL = env_AHDSR(_t,tL,1.3,0.,.1,1.,.4)*s_atan((1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))))
      +.4*env_AHDSR(_t,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)));
                        amaysynR = env_AHDSR(_t2,tL,1.3,0.,.1,1.,.4)*s_atan((1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))))
      +.4*env_AHDSR(_t2,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)));
                    }
                    else if(syn == 12){
                        amaysynL = env_AHDSR(_t,tL,.002,0.,.15,.25,.13)*bandpassBPsaw1(_t,f,tL,vel,(2000.+(1500.*_sin(.25*B))),10.,100.);
                        amaysynR = env_AHDSR(_t2,tL,.002,0.,.15,.25,.13)*bandpassBPsaw1(_t2,f,tL,vel,(2000.+(1500.*_sin(.25*B))),10.,100.);
                    }
                    else if(syn == 32){
                        amaysynL = ((QFM((_t-0.0*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-0.0*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-4.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-4.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t-8.0e-03*(1.+3.*_sin(.1*_t))),f,0.,.00787*71.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t-8.0e-03*(1.+3.*_sin(.1*_t))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.))*env_AHDSR(_t,tL,.082,.001,.062,.521,.153));
                        amaysynR = ((QFM((_t2-0.0*(1.+3.*_sin(.1*_t2))),f,0.,.00787*71.,.00787*env_AHDSR((_t2-0.0*(1.+3.*_sin(.1*_t2))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t2-0.0*(1.+3.*_sin(.1*_t2))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t2-0.0*(1.+3.*_sin(.1*_t2))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t2-4.0e-03*(1.+3.*_sin(.1*_t2))),f,0.,.00787*71.,.00787*env_AHDSR((_t2-4.0e-03*(1.+3.*_sin(.1*_t2))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t2-4.0e-03*(1.+3.*_sin(.1*_t2))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t2-4.0e-03*(1.+3.*_sin(.1*_t2))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.)+QFM((_t2-8.0e-03*(1.+3.*_sin(.1*_t2))),f,0.,.00787*71.,.00787*env_AHDSR((_t2-8.0e-03*(1.+3.*_sin(.1*_t2))),tL,.18,.165,.229,.056,0.)*51.,.00787*env_AHDSR((_t2-8.0e-03*(1.+3.*_sin(.1*_t2))),tL,.059,.205,.04,.098,0.)*5.,.00787*env_AHDSR((_t2-8.0e-03*(1.+3.*_sin(.1*_t2))),tL,.108,.165,.094,.113,0.)*44.,.5,1.,1.001,1.,.00787*103.,.00787*20.,.00787*93.,.00787*92.,4.))*env_AHDSR(_t2,tL,.082,.001,.062,.521,.153));
                    }
                    else if(syn == 33){
                        amaysynL = (MADD(_t,f,0.,32,1,-1.,69.85,2.23,.47,2.79,.01,(.4+(.25*_sin(3.*B))),0)*(1.0*(s_atan(GAC(((_t-0.0)-0.0*(1.+.5*_sin(.5*(_t-0.0)))),0.,-.913,.201,.763,.773,.019,.043,-.704)*MADD(((_t-0.0)-0.0*(1.+.5*_sin(.5*(_t-0.0)))),f,0.,32,1,-1.,69.85,2.23,.47,2.79,.01,(.4+(.25*_sin(3.*B))),0))))+.4*env_AHDSR(_t,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(1.02*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.4,69.85,2.23,.47,2.79,0.,.39,0),.05*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))));
                        amaysynR = (MADD(_t2,f,0.,32,1,-1.,69.85,2.23,.47,2.79,.01,(.4+(.25*_sin(3.*B))),0)*(1.0*(s_atan(GAC(((_t2-0.0)-0.0*(1.+.5*_sin(.5*(_t2-0.0)))),0.,-.913,.201,.763,.773,.019,.043,-.704)*MADD(((_t2-0.0)-0.0*(1.+.5*_sin(.5*(_t2-0.0)))),f,0.,32,1,-1.,69.85,2.23,.47,2.79,.01,(.4+(.25*_sin(3.*B))),0))))+.4*env_AHDSR(_t2,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(1.02*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.4,69.85,2.23,.47,2.79,0.,.39,0),.05*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))));
                    }
                    else if(syn == 34){
                        amaysynL = (env_AHDSR(_t,tL,1.3,0.,.1,1.,.4)*s_atan((1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))))+.4*env_AHDSR(_t,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(2.83*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))));
                        amaysynR = (env_AHDSR(_t2,tL,1.3,0.,.1,1.,.4)*s_atan((1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9)))+(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))))+.4*env_AHDSR(_t2,tL,1.3,0.,.1,1.,.4)*(1.0*s_atan(2.83*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4)*supershape(MADD((_t2-0.0),f,0.,200,1,-.59,171.92,2.27,4.68,2.91,.02,.88,0),.08*env_AHDSR((_t2-0.0),tL,1.3,0.,.1,1.,.4),0.,.8,.05,1.,.9))));
                    }
                    
                    sL += amtL * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynL);
                    sR += amtR * s_atan(trk_norm(trk) * clamp(env,0.,1.) * amaysynR);
                }
            }
        }
    }
    return vec2(s_atan(sidechain * sL + dL), s_atan(sidechain * sR + dR));
}

vec2 mainSound(float t)
{
    return mainSynth(t);
}

void main()
{
   float t = (iBlockOffset + (gl_FragCoord.x - .5) + (gl_FragCoord.y - .5)*iTexSize)/iSampleRate;
   vec2 y = mainSound( t );
   vec2 v  = floor((0.5+0.5*y)*65535.0);
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = floor(v/256.0)/255.0;
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}
