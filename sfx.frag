#version 130
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0., 0.01, x); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _sq(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq(float a) { return clip(50.*_sin(a)); }
float _psq(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float quant(float a,float div,float invdiv) { return floor(div*a+.5)*invdiv; }
float freqC1(float note){ return 32.7 * pow(2.,note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return round(sin(.5*PI*float(n))); }
float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

const float BPM = 25.;
const float BPS = BPM/60.;
const float SPB = 60./BPM;

const float Fsample = 44100.;
const float Tsample = 1./Fsample;

const float filterthreshold = 1e-3;

//TEXCODE




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

float reverbFsaw3_IIR(float time, float f, float tL, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4)
{
    int imax = int(log(filterthreshold)/log(IIRgain));
    float delay[4] = float[4](IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    
    float sum = 0.;
    
    // 4 IIR comb filters
    for(int d=0; d<8; d++)
    {
        float fac = 1.;
        
        for(int i=0; i<imax; i++)
        {
            float _TIME = time - float(i)*delay[d] * (.8 + .4*pseudorandom(sum));
            sum += fac*(theta(_TIME*SPB)*exp(-8.*_TIME*SPB)*((.5+(.5*_psq(8.*_TIME*SPB)))*(0.+(1.*(2.*fract(f*_TIME+0.)-1.)))));
            fac *= -IIRgain;
        }
    }
    return .25*sum;
}

float reverbFsaw3_AP1(float time, float f, float tL, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1)
{
    // first allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbFsaw3_IIR(_TIME, f, tL, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
    float fac = 1. - APgain * APgain;
    
    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));
    
    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel1 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbFsaw3_IIR(_TIME, f, tL, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}

float reverbFsaw3(float time, float f, float tL, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1, float APdel2)
{   // // based on this Schroeder Reverb from Paul Wittschen: http://www.paulwittschen.com/files/schroeder_paper.pdf
    // todo: add some noise...
    // second allpass delay line
    float _TIME = time;
    float sum = -APgain * reverbFsaw3_AP1(_TIME, f, tL, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
    float fac = 1. - APgain * APgain;

    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));

    for(int i=0; i<imax; i++)
    {
        _TIME -= APdel2 * (.9 + 0.2*pseudorandom(time));
        sum += fac * reverbFsaw3_AP1(_TIME, f, tL, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);
        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));
    }
    return sum;        
}




float AMAYSYN(float t, float B, float Bon, float Boff, float note, int Bsyn, float Brel)
{
    float Bprog = B-Bon;
    float Bproc = Bprog/(Boff-Bon);
    float L = Boff-Bon;
    float tL = SPB*L;
    float _t = SPB*(B-Bon);
    float f = freqC1(note);
	float vel = 1.; //implement later

    float env = theta(B-Bon) * (1. - smoothstep(Boff, Boff+Brel, B));
	float s = _sin(t*f);

	if(Bsyn == 0){}
    else if(Bsyn == 2){
      s = theta(Bprog)*exp(-16.*mod(Bprog,.125))*theta(Bprog)*exp(-1.5*Bprog)*(s_atan((0.+(1.*(2.*fract(f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.01)*f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.033)*f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.04)*f*t+0.)-1.))))+.6*s_atan((0.+(1.*(2.*fract(.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.-.05)*.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.+.03)*.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.+.02)*.5*f*t+.01)-1.)))));
    }
    else if(Bsyn == 3){
      s = (0.+(1.*_sin(f*t)))
      +-.1*GAC(t,0.,1.,2.,-.5,3.,2.,2.,-.25)*(0.+(1.*_sin(f*t)))
      +.1*GAC(t,0.,1.,2.,-.5,3.,2.,2.,-.25)*supershape((0.+(1.*_sin(f*t))),1.,.01,.7,.1,.6,.8);
    }
    else if(Bsyn == 6){
      s = reverbFsaw3(_t,f,tL,.1,.00297,.00371,.00411,.00437,.3,.00017,.0005);
env = theta(B-Bon)*pow(1.-smoothstep(Boff, Boff+Brel, B),.005);
    }
    
    
	return clamp(env,0.,1.) * s_atan(s);
}


uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;
uniform int iTexS;
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
float rfloat(float off)
{
    float d = rshort(off);
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

#define NTRK 4.
#define NMOD 24.
#define NPTN 5.
#define NNOT 166.

int trk_sep(int index){return int(rfloat(float(index)));}
int trk_syn(int index){return int(rfloat(NTRK+1.+float(index)));}
float trk_norm(int index){return rfloat(2.*NTRK+1.+float(index));}
float trk_rel(int index){return rfloat(3.*NTRK+1.+float(index));}
float mod_on(int index){return rfloat(4.*NTRK+1.+float(index));}
float mod_off(int index){return rfloat(4.*NTRK+NMOD+1.+float(index));}
int mod_ptn(int index){return int(rfloat(4.*NTRK+2.*NMOD+1.+float(index)));}
float mod_transp(int index){return rfloat(4.*NTRK+3.*NMOD+1.+float(index));}
int ptn_sep(int index){return int(rfloat(4.*NTRK+4.*NMOD+1.+float(index)));}
float note_on(int index){return rfloat(4.*NTRK+4.*NMOD+NPTN+2.+float(index));}
float note_off(int index){return rfloat(4.*NTRK+4.*NMOD+NPTN+NNOT+2.+float(index));}
float note_pitch(int index){return rfloat(4.*NTRK+4.*NMOD+NPTN+2.*NNOT+2.+float(index));}
float note_vel(int index){return rfloat(4.*NTRK+4.*NMOD+NPTN+3.*NNOT+2.+float(index));}

float mainSynth(float time)
{
    float max_mod_off = 26.;
    int drum_index = 15;
    float drum_synths = 2.;
    
    
    float r = 0.;
    float d = 0.;

    // mod for looping
    float BT = mod(BPS * time, max_mod_off);
    if(BT > max_mod_off) return r;
    time = SPB * BT;

    float r_sidechain = 1.;

    float Bon = 0.;
    float Boff = 0.;

    int iNTRK = int(NTRK);
    for(int trk = 0; trk < iNTRK; trk++)
    {
        int tsep = trk_sep(trk);
        int tlen = trk_sep(trk+1) - tsep;

        int _modU = tlen-1;
        for(int i=0; i<tlen-1; i++) if(BT < mod_on(tsep + i)) {_modU = i; break;}
               
        int _modL = tlen-1;
        for(int i=0; i<tlen-1; i++) if(BT < mod_off(tsep + i) + trk_rel(trk)) {_modL = i; break;}
       
        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            float B = BT - mod_on(tsep + _mod);

            int ptn = mod_ptn(tsep + _mod);
            int psep = ptn_sep(ptn);
            int plen = ptn_sep(ptn+1) - psep;
            
            int _noteU = plen-1;
            for(int i=0; i<plen-1; i++) if(B < note_on(psep + i + 1) + trk_rel(trk)) {_noteU = i; break;}

            int _noteL = plen-1;
            for(int i=0; i<plen-1; i++) if(B <= note_off(psep + i ) + trk_rel(trk)) {_noteL = i; break;}
           
            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                Bon    = note_on(psep + _note);
                Boff   = note_off(psep + _note);

                if(trk_syn(trk) == drum_index)
                {
                    int Bdrum = int(mod(note_pitch(psep + _note), drum_synths));
                    float Bvel = note_vel(psep + _note) * pow(2., mod_transp(tsep + _mod)/6.);

                    //0 is for sidechaining - am I doing this right?
                    if(Bdrum == 0)
                        r_sidechain = 1. - smoothstep(Bon,Bon+1e-4,B) + smoothstep(Bon,Boff,B);
                    else
                        d += trk_norm(trk) * AMAYSYN(time, B, Bon, Boff, Bvel, -Bdrum, trk_rel(trk));
                }
                else
                {
                    r += trk_norm(trk) * AMAYSYN(time, B, Bon, Boff, note_pitch(psep+_note) + mod_transp(tsep+_mod), trk_syn(trk), trk_rel(trk));
                }
            }
        }
    }

    return s_atan(s_atan(r_sidechain * r + d));
}

vec2 mainSound(float t)
{
    //enhance the stereo feel
    float stereo_delay = 2e-4;
      
    return vec2(mainSynth(t), mainSynth(t-stereo_delay));
}
