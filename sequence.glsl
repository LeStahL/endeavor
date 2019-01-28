#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0., 0.01, x); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin_(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq(float a) { return clip(50.*_sin(a)); }
float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
////float quant(float a,float div) { return floor(div*a+.5)/div; }
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

const float sequence_texture[784] = float[784](0.,6.,9.,21.,24.,2.,8.,-1.,1.,.64990234375,.35009765625,0.,.239990234375,0.,0.,.7998046875,.0999755859375,0.,4.,8.,12.,16.,20.,0.,4.,8.,0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,12.,16.,20.,4.,8.,12.,16.,20.,24.,4.,8.,12.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,16.,20.,24.,0.,0.,0.,3.,3.,3.,1.,1.,1.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,4.,4.,4.,0.,-2.,0.,0.,0.,0.,0.,-2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,12.,12.,12.,0.,1.,65.,113.,124.,166.,0.,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,2.0625,2.125,2.1875,2.25,2.3125,2.375,2.4375,2.5,2.5625,2.625,2.6875,2.75,2.8125,2.875,2.9375,3.,3.0625,3.125,3.1875,3.25,3.3125,3.375,3.4375,3.5,3.5625,3.625,3.6875,3.75,3.8125,3.875,3.9375,0.,.03125,.09375,.125,.15625,.21875,.25,.28125,.34375,.375,.40625,.46875,.5,.53125,.59375,.625,.65625,.71875,.75,.78125,.84375,.875,.90625,.96875,1.,1.03125,1.09375,1.125,1.15625,1.21875,1.25,1.28125,1.34375,1.375,1.40625,1.46875,1.5,1.53125,1.59375,1.625,1.65625,1.71875,1.75,1.78125,1.84375,1.875,1.90625,1.96875,0.,.5,1.,1.25,1.5,2.,2.25,2.5,3.,3.25,3.5,0.,0.,0.,0.,.5,.5,.5,.5,.875,1.,1.,1.,1.,1.25,1.25,1.25,1.5,1.5,1.5,1.5,2.,2.,2.,2.,2.25,2.25,2.25,2.5,2.5,2.5,2.5,3.,3.,3.,3.,3.25,3.25,3.25,3.5,3.5,3.5,3.5,3.625,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,2.0625,2.125,2.1875,2.25,2.3125,2.375,2.4375,2.5,2.5625,2.625,2.6875,2.75,2.8125,2.875,2.9375,3.,3.0625,3.125,3.1875,3.25,3.3125,3.375,3.4375,3.5,3.5625,3.625,3.6875,3.75,3.8125,3.875,3.9375,4.,.125,.0625,.125,.25,.1875,.25,.375,.3125,.375,.5,.4375,.5,.625,.5625,.625,.75,.6875,.75,.875,.8125,.875,1.,.9375,1.,1.125,1.0625,1.125,1.25,1.1875,1.25,1.375,1.3125,1.375,1.5,1.4375,1.5,1.625,1.5625,1.625,1.75,1.6875,1.75,1.875,1.8125,1.875,2.,1.9375,2.,.5,1.,1.25,1.5,2.,2.25,2.5,3.,3.25,3.5,4.,.5,.5,.5,.5,1.,1.,.875,1.,1.,1.25,1.25,1.5,1.25,1.5,1.5,1.5,2.,2.,2.,2.,2.25,2.25,2.5,2.25,2.5,2.5,2.5,3.,3.,3.,3.,3.25,3.25,3.5,3.25,3.5,3.5,3.5,4.,4.,4.,4.,13.,37.,49.,54.,48.,49.,54.,48.,49.,42.,49.,52.,47.,50.,49.,45.,47.,42.,49.,54.,48.,49.,54.,48.,49.,42.,49.,57.,54.,61.,57.,66.,53.,49.,53.,42.,54.,53.,54.,47.,57.,56.,57.,49.,54.,61.,57.,66.,69.,54.,49.,54.,48.,49.,54.,48.,49.,42.,49.,52.,47.,54.,49.,57.,44.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,27.,28.,28.,21.,19.,21.,19.,26.,21.,25.,30.,21.,18.,26.,33.,40.,57.,21.,31.,38.,55.,19.,57.,40.,33.,61.,21.,31.,38.,19.,45.,38.,66.,26.,33.,40.,69.,21.,37.,44.,25.,49.,42.,68.,18.,33.,40.,61.,21.,30.,37.,18.,38.,45.,66.,19.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.);





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
float bandpassBPsaw1(float time, float f, float tL, float fcenter, float bw, float M)
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
        y += w*h*(0.+(1.*(2.*fract(f*_TIME+0.)-1.)));
    }
    
    return s_atan(M*M*y); // I DO NOT CARE ANYMORE
}




float AMAYSYN(float t, float B, float Bon, float Boff, float note, int Bsyn, float Bvel, float Brel)
{
    float Bprog = B-Bon;
    float Bproc = Bprog/(Boff-Bon);
    float L = Boff-Bon;
    float tL = SPB*L;
    float _t = SPB*(B-Bon);
    float f = freqC1(note);

    float env = theta(B-Bon) * (1. - smoothstep(Boff, Boff+Brel, B));
	float s = _sin(t*f);

	if(Bsyn == 0){}
    else if(Bsyn == 1){
      s = env_AHDSR(_t,tL,.01,1.,1.,0.,0.)*(0.+(1.*_sin_(f*t,(.5+(.5*(2.*fract(2.*B+0.)-1.)))*(0.+(1.*_sin_(.999*f*t,.35*(0.+(1.*_sin(f*t)))))))))
      +env_AHDSR(_t,tL,.01,1.,1.,0.,0.)*.4*(0.+(1.*_sin(.5*f*t)))
      +env_AHDSR(_t,tL,.01,1.,1.,0.,0.)*.4*(0.+(1.*_sin(.501*f*t)));
    }
    else if(Bsyn == 2){
      s = theta(Bprog)*exp(-16.*mod(Bprog,.125))*theta(Bprog)*exp(-1.5*Bprog)*(s_atan((0.+(1.*(2.*fract(f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.01)*f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.033)*f*t+0.)-1.)))+(0.+(1.*(2.*fract((1.-.04)*f*t+0.)-1.))))+.6*s_atan((0.+(1.*(2.*fract(.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.-.05)*.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.+.03)*.5*f*t+.01)-1.)))+(0.+(1.*(2.*fract((1.+.02)*.5*f*t+.01)-1.)))));
    }
    else if(Bsyn == 8){
      s = bandpassBPsaw1(_t,f,tL,(2000.+(1500.*_sin(.25*B))),10.,100.);
    }
    
    
	return clamp(env,0.,1.) * s_atan(s);
}


float rfloat(int off){return sequence_texture[off];}

#define NTRK 4
#define NMOD 24
#define NPTN 5
#define NNOT 166

int trk_sep(int index)      {return int(rfloat(index));}
int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}
float mod_on(int index)     {return     rfloat(index+1+4*NTRK);}
float mod_off(int index)    {return     rfloat(index+1+4*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloat(index+1+4*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloat(index+1+4*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloat(index+1+4*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloat(index+2+4*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloat(index+2+4*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloat(index+2+4*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_vel(int index)   {return     rfloat(index+2+4*NTRK+4*NMOD+NPTN+3*NNOT);}

float mainSynth(float time)
{
    float max_mod_off = 24.1;
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

    for(int trk = 0; trk < NTRK; trk++)
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

                    if(Bdrum == 0) // "sidechaining"
                        r_sidechain = 1. - smoothstep(Bon,Bon+1e-4,B) + smoothstep(Bon,Boff,B);
                    else
                        d += trk_norm(trk) * AMAYSYN(time, B, Bon, Boff, mod_transp(tsep + _mod),
                                                     -Bdrum, note_vel(psep + _note), trk_rel(trk));
                }
                else
                {
                    r += trk_norm(trk) * AMAYSYN(time, B, Bon, Boff, note_pitch(psep+_note) + mod_transp(tsep+_mod),
                                                 trk_syn(trk), note_vel(psep + _note), trk_rel(trk));
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
