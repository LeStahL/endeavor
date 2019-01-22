/*  Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
 *  Copyright (C) 2017  QM <TODO>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
 
#version 130

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;
uniform int iTexS;
#define PI radians(180.)
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0., 0.01, x); }
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _sin(float a, float p) { return sin(2. * PI * mod(a,1.) + p); }
float _unisin(float a,float b) { return (.5*_sin(a) + .5*_sin((1.+b)*a)); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _sq(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _psq(float a) { return clip(50.*_sin(a)); }
float _psq(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } 
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float _saw(float a) { return (2.*fract(a) - 1.); }
float quant(float a,float div,float invdiv) { return floor(div*a+.5)*invdiv; }
float quanti(float a,float div) { return floor(div*a+.5)/div; }
float freqC1(float note){ return 32.7 * pow(2.,note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }
float minus1hochNminus1halbe(int n) { return round(sin(.5*PI*float(n))); }

#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d

const float BPM = 15.;
const float BPS = BPM/60.;
const float SPB = 60./BPM;

const float Fsample = 44100.; // I think?
const float Tsample = 2.267573696e-5;

float doubleslope(float t, float a, float d, float s)
{
    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);
}

float env_AHD(float t, float a, float h, float d)
{
    return t<a ? t/a : t<a+h ? 1. : t < a+h+d ? 1.+(1.-t)*(a+h)/d : 0.;
}

float env_ADSR(float x, float L, float A, float D, float S, float R)
{
    float att = x/A;
    float dec = 1. - (1.-S)*(x-A)/D;
    float rel = (x <= L-R) ? 1. : (L-x)/R;
    return x<A ? att : x<A+D ? dec : x<= L-R ? S : x<=L ? (L-x)/R : 0.;
}

float env_ADSRexp(float x, float L, float A, float D, float S, float R)
{
    float att = pow(x/A,8.);
    float dec = S + (1.-S) * exp(-(x-A)/D);
    float rel = (x <= L-R) ? 1. : pow((L-x)/R,4.);
    return (x < A ? att : dec) * rel;    
}

float s_atan(float a) { return 2./PI * atan(a); }
float s_crzy(float amp) { return clamp( s_atan(amp) - 0.1*cos(0.9*amp*exp(amp)), -1., 1.); }
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
    return t<0. ? 0. : a + b*t + c*t*t + d*sin(e*t) + f*exp(-g*t);
}

float comp_SAW(int N, float inv_N) {return inv_N * minus1hochN(N);}
float comp_TRI(int N, float inv_N) {return N % 2 == 0 ? 0. : inv_N * inv_N * minus1hochNminus1halbe(N);}
float comp_SQU(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}
float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (minus1hochN(N)*_sin(PW*float(N)+.25) - 1.);}

float MACESQ(float t, float f, float phase, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)
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


float TRISQ(float t, float f, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
{
    float ret = 0.;
   
    int Ninc = 1; // try this: leaving out harmonics...
   
    for(int N=0; N<=MAXN; N+=Ninc)
    {
        float mode     = 2.*float(N) + 1.;
        float inv_mode = 1./mode;         // avoid division? save table of Nmax <= 20 in some array or whatever
        float comp_TRI = (N % 2 == 1 ? -1. : 1.) * inv_mode*inv_mode;
        float comp_SQU = inv_mode;
        float filter_N = pow(1. + pow(float(N) * INR,2.*NDECAY),-.5) + RES * exp(-pow(float(N)*INR*RES_Q,2.));

        ret += (MIX * comp_TRI + (1.-MIX) * comp_SQU) * filter_N * _sin(mode * f * t);
    }
   
    return ret;
}

// CHEERS TO metabog https://www.shadertoy.com/view/XljSD3 - thanks for letting me steal
float resolpsomesaw1(float time, float f, float tL, float fa, float reso)
{
    int maxTaps = 128;
    fa = sqrt(fa*Tsample);
    float c = pow(0.5, (128.0-fa*128.0)  / 16.0);
    float r = pow(0.5, (reso*128.0+24.0) / 16.0);
    
    float v0 = 0.;
    float v1 = 0.;
    
    for(int i = 0; i < maxTaps; i++)
    {
          float _TIME = time - float(maxTaps-i)/Fsample;
          float Bprog = _TIME * BPS; //might need that
          float inp = (2.*fract(f*_TIME+0.)-1.);
          v0 =  (1.0-r*c)*v0  -  (c)*v1  + (c)*inp;
          v1 =  (1.0-r*c)*v1  +  (c)*v0;
    }
    return v1;
}
float hardkick(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
   
    float f   = 60. + 150. * smoothstep(-0.3, 0., -t);
    float env = smoothstep(0.,0.01,t) * smoothstep(-0.1, 0.2, 0.3 - t);
   
    float kick_body = env * .1*TRISQ(t, f, 100, 1., 1., .1, 16., 10.); // more heavy bass drum: increase reso parameters?
  
    kick_body += .7 * (smoothstep(0.,0.01,t) * smoothstep(-0.2, 0.2, 0.3 - t)) * _sin(t*f*.5);

    float kick_click = 1.5*step(t,0.05) * _sin(t*5000. * _saw(t*1000.));
   
    kick_click = s_atan(40.*(1.-exp(-1000.*t))*exp(-80.*t) * _sin((1200.-1000.*sin(1000.*t*sin(30.*t)))*t));
   
    float kick_blobb = s_crzy(10.*(1.-exp(-1000.*t))*exp(-30.*t) * _sin((300.-300.*t)*t));
   
    return vel * 2.*clamp(kick_body + kick_blobb + kick_click,-1.5,1.5);
}
// CHEERS TO metabog https://www.shadertoy.com/view/XljSD3 - thanks for letting me steal
float resolpsaw2D(float time, float f, float tL, float fa, float reso)
{
    int maxTaps = 128;
    fa = sqrt(fa*Tsample);
    float c = pow(0.5, (128.0-fa*128.0)  / 16.0);
    float r = pow(0.5, (reso*128.0+24.0) / 16.0);
    
    float v0 = 0.;
    float v1 = 0.;
    
    for(int i = 0; i < maxTaps; i++)
    {
          float _TIME = time - float(maxTaps-i)/Fsample;
          float Bprog = _TIME * BPS; //might need that
          float inp = s_atan((2.*fract((f+.3*_sin(5.*Bprog)*env_ADSR(_TIME,tL,.2,.3,.8,.2))*_TIME+0.)-1.)+(2.*fract((1.-.01)*(f+.3*_sin(5.*Bprog)*env_ADSR(_TIME,tL,.2,.3,.8,.2))*_TIME+0.)-1.)+(2.*fract((1.-.011)*(f+.3*_sin(5.*Bprog)*env_ADSR(_TIME,tL,.2,.3,.8,.2))*_TIME+0.)-1.)+(2.*fract((1.+.02)*(f+.3*_sin(5.*Bprog)*env_ADSR(_TIME,tL,.2,.3,.8,.2))*_TIME+0.)-1.));
          v0 =  (1.0-r*c)*v0  -  (c)*v1  + (c)*inp;
          v1 =  (1.0-r*c)*v1  +  (c)*v0;
    }
    return v1;
}
// CHEERS TO metabog https://www.shadertoy.com/view/XljSD3 - thanks for letting me steal
float resolpA1oscmixW(float time, float f, float tL, float fa, float reso)
{
    int maxTaps = 128;
    fa = sqrt(fa*Tsample);
    float c = pow(0.5, (128.0-fa*128.0)  / 16.0);
    float r = pow(0.5, (reso*128.0+24.0) / 16.0);
    
    float v0 = 0.;
    float v1 = 0.;
    
    for(int i = 0; i < maxTaps; i++)
    {
          float _TIME = time - float(maxTaps-i)/Fsample;
          float Bprog = _TIME * BPS; //might need that
          float inp = supershape((s_atan(_sq(.25*f*_TIME,.2*(2.*fract(2.*f*_TIME+.4*_tri(.5*f*_TIME+0.))-1.))+_sq((1.-.004)*.25*f*_TIME,.2*(2.*fract(2.*f*_TIME+.4*_tri(.5*f*_TIME+0.))-1.)))+.8*(2.*fract(2.*f*_TIME+.4*_tri(.5*f*_TIME+0.))-1.)),(2.*fract(2.*f*_TIME+.4*_tri(.5*f*_TIME+0.))-1.),.1,.3,.3,.8,.8);
          v0 =  (1.0-r*c)*v0  -  (c)*v1  + (c)*inp;
          v1 =  (1.0-r*c)*v1  +  (c)*v0;
    }
    return v1;
}
// CHEERS TO metabog https://www.shadertoy.com/view/XljSD3 - thanks for letting me steal
float resolpA24_mix(float time, float f, float tL, float fa, float reso)
{
    int maxTaps = 128;
    fa = sqrt(fa*Tsample);
    float c = pow(0.5, (128.0-fa*128.0)  / 16.0);
    float r = pow(0.5, (reso*128.0+24.0) / 16.0);
    
    float v0 = 0.;
    float v1 = 0.;
    
    for(int i = 0; i < maxTaps; i++)
    {
          float _TIME = time - float(maxTaps-i)/Fsample;
          float Bprog = _TIME * BPS; //might need that
          float inp = (.7*clip(2.5*(fract(16.*Bprog+0.)+0.))*(2.*fract(.99*f*_TIME+0.)-1.)+.7*clip(2.5*(fract(16.*Bprog+0.)+0.))*_sq(.5*f*_TIME,.2)+.7*clip(2.5*(fract(16.*Bprog+0.)+0.))*_sin(.48*f*_TIME,.25)+-.35);
          v0 =  (1.0-r*c)*v0  -  (c)*v1  + (c)*inp;
          v1 =  (1.0-r*c)*v1  +  (c)*v0;
    }
    return v1;
}


float bitexplosion(float time, float B, int dmaxN, float fvar, float B2amt, float var1, float var2, float var3, float decvar)
{
    float snd = 0.;
    float B2 = mod(B,2.);
    float f = 60.*fvar;
	float dt = var1 * 2.*PI/15. * B/sqrt(10.*var2-.5*var3*B);
    int maxN = 10 + dmaxN;
    for(int i=0; i<2*maxN+1; i++)
    {
        float t = time + float(i - maxN)*dt;
        snd += _sin(f*t + .5*(1.+B2amt*B2)*_sin(.5*f*t));
    }
    float env = exp(-2.*decvar*B);
    return atan(snd * env);
}

float AMAYSYN(float t, float B, float Bon, float Boff, float note, int Bsyn)
{
    float Bprog = B-Bon;
    float Bproc = Bprog/(Boff-Bon);
    float L = Boff-Bon;
    float tL = SPB*L;
    float _t = SPB*(B-Bon);
    float f = freqC1(note);
	float vel = 1.;

    float env = theta(B-Bon) * theta(Boff-B);
	float s = _sin(t*f);

	if(Bsyn == 0){}
    else if(Bsyn == 2){
      s = theta(Bprog)*exp(-16.*mod(Bprog,.125))*theta(Bprog)*exp(-1.5*Bprog)*(s_atan((2.*fract(f*t+0.)-1.)+(2.*fract((1.-.01)*f*t+0.)-1.)+(2.*fract((1.-.033)*f*t+0.)-1.)+(2.*fract((1.-.04)*f*t+0.)-1.))+.6*s_atan((2.*fract(.5*f*t+.01)-1.)+(2.*fract((1.-.05)*.5*f*t+.01)-1.)+(2.*fract((1.+.03)*.5*f*t+.01)-1.)+(2.*fract((1.+.02)*.5*f*t+.01)-1.)));}	
    else if(Bsyn == 5){
      s = env_ADSR(_t,tL,.2,.3,.8,.2)*resolpsaw2D(_t,f,tL,300.*env_ADSR(Bprog,L,.5,.5,.4,0.),0.);}
    else if(Bsyn == 6){
      s = env_AHD(_t,0.,.2,.2)*resolpA1oscmixW(_t,f,tL,800.*env_AHD(_t,.0001,.0001,.4),.1)*5.;}
    else if(Bsyn == 12){
      s = _sin(4.*f*t,(.2+(.5*_sin(.21*Bprog)))*_sin(.25*f*t));}
    else if(Bsyn == 13){
      s = env_ADSRexp(Bprog,L,.75,.5,.1,.8)*floor(8.*MACESQ(floor(128.*_t+.5)/128.,f,0.,200,1,-1.,(500.+(4500.*env_ADSRexp(Bprog,L,.5,.5,.1,10.))),20.,0.,0.,.018,0.,0)+.5)*1.2e-01;}
    else if(Bsyn == 18){
      s = (200.+(300.*clip(4.*(fract(-16.*Bprog+0.)+.5))))*env_ADSR(_t,tL,0.,.2,.2,.5)*resolpA24_mix(_t,f,tL,(200.+.3*f),.3);}
    
    else if(Bsyn == -1){
          s = hardkick(_t,0.,1.);
    }
//       s = s_atan(vel*smoothstep(0.,.1,_t)*smoothstep(.1+.3,.3,_t)*(clip(10.*_tri((71.+(133.7-71.)*smoothstep(-.1, 0.,-_t))*_t))+_sin(.5*(71.+(133.7-71.)*smoothstep(-.1, 0.,-_t))*_t)))+1.2*step(_t,.05)*_sin(5000.*_t*.8*_saw(1000.*_t*.8));}
    else if(Bsyn == -2){
          s = hardkick(_t,0.,1.);
}
//       s = 3.*s_atan(vel*smoothstep(0.,.015,_t)*smoothstep(.1+.15,.15,_t)*MACESQ(_t,(50.+(200.-50.)*smoothstep(-.12, 0.,-_t)),5.,10,1,.8,1.,1.,1.,.1,.1,0.,1) + .4*.5*step(_t,.03)*_sin(_t*1100.*1.*_saw(_t*800.*1.)) + .4*(1.-exp(-1000.*_t))*exp(-40.*_t)*_sin((400.-200.*_t)*_t*_sin(1.*(50.+(200.-50.)*smoothstep(-.12, 0.,-_t))*_t)));}
    else if(Bsyn == -3){
      s = hardkick(_t,0.,1.);
      //2.*s_atan(vel*smoothstep(0.,.01,_t)*smoothstep(.3+.1,.1,_t)*MACESQ(_t,(60.+(210.-60.)*smoothstep(-.3, 0.,-_t)),5.,10,1,.8,1.,1.,1.,.1,0.,0.,1) + 1.5*.5*step(_t,.05)*_sin(_t*1100.*5.*_saw(_t*800.*5.)) + 1.5*(1.-exp(-1000.*_t))*exp(-40.*_t)*_sin((400.-200.*_t)*_t*_sin(1.*(60.+(210.-60.)*smoothstep(-.3, 0.,-_t))*_t)));}
      }
    else if(Bsyn == -4){
      s = .7*vel*fract(sin(t*100.*.9)*50000.*.9)*doubleslope(_t,.03,.15,.15);}
    else if(Bsyn == -5){
      s = vel*bitexplosion(t, Bprog, 1,2.,2.,1.5,2.,1.,1.);}
    else if(Bsyn == -6){
      s = .4*(.6+.25*_psq(4.*B,0.))*vel*fract(sin(t*100.*.3)*50000.*2.)*doubleslope(_t,0.,.05,0.);}
    else if(Bsyn == -7){
      s = vel*clamp(1.6*_tri(_t*(350.+(6000.-800.)*smoothstep(-.01,0.,-_t)+(800.-350.)*smoothstep(-.01-.01,-.01,-_t)))*smoothstep(-.1,-.01-.01,-_t) + .7*fract(sin(t*90.)*4.5e4)*doubleslope(_t,.05,.3,.3),-1., 1.)*doubleslope(_t,0.,.25,.3);}
    
	return clamp(env,0.,1.) * s_atan(s);
}

float BA8(float x, int pattern)
{
    x = mod(x,1.);
    float ret = 0.;
	for(int b = 0; b < 8; b++)
    	if ((pattern & (1<<b)) > 0) ret += step(x,float(7-b)/8.);
    return ret * .125;
}

float mainSynth(float time)
{
    int NO_trks = 8;
    int trk_sep[9] = int[9](0,1,5,7,15,19,25,27,39);
    int trk_syn[8] = int[8](19,6,5,6,2,2,12,19);
    float trk_norm[8] = float[8](.5,.3,.4,.3,.1,.1,.1,0.5);    
    float trk_rel[8] = float[8](0.,0.,.8,0.,.5,.5,0.,0.);
    float mod_on[39] = float[39](4.,0.,1.,2.,3.,0.,12.,0.,1.,2.,3.,16.,17.,18.,19.,8.,9.,10.,11.,6.,7.,8.,9.,10.,11.,12.,14.,4.,5.,6.,7.,12.,13.,14.,15.,16.,17.,18.,19.);
    float mod_off[39] = float[39](8.,1.,2.,3.,4.,4.,16.,1.,2.,3.,4.,17.,18.,19.,20.,9.,10.,11.,13.,7.,8.,9.,10.,11.,12.,14.,16.,5.,6.,7.,8.,13.,14.,15.,16.,17.,18.,19.,20.);
    int mod_ptn[39] = int[39](7,0,0,0,1,2,11,0,0,0,1,0,0,0,1,5,5,5,8,0,1,0,0,0,1,10,10,3,3,3,3,9,9,9,9,12,12,12,12);
    float mod_transp[39] = float[39](0.,12.,12.,12.,12.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-12.,-12.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    float max_mod_off = 22.;

    int drum_index = 19;
    float drum_synths = 8.;
    int NO_ptns = 13;
    int ptn_sep[14] = int[14](0,14,38,40,72,72,105,105,106,129,167,211,223,251);
    float note_on[251] = float[251](0.,.0625,.1875,.25,.3125,.34375,.4375,.5625,.6875,.75,.8125,.84375,.90625,.9375,0.,.0625,.09375,.125,.125,.1875,.25,.25,.3125,.34375,.34375,.4375,.5625,.59375,.625,.625,.6875,.71875,.75,.84375,.875,.90625,.9375,.96875,0.,0.,0.,.03125,.0625,.09375,.125,.15625,.1875,.21875,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.75,.78125,.8125,.84375,.875,.90625,.9375,.96875,0.,.03125,.03125,.0625,.09375,.125,.15625,.1875,.21875,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.75,.78125,.8125,.84375,.875,.90625,.9375,.96875,0.,0.,.03125,.0625,.09375,.125,.1875,.25,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.6875,.75,.8125,.875,.9375,0.,.03125,.0625,.09375,.125,.15625,.1875,.21875,.21875,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.71875,.75,.78125,.8125,.84375,.875,.90625,.90625,.9375,.9375,.9375,.96875,.96875,0.,.03125,.0625,.09375,.15625,.21875,.28125,.34375,.375,.4375,.5,.5,.53125,.5625,.59375,.65625,.71875,.78125,.84375,.875,.9375,.96875,1.,1.03125,1.0625,1.09375,1.15625,1.21875,1.28125,1.34375,1.375,1.4375,1.46875,1.5,1.53125,1.5625,1.59375,1.65625,1.71875,1.78125,1.84375,1.875,1.9375,1.96875,0.,.5,.875,1.,1.5,2.,2.5,2.875,3.,3.25,3.375,3.5,0.,.03125,.0625,.09375,.1875,.21875,.25,.28125,.3125,.34375,.34375,.4375,.46875,.5,.53125,.5625,.59375,.625,.6875,.71875,.75,.78125,.8125,.84375,.84375,.90625,.9375,.9375);
    float note_off[251] = float[251](.0625,.1875,.25,.3125,.34375,.4375,.5625,.6875,.75,.8125,.84375,.90625,.9375,1.,.0625,.09375,.125,.1875,.1875,.25,.3125,.3125,.34375,.4375,.5,.5625,.59375,.625,.6875,.6875,.71875,1.,.8125,.875,.90625,.9375,.96875,1.,4.,4.,.03125,.0625,.09375,.125,.15625,.1875,.21875,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.75,.78125,.8125,.84375,.875,.90625,.9375,.96875,1.,.03125,.0625,.0625,.09375,.125,.15625,.1875,.21875,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.75,.78125,.8125,.84375,.875,.90625,.9375,.96875,1.,4.,.03125,.0625,.09375,.125,.1875,.25,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.6875,.75,.9375,1.1875,1.25,1.4375,.03125,.0625,.09375,.125,.15625,.1875,.21875,.25,.25,.28125,.3125,.34375,.375,.40625,.4375,.46875,.5,.53125,.5625,.59375,.625,.65625,.6875,.71875,.75,.75,.78125,.8125,.84375,.875,.90625,.9375,.9375,.96875,.96875,.96875,1.,1.,.0625,.09375,.125,.15625,.21875,.28125,.34375,.40625,.4375,.5,.53125,.5625,.59375,.625,.65625,.71875,.78125,.84375,.90625,.9375,1.,1.,1.0625,1.09375,1.125,1.15625,1.21875,1.28125,1.34375,1.40625,1.4375,1.5,1.5,1.5625,1.59375,1.625,1.65625,1.71875,1.78125,1.84375,1.90625,1.9375,2.,2.,.46875,.84375,1.,1.46875,1.96875,2.46875,2.84375,3.,3.375,4.,3.5,4.,.03125,.0625,.09375,.125,.21875,.25,.28125,.3125,.34375,.375,.375,.46875,.5,.53125,.5625,.59375,.625,.65625,.71875,.75,.78125,.8125,.84375,.875,.875,.9375,.96875,.96875);
    float note_pitch[251] = float[251](21.,21.,21.,21.,21.,21.,21.,21.,21.,21.,21.,24.,21.,26.,21.,21.,21.,28.,45.,21.,43.,24.,21.,26.,42.,21.,21.,21.,40.,29.,21.,24.,38.,36.,38.,40.,45.,48.,21.,33.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,30.,40.,45.,40.,48.,40.,45.,57.,52.,48.,60.,59.,55.,50.,47.,48.,52.,47.,60.,55.,48.,40.,45.,57.,52.,48.,60.,59.,55.,64.,59.,60.,50.,52.,29.,40.,45.,48.,40.,45.,52.,60.,55.,50.,47.,48.,52.,47.,60.,55.,48.,40.,45.,52.,55.,57.,59.,55.,27.,30.,27.,30.,31.,30.,30.,27.,30.,27.,30.,27.,30.,31.,30.,27.,30.,27.,30.,27.,30.,31.,30.,30.,27.,30.,27.,30.,27.,30.,31.,27.,30.,31.,30.,27.,30.,27.,33.,21.,33.,21.,21.,21.,33.,21.,33.,28.,45.,31.,19.,31.,19.,19.,19.,31.,19.,31.,26.,43.,33.,21.,33.,21.,21.,21.,30.,18.,30.,25.,42.,26.,14.,26.,14.,14.,14.,26.,14.,26.,21.,38.,45.,43.,45.,49.,54.,45.,43.,45.,38.,21.,40.,33.,27.,30.,27.,30.,27.,30.,27.,30.,27.,30.,27.,27.,30.,31.,30.,27.,30.,31.,27.,30.,27.,30.,27.,31.,30.,27.,31.,30.);
    float note_vel[251] = float[251](1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.);
    
    float r = 0.;
    float d = 0.;

    // mod for looping
    float BT = mod(BPS * time, max_mod_off);
    if(BT > max_mod_off) return r;
    time = SPB * BT;

    float r_sidechain = 1.;

    float Bon = 0.;
    float Boff = 0.;

    for(int trk = 0; trk < NO_trks; trk++)
    {
        int TLEN = trk_sep[trk+1] - trk_sep[trk];

        int _modU = TLEN-1;
        for(int i=0; i<TLEN-1; i++) if(BT < mod_on[(trk_sep[trk]+i)]) {_modU = i; break;}
               
        int _modL = TLEN-1;
        for(int i=0; i<TLEN-1; i++) if(BT < mod_off[(trk_sep[trk]+i)] + trk_rel[trk]) {_modL = i; break;}
       
        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            float B = BT - mod_on[trk_sep[trk]+_mod];

            int ptn = mod_ptn[trk_sep[trk]+_mod];
            int PLEN = ptn_sep[ptn+1] - ptn_sep[ptn];
           
            int _noteU = PLEN-1;
            for(int i=0; i<PLEN-1; i++) if(B < note_on[(ptn_sep[ptn]+i+1)]) {_noteU = i; break;}

            int _noteL = PLEN-1;
            for(int i=0; i<PLEN-1; i++) if(B <= note_off[(ptn_sep[ptn]+i)] + trk_rel[trk]) {_noteL = i; break;}
           
            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                Bon    = note_on[(ptn_sep[ptn]+_note)];
                Boff   = note_off[(ptn_sep[ptn]+_note)] + trk_rel[trk];

                float anticlick = 1.-exp(-1000.*(B-Bon)); //multiply this elsewhere?

                if(trk_syn[trk] == drum_index)
                {
                    int Bdrum = int(mod(note_pitch[ptn_sep[ptn]+_note], drum_synths));
                    float Bvel = note_vel[(ptn_sep[ptn]+_note)] * pow(2.,mod_transp[trk_sep[trk]+_mod]/6.);

                    //0 is for sidechaining - am I doing this right?
                    if(Bdrum == 0)
                        r_sidechain = anticlick - .999 * theta(B-Bon) * smoothstep(Boff,Bon,B);
                    else
                        d += trk_norm[trk] * AMAYSYN(time, B, Bon, Boff, Bvel, -Bdrum);
                }
                else
                {
                    r += trk_norm[trk] * AMAYSYN(time, B, Bon, Boff,
                                                   note_pitch[(ptn_sep[ptn]+_note)] + 12. + mod_transp[trk_sep[trk]+_mod], trk_syn[trk]);
                }
            }
        }
    }

    return s_atan(s_atan(r_sidechain * r + d));
//    return sign(snd) * sqrt(abs(snd)); // eine von Matzes "besseren" Ideen
}

vec2 mainSound(float t)
{
    //enhance the stereo feel
    float stereo_delay = 2e-4;
      
    return vec2(mainSynth(t), mainSynth(t-stereo_delay));
}


void main() 
{
   // compute time `t` based on the pixel we're about to write
   // the 512.0 means the texture is 512 pixels across so it's
   // using a 2 dimensional texture, 512 samples per row
   float t = (iBlockOffset + (gl_FragCoord.x-0.5) + (gl_FragCoord.y-0.5)*float(iTexS))/iSampleRate;
    
//    t = mod(t, 4.5);
    
   // Get the 2 values for left and right channels
   vec2 y = iVolume * mainSound( t );

   // convert them from -1 to 1 to 0 to 65536
   vec2 v  = floor((0.5+0.5*y)*65535.0);

   // separate them into low and high bytes
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = mod(floor(v/256.0), 256.)/255.0;

   // write them out where 
   // RED   = channel 0 low byte
   // GREEN = channel 0 high byte
   // BLUE  = channel 1 low byte
   // ALPHA = channel 1 high byte
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}
