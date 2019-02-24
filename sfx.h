/* Generated with shader-compressor by NR4/Team210. */
#ifndef SFX_H
#define SFX_H
const char * sfx_frag =
"#version 130\n"
"#define PI radians(180.)\n"
"float clip(float a) { return clamp(a,-1.,1.); }\n"
"float theta(float x) { return clamp(x*1e4,0.,1.); }\n"
"float _sin(float a) { return sin(2. * PI * mod(a,1.)); }\n"
"float _sq(float a) { return sign(2.*fract(a) - 1.); }\n"
"float _sq_(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }\n"
"float _psq(float a) { return clip(50.*_sin(a)); }\n"
"float _psq_(float a, float pwm) { return clip(50.*(_sin(a) - pwm)); } \n"
"float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }\n"
"float _saw(float a) { return (2.*fract(a) - 1.); }\n"
"float freqC1(float note){ return 32.7 * pow(2., note/12.); }\n"
"float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }\n"
"float minus1hochNminus1halbe(int n) { return round(sin(.5*PI*float(n))); }\n"
"float pseudorandom(float x) { return fract(sin(dot(vec2(x),vec2(12.9898,78.233))) * 43758.5453); }\n"
"\n"
"#define pat4(a,b,c,d,x) mod(x,1.)<.25 ? a : mod(x,1.)<.5 ? b : mod(x,1.) < .75 ? c : d\n"
"\n"
"const float BPM = 35.;\n"
"const float BPS = BPM/60.;\n"
"const float SPB = 60./BPM;\n"
"\n"
"const float Fsample = 44100.;\n"
"const float Tsample = 1./Fsample;\n"
"\n"
"const float filterthreshold = 1e-3;\n"
"\n"
"//TEXCODE\n"
"\n"
"float doubleslope(float t, float a, float d, float s)\n"
"{\n"
"    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);\n"
"}\n"
"\n"
"\n"
"\n"
"// One-dimensional value noise from https://www.shadertoy.com/view/wdj3D1 (NR4)\n"
"\n"
"float lpnoise(float t, float fq) // kudos to Dmitry Andreev - and'2014!\n"
"{\n"
"    t *= fq;\n"
"\n"
"    float tt = fract(t);\n"
"    float tn = t - tt;\n"
"    tt = smoothstep(0.0, 1.0, tt);\n"
"\n"
"    // does pseudorandom(...) somehow equal hash22 noise?\n"
"    float n0 = pseudorandom(floor(tn + 0.0) / fq);\n"
"    float n1 = pseudorandom(floor(tn + 1.0) / fq);\n"
"\n"
"    return mix(n0, n1, tt);\n"
"}\n"
"\n"
"\n"
"float env_AHDSR(float x, float L, float A, float H, float D, float S, float R)\n"
"{\n"
"    float att = x/A;\n"
"    float dec = 1. - (1.-S)*(x-H-A)/D;\n"
"    float rel = (x <= L-R) ? 1. : (L-x)/R;\n"
"    return (x<A ? att : x<A+H ? 1 : x<A+H+D ? dec : x<=L-R ? S : x<=L ? (L-x)/R : 0.);\n"
"}\n"
"\n"
"\n"
"float s_atan(float a) { return 2./PI * atan(a); }\n"
"float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } \n"
"\n"
"float supershape(float s, float amt, float A, float B, float C, float D, float E)\n"
"{\n"
"    float w;\n"
"    float m = sign(s);\n"
"    s = abs(s);\n"
"\n"
"    if(s<A) w = B * smoothstep(0.,A,s);\n"
"    else if(s<C) w = C + (B-C) * smoothstep(C,A,s);\n"
"    else if(s<=D) w = s;\n"
"    else if(s<=1.)\n"
"    {\n"
"        float _s = (s-D)/(1.-D);\n"
"        w = D + (E-D) * (1.5*_s*(1.-.33*_s*_s));\n"
"    }\n"
"    else return 1.;\n"
"    \n"
"    return m*mix(s,w,amt);\n"
"}\n"
"\n"
"\n"
"float comp_SAW(int N, float inv_N) {return inv_N * minus1hochN(N);}\n"
"float comp_TRI(int N, float inv_N) {return N % 2 == 0 ? 0. : inv_N * inv_N * minus1hochNminus1halbe(N);}\n"
"float comp_SQU(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (1. - minus1hochNminus1halbe(N))*_sin(PW);}\n"
"float comp_HAE(int N, float inv_N, float PW) {return N % 2 == 0 ? 0. : inv_N * (minus1hochN(N)*_sin(PW*float(N)+.25) - 1.);}\n"
"\n"
"float MADD(float t, float f, float phase, int NMAX, int NINC, float MIX, float CO, float NDECAY, float RES, float RES_Q, float DET, float PW, int keyF)\n"
"{\n"
"    float ret = 0.;\n"
"    float INR = keyF==1 ? 1./CO : f/CO;\n"
"    float IRESQ = keyF==1 ? 1./RES_Q : 1./(RES_Q*f);\n"
"    \n"
"    float p = f*t + phase;\n"
"    for(int N=1; N<=NMAX; N+=NINC)\n"
"    {\n"
"        float float_N = float(N);\n"
"        float inv_N = 1./float_N;\n"
"        float comp_mix = MIX < 0. ? (MIX+1.) * comp_TRI(N,inv_N)    +  (-MIX)  * comp_SAW(N,inv_N)\n"
"                       : MIX < 1. ?   MIX    * comp_TRI(N,inv_N)    + (1.-MIX) * comp_SQU(N,inv_N,PW)\n"
"                                  : (MIX-1.) * comp_HAE(N,inv_N,PW) + (2.-MIX) * comp_SQU(N,inv_N,PW);\n"
"\n"
"        float filter_N = pow(1. + pow(float_N*INR,NDECAY),-.5) + RES * exp(-pow((float_N*f-CO)*IRESQ,2.));\n"
"        \n"
"        if(abs(filter_N*comp_mix) < 1e-6) break; //or is it wise to break already?\n"
"        \n"
"        ret += comp_mix * filter_N * (_sin(float_N * p) + _sin(float_N * p * (1.+DET)));\n"
"    }\n"
"    return s_atan(ret);\n"
"}\n"
"\n"
"float QFM_FB(float PH, float FB) // my guessing of feedback coefficients, FB>0 'saw', FB<0 'sq'\n"
"{\n"
"    if(FB > 0.) return abs(FB) * .8*_sin(PH + .35*_sin(PH));\n"
"    else return abs(FB) * _sin(PH + .5*PI);\n"
"}\n"
"\n"
"float QFM(float t, float f, float phase, float LV1, float LV2, float LV3, float LV4, float FR1, float FR2, float FR3, float FR4, float FB1, float FB2, float FB3, float FB4, float ALGO)\n"
"{\n"
"    int iALGO = int(ALGO);\n"
"    float PH1 = FR1 * f * t + phase;\n"
"    float PH2 = FR2 * f * t + phase;\n"
"    float PH3 = FR3 * f * t + phase;\n"
"    float PH4 = FR4 * f * t + phase;\n"
"    \n"
"    float LINK41 = 0., LINK42 = 0., LINK43 = 0., LINK32 = 0., LINK31 = 0., LINK21 = 0.; \n"
"    if(iALGO == 1)       {LINK43 = 1.; LINK32 = 1.; LINK21 = 1.;}\n"
"    else if(iALGO == 2)  {LINK42 = 1.; LINK32 = 1.; LINK21 = 1.;}    \n"
"    else if(iALGO == 3)  {LINK41 = 1.; LINK32 = 1.; LINK21 = 1.;}\n"
"    else if(iALGO == 4)  {LINK42 = 1.; LINK43 = 1.; LINK31 = 1.; LINK21 = 1.;}\n"
"    else if(iALGO == 5)  {LINK41 = 1.; LINK31 = 1.; LINK21 = 1.;}\n"
"    else if(iALGO == 6)  {LINK43 = 1.; LINK32 = 1.;}\n"
"    else if(iALGO == 7)  {LINK43 = 1.; LINK32 = 1.; LINK31 = 1.;}\n"
"    else if(iALGO == 8)  {LINK21 = 1.; LINK43 = 1.;}\n"
"    else if(iALGO == 9)  {LINK43 = 1.; LINK42 = 1.; LINK41 = 1.;}\n"
"    else if(iALGO == 10) {LINK43 = 1.; LINK42 = 1.;}\n"
"    else if(iALGO == 11) {LINK43 = 1.;}\n"
"\n"
"    float OP4 = LV4 * _sin(PH4 + QFM_FB(PH4, FB4));\n"
"    float OP3 = LV3 * _sin(PH3 + QFM_FB(PH3, FB3) + LINK43*OP4);\n"
"    float OP2 = LV2 * _sin(PH2 + QFM_FB(PH2, FB2) + LINK42*OP4 + LINK32*OP3);\n"
"    float OP1 = LV1 * _sin(PH1 + QFM_FB(PH1, FB1) + LINK41*OP4 + LINK31*OP3 + LINK32*OP2);\n"
"    \n"
"    float sum = OP1;\n"
"    if(LINK21 > 0.) sum += OP2;\n"
"    if(LINK31 + LINK32 > 0.) sum += OP3;\n"
"    if(LINK41 + LINK42 + LINK43 > 0.) sum += OP4;\n"
"    \n"
"    return s_atan(sum);\n"
"}\n"
"\n"
"float reverbFsaw3_IIR(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4)\n"
"{\n"
"    int imax = int(log(filterthreshold)/log(IIRgain));\n"
"    float delay[4] = float[4](IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"    \n"
"    float sum = 0.;\n"
"    \n"
"    // 4 IIR comb filters\n"
"    for(int d=0; d<8; d++)\n"
"    {\n"
"        float fac = 1.;\n"
"        \n"
"        for(int i=0; i<imax; i++)\n"
"        {\n"
"            float _TIME = time - float(i)*delay[d] * (.8 + .4*pseudorandom(sum));\n"
"            sum += fac*(theta(_TIME*SPB)*exp(-8.*_TIME*SPB)*((.5+(.5*_psq(8.*_TIME*SPB)))*(2.*fract(f*_TIME+0.)-1.)));\n"
"            fac *= -IIRgain;\n"
"        }\n"
"    }\n"
"    return .25*sum;\n"
"}\n"
"\n"
"float reverbFsaw3_AP1(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1)\n"
"{\n"
"    // first allpass delay line\n"
"    float _TIME = time;\n"
"    float sum = -APgain * reverbFsaw3_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"    float fac = 1. - APgain * APgain;\n"
"    \n"
"    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));\n"
"    \n"
"    for(int i=0; i<imax; i++)\n"
"    {\n"
"        _TIME -= APdel1 * (.9 + 0.2*pseudorandom(time));\n"
"        sum += fac * reverbFsaw3_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));\n"
"    }\n"
"    return sum;        \n"
"}\n"
"\n"
"float reverbFsaw3(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1, float APdel2)\n"
"{   // // based on this Schroeder Reverb from Paul Wittschen: http://www.paulwittschen.com/files/schroeder_paper.pdf\n"
"    // todo: add some noise...\n"
"    // second allpass delay line\n"
"    float _TIME = time;\n"
"    float sum = -APgain * reverbFsaw3_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);\n"
"    float fac = 1. - APgain * APgain;\n"
"\n"
"    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));\n"
"\n"
"    for(int i=0; i<imax; i++)\n"
"    {\n"
"        _TIME -= APdel2 * (.9 + 0.2*pseudorandom(time));\n"
"        sum += fac * reverbFsaw3_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);\n"
"        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));\n"
"    }\n"
"    return sum;        \n"
"}\n"
"float bandpassBPsaw1(float time, float f, float tL, float vel, float fcenter, float bw, float M)\n"
"{\n"
"    float y = 0.;\n"
"        \n"
"    float facM = 2.*PI/M;\n"
"    float facL = 2.*PI*Tsample * (fcenter - bw);\n"
"    float facH = 2.*PI*Tsample * (fcenter + bw);\n"
"    \n"
"    if(facL < 0.) facL = 0.;\n"
"    if(facH > PI) facH = PI;\n"
"    \n"
"    float _TIME, mm, w, h;\n"
"    \n"
"    M--;\n"
"    for(float m=1.; m<=M; m++)\n"
"    {\n"
"        mm = m - .5*M;\n"
"        w = .42 - .5 * cos(mm*facM) - .08 * cos(2.*mm*facM);\n"
"        h = 1./(PI*mm) * (sin(mm*facH) - sin(mm*facL));\n"
"        \n"
"        _TIME = time - m*Tsample;\n"
"        y += w*h*(2.*fract(f*_TIME+0.)-1.);\n"
"    }\n"
"    \n"
"    return s_atan(M*M*y); // I DO NOT CARE ANYMORE\n"
"}\n"
"float avglpBDbody3f(float time, float f, float tL, float vel, float N)\n"
"{    \n"
"    int iN = int(N);\n"
"\n"
"    float _TIME = time;\n"
"    float avg = 0.;\n"
"    \n"
"    for(int i = 0; i < iN; i++)\n"
"    {\n"
"          _TIME = time - float(i)*Tsample;\n"
"          avg += s_atan(smoothstep(0.,.01,_TIME)*smoothstep(.3+.1,.1,_TIME)*MADD(_TIME,(60.+(150.-60.)*smoothstep(-.2, 0.,-_TIME)),5.,10,1,.8,1.,1.,1.,.1,0.,0.,1) + 1.5*.5*step(_TIME,.05)*_sin(_TIME*1100.*5.*_saw(_TIME*800.*5.)) + 1.5*(1.-exp(-1000.*_TIME))*exp(-40.*_TIME)*_sin((400.-200.*_TIME)*_TIME*_sin(1.*(60.+(150.-60.)*smoothstep(-.2, 0.,-_TIME))*_TIME))) / N;\n"
"    }\n"
"    return avg;\n"
"}\n"
"float avglpBDbody3ff(float time, float f, float tL, float vel, float N)\n"
"{    \n"
"    int iN = int(N);\n"
"\n"
"    float _TIME = time;\n"
"    float avg = 0.;\n"
"    \n"
"    for(int i = 0; i < iN; i++)\n"
"    {\n"
"          _TIME = time - float(i)*Tsample;\n"
"          avg += avglpBDbody3f(_TIME,f,tL,vel,2.) / N;\n"
"    }\n"
"    return avg;\n"
"}\n"
"float reverbsnrrev_IIR(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4)\n"
"{\n"
"    int imax = int(log(filterthreshold)/log(IIRgain));\n"
"    float delay[4] = float[4](IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"    \n"
"    float sum = 0.;\n"
"    \n"
"    // 4 IIR comb filters\n"
"    for(int d=0; d<8; d++)\n"
"    {\n"
"        float fac = 1.;\n"
"        \n"
"        for(int i=0; i<imax; i++)\n"
"        {\n"
"            float _TIME = time - float(i)*delay[d] * (.8 + .4*pseudorandom(sum));\n"
"            sum += fac*clip((1.+1.6)*(_tri(_TIME*(350.+(6000.-800.)*smoothstep(-.01,0.,-_TIME)+(800.-350.)*smoothstep(-.01-.01,-.01,-_TIME)))*smoothstep(-.1,-.01-.01,-_TIME) + .7*fract(sin(_TIME*90.)*4.5e4)*doubleslope(_TIME,.05,.3,.3),-1., 1.)*doubleslope(_TIME,0.,.25,.3));\n"
"            fac *= -IIRgain;\n"
"        }\n"
"    }\n"
"    return .25*sum;\n"
"}\n"
"\n"
"float reverbsnrrev_AP1(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1)\n"
"{\n"
"    // first allpass delay line\n"
"    float _TIME = time;\n"
"    float sum = -APgain * reverbsnrrev_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"    float fac = 1. - APgain * APgain;\n"
"    \n"
"    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));\n"
"    \n"
"    for(int i=0; i<imax; i++)\n"
"    {\n"
"        _TIME -= APdel1 * (.9 + 0.2*pseudorandom(time));\n"
"        sum += fac * reverbsnrrev_IIR(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4);\n"
"        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));\n"
"    }\n"
"    return sum;        \n"
"}\n"
"\n"
"float reverbsnrrev(float time, float f, float tL, float vel, float IIRgain, float IIRdel1, float IIRdel2, float IIRdel3, float IIRdel4, float APgain, float APdel1, float APdel2)\n"
"{   // // based on this Schroeder Reverb from Paul Wittschen: http://www.paulwittschen.com/files/schroeder_paper.pdf\n"
"    // todo: add some noise...\n"
"    // second allpass delay line\n"
"    float _TIME = time;\n"
"    float sum = -APgain * reverbsnrrev_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);\n"
"    float fac = 1. - APgain * APgain;\n"
"\n"
"    int imax = 1 + int((log(filterthreshold)-log(fac))/log(APgain));\n"
"\n"
"    for(int i=0; i<imax; i++)\n"
"    {\n"
"        _TIME -= APdel2 * (.9 + 0.2*pseudorandom(time));\n"
"        sum += fac * reverbsnrrev_AP1(_TIME, f, tL, vel, IIRgain, IIRdel1, IIRdel2, IIRdel3, IIRdel4, APgain, APdel1);\n"
"        fac *= APgain * (1. + 0.01*pseudorandom(_TIME));\n"
"    }\n"
"    return sum;        \n"
"}\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"uniform float iBlockOffset;\n"
"uniform float iSampleRate;\n"
"uniform float iVolume;\n"
"uniform int iTexSize;\n"
"uniform sampler2D iSequence;\n"
"uniform float iSequenceWidth;\n"
"\n"
"// Read short value from texture at index off\n"
"float rshort(float off)\n"
"{\n"
"    float hilo = mod(off, 2.);\n"
"    off *= .5;\n"
"    vec2 ind = (vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth))+.05)/iSequenceWidth;\n"
"    vec4 block = texture(iSequence, ind);\n"
"    vec2 data = mix(block.rg, block.ba, hilo);\n"
"    return round(dot(vec2(255., 65280.), data));\n"
"}\n"
"\n"
"// Read float value from texture at index off\n"
"float rfloat(int off)\n"
"{\n"
"    float d = rshort(float(off));\n"
"    float sign = floor(d/32768.),\n"
"        exponent = floor(d/1024.-sign*32.),\n"
"        significand = d-sign*32768.-exponent*1024.;\n"
"\n"
"    if(exponent == 0.)\n"
"         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;\n"
"    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);\n"
"}\n"
"\n"
"#define NTRK 4\n"
"#define NMOD 19\n"
"#define NPTN 5\n"
"#define NNOT 62\n"
"\n"
"int trk_sep(int index)      {return int(rfloat(index));}\n"
"int trk_syn(int index)      {return int(rfloat(index+1+1*NTRK));}\n"
"float trk_norm(int index)   {return     rfloat(index+1+2*NTRK);}\n"
"float trk_rel(int index)    {return     rfloat(index+1+3*NTRK);}\n"
"float trk_slide(int index)  {return     rfloat(index+1+4*NTRK);}\n"
"float mod_on(int index)     {return     rfloat(index+1+5*NTRK);}\n"
"float mod_off(int index)    {return     rfloat(index+1+5*NTRK+1*NMOD);}\n"
"int mod_ptn(int index)      {return int(rfloat(index+1+5*NTRK+2*NMOD));}\n"
"float mod_transp(int index) {return     rfloat(index+1+5*NTRK+3*NMOD);}\n"
"int ptn_sep(int index)      {return int(rfloat(index+1+5*NTRK+4*NMOD));}\n"
"float note_on(int index)    {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN);}\n"
"float note_off(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+1*NNOT);}\n"
"float note_pitch(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+2*NNOT);}\n"
"float note_vel(int index)   {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+3*NNOT);}\n"
"float note_slide(int index) {return     rfloat(index+2+5*NTRK+4*NMOD+NPTN+4*NNOT);}\n"
"\n"
"float mainSynth(float time)\n"
"{\n"
"    float max_mod_off = 12.;\n"
"    int drum_index = 25;\n"
"    \n"
"    float r = 0.;\n"
"    float d = 0.;\n"
"\n"
"    // mod for looping\n"
"    float BT = mod(BPS * time, max_mod_off);\n"
"    if(BT > max_mod_off) return r;\n"
"    time = SPB * BT;\n"
"\n"
"    float r_sidechain = 1.;\n"
"    float amaysyn;\n"
"\n"
"    float B, Bon, Boff, Bprog, Bproc, L, tL, _t, vel, Bsyn, trel, f;\n"
"    int tsep0, tsep1, _modU, _modL, ptn, psep0, psep1, _noteU, _noteL, Bdrum;\n"
"\n"
"    for(int trk = 0; trk < NTRK; trk++)\n"
"    {\n"
"        tsep0 = trk_sep(trk);\n"
"        tsep1 = trk_sep(trk + 1);\n"
"        trel  = trk_rel(trk);\n"
" \n"
"        for(_modU = tsep0; (_modU < tsep1 - 1) && (BT > mod_on(_modU + 1)); _modU++);             \n"
"        for(_modL = tsep0; (_modL < tsep1 - 1) && (BT >= mod_off(_modL) + trel); _modL++);\n"
"\n"
"        for(int _mod = _modL; _mod <= _modU; _mod++)\n"
"        {\n"
"            B = BT - mod_on(_mod);\n"
"\n"
"            ptn   = mod_ptn(_mod);\n"
"            psep0 = ptn_sep(ptn);\n"
"            psep1 = ptn_sep(ptn + 1);\n"
"                         \n"
"            for(_noteU = psep0; (_noteU < psep1 - 1) && (B > note_on(_noteU + 1)); _noteU++);\n"
"            for(_noteL = psep0; (_noteL < psep1 - 1) && (B >= note_off(_noteL) + trel); _noteL++);\n"
"           \n"
"            for(int _note = _noteL; _note <= _noteU; _note++)\n"
"            {\n"
"                amaysyn = 0.;\n"
"                Bon     = note_on(_note);\n"
"                Boff    = note_off(_note);\n"
"                L       = Boff - Bon;\n"
"                tL      = L * SPB;\n"
"                Bprog   = B - Bon;\n"
"                Bproc   = Bprog / L;\n"
"                _t      = Bprog * SPB;\n"
"                vel     = note_vel(_note);\n"
"                Bsyn    = trk_syn(trk);\n"
"//                time = tL; // does every place use the right time / tL ??\n"
"\n"
"                Boff   += trel;\n"
"\n"
"                if(Bsyn == drum_index)\n"
"                {\n"
"                    amaysyn = 0.;\n"
"                    Bdrum = int(note_pitch(_note));\n"
"                    if(Bdrum == 0) { r_sidechain = min(r_sidechain, 1. - min(1e4 * Bprog,1.) + pow(Bproc,8.)); }\n"
"                    else if(Bdrum == 3){\n"
"                        amaysyn = 4.*avglpBDbody3ff(_t,f,tL,vel,2.);\n"
"                    }\n"
"                    \n"
"                    d += trk_norm(trk) * amaysyn * theta(Bprog);\n"
"                }\n"
"                else\n"
"                {\n"
"                    f = freqC1(note_pitch(_note) + mod_transp(_mod));\n"
"\n"
"                    if(false && abs(note_slide(_note))>1e-2) // THIS IS SLIDEY BIZ - BUT SHITTY AS CRAPPY FUCK RIGHT NOW!\n"
"                    {\n"
"                        float Bslide = trk_slide(trk);\n"
"                        float fac    = note_slide(_note) * log(2.)/36.;\n"
"                        float help   = 1. - min(Bprog, Bslide)/(Bslide);\n"
"                        float corr   = ((Bslide*SPB) * (1. + fac - help*(1.+fac*help*help)));\n"
"\n"
"                        //LIN PHASE SLIDEY\n"
"                        fac = note_slide(_note) * log(2.)/12.;\n"
"                        corr = exp(fac) * (Bslide*SPB)/fac * (1. - exp(-fac * (1.-help)));\n"
"\n"
"                        f = (B <= Bslide) ? f * corr / _t : f * (1. + corr / _t); \n"
"                    }\n"
"\n"
"                    float env = theta(Bprog) * (1. - smoothstep(Boff-trel, Boff, B));\n"
"                    if(Bsyn == 0){amaysyn = _sin(f*_t);}\n"
"                    else if(Bsyn == 4){\n"
"                        amaysyn = .8*env_AHDSR(_t,tL,.001,0.,.1,1.,.3)*(supershape(clip(1.6*QFM((_t-0.0*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)\n"
"      +supershape(clip(1.6*QFM((_t-2.0e-03*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8)\n"
"      +supershape(clip(1.6*QFM((_t-4.0e-03*(1.+2.*_sin(.15*_t))),f,0.,.00787*127.*pow(vel,12.*7.87e-3),.00787*112.*pow(vel,63.*7.87e-3),.00787*127.*pow(vel,26.*7.87e-3),.00787*96.*pow(vel,120.*7.87e-3),.5,1.,1.5,1.,.00787*0.,.00787*0.,.00787*0.,.00787*50.,8.)),.3,.2,.8,.4,.8,.8));\n"
"                    }\n"
"                    else if(Bsyn == 10){\n"
"                        amaysyn = env_AHDSR(_t,tL,.002,0.,.15,.25,.13)*bandpassBPsaw1(_t,f,tL,vel,(2000.+(1500.*_sin(.25*B))),10.,100.);\n"
"                    }\n"
"                    \n"
"                    r += trk_norm(trk) * clamp(env,0.,1.) * amaysyn;\n"
"                }\n"
"            }\n"
"        }\n"
"    }\n"
"    return s_atan(r_sidechain * r + d);\n"
"}\n"
"\n"
"vec2 mainSound(float t)\n"
"{\n"
"    //enhance the stereo feel\n"
"    float stereo_delay = 2e-4;\n"
"\n"
"    return vec2(mainSynth(t), mainSynth(t-stereo_delay));\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"   float t = (iBlockOffset + (gl_FragCoord.x - .5) + (gl_FragCoord.y - .5)*iTexSize)/iSampleRate;\n"
"   vec2 y = mainSound( t );\n"
"   vec2 v  = floor((0.5+0.5*y)*65535.0);\n"
"   vec2 vl = mod(v,256.0)/255.0;\n"
"   vec2 vh = floor(v/256.0)/255.0;\n"
"   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);\n"
"}\n"
"\n"
;
#endif
