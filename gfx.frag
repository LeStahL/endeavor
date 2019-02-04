/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#version 130

float iScale, iNBeats = 0.;
uniform float iTime;
uniform vec2 iResolution;
uniform sampler2D iFont;
uniform float iFontWidth;
uniform sampler2D iSequence;
uniform float iSequenceWidth;
uniform float iExecutableSize;

// Global constants
const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);
float a; // Aspect ratio
#define FSAA 2 // Antialiasing

// Global variables
// vec3 col = c.yyy;

// Read short value from texture at index off
float rshort(float off)
{
    // Parity of offset determines which byte is required.
    float hilo = mod(off, 2.);
    // Find the pixel offset your data is in (2 unsigned shorts per pixel).
    off *= .5;
    // - Determine texture coordinates.
    //     offset = i*iFontWidth+j for (i,j) in [0,iFontWidth]^2
    //     floor(offset/iFontWidth) = floor((i*iFontwidth+j)/iFontwidth)
    //                              = floor(i)+floor(j/iFontWidth) = i
    //     mod(offset, iFontWidth) = mod(i*iFontWidth + j, iFontWidth) = j
    // - For texture coordinates (i,j) has to be rescaled to [0,1].
    // - Also we need to add an extra small offset to the texture coordinate
    //   in order to always "hit" the right pixel. Pixel width is
    //     1./iFontWidth.
    //   Half of it is in the center of the pixel.
    vec2 ind = (vec2(mod(off, iFontWidth), floor(off/iFontWidth))+.05)/iFontWidth;
    // Get 4 bytes of data from the texture
    vec4 block = texture(iFont, ind);
    // Select the appropriate word
    vec2 data = mix(block.rg, block.ba, hilo);
    // Convert bytes to unsigned short. The lower bytes operate on 255,
    // the higher bytes operate on 65280, which is the maximum range 
    // of 65535 minus the lower 255.
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloat(float off)
{
    // Convert the bytes to unsigned short as first step.
    float d = rshort(off);
    
    // Convert bytes to IEEE 754 float16. That is
    // 1 sign bit, 5 bit exponent, 11 bit mantissa.
    // Also it has a weird conversion rule that is not evident at all.
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    // Return full float16
    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

// Read short value from texture at index off
float rshorts(float off)
{
    float hilo = mod(off, 2.);
    off *= .5;
    vec2 ind = (vec2(mod(off, iSequenceWidth), floor(off/iSequenceWidth))+.05)/iSequenceWidth;
    vec4 block = texture(iSequence, ind);
    vec2 data = mix(block.rg, block.ba, hilo);
    return round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
float rfloats(int off)
{
    float d = rshorts(float(off));
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    if(exponent == 0.)
         return mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    return mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
}

// TODO: COPY THIS FROM SFX SHADER TO ACHIEVE SYNC
#define NTRK 2
#define NMOD 5
#define NPTN 2
#define NNOT 32

int trk_sep(int index)      {return int(rfloats(index));}
int trk_syn(int index)      {return int(rfloats(index+1+1*NTRK));}
float trk_norm(int index)   {return     rfloats(index+1+2*NTRK);}
float trk_rel(int index)    {return     rfloats(index+1+3*NTRK);}
float mod_on(int index)     {return     rfloats(index+1+4*NTRK);}
float mod_off(int index)    {return     rfloats(index+1+4*NTRK+1*NMOD);}
int mod_ptn(int index)      {return int(rfloats(index+1+4*NTRK+2*NMOD));}
float mod_transp(int index) {return     rfloats(index+1+4*NTRK+3*NMOD);}
int ptn_sep(int index)      {return int(rfloats(index+1+4*NTRK+4*NMOD));}
float note_on(int index)    {return     rfloats(index+2+4*NTRK+4*NMOD+NPTN);}
float note_off(int index)   {return     rfloats(index+2+4*NTRK+4*NMOD+NPTN+1*NNOT);}
float note_pitch(int index) {return     rfloats(index+2+4*NTRK+4*NMOD+NPTN+2*NNOT);}
float note_vel(int index)   {return     rfloats(index+2+4*NTRK+4*NMOD+NPTN+3*NNOT);}

const float BPM = 35.;
const float BPS = BPM/60.;
const float SPB = 60./BPM;

// Extract drum signal
float scale(float t)
{
    float max_mod_off = 4.;
    int drum_index = 24;
    float d = 0.;

    // mod for looping
    float BT = mod(BPS * t, max_mod_off);
    float Bon = 0.;
    float Boff = 0.;
    
    for(int trk = 0; trk < NTRK; trk++)
    {
        if(trk_syn(trk) != drum_index) continue;
        int tsep = trk_sep(trk);
        int tlen = trk_sep(trk+1) - tsep;

        int _modU = tlen-1;
        for(int i=0; i<tlen-1; i++) if(BT < mod_on(tsep + i + 1)) {_modU = i; break;}
               
        int _modL = tlen-1;
        for(int i=0; i<tlen-1; i++) if(BT < mod_off(tsep + i) + trk_rel(trk)) {_modL = i; break;}
        
        for(int _mod = _modL; _mod <= _modU; _mod++)
        {
            float B = BT - mod_on(tsep + _mod);

            int ptn = mod_ptn(tsep + _mod);
            int psep = ptn_sep(ptn);
            int plen = ptn_sep(ptn+1) - psep;
            
            int _noteU = plen;
            for(int i=0; i<plen; i++) if(B < note_on(psep + i + 1)) {_noteU = i; break;}

            int _noteL = plen;
            for(int i=0; i<plen; i++) if(B <= note_off(psep + i ) + trk_rel(trk)) {_noteL = i; break;}
           
            for(int _note = _noteL; _note <= _noteU; _note++)
            {
                Bon    = note_on(psep + _note);
                Boff   = note_off(psep + _note);

                int Bdrum = int(note_pitch(psep + _note));
                if(Bdrum != 0)
                {
                    float env = (1.-smoothstep(Bon, Bon+.1, B));
                    d = max(d, env);
                    iNBeats += 1.;
                }
            }
            
            return d;
        }
    }
    return 0.;
}

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

// One-dimensional perlin noise
float snoise_1d(float t)
{
    float i = floor(t);
    t = fract(t);
    t = ((6.*t-15.)*t+10.)*t*t*t;
    return mix(-1.+2.*rand(i*c.xx), -1.+2.*rand((i+1.)*c.xx), t);
}

// Two-dimensional perlin noise
float snoise_2d(vec2 t)
{
    vec2 i = floor(t);
    t = fract(t);
    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise
    t = smoothstep(c.yy, c.xx, t); // TODO: add this for faster value noise
    vec2 v1 = vec2(rand(i), rand(i+c.xy)),
        v2 = vec2(rand(i+c.yx), rand(i+c.xx));
    v1 = c.zz+2.*mix(v1, v2, t.y);
    return mix(v1.x, v1.y, t.x);
}

// Multi-frequency simplex noise
float mfsnoise_2d(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    float n = 0.;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*snoise_2d(f*x) + sum;
        a = a*phi;
        n += 1.;
    }
    sum *= (1.-phi)/(1.-pow(phi, n));
    return sum;
}

// 3D rotational matrix
mat3 rot(vec3 p)
{
    return mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// add object to scene
vec2 add(vec2 sda, vec2 sdb)
{
    return mix(sda, sdb, step(sdb.x, sda.x));
}

vec2 sub(vec2 sda, vec2 sdb)
{
    return mix(-sda, sdb, step(sda.x, sdb.x));
}

// Distance to line segment
float lineseg(vec2 x, vec2 p1, vec2 p2)
{
    vec2 d = p2-p1;
    return length(x-mix(p1, p2, clamp(dot(x-p1, d)/dot(d,d),0.,1.)));
}

float lineseg(vec3 x, vec3 p1, vec3 p2)
{
    vec3 d = p2-p1;
    return length(x-mix(p1, p2, clamp(dot(x-p1, d)/dot(d,d),0.,1.)));
}

// Distance to circle
float circle(vec2 x, float r)
{
    return abs(length(x)-r);
}

// Distance to circle segment
float circlesegment(vec2 x, float r, float p0, float p1)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))
    	return abs(length(x)-r);
    return min(
        length(x-r*vec2(cos(p0), sin(p0))),
        length(x-r*vec2(cos(p1), sin(p1)))
        );
}

// compute distance to regular polygon
float dpoly_min(vec2 x, float N, float R)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    return R-length(x)*cos(t)/cos(.5*d);
}

// 2D box
float box(vec2 x, vec2 b)
{
  vec2 d = abs(x) - b;
  return length(max(d,c.yy)) + min(max(d.x,d.y),0.);
}

// Get glyph data from texture
float dglyph(vec2 x, float ordinal, float size)
{
    float dis = box(x, 2.*size*c.xx);
    if(dis > 0.)
        return dis+.5*size;

    // Find glyph offset in glyph index
    float nglyphs = rfloat(1.),
        offset = 0;
        
    for(float i=0.; i<nglyphs; i+=1.)
    {
        float ord = floor(rfloat(2.+2.*i));
        if(ord == ordinal)
        {
            offset = floor(rfloat(2.+2.*i+1.));
            break;
        }
    }
    
    if(offset == 0.) return 1.;
    
    // Get distance from glyph data
    float d = 1.;
    
    // Lines
    float nlines = floor(rfloat(offset));
    offset += 1.;
    for(float i=0.; i<nlines; i+=1.)
    {
        float x1 = rfloat(offset);
        offset += 1.;
        float y1 = rfloat(offset);
        offset += 1.;
        float x2 = rfloat(offset);
        offset += 1.;
        float y2 = rfloat(offset);
        offset += 1.;
        d = min(d, lineseg(x, size*vec2(x1,y1), size*vec2(x2, y2)));
    }
    
    // Circles
    float ncircles = floor(rfloat(offset));
    offset += 1.;
    for(float i=0.; i<ncircles; i+=1.)
    {
        float xc = rfloat(offset);
        offset += 1.;
        float yc = rfloat(offset);
        offset += 1.;
        float r = rfloat(offset);
        offset += 1.;
        d = min(d, circle(x-size*vec2(xc, yc), size*r));
    }
    
    // Circle segments
    float nsegments = floor(rfloat(offset));
    offset += 1.;
    for(float i=0.; i<nsegments; i+=1.)
    {
        float xc = rfloat(offset);
        offset += 1.;
        float yc = rfloat(offset);
        offset += 1.;
        float r = rfloat(offset);
        offset += 1.;
        float phi0 = rfloat(offset);
        offset += 1.;
        float phi1 = rfloat(offset);
        offset += 1.;
        d = min(d, circlesegment(x-size*vec2(xc,yc), size*r, phi0, phi1));
    }
    
    if(nlines+ncircles+nsegments == 0.)
        return dis;
    
    return d;
}

// Get distance to string from database
float dstring(vec2 x, float ordinal, float size)
{
    // Get string database offset
    float stroff0 = floor(rfloat(0.));
    
    // Return 1 if wrong ordinal is supplied
    float nstrings = floor(rfloat(stroff0));
    if(ordinal >= nstrings)
        return 1.;
        
    // Get offset and length of string from string database index
    float stroff = floor(rfloat(stroff0+1.+2.*ordinal));
    float len = floor(rfloat(stroff0+2.+2.*ordinal));
    
    /* Slower code
    float d = 1.;
    for(float i=0.; i<len; i+=1.)
        d = min(d, dglyph(x-2.1*i*size*c.xy,floor(rfloat(0.+stroff+i)), .8*size));
    return d;
    */
    
    // Draw glyphs
    vec2 dx = mod(x-size, 2.*size)-size, 
        ind = ceil((x-dx+size)/2./size);
    
    // Bounding box
    float bound = box(x-size*(len-3.)*c.xy, vec2(size*len, 1.*size));
    if(bound > 0.)
        return bound+.5*size;
    return dglyph(dx, floor(rfloat(stroff+ind.x)), .7*size);
}

// distance to a floating point number string
// for debugging stuff while shader is loaded
float dfloat(vec2 x, float num, float size)
{
    float d = 1., index = 0.;
    
    // Determine sign and output it if present
    float sign = sign(num), exp = 0.;
    if(sign<0.)
    {
        d = min(d, dglyph(x, 45., .7*size));
        index += 1.;
        num *= -1.;
    }
    
    // The first power of ten that floors num to anything not zero is the exponent
    for(exp = -15.; exp < 15.; exp += 1.)
        if(floor(num*pow(10.,exp)) != 0.)
            break;
    exp *= -1.;
    // Determine the significand and output it
    for(float i = exp; i >= exp-5.; i -= 1.)
    {
        float po = pow(10.,i);
        float ca = floor(num/po);
        num -= ca*po;
        
        d = min(d, dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size));
        index += 1.;
        if(i == exp) // decimal point
        {
            d = min(d, dglyph(x-2.*index*size*c.xy, 46., .7*size));
            index += 1.;
        }
    }
    
    // Output the exponent
    d = min(d, dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 101., .7*size));
    index += 1.;
    if(exp < 0.) // Sign
    {
        d = min(d, dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 45., .7*size));
        index += 1.;
        exp *= -1.;
    }
    float ca = floor(exp/10.);
    d = min(d, dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size));
    index += 1.;
    ca = floor(exp-10.*ca);
    d = min(d, dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size));
    index += 1.;
    
    return d;
}

// Distance to 210 logo
float logo(vec2 x, float r)
{
    return min(
        min(circle(x+r*c.zy, r), lineseg(x,r*c.yz, r*c.yx)),
        circlesegment(x+r*c.xy, r, -.5*pi, .5*pi)
    );
}

// Distance to stroke for any object
float stroke(float d, float w)
{
    return abs(d)-w;
}

// Distance to hexagon pattern
vec2 ind;
float hexagon( vec2 p ) 
{
    vec2 q = vec2( p.x*1.2, p.y + p.x*0.6 );
    
    vec2 pi = floor(q);
    vec2 pf = fract(q);

    float v = mod(pi.x + pi.y, 3.);

    float ca = step(1.,v);
    float cb = step(2.,v);
    vec2  ma = step(pf.xy,pf.yx);
    
    ind = pi + ca - cb*ma;
    
    return dot( ma, 1.0-pf.yx + ca*(pf.x+pf.y-1.0) + cb*(pf.yx-2.0*pf.xy) );
}

// extrusion
float zextrude(float z, float d2d, float h)
{
    vec2 w = vec2(-d2d, abs(z)-.5*h);
    return length(max(w,0.));
}

float box( vec3 x, vec3 b )
{
    return length(max(abs(x) - b,0.));
}

vec2 inset(vec3 x)
{
    float rs = 1.9;
    return vec2(-.11+min(x.y+.4, abs(length(x)-rs)), 2.);
}

// Hangar scene
vec2 scene(vec3 x) 
{
    // Start with floor (floor material: 1)
    vec2 sdf = vec2(x.y+.4/*+.01*snoise_2d(2.*x.xz-iTime)+.01*snoise_2d(4.1*x.xz-iTime*c.yx)*/, 1.);
        
    // Add glass sphere (glass material: 2)
    float rs = 1.9;
    //sdf = add(sdf, vec2(stroke(length(x)-rs,.05), 2.));

    // Add skydome
    sdf = add(sdf, vec2(abs(length(x)-2.*rs), 0.));
    
    // Add hexagonal windows to glass sphere (ceil material: 3)
    vec2 pt = vec2(atan(x.x,x.y+.4), -acos(x.z/length(x+.4*c.yxy)));
    float d = stroke(zextrude(length(x)-rs,-stroke(hexagon(vec2(5.,10.)*pt), .1),.1), .05);
    sdf = add(sdf, vec2(d, 3.));
    
    // Make some of the windows closed.
    if(rand(ind) < .5)
    {
        float d = stroke(zextrude(length(x)-rs,stroke(hexagon(vec2(5.,10.)*pt), .1),.1), .01);
        sdf = add(sdf, vec2(d, 2.));
    }
    
    // Add floor panel below windows, material: 4
    d = stroke(zextrude(x.y+.4, -stroke(length(x.xz)-rs,.2),.1),.05);
    sdf = add(sdf, vec2(d, 4.));
    
    // Add lamps
    // TODO: circle of lamps that look like cups with spheres in them
    /*
    vec3 z = x+.7*rs*c.yyx-.6*c.yxy;
    float rl = .2;
    d = length(z)-rl;
    sdf = add(sdf, vec2(d, 6.));
    */
    // Add piano
    
    // Add guard objects for debugging
    float dr = .1;
    vec3 y = mod(x,dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    sdf.x = min(sdf.x, guard);
    
    return sdf;
}

// Greetings scene
vec2 greetings(vec3 x)
{
    vec2 sdf = c.xy;
    
    return sdf;
}

// graph traversal for text effect
vec2 textpre(vec3 x)
{
    vec2 sdf =  vec2(x.z, 7.);
    float structure = stroke(logo(x.xy+.3*c.xy,.6),.25);
    float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
    if(structure < 0. && blend >= 1.e-3)
    {
        float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
        sdf = vec2(stroke(zextrude(x.z, 1.5*x.z-stroke(logo(x.xy+.3*c.xy,.6),.25), 1.*blend*clamp(1.-exp(-(x.x-34.)-8.*iTime), 0., .5)), .05*blend), 7.);
    }
    sdf.x = abs(sdf.x)-.3;
    return sdf;
}

// 3D Effect on text in intro (210 logo)
vec2 texteffect(vec3 x)
{
    // Start with z=0 plane
    vec2 sdf = vec2(x.z, 7.);
    float hex = hexagon(18.*x.xy);
    
    // compute hexagon indices in cartesian coordinates
    vec2 cind = ind/18.;
    cind = vec2(cind.x/1.2, cind.y);
    cind = vec2(cind.x, cind.y-cind.x*.6);
    
    // build up team210 logo (t < 12.)
    float structure = stroke(logo(cind+.3*c.xy,.6),.25);
    float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
    if(structure < 0. && blend >= 1.e-3)
    {
        float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
        sdf = vec2(stroke(zextrude(x.z, 2.*x.z-stroke(logo(cind.xy+.3*c.xy,.6),.25), (.5+.5*snoise_2d(24.*cind.xy-iTime))*blend*clamp(1.-exp(-(ind.x-34.)-8.*iTime), 0., 1.)), .05*blend), 7.);
    }
    
    // Add guard objects for debugging
    float dr = .03;
    vec3 y = mod(x,dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    sdf.x = min(sdf.x, guard);
    
    return sdf;
}

vec2 texteffect2(vec3 x) // text effect for endeavor text (bounce with rhythm
{
    vec2 sdf = vec2(x.z, 7.);
    float hex = hexagon(18.*x.xy);
    
    // compute hexagon indices in cartesian coordinates
    vec2 cind = ind/18.;
    cind = vec2(cind.x/1.2, cind.y);
    cind = vec2(cind.x, cind.y-cind.x*.6);
    
    // build up endeavour text
    // Show demo name: "Endeavor" (t < 25.)
    float endeavor = dstring(cind+2.*(-6.+1.2*iTime-1.2*14.)*c.xy, 0., .8);
    endeavor = stroke(endeavor, .2);
    float structure = mix(0., endeavor, clamp(.25*(iTime-14.), 0., 1.));
    float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
    if(structure < 0. && blend >= 1.e-3)
    {
        float blend = smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
        sdf = vec2(stroke(zextrude(x.z, 2.*x.z-stroke(logo(cind.xy+.3*c.xy,.6),.25), (.5+.5*snoise_2d(24.*cind.xy-iTime))*blend*clamp(1.-exp(-(ind.x-34.)-8.*iTime), 0., 1.)), .05*blend), 7.);
    }
//     float endeavor = dstring(cind+2.*(-6.+1.2*iTime-1.2*14.)*c.xy, 0., .8);
//     endeavor = stroke(endeavor, .2);
//     structure = mix(structure, endeavor, clamp(.25*(iTime-14.), 0., 1.));
    
    return sdf;
}

vec3 post1(vec2 uv, vec3 col)
{
    if(uv.y < .8) 
    {
        // scanlines
        col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
        col = clamp(col,c.yyy,c.xxx);
        return col;
    }
    
    // Preparations
    vec3 blu = vec3(.2, .68, 1.);
    float px = 1.5/iResolution.y;
    
    // 210 logo
    float dt0 = logo(uv-2.*vec2(-.45*a,.45),.04);
    dt0 = stroke(dt0, .01);
    col = mix(col, mix(col, blu, .5), smoothstep(px, -px ,dt0));
    dt0 = stroke(dt0, .0025);
    col = mix(col, blu, smoothstep(px, -px ,dt0));
    
    // bounding box for time display
    dt0 = stroke(lineseg(uv-2.*vec2(-.45*a, .45)-.2*c.xy, c.yy, 1.*c.xy), .05);
    col = mix(col, mix(col, blu, .5), smoothstep(px, -px, dt0));
    float dt1 = stroke(dt0, .0025);
    col = mix(col, blu, smoothstep(px, -px, dt1));
    
    // "elapsed:" text with time display
    float dta = dstring(uv-2.*vec2(-.45*a,.45)-.3*c.xy,1., .025);
    dta = min(dta, dfloat(uv-2.*vec2(-.45*a,.45)-.7*c.xy, iTime, .025));
    dta = stroke(dta, .0025);
    col = mix(col, clamp(1.*blu, 0., 1.), smoothstep(px, -px, dta));
    
    // bounding box for size display
    dt0 = stroke(lineseg(uv-2.*vec2(-.45*a, .45)-1.4*c.xy, c.yy, 1.*c.xy), .05);
    col = mix(col, mix(col, blu, .5), smoothstep(px, -px, dt0));
    dt1 = stroke(dt0, .0025);
    col = mix(col, blu, smoothstep(px, -px, dt1));
    
    // "size:" text with size
    dta = dstring(uv-2.*vec2(-.45*a,.45)-1.5*c.xy,2., .025);
    dta = min(dta, dfloat(uv-2.*vec2(-.45*a,.45)-1.7*c.xy, iExecutableSize, .025));
    dta = stroke(dta, .0025);
    col = mix(col, clamp(1.*blu, 0., 1.), smoothstep(px, -px, dta));
    
    // scanlines
    col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);
    col = clamp(col,c.yyy,c.xxx);
    
    return col;
}

//performs raymarching
//scene: name of the scene function
//xc: 	 name of the coordinate variable
//ro:	 name of the ray origin variable
//d:	 name of the distance variable
//dir:	 name of the direction variable
//s:	 name of the scenestruct variable
//N:	 number of iterations used
//eps:	 exit criterion
//flag:  name of the flag to set if raymarching succeeded
#define raymarch(scene, xc, ro, d, dir, s, N, eps, flag) \
    flag = false;\
    for(int ia=0; ia<N; ++ia)\
    {\
        xc = ro + d*dir;\
        s = scene(xc);\
        if(s.x < eps)\
        {\
            flag = true;\
            break;\
        }\
        d += s.x;\
    }

//computes normal with finite differences
//scene: name of the scene function
//n:	 name of the normal variable
//eps:	 precision of the computation
//xc:	 location of normal evaluation
#define calcnormal(scene, n, eps, xc) \
    {\
        float ss = scene(xc).x;\
        n = normalize(vec3(scene(xc+eps*c.xyy).x-ss,\
                        scene(xc+eps*c.yxy).x-ss,\
                        scene(xc+eps*c.yyx).x-ss));\
    }

//camera setup
//camera: camera function with camera(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
//ro:	  name of the ray origin variable
//r:	  name of the right variable
//u:	  name of the up variable
//t:	  name of the target variable
//uv:	  fragment coordinate
//dir:	  name of the dir variable
#define camerasetup(camera, ro, r, u, t, uv, dir) \
    {\
        camera(ro, r, u, t);\
        t += uv.x*r+uv.y*u;\
        dir = normalize(t-ro);\
    }

        //uv += .02*vec2(snoise_2d(uv-iTime+2.),snoise_2d(uv-iTime+3.));\
//post processing: 210 logo and trendy display lines
//col: output color
//uv:  fragment coordinate
#define post(color, uv) \
    {\
        color = post1(uv, color);\
    }
    
//camera for scene 1
void camera1(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
{
    ro = c.yyx-1.*c.yyx+.1*(iTime-28.)*c.yyx;
    r = c.xyy;
    u = c.yxy;
    t = c.yyy-1.*c.yyx+.1*(iTime-28.)*c.yyx;
}

// static camera
void camera0(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
{
    float blend = 0.;// smoothstep(2., 6., iTime)*(1.-smoothstep(6.,12.,iTime));
    ro = c.yyx-.5*c.yxy*blend;
    r = c.xyy;
    u = c.yxy+.5*c.yyx*blend;
    t = .1*c.yyx;
}

vec3 stdcolor(vec2 x)
{
    return 0.5 + 0.5*cos(iTime+x.xyx+vec3(0,2,4));
}

float star(vec2 x, float r0)
{
    return 1.-smoothstep(.5*r0, r0, length(x));
}

vec3 background(vec2 x)
{
    // Add sky gradient
    vec3 col = mix(vec3(1., .56, .44), vec3(1.,1.,.87), abs(x.y));
    
    // Add sun
    float d = length(vec2(x.x, abs(x.y+.15))-.3*c.yx)-.15;
    col = mix(col, c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
    
    // Add clouds
    float da = .5*snoise_2d(2.*vec2(x.x, abs(x.y+.15))-.4*iTime);
    float dx = da+mfsnoise_2d(vec2(x.x-.2*iTime, abs(x.y+.15)), 1.e0, 1.e3, .55);
    col = mix(col, .8*vec3(1.,.7,.57), clamp(2.5+dx, 0., 1.));
    col = mix(col, .9*vec3(1.,.7,.57), clamp(2.3+dx, 0., 1.));
    col = mix(col, vec3(1.,.7,.57), clamp(2.1+dx, 0., 1.));
    
    // And more clouds
    da = .5*snoise_2d(2.*vec2(x.x, abs(x.y+.15))-.4*iTime-15.);
    dx = da+mfsnoise_2d(vec2(x.x-.1*iTime-15., abs(x.y+.15)), 1.e0, 1.e3, .55);
    col = mix(col, .8*vec3(1.,1.,.87), clamp(2.5+dx, 0., 1.));
    col = mix(col, .9*vec3(1.,1.,.87), clamp(2.3+dx, 0., 1.));
    col = mix(col, vec3(1.,1.,.87), clamp(1.6+dx, 0., 1.));
    
    col = mix(col, c.xxx, .9*smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
    
    
    return col;
}

// Initial intro
vec3 background2(vec2 uv)
{
    // hexagonal grid
    float d = stroke(-hexagon(18.*uv), .1);
    
    // compute hexagon indices in cartesian coordinates
    vec2 cind = ind/18.;
    cind = vec2(cind.x/1.2, cind.y);
    cind = vec2(cind.x, cind.y-cind.x*.6);
    
    // build up team210 logo (t < 12.)
    float structure = exp(-(ind.x-34.)-8.*iTime)+stroke(logo(cind+.3*c.xy,.6),.25);
    
    // blend to gray (t < 16.)
    structure = mix(structure, hexagon(18.*uv), clamp(.25*(iTime-12.), 0., 1.));
    
    vec2 dind = cind;
    
    // Show demo name: "Endeavor" (t < 25.)
    float endeavor = dstring(cind+2.*(-6.+1.2*iTime-1.2*14.)*c.xy, 0., .8);
    endeavor = stroke(endeavor, .2);
    structure = mix(structure, endeavor, clamp(.25*(iTime-14.), 0., 1.));
    
    // blend hexagons smaller (t < 27.)
    d = mix(d, stroke(-hexagon(8.*uv), .1), clamp(.5*(iTime-24.), 0., 1.));
    structure = mix(structure, hexagon(8.*uv), clamp(.5*(iTime-24.), 0., 1.));
    dind = ind/8.;
    dind = vec2(dind.x/1.2, dind.y);
    dind = vec2(dind.x, dind.y-dind.x*.6);
    cind = mix(cind, dind,  clamp(.5*(iTime-24.), 0., 1.));
    
    // make background change the color with time
    vec2 dt = vec2(snoise_2d(5.*cind+2.),snoise_2d(5.*cind+3.));
    float m = (.5+.5*snoise_2d(50.*cind)
        + mix(-1.,clamp(.5+.5*snoise_2d(.5*(cind)-dt-2.*iTime*c.xx),0.,1.), clamp(.125*(iTime-1.),0.,1.)));
    vec3 c1 = mix(c.yyy, c.yyy,m)*smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d);
        c1 = mix(c1, mix(c.yyy, vec3(1.,0.27,0.),m), smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, stroke(structure,.05)))*smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d);
    c1 = clamp(c1, 0., 1.);
    
    // grayscale everything outside the structure
    if(structure > 0.)
        c1 = mix(.3*length(c1)*c.xxx/sqrt(3.), c1, clamp(.5*(iTime-24.), 0., 1.));
    
    // blend to black at the end
    c1 = mix(c1, c.yyy, clamp(iTime-27., 0., 1.));
    
    return clamp(c1,0.,1.);
}

vec3 color(float rev, float ln, float mat, vec2 uv, vec3 x)
{
    if(mat == 7.)
        return background2(x.xy);
    if(mat == 6.)
        return clamp(.7*c.xxx + .7*c.xxy*(ln) + c.xxx*abs(pow(rev,8.)), 0., 1.);
    if(mat == 2.)
    {
        vec3 col = .1*c.xyy + .3*c.xyy * abs(ln) + .8*c.xxy * abs(pow(rev,8.));
        vec2 pt = vec2(atan(x.x,x.y+.4), -acos(x.z/length(x+.4*c.yxy)));
        float d = stroke(-hexagon(6.*vec2(5.,10.)*pt), .1);
        float m = rand(ind/*+floor(2.*iTime)*/);
        col = mix(col, mix(col, vec3(1.,0.27,0.),m), smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d));
        return col;
    }
    return .1*c.xyy + .3*c.xyy * abs(ln) + .8*c.xxy * abs(pow(rev,8.));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
//     iScale = scale(iTime+.1); //TODO: ADD THIS! IT SYNCS
    a = iResolution.x/iResolution.y;
    vec3 ro, r, u, t, x, dir;
    vec2 s, uv;
    
    float d = 0.;
    bool hit;
    
    vec3 col = c.yyy;
    
    // Antialiasing
#if FSAA!=1
    for(int jii=0; jii<FSAA; ++jii)
    for(int jij=0; jij<FSAA; ++jij)
    {
    vec2 o = vec2(float(jii),float(jij)) / float(FSAA) - .5;
    uv = (-iResolution.xy + 2.*(fragCoord+o))/iResolution.y;
#else 
    uv = (-iResolution.xy + 2.*fragCoord)/iResolution.y;
#endif
//     if(iTime < 1000.)
//     {
// //         vec3 c1 = texture(iSequence, uv).rgb;
//         vec3 c1 = c.yyy;
// //         float st = scale(iTime);
//         iScale = scale(iTime-uv.x);
//         if(uv.y > 0.)  
//             c1 += step(uv.y, iNBeats/16.);
//         c1 = mix(c1, c.xxy, step(abs(uv.x), .01));
//         col += c1;
//     }
//     if(iTime < 1000.)
//     {
//         float d = dglyph(uv, 57., .1);
//         //float d = dstring(uv-.1, 1., .05);
//         float d = dfloat(uv, rfloats(3.), .05);
//         if(d == 1.)col += c.yxy;
//         else
//         {
//             d = stroke(d, .01);
//             col +=  mix(c.yyy, c.xyy, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d));
//         }
//     }
//     else
    if(iTime < 28.) // "Enter the Logic Farm" logo/title, t < 31.
    {
        vec3 c1 = c.yyy;
        
        camerasetup(camera0, ro, r, u, t, uv, dir);
//         d = (.8-ro.z)/dir.z;
        raymarch(textpre, x, ro, d, dir, s, 100, 2.e-5, hit);
        if(hit) hit = false;
        else d = -ro.z/dir.z;
        raymarch(texteffect, x, ro, d, dir, s, 200, 2.e-5, hit);
        
        if(hit)
        {
            vec3 n;
            calcnormal(scene, n, 2.e-4, x);

            float rs = 1.9;
            vec3 l = x+1.*c.yyx,
                //l = -1.*c.yxy+1.5*c.yyx, 
                re = normalize(reflect(-l,n)), v = normalize(x-ro);
            float rev = (dot(re,v)), ln = (dot(l,n));

            c1 = color(rev, ln, s.y, uv, x);
        }
        else c1 = background2((ro-ro.z/dir.z*dir).xy);
        
        col += c1;
    }
    else if(iTime < 10000.)
    {
        vec3 c1 = c.yyy;
        camerasetup(camera1, ro, r, u, t, uv, dir);
        d = 0.;
        raymarch(inset, x, ro, d, dir, s, 40, 1.e-4, hit);
        if(hit) hit = false;
        raymarch(scene, x, ro, d, dir, s, 300, 1.e-4, hit);
        
        if(hit)
        {
            vec3 n;
            calcnormal(scene, n, 2.e-4, x);

            float rs = 1.9;
            vec3 //l = .7*rs*c.yyx-.6*c.yxy,
                l = -1.*c.yxy+1.5*c.yyx, 
                re = normalize(reflect(-l,n)), v = normalize(x-ro);
            float rev = (dot(re,v)), ln = (dot(l,n));

            c1 = color(rev, ln, s.y, uv, x);

//             // Soft shadows
//             /*
//             vec3 ddir = normalize(x-l), xx;
//             vec2 ss;
//             float dd;
//             raymarch(sscene, xx, l, dd, ddir, ss, 450, 1.e-4, hit);
//             if(dd<1.e-4)
//                 col = mix(col, c.yyy, .5);
//             */
// 
//             // Reflections
            if(s.y == 1.)
            {
                for(float k = .7; k >= .7; k -= .1)
                {
                    dir = (reflect(dir, n));
                    d = 2.e-4;
                    ro = x;
                    raymarch(inset, x, ro, d, dir, s, 20, 1.e-4, hit);
                    raymarch(scene, x, ro, d, dir, s, 300, 1.e-4, hit);
                    if(hit)
                    {
                        vec3 n2;
                        calcnormal(scene, n2, 2.e-4, x);
//                         l = -1.*c.yxy+1.5*c.yyx;
                        re = normalize(reflect(-l,n)); 
                        v = normalize(x-ro);
                        rev = abs(dot(re,v));
                        ln = abs(dot(l,n));

                        c1 = mix(c1, color(rev, ln, s.y, uv, x), k);
                    }
                    else c1 = mix(c1,background(uv), k);
                }
            }
        }
        else c1 = background(uv);

        // lens flare
        for(float k=0.; k<8.; k+=1.)
        {
            vec2 dx = .15*vec2(-1.+2.*rand(c.xx+k), -1.+2.*rand(c.xx+k+1.));
            vec3 cx = c.xxx-.2*vec3(rand(c.xx+k+2.), rand(c.xx+k+3.), rand(c.xx+k+4.));
            float sx = .05+.05*rand(c.xx+k+5.);
            float da = dpoly_min(uv-.15*c.yx+dx, 6., sx);
            c1 = mix(c1, mix(c1,cx, .5), smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, da));
        }
        
        // fog
        //c1 = mix(c1, mix(vec3(1., .56, .44), vec3(1.,1.,.87), abs(uv.y)), tanh(.*d));
        
        col += mix(mix(c.yyy, c1, clamp(iTime-29., 0., 1.)), c.yyy, clamp(iTime-44., 0., 1.));
    }
    
    
        
#if FSAA!=1
    }
    col /= float(FSAA*FSAA);
#else
#endif
    
    // Post-process
    post(col, uv);
    
    // Set the fragment color
    fragColor = vec4(col, 1.);    
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
