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

uniform float iNBeats;
uniform float iScale;
uniform float iTime;
uniform vec2 iResolution;
uniform sampler2D iFont;
uniform float iFontWidth;

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

// Hash function
float rand(vec2 x)
{
    return fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

/* Simplex noise -
Copyright (C) 2011 by Ashima Arts (Simplex noise)
Copyright (C) 2011-2016 by Stefan Gustavson (Classic noise and others)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
vec3 taylorInvSqrt(vec3 r) 
{     
    return 1.79284291400159-0.85373472095314*r; 
}

vec3 permute(vec3 x)
{
    return mod((x*34.+1.)*x, 289.);
}

float snoise(vec2 P) 
{     
    const vec2 C = vec2 (0.211324865405187134, 0.366025403784438597);  
    vec2 i = floor(P+dot(P, C.yy)) ; 
    vec2 x0 = P-i+dot(i, C.xx) ; 
    // Other  corners 
    vec2 i1 ; 
    i1.x = step ( x0.y , x0.x ) ;  //  1.0  i f  x0 . x > x0 . y ,  e l s e  0.0 
    i1.y = 1.0 - i1.x ; 
    // x1 = x0 − i1 + 1.0 ∗ C. xx ;  x2 = x0 − 1.0 + 2.0 ∗ C. xx ; 
    vec4 x12 = x0.xyxy + vec4 ( C.xx , C.xx * 2.0 - 1.0) ; 
    x12.xy -= i1 ; 
    //  Permutations 
    i = mod( i ,  289.0) ;  // Avoid  truncation  in  polynomial  evaluation 
    vec3 p = permute ( permute ( i.y + vec3 (0.0 , i1.y ,  1.0  ) ) + i.x + vec3 (0.0 , i1.x ,  1.0  ) ) ; 
    //  Circularly  symmetric  blending  kernel
    vec3 m = max(0.5 - vec3 ( dot ( x0 , x0 ) ,  dot ( x12.xy , x12.xy ) , dot ( x12.zw , x12.zw ) ) ,  0.0) ; 
    m = m * m ; 
    m = m * m ; 
    //  Gradients  from 41  points  on a  line ,  mapped onto a diamond 
    vec3 x = fract ( p * (1.0  /  41.0) ) * 2.0 - 1.0  ; 
    vec3 gy = abs ( x ) - 0.5  ; 
    vec3 ox = floor ( x + 0.5) ;  // round (x)  i s  a GLSL 1.30  feature 
    vec3 gx = x - ox ; //  Normalise  gradients  i m p l i c i t l y  by  s c a l i n g m 
    m *= taylorInvSqrt ( gx * gx + gy * gy ) ; // Compute  f i n a l  noise  value  at P 
    vec3 g ; 
    g.x = gx.x * x0.x + gy.x * x0.y ; 
    g.yz = gx.yz * x12.xz + gy.yz * x12.yw ; 
    //  Scale  output  to  span  range  [ − 1 ,1] 
    //  ( s c a l i n g  f a c t o r  determined by  experiments ) 
    return  -1.+2.*(130.0 * dot ( m , g ) ) ; 
}
/* End of Simplex Noise */

// Multi-frequency simplex noise
float mfsnoise(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*snoise(f*x) + sum;
        a = a*phi;
    }
    
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
    return length(x)-r;
}

// Distance to circle segment
float circlesegment(vec2 x, float r, float p0, float p1)
{
    float p = mod(atan(x.y, x.x), 2.*pi);
    float pa = clamp(p, p0, p1);
    float d = 1.;
    if(pa != p)
    {
        d = min(
            length(x-r*vec2(cos(p0), sin(p0))),
            length(x-r*vec2(cos(p1), sin(p1)))
            );
    }
    else
	    d = length(x-r*vec2(cos(pa), sin(pa)));
       
    return d;
}

// Get glyph data from texture
float dglyph(vec2 x, float ordinal, float size)
{
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

// compute distance to regular polygon
float dpoly_min(vec2 x, float N, float R)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    return R-length(x)*cos(t)/cos(.5*d);
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

vec2 scene(vec3 x)
{
    // Start with floor (floor material: 1)
    vec2 sdf = vec2(x.y+.4/*+.01*snoise(2.*x.xz-iTime)+.01*snoise(4.1*x.xz-iTime*c.yx)*/, 1.);
        
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
    for(int i=0; i<N; ++i)\
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

//post processing: 210 logo and trendy display lines
//col: output color
//uv:  fragment coordinate
#define post(color, uv) \
    {\
        col = mix(clamp(col,c.yyy,c.xxx), c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, stroke(logo(uv-2.*vec2(-.45*a,.45),.04),.01)));\
        col += vec3(0., 0.05, 0.1)*sin(uv.y*1050.+ 5.*iTime);\
    }
    
//camera for scene 1
void camera1(out vec3 ro, out vec3 r, out vec3 u, out vec3 t)
{
    ro = c.yyx-1.*c.yyx+.1*(iTime-28.)*c.yyx;
    r = c.xyy;
    u = c.yxy;
    t = c.yyy-1.*c.yyx+.1*(iTime-28.)*c.yyx;
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
    float da = .5*snoise(2.*vec2(x.x, abs(x.y+.15))-.4*iTime);
    float dx = da+mfsnoise(vec2(x.x-.2*iTime, abs(x.y+.15)), 1.e0, 1.e3, .55);
    col = mix(col, .8*vec3(1.,.7,.57), clamp(2.5+dx, 0., 1.));
    col = mix(col, .9*vec3(1.,.7,.57), clamp(2.3+dx, 0., 1.));
    col = mix(col, vec3(1.,.7,.57), clamp(2.1+dx, 0., 1.));
    
    // And more clouds
    da = .5*snoise(2.*vec2(x.x, abs(x.y+.15))-.4*iTime-15.);
    dx = da+mfsnoise(vec2(x.x-.1*iTime-15., abs(x.y+.15)), 1.e0, 1.e3, .55);
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
    
    // blend hexagons smaller (t < 18.)
    d = mix(d, stroke(-hexagon(38.*uv), .1), clamp(.5*(iTime-16.), 0., 1.));
    structure = mix(structure, hexagon(38.*uv), clamp(.5*(iTime-16.), 0., 1.));
    vec2 dind = ind/38.;
    dind = vec2(dind.x/1.2, dind.y);
    dind = vec2(dind.x, dind.y-dind.x*.6);
    cind = mix(cind, dind,  clamp(.5*(iTime-16.), 0., 1.));
    
    // Show demo name: "Endeavor" (t < 25.)
    float dp = 6.*pi/8.;
    // e
    float endeavor = circlesegment(cind*c.zx-vec2(1.2, .5), .3, -dp, dp);
    endeavor = min(endeavor, lineseg(cind, vec2(-1.1,.5), vec2(-1.3,.5)));
    // n
    endeavor = min(endeavor, circlesegment(cind.yx-vec2(-.4,.5).yx, .3, -dp, dp));
    // D
    endeavor = min(endeavor, circlesegment(cind-vec2(.4,.5), .3, -dp, dp));
    endeavor = min(endeavor, lineseg(cind, vec2(.3, .45), vec2(.3,.55)));
    // e
    endeavor = min(endeavor, circlesegment(cind*c.zx-vec2(-1.3, .5), .3, -dp, dp));
    endeavor = min(endeavor, lineseg(cind, vec2(1.2, .5), vec2(1.4,.5)));
    // A
    endeavor = min(endeavor, circlesegment(cind.yx-vec2(-1.2,-.5).yx, .3, -dp, dp));
    endeavor = min(endeavor, lineseg(cind, vec2(-1.1,-.5), vec2(-1.3,-.5)));
    // v
    endeavor = min(endeavor, circlesegment((cind.yx-vec2(-.3,-.5).yx)*c.zx, .3, -dp, dp));
    // o
    endeavor = min(endeavor, circle(cind-vec2(.6,-.5), .3));
    // r
    endeavor = min(endeavor, circlesegment(cind*c.zx-vec2(-1.5, -.5), .3, 0., .66*dp));
    endeavor = min(endeavor, lineseg(cind, vec2(1.2, -.5), vec2(1.2,-.8)));
    // stroke
    endeavor = stroke(endeavor, .13);
    //endeavor *= clamp(.25*(iTime-18.), 0., 1.)*exp(-(ind.x-34.)-8.*(iTime-18.));
    structure = mix(structure, endeavor, clamp(.25*(iTime-14.), 0., 1.));
    
    // blend hexagons smaller (t < 27.)
    d = mix(d, stroke(-hexagon(8.*uv), .1), clamp(.5*(iTime-24.), 0., 1.));
    structure = mix(structure, hexagon(8.*uv), clamp(.5*(iTime-24.), 0., 1.));
    dind = ind/8.;
    dind = vec2(dind.x/1.2, dind.y);
    dind = vec2(dind.x, dind.y-dind.x*.6);
    cind = mix(cind, dind,  clamp(.5*(iTime-24.), 0., 1.));
    
    // make background change the color with time
    vec2 dt = vec2(snoise(cind+2.),snoise(cind+3.));
    float m = 1.5+.5*snoise(10.*cind)
        + mix(-2.,clamp(.5+.5*snoise(.05*(cind)-dt-iTime*c.xx),0.,1.), clamp(.125*(iTime-1.),0.,1.));
    vec3 c1 = mix(c.yyy, .8*c.yyy,m)*smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d);
        c1 = mix(c1, mix(c.yyy, vec3(1.,0.27,0.),m), smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, stroke(structure,.05)))*smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d);
    c1 = clamp(c1, 0., 1.);
    
    // grayscale everything outside the structure
    if(structure > 0.)
        c1 = mix(.7*length(c1)*c.xxx/sqrt(3.), c1, clamp(.5*(iTime-24.), 0., 1.));
    
    // blend to black at the end
    c1 = mix(c1, c.yyy, clamp(iTime-27., 0., 1.));
    
    return clamp(c1,0.,1.);
}

vec3 color(float rev, float ln, float mat, vec2 uv, vec3 x)
{
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
    a = iResolution.x/iResolution.y;
    vec3 ro, r, u, t, x, dir;
    vec2 s, uv;
    
    float d = 0.;
    bool hit;
    
    vec3 col = c.yyy;
    
    // Antialiasing
#if FSAA!=1
    for(int i=0; i<FSAA; ++i)
    for(int j=0; j<FSAA; ++j)
    {
    vec2 o = vec2(float(i),float(j)) / float(FSAA) - .5;
    uv = (-iResolution.xy + 2.*(fragCoord+o))/iResolution.y;
#else 
    uv = (-iResolution.xy + 2.*fragCoord)/iResolution.y;
    
#endif
    if(iTime < 1000.)
    {
        float d = dglyph(uv, 99., .5);
        if(d == 1.)col += c.yxy;
        else
        {
            d = stroke(d, .01);
            col +=  mix(c.yyy, c.xyy, smoothstep(-1.5/iResolution.y, 1.5/iResolution.y, d));
        }
    }
    else
    if(iTime < 28.) // "Enter the Logic Farm" logo/title, t < 31.
    {
        col += background2(uv);
    }
    else if(iTime < 10000.)
    {
        vec3 c1 = c.yyy;
    
        camerasetup(camera1, ro, r, u, t, uv, dir);
        raymarch(inset, x, ro, d, dir, s, 40, 1.e-4, hit);
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
        for(float i=0.; i<8.; i+=1.)
        {
            vec2 dx = .15*vec2(-1.+2.*rand(c.xx+i), -1.+2.*rand(c.xx+i+1.));
            vec3 cx = c.xxx-.2*vec3(rand(c.xx+i+2.), rand(c.xx+i+3.), rand(c.xx+i+4.));
            float sx = .05+.05*rand(c.xx+i+5.);
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
