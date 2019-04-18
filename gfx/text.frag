/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
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

uniform float iFontWidth, iExecutableSize, iTime;
uniform vec2 iResolution;
uniform sampler2D iChannel0, iFont;

out vec4 gl_FragColor;

// Global constants
const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);
float a; // Aspect ratio

// Hash function
void rand(in vec2 x, out float num)
{
    num = fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);
}

// Arbitrary-frequency 2D noise
void lfnoise(in vec2 t, out float num)
{
    vec2 i = floor(t);
    t = fract(t);
    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise
    t = smoothstep(c.yy, c.xx, t); // TODO: add this for faster value noise
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    num = mix(v1.x, v1.y, t.x);
}

// Read short value from texture at index off
void rshort(in float off, out float val)
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
    val = round(dot(vec2(255., 65280.), data));
}

// Read float value from texture at index off
void rfloat(in float off, out float val)
{
    // Convert the bytes to unsigned short as first step.
    float d;
    rshort(off, d);
    
    // Convert bytes to IEEE 754 float16. That is
    // 1 sign bit, 5 bit exponent, 11 bit mantissa.
    // Also it has a weird conversion rule that is not evident at all.
    float sign = floor(d/32768.),
        exponent = floor(d/1024.-sign*32.),
        significand = d-sign*32768.-exponent*1024.;

    // Return full float16
    if(exponent == 0.)
    {
        val = mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    }
    else
    {
        val = mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
    }
}


// 2D box
void box(in vec2 x, in vec2 b, out float dst)
{
    vec2 d = abs(x) - b;
    dst = length(max(d,c.yy)) + min(max(d.x,d.y),0.);
}

// Distance to circle
void circle(in vec2 x, out float d)
{
    d = abs(length(x)-1.0);
}

// Distance to line segment
void lineseg(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// 2D rhomboid
void rhomboid(in vec2 x, in vec2 b, in float tilt, out float dst)
{
    x.x -= tilt/2./b.y*x.y;
    box(x,b,dst);
}

// Distance to hexagon pattern
void dhexagonpattern(in vec2 p, out float d, out vec2 ind) 
{
    vec2 q = vec2( p.x*1.2, p.y + p.x*0.6 );
    
    vec2 pi = floor(q);
    vec2 pf = fract(q);

    float v = mod(pi.x + pi.y, 3.0);

    float ca = step(1.,v);
    float cb = step(2.,v);
    vec2  ma = step(pf.xy,pf.yx);
    
    d = dot( ma, 1.0-pf.yx + ca*(pf.x+pf.y-1.0) + cb*(pf.yx-2.0*pf.xy) );
    ind = pi + ca - cb*ma;
    ind = vec2(ind.x/1.2, ind.y);
    ind = vec2(ind.x, ind.y-ind.x*.6);
}

// Distance to circle segment
void circlesegment(in vec2 x, in float r, in float p0, in float p1, out float d)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))
    {
        d = abs(length(x)-r);
        return;
    }
    d = min(
        length(x-r*vec2(cos(p0), sin(p0))),
        length(x-r*vec2(cos(p1), sin(p1)))
        );
}

// Compute distance to stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0) - s;
}

// Get glyph data from texture
void dglyph(in vec2 x, in float ordinal, in float size, out float dst)
{
    float dis;
    box(x, 2.*size*c.xx, dis);
    if(dis > 0.)
    {
        dst = dis+.5*size;
        return;
    }

    // Find glyph offset in glyph index
    float nglyphs, offset = 0;
    rfloat(1., nglyphs);
        
    for(float i=0.; i<nglyphs; i+=1.)
    {
        float ord;
        rfloat(2.+2.*i, ord);
        ord = floor(ord);
        
        if(ord == ordinal)
        {
            rfloat(2.+2.*i+1., offset);
            offset = floor(offset);
            break;
        }
    }
    
    if(offset == 0.) 
    {
        dst = 1.;
        return;
    }
    
    // Get distance from glyph data
    float d = 1.;
    
    // Lines
    float nlines;
    rfloat(offset, nlines);
    nlines = floor(nlines);
    offset += 1.;
    for(float i=0.; i<nlines; i+=1.)
    {
        float x1;
        rfloat(offset, x1);
        offset += 1.;
        float y1;
        rfloat(offset, y1);
        offset += 1.;
        float x2;
        rfloat(offset, x2);
        offset += 1.;
        float y2;
        rfloat(offset, y2);
        offset += 1.;
        float da;
        lineseg(x, size*vec2(x1,y1), size*vec2(x2, y2), da);
        d = min(d,da);
    }
    
    // Circles
    float ncircles;
    rfloat(offset, ncircles);
    ncircles = floor(ncircles);
    offset += 1.;
    for(float i=0.; i<ncircles; i+=1.)
    {
        float xc;
        rfloat(offset, xc);
        offset += 1.;
        float yc;
        rfloat(offset, yc);
        offset += 1.;
        float r;
        rfloat(offset, r);
        offset += 1.;
        float da;
        circle( (x-size*vec2(xc, yc))/size/r,da);
        d = min(d, da*size*r);
    }
    
    // Circle segments
    float nsegments;
    rfloat(offset, nsegments);
    nsegments = floor(nsegments);
    offset += 1.;
    for(float i=0.; i<nsegments; i+=1.)
    {
        float xc;
        rfloat(offset, xc);
        offset += 1.;
        float yc;
        rfloat(offset, yc);
        offset += 1.;
        float r;
        rfloat(offset, r);
        offset += 1.;
        float phi0;
        rfloat(offset, phi0);
        offset += 1.;
        float phi1;
        rfloat(offset, phi1);
        offset += 1.;
        float da;
        circlesegment(x-size*vec2(xc,yc), size*r, phi0, phi1, da);
        d = min(d, da);
    }
    
    if(nlines+ncircles+nsegments == 0.)
        dst = dis;
    else dst = d;
}

// Get distance to string from database
void dstring(in vec2 x, in float ordinal, in float size, out float dst)
{
    // Get string database offset
    float stroff0;
    rfloat(0., stroff0);
    stroff0 = floor(stroff0);
    
    // Return 1 if wrong ordinal is supplied
    float nstrings;
    rfloat(stroff0, nstrings);
    nstrings = floor(nstrings);
    if(ordinal >= nstrings)
    {
        dst = 1.;
        return;
    }
    
    // Get offset and length of string from string database index
    float stroff;
    rfloat(stroff0+1.+2.*ordinal, stroff);
    stroff = floor(stroff);
    float len;
    rfloat(stroff0+2.+2.*ordinal, len);
    len = floor(len);
    
    // Draw glyphs
    vec2 dx = mod(x-size, 2.*size)-size, 
        ind = ceil((x-dx+size)/2./size);
    
    // Bounding box
    float bound;
    box(x-size*(len-3.)*c.xy, vec2(size*len, 1.*size), bound);
    if(bound > 0.)
    {
        dst = bound+.5*size;
        return;
    }
    
    float da;
    rfloat(stroff+ind.x, da);
    da = floor(da);
    dglyph(dx, da, .7*size, dst);
}

// distance to a floating point number string
// for debugging stuff while shader is loaded
void dfloat(in vec2 x, in float num, in float size, out float dst)
{
    float d = 1., index = 0.;
    
    // Determine sign and output it if present
    float sign = sign(num), exp = 0.;
    if(sign<0.)
    {
        float da;
        dglyph(x, 45., .7*size, da);
        d = min(d, da);
        index += 1.;
        num *= -1.;
    }
    
    // The first power of ten that floors num to anything not zero is the exponent
    for(exp = -15.; exp < 15.; exp += 1.)
        if(floor(num*pow(10.,exp)) != 0.)
            break;
    exp *= -1.;
    // Determine the significand and output it
    for(float i = exp; i >= max(exp-5.,-33); i -= 1.)
    {
        float po = pow(10.,i);
        float ca = floor(num/po);
        num -= ca*po;
        
        float da;
        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);
        d = min(d, da);
        index += 1.;
        if(i == exp) // decimal point
        {
            dglyph(x-2.*index*size*c.xy, 46., .7*size, da);
            d = min(d, da);
            index += 1.;
        }
    }
    
    // Output the exponent
    float db;
    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 101., .7*size, db);
    d = min(d, db);
    index += 1.;
    if(exp < 0.) // Sign
    {
        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 45., .7*size,db);
        d = min(d, db);
        index += 1.;
        exp *= -1.;
    }
    float ca = floor(exp/10.);
    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);
    d = min(d, db);
    index += 1.;
    ca = floor(exp-10.*ca);
    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);
    d = min(d, db);
    index += 1.;
    
    dst = d;
}

// Add scene contents
void add(in vec4 src1, in vec4 src2, out vec4 dst)
{
    dst = mix(src1, src2, smoothstep(0., 1.5/iResolution.y, -src2.x));
}

void blendadd(in vec4 src1, in vec4 src2, in float tlo, in float thi, out vec4 dst)
{
    vec4 added;
    add(src1, src2, added);
    dst = mix(src1, added, smoothstep(tlo-.5,tlo+.5,iTime)*(1.-smoothstep(thi-.5,thi+.5,iTime)));
}

// UI Window Control
void window(in vec2 x, in vec2 size, in vec3 bg, in float title_index, out vec4 col)
{
    size.x *= .5;
    col = vec4(1., bg);
    
    const float cellsize = .015, bordersize = .005;
    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),
        bordercolor = vec3(1.00,0.71,0.02);
    vec4 c2 = vec4(1., titlecolor);
    
    float dhx, dhy;
    vec2 ind;
    dhexagonpattern(72.*x,  dhx, ind);
    stroke(dhx, .1, dhx);
    lfnoise(ind-iTime, dhy);
    
    // Window background
    box(x+.5*size*c.yx,size*vec2(1.,.5),c2.x);
    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/size.y), .5+.5*dhy*step(0.,dhx));
    add(col, c2, col);
    
    // Title bar
    c2.gba = titlecolor;
    rhomboid(x+.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    rhomboid(x, vec2(.65*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    rhomboid(x-.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    
    // Border of title bar
    c2 = vec4(1., bordercolor);
    stroke(col.x,bordersize,c2.x);
    add(col,c2,col);
    
    // Window Border
    lineseg(x, -.9*size.x*c.xy, -size.x*c.xy, c2.x);
    float d;
    lineseg(x, -size.x*c.xy, -size, d);
    c2.x = min(c2.x, d);
    lineseg(x, -size, size*c.xz, d);
    c2.x = min(c2.x, d);
    lineseg(x, size*c.xz, size*c.xy, d);
    c2.x = min(c2.x, d);
    lineseg(x, .9*size.x*c.xy, size.x*c.xy, d);
    c2.x = min(c2.x, d);
    stroke(c2.x,.25*bordersize,c2.x);
    add(col, c2, col);
}

void progressbar(in vec2 x, in float width, in float progress, out vec4 col)
{
    const float cellsize = .015, bordersize = .005;
    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),
        bordercolor = vec3(1.00,0.71,0.02), bg = c.yyy;
    vec4 c2 = vec4(1., titlecolor);
    
    // Window background
    box(x+.5*width*c.yx,width*c.xy,c2.x);
    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/cellsize), .5);
    add(col, c2, col);
    
    // Bar background
    c2.gba = titlecolor;
    rhomboid(x, vec2(.5*width,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    
    // Border
    c2.gba = bordercolor;
    stroke(c2.x,.5*bordersize,c2.x);
    add(col, c2, col);
    
    // Progress
    float wc = width/cellsize;
    x.x -= .5*x.y;
    vec2 y = vec2(mod(x.x, 1.2*cellsize)-.6*cellsize, x.y),
        index = (x-y)/.6/cellsize;
    if(abs(index.x) < .8*wc && -index.x > .8*wc*(1.-2.*progress))
    {
        box(y, vec2(.5*cellsize, .8*cellsize), c2.x);
        add(col, c2, col);
    }
}

// Revision logo of width 1.
void drevision(in vec2 x, in float r, out float dst)
{
    float l = length(x),
        p = atan(x.y,x.x),
	    d = abs(l-r*.07)-.02, 
        k1 = abs(l-r*.16)-.03,
        k2 = abs(l-r*.21)-.02, 
        k3 = abs(l-r*.35)-.03,
        k4 = abs(l-r*.45)-.02;
    vec4 n1;
    lfnoise(1.*c.xx-3.*iTime, n1.x);
    lfnoise(2.*c.xx-2.4*iTime, n1.y);
    lfnoise(3.*c.xx-2.9*iTime, n1.z);
    lfnoise(4.*c.xx-3.1*iTime, n1.w);
    n1 = mix(n1,c.yyyy, clamp((iTime-24.)/2.,0.,1.));
    d = min(d, mix(d, abs(l-.11)-.03, step(p, -1.71)*step(-2.73, p)));
    d = min(d, mix(d, k1, step(p+n1.x, 3.08)*step(2.82,p)));
    d = min(d, mix(d, k1, step(p+n1.x, 1.47)*step(.81,p)));
    d = min(d, mix(d, k1, step(p+n1.x, -.43)*step(-1.19,p)));
    d = min(d, mix(d, k2, step(p+n1.y, -2.88)*step(-pi,p)));
    d = min(d, mix(d, k2, step(p+n1.y, pi)*step(2.38,p)));
    d = min(d, mix(d, k2, step(p+n1.y, 2.1)*step(.51,p)));
    d = min(d, mix(d, k2, step(p+n1.y, .3)*step(-1.6,p)));
    d = min(d, abs(l-.24)-.02);
    d = min(d, mix(d, k3, step(p+n1.z, -2.18)*step(-pi, p)));
    d = min(d, mix(d, k3, step(p+n1.z, -1.23)*step(-1.7, p)));
    d = min(d, mix(d, k3, step(p+n1.z, -.58)*step(-.78, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 0.)*step(-.29, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 1.25)*step(1.06, p)));
    d = min(d, mix(d, k3, step(p+n1.z, 1.99)*step(.5*pi, p)));
    d = min(d, abs(l-.41)-.03);
    d = min(d, mix(d, k4, step(p+n1.w, 1.04)*step(.04, p)));
    d = min(d, mix(d, k4, step(p+n1.w, -2.2)*step(-2.34, p)));
    
    dst = d-.005;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    
    float d;
    
    vec4 old = vec4(-1.,texture(iChannel0, fragCoord/iResolution.xy).rgb), new = old; // Scene
    
    // Display time
    vec4 b = vec4(1., vec3(0.99,0.64,0.02)), bd = vec4(1., .5*vec3(0.99,0.64,0.02));
    box(uv-vec2(-.48,.45)-.03*sin(iTime)*c.xy, vec2(.2,.02), b.x);    
    stroke(b.x, .001, bd.x);
    add(b, bd, b);
    box(uv-vec2(-.08,.45)-.03*sin(iTime)*c.xy, vec2(.2,.02), bd.x);
    bd.gba = vec3(0.60,0.06,0.00);
    add(b, bd, b);
    stroke(bd.x, .001, bd.x);
    add(b, bd, b);
    dfloat(uv-vec2(-.63,.45)-.03*sin(iTime)*c.xy, iTime, .018, bd.x);
    stroke(bd.x, .004, bd.x);
    add(b, bd, b);
    //dfloat(uv-vec2(-.23,.45)-.03*sin(iTime)*c.xy, iExecutableSize, .018, bd.x);
    dstring(uv-vec2(-.225,.45)-.03*sin(iTime)*c.xy, 6., .018, bd.x);
    stroke(bd.x, .004, bd.x);
    bd.gba = vec3(0.99,0.64,0.02);
    add(b, bd, b);
    b.gba = mix(old.gba, b.gba, .8);
    
    blendadd(old, b, 5., 999., old);
    
    if(iTime < 15.)
    {
        dstring(uv+.6*c.xy, 3., .05, d); // Team210 present
        stroke(d, .01, d);
        new = vec4(d, mix(old.gba, c.xxx, .6));
        blendadd(old,new,5.,13.,new);
        
        dstring(uv+.6*c.xy+.1*c.yx, 4., .03, d); // A production made of joy
        stroke(d, .005, d);
        old = vec4(d, mix(old.gba, c.xxx, .6));
        blendadd(new,old,7.,13.,new);
    }
    else if(iTime < 37.)
    {
        vec4 c0, c2;
        
        dstring(uv+vec2(.55,.4), 0., .09, d); // Endeavour
        stroke(d, .018, d);
        new = vec4(d, mix(old.gba, c.xxx, .6));
        blendadd(old,new,16.,24.,new);
        
        // Add revision loading logo
        float d;
        window(uv-.5*c.xy, vec2(.5,.5), new.gba, 0., c2);
        add(c0, c2, c0);
        drevision(3.*uv-3.*vec2(.5,-.25), 1., d);
        vec3 revcol = vec3(1.00,0.71,0.02);
        //mix(new.gba, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*uv.y), .5);
        add(c0, vec4(d, revcol), c0);
        blendadd(new, c0, 24.,28., new);
        
        window(uv+.4*c.xy, vec2(.6,.4), new.gba, 0., c2);
        add(c0, c2, c0);
        progressbar(uv+vec2(.4,.3), .35, .2+.2*sin(iTime), c2);
        add(c0, c2, c0);
        dstring(uv+vec2(.6,.05), 72., .014, d); // Booze barrel
        stroke(d, .003, d);
        c2 = vec4(d, vec3(1.00,0.71,0.02));
        add(c0,c2,c0);
        progressbar(uv+vec2(.4,.2), .35, .5+.5*sin(iTime), c2);
        add(c0, c2, c0);
        dstring(uv+vec2(.6,.15), 73., .014, d); // Electric Energy
        stroke(d, .003, d);
        c2 = vec4(d, vec3(1.00,0.71,0.02));
        add(c0,c2,c0);
        progressbar(uv+vec2(.4,.1), .35, .9+.1*sin(iTime), c2);
        add(c0, c2, c0);
        dstring(uv+vec2(.6,.25), 74., .014, d); // Atomic Diesel reserves
        stroke(d, .003, d);
        c2 = vec4(d, vec3(1.00,0.71,0.02));
        add(c0,c2,c0);
        blendadd(new,c0,28.,35.,new);
        mix(new,c.xyyy, clamp(iTime-36., 0., 1.));
    }
    else if(iTime < 130.)
    {
        vec4 bd,bda;
        box(uv-vec2(-.08,.25), vec2(.42,.03), bd.x);
        bd.gba = mix(old.gba, vec3(0.84,0.18,0.53),.8);
        dstring(uv-vec2(-.42,.25), 70., .020, d); // Dont forget to take
        stroke(d, .0045, d);
        bda = vec4(d, mix(new.gba, c.xxx, .9));
        add(bd,bda,bd);
        blendadd(old, bd, 118.,122., new);
        
        box(uv-vec2(.1,.15), vec2(.34,.03), bd.x);
        bd.gba = mix(new.gba, vec3(0.84,0.18,0.53),.8);
        dstring(uv-vec2(-.08,.15), 71., .020, d); // your medicine
        stroke(d, .0045, d);
        bda = vec4(d, mix(new.gba, c.xxx, .9));
        add(bd,bda,bd);
        blendadd(new, bd, 119.,122., new);
    }
    else if(iTime < 140.)
    {
        vec4 bd,bda;
        box(uv-vec2(-.08,.25)-.03*sin(iTime+.1)*c.xy, vec2(.42,.03), bd.x);
        bd.gba = mix(old.gba, vec3(0.73,0.90,0.22),.8);
        dstring(uv-vec2(-.42,.25)-.03*sin(iTime+.1)*c.xy, 67., .020, d); // Two times two-ten is
        stroke(d, .004, d);
        bda = vec4(d, mix(new.gba, vec3(0.23,0.27,0.16), .6));
        add(bd,bda,bd);
        blendadd(old, bd, 132.,138., new);
        
        box(uv-vec2(-.08,.15)-.03*sin(iTime+.4)*c.xy, vec2(.18,.03), bd.x);
        bd.gba = mix(new.gba, vec3(0.73,0.90,0.22),.8);
        dstring(uv-vec2(-.18,.15)-.03*sin(iTime+.4)*c.xy, 68., .020, d); // today is
        stroke(d, .004, d);
        bda = vec4(d, mix(new.gba, vec3(0.23,0.27,0.16), .6));
        add(bd,bda,bd);
        blendadd(new, bd, 133.,138., new);
        
        box(uv-vec2(-.08,-.05)-.03*sin(iTime+.9)*c.xy, vec2(.9,.08), bd.x);
        bd.gba = mix(new.gba,c.xxy,.8);
        dstring(uv-vec2(-.64,-.05)-.03*sin(iTime+.9)*c.xy, 69., .07, d); // four-twenty
        stroke(d, .016, d);
        bda = vec4(d, mix(new.gba, vec3(0.23,0.27,0.16), .6));
        add(bd,bda,bd);
        blendadd(new, bd, 134.,138., new);
        
        box(uv-vec2(+.14,-.17)-.03*sin(iTime+.13)*c.xy, vec2(.42,.03), bd.x);
        bd.gba = mix(new.gba, vec3(0.73,0.90,0.22),.8);
        dstring(uv-vec2(-.18,-.17)-.03*sin(iTime+.13)*c.xy, 75., .020, d); // whatever that means
        stroke(d, .004, d);
        bda = vec4(d, mix(new.gba, vec3(0.23,0.27,0.16), .6));
        add(bd,bda,bd);
        blendadd(new, bd, 135.,138., new);
    }
    else new = old;
    
    fragColor = vec4(new.gba, 1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
