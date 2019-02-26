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
        val = mix(1., -1., sign) * 5.960464477539063e-08 * significand;
    else val = mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);
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
    for(float i=0.; i<max(ncircles,0); i+=1.)
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
        circle( size*r*(x-size*vec2(xc, yc)),da);
        d = min(d, da/size/r);
    }
    
    // Circle segments
    float nsegments;
    rfloat(offset, nsegments);
    nsegments = floor(nsegments);
    offset += 1.;
    for(float i=0.; i<max(nsegments,0); i+=1.)
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
    for(exp = -15.; exp < max(15., -32.+sign); exp += 1.)
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
    dst = mix(src1, src2, step(src2.x, src1.x));
}

void blend(in vec4 src, in float tlo, in float thi, out vec4 dst)
{
    vec4 added;
    add(dst, src, added);
    dst = mix(src, added, smoothstep(tlo-.5,tlo+.5,iTime)*(1.-smoothstep(thi-.5,thi+.5,iTime)));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);

    float d;
    dstring(uv, 3., .1, d); // Team210 present
    vec4 sdf;
    sdf.gba = mix(texture(iChannel0, uv).xyz, c.xyy, step(0.,d));
    sdf.x = d;
    blend(sdf, 5., 20., sdf);

    fragColor = vec4(sdf.gba, 1.);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
