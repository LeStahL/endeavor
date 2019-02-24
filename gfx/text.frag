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

uniform float iFontWidth, iExecutableSize;
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

// Get glyph data from texture
void dglyph(in vec2 x, in float ordinal, in float size, out float dst)
{
    float dis;
    box(x, 2.*size*c.xx), dis;
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
    rfloat(offset, nlines)
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
        //TODO: continue here
        d = min(d, lineseg(x, size*vec2(x1,y1), size*vec2(x2, y2)));
    }
    
    // Circles
    float ncircles = floor(rfloat(offset));
    offset += 1.;
    for(float i=0.; i<max(ncircles,0); i+=1.)
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
    for(float i=0.; i<max(nsegments,0); i+=1.)
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
    {
        return 1.;
    }
    
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
    {
        return bound+.5*size;
    }
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec4 col = vec4(0.);
    float bound = sqrt(iFSAA)-1.;
   	for(float i = -.5*bound; i<=.5*bound; i+=1.)
        for(float j=-.5*bound; j<=.5*bound; j+=1.)
     		col += texture(iChannel0, fragCoord/iResolution.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy);
    col /= iFSAA;
    fragColor = col;
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
