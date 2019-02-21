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

uniform float iTime, iProgress;
uniform vec2 iResolution;

out vec4 gl_FragColor;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

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

// Multi-frequency 2D noise
void mfnoise(in vec2 x, in float fmin, in float fmax, in float alpha, out float num)
{
    num = 0.;
    float a = 1., nf = 0., buf;
    for(float f = fmin; f<fmax; f = f*2.)
    {
        lfnoise(f*x, buf);
        num += a*buf;
        a *= alpha;
        nf += 1.;
    }
    num *= (1.-alpha)/(1.-pow(alpha, nf));
}

// Compute distance to regular polygon
void dpolygon(in vec2 x, in float N, out float d)
{
    d = 2.0*pi/N;
    float t = mod(acos(x.x/length(x)), d)-0.5*d;
    d = 0.5-length(x)*cos(t)/cos(0.5*d);
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

// Compute distance to stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0) - s;
}

// Apply distance function to colors
void render(in vec3 c1, in vec3 c2, in float d, out vec3 col)
{
    col = mix(c1, c2, step(0.0,d));
}

// Image function
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    a = iResolution.x/iResolution.y;
    
    float d, d0, d1, rn, p, ln;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), ind;
    vec3 col, col1, col2;
    
    uv *= 2.0;
    
    // Pattern
    dhexagonpattern((30.0)*uv, d, ind);
    stroke(d, 0.2, d);
    d = min(d, d0);
    
    // Make colors
    rand(ind, rn);
    col1 = mix((.85+.15*rn)*vec3(1.0, 0., 0.), vec3(1.0, 1.0, .2), clamp(rn, 0.0, 1.0));
    col = col1;
    
    // Big hexagon
    dpolygon(2.0*ind/66.0, 6.0, d0);
    stroke(d0, 0.2, d0);
    d0 = -d0;
    
    // Fill hexagon with colorful small hexagons that scale with progress
    mfnoise(ind, .5,100., .35, ln);
    ln = .5+.5*ln;
    p = atan(ind.y, ind.x);
    col2 = .2*col;//length(col)/length(c.xxx)*c.xxx;
    col = mix(c.yyy, col2, ln);
    col = mix(col, clamp(3.*col1,0.,1.), step(p, -pi+2.*pi*iProgress));
    render(c.yyy, col, d, col);
    render(c.yyy, col, d0, col);
    
    // Big hexagon border
    stroke(d0+.1, 0.05, d0);
    render(col, mix(c.yyy, mix(col1, col2, iProgress), .75+.25*ln), -d0, col);
    stroke(d, 0.1, d);
    render(col, c.yyy, d, col);
    
    // Add lightning
    if(d0 > .1 && length(ind/66.) > .4)
    {
        float p0 = mod(p, pi/3.)-pi/6., pi = (p-p0)*3./pi;
        float d3, off;
        rand(pi*c.xx, off);
        lfnoise(.5*length(ind)*c.xx+50.*off-5.*iTime, d3);
        d3 = p0-.15*d3;
        stroke(d3, .1, d3);
        render(col, .7*col1, -d3, col);
        stroke(d3, .05, d3);
        render(col, .2*col1, -d3, col);
    }
    render(col, c.yyy, d, col);
    
    // Small hexagons
    for(float i=0.; i<30.; i+=1.)
    {
        float size;
        vec2 velocity;
        rand(i*c.xx, size);
        size = .1+.5*size;
        rand((i+1.3)*c.xy, velocity.x);
        rand((i+2.1)*c.yx, velocity.y);
        velocity = -1.+2.*velocity;
        
        vec2 location = -vec2(a,1.)+vec2(a,2.)*mod(velocity*.5*iTime, vec2(a,1.));
        
        dpolygon((uv-location)/size, 6.0, d0);
        stroke(d0, 0.01, d);
        d = -d;
        render(col, mix(col, mix(c.xyy, c.xxy, .2+.5*velocity.x), .2), d0, col);
        render(col, mix(c.xyy, c.xxy, clamp(.3+.5*velocity.x,0.,1.)), d, col);
    }
    
    // Output to screen
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
