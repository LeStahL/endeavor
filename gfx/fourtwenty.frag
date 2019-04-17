/*
 * Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#version 130

uniform float iTime, iProgress;
uniform vec2 iResolution;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

// Hash function
void rand(in vec2 x, out float num)
{
    x += 400.;
    num = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
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

// compute distance to regular star
void dstar(in vec2 x, in float N, in vec2 R, out float dst)
{
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
    dst = dot(x-p1,ff)/length(ff);
}

// 2D box
void dbox(in vec2 x, in vec2 b, out float d)
{
	vec2 da = abs(x)-b;
	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

// Distance to circle
void dcircle(in vec2 x, in float R, out float d)
{
    d = length(x)-R;
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Distance to line segment
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

void colorize(in vec3 old, in float ind, in vec2 uv, out vec3 col)
{
    col = old;
    float d=1., da, R = .5, dd = 1.;
    vec2 p0 = -.3*c.yx, dis;
    
    lfnoise(2.*uv.x*c.xx, dis.x);
    lfnoise(2.*uv.y*c.xx, dis.y);
    
    // Draw Leaves
    for(float phi=-pi/2.+pi/5.; phi<3.*pi/2.; phi += pi/5.)
    {
        vec2 p = R*vec2(cos(phi), sin(phi));
        dlinesegment(uv-.1*dis,p0,p,da);
        d = min(d,da);
        
        vec2 dad = p-p0;
        float dpd = dot(uv-p0, dad)/dot(dad,dad),
        	arg = 2.*pi*dpd*36.*length(dad),
            q = .3;
        stroke(da, (.08*length(dad)+.01-.005*sin(arg)/(1.+q*q-2.*q*cos(arg)))*smoothstep(0.,.5,dpd)*(1.-smoothstep(.5,1.,dpd)), da);
        dd = min(dd, da);
    }
    
    float scale;
    rand(ind*c.xx, scale);
    
    col = mix(col, (.5+scale)*vec3(0.24,0.54,0.01), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, dd)); 
    stroke(dd,.002, dd);
    col = mix(col, (1.-scale)*vec3(0.16,0.40,0.05), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, dd)); 

    stroke(d,.002, d);
    col = mix(col, (1.-scale)*vec3(0.16,0.40,0.05), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d)); 
    
    col = clamp(col, 0.,1.);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = c.yyy;
    
    for(int i=0; i<80; ++i)
    {
        float phi;
        lfnoise(float(i)*c.xx-iTime, phi);
        mat2 RR = mat2(cos(phi),-sin(phi),sin(phi), cos(phi));
        float scale;
        rand(float(i)*c.xx+.1, scale);
        vec2 delta;
        rand(float(i)*c.xx+.3, delta.x);
        delta.x *= a;
        rand(float(i)*c.xx+.5, delta.y);
        delta -= .5*vec2(a,1.);
        colorize(col, float(i), (1.+3.*scale)*RR*(uv-delta), col);
    }
    col = clamp(col, 0.,1.);
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
