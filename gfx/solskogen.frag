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

uniform float iTime;
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

// Distance to line segment
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Arbitrary-frequency 2D sharp noise
void lfsharp(in vec2 t, out float num)
{
    vec2 i = floor(t);
    t = fract(t);
    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise
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
    //R.x = mix(0.,R.x, (10.*x.y));
    x.y += .9*abs(x.x);
    vec2 y = x;
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
  	dst = max(y.y,-dot(x-p1,ff)/length(ff));
}

// compute distance to regular star
void dstarhalf(in vec2 x, in float N, in vec2 R, out float dst)
{
    //R.x = mix(0.,R.x, (10.*x.y));
    x.y += .9*abs(x.x);
    vec2 y = x;
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
  	dst = max(y.y,-dot(x-p1,ff)/length(ff));
    dst = max(dst, y.x);
}

void stroke(in float d, in float h, out float dst)
{
    dst = abs(d)-h;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    float da =  -clamp(uv.x,-.5*a,0.)/.5/a;
    vec3 col = c.yyy,
        ccenter =  mix(vec3(0.79,0.33,0.61), vec3(0.77,0.30,0.33),da);
    
    if(uv.y < -.25);
    else if(uv.y < .2)
	    col = mix(vec3(0.37,0.23,0.53), ccenter, clamp((uv.y+mix(.25, .45,  da))/mix(.45, .65, da),0.,1.));
    else
        col = mix(ccenter, vec3(0.36,0.23,0.52), clamp((uv.y-.2)/.3,0.,1.));
    
    vec3 mc = mix(vec3(0.27,0.17,0.45),vec3(0.25,0.10,0.30),1.-clamp((uv.y+.5)/.7, 0., 1.));
    float mountains, smountains;
    lfnoise(2.*(uv.x*c.xx), mountains);
    lfsharp(22.*abs(uv.x*c.xx), smountains);
    float d = -.05+uv.y-mix(-.25, .05, abs(uv.x-.1)/.5/a) + .1*mix(mountains, .5*smountains, smoothstep(.4,.5,clamp(abs(uv.x)/.5/a,0.,1.)))*uv.x;
    col = mix(col,mc, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
    
    if(uv.y < -.25)
    {
        uv.x -= .1;
        vec3 ca = mix( vec3(0.11,0.06,0.20),vec3(0.44,0.14,0.36), clamp((uv.y+.5)/.3, 0.,1.));
        vec3 sea;
        lfnoise(5.*uv.y*c.xx, sea.x);
        lfnoise(25.*uv.y*c.xx+122.1*step(sign(uv.x),0.), d);
        sea.x = .5*(.3*sea.x+1.4*d);
        sea = .5*c.xxx+.5*sea;
        sea *= clamp((1.-smoothstep(-.375,-.15,uv.y)),0.,1.);
        stroke(uv.x, 2.*sea.x, d);
        col = mix(col, ca, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
        
        if(d < 0.)
        {
            vec2 delta;
            lfnoise(43.*uv-iTime-133., delta.x);
            lfnoise(46.*uv-iTime-333., delta.y);
            uv += .005*delta;
        }
        
        ca = mix( vec3(0.47,0.15,0.36), vec3(0.68,0.19,0.43), clamp((uv.y+.5)/.3, 0.,1.));
        lfnoise(5.*uv.y*c.xx, sea.x);
        lfnoise(25.*uv.y*c.xx+122.1*step(sign(uv.x),0.)-.24*iTime, d);
        sea.x = .5*(.01*sea.x+.99*d);
        sea = .5*c.xxx+.5*sea;
        sea *= clamp(cos(smoothstep(-.5,-.375,uv.y))*(1.-smoothstep(-.375,-.2,uv.y)),0.,1.);
        stroke(uv.x, 1.2*sea.x, d);
        col = mix(col, ca, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
        
        ca = vec3(0.82,0.26,0.55);
        lfnoise(12.*uv.y*c.xx, sea.x);
        lfnoise(45.*uv.y*c.xx+122.1*step(sign(uv.x),0.)-.2*iTime, d);
        sea.x = -.5+.5*(1.9*abs(sea.x)+.99*d);
        sea = .5*c.xxx+.5*sea;
        sea *= clamp(cos(smoothstep(-.5,-.175,uv.y))*(1.-smoothstep(-.375,-.2,uv.y)),0.,1.);
        stroke(abs(uv.x), .4*sea.x, d);
        col = mix(col, ca, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
        uv.x += .1;
    }

    vec3 dc;
    rand(uv-1.e-1*iTime*c.xx+.13, dc.x);
    rand(uv-1.e-1*iTime*c.xx+.33, dc.y);
    rand(uv-1.e-1*iTime*c.xx+.53, dc.z);
    
    uv.y += .2;
    uv.y -= .3*sign(uv.x);
    uv.x  = abs(uv.x)-.75;
    
    // Trees
    float tree;
    vec2 delta;
	lfnoise(13.*uv-133.-iTime, delta.x);
    lfnoise(16.*uv-333.-iTime-.1, delta.y);
    for(float i=9.; i>=0.; i -= 1.)
    {
        vec3 cc = mix(vec3(0.03,0.04,0.24), vec3(0.03,0.04,0.12), i/10.),
            cc1 = mix(vec3(0.06,0.10,0.31), vec3(0.06,0.06,0.18), i/10.);
        vec3 rii;
        rand(i*c.xx, rii.x);
		rand(i*c.xx+.3, rii.x);
        rand(i*c.xx+.6, rii.x);
        uv.x += .01*rii.y;
        
    	dstar(uv-.005*delta+(-.5+i*.1)*c.yx, 12.+4.*rii.x, vec2(.105+.03*i+.03*rii.y,.18+.03*i+.03*rii.y), tree);
    	col = mix(col, cc, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, tree));
        
        dstar(uv-.005*delta+(-.5+i*.1)*c.yx, 12.+4.*rii.x, vec2(.105+.03*i-.08+.03*rii.y,.18+.03*i-.08+.03*rii.y), tree);
    	col = mix(col, cc1, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, tree));
        uv.x -= .01*rii.y;
    }
    
    col = col + .1*(-.5+dc);
    
    col = clamp(col, 0., 1.);
    col = mix(c.yyy, col, smoothstep(0.,.5,iTime)*(1.-smoothstep(9.5,10.,iTime)));
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
