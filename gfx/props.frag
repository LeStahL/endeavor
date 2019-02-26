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

uniform float iTime;
uniform vec2 iResolution;

out vec4 gl_FragColor;

// Global constants
const float pi = acos(-1.);
const vec3 c = vec3(1.0, 0.0, -1.0);
float a = 1.0;

// 2D box
void dbox(in vec2 x, in vec2 b, out float d)
{
	vec2 da = abs(x)-b;
	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

// Compute distance to regular polygon
void dpolygon(in vec2 x, in float N, out float d)
{
    d = 2.0*pi/N;
    float t = mod(acos(x.x/length(x)), d)-0.5*d;
    d = -0.5+length(x)*cos(t)/cos(0.5*d);
}

// Distance to circle
void dcircle(in vec2 x, out float d)
{
    d = length(x)-1.0;
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

// 2D rotational matrix
void rot(in float phi, out mat2 m)
{
    vec2 cs = vec2(cos(phi), sin(phi));
    m = mat2(cs.x, -cs.y, cs.y, cs.x);
}

// Distance to pig ear
void dear(in vec2 x, out float d)
{
    d = abs(2.*x.y)
        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))
        -.5*min(-abs(x.x),.01);
}

// Distance to spacepigs logo in hexagon
void dspacepigs(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Head
    dcircle(2.5*x,d0);
    d0 /= 2.5;
    
    // Ears
    dear(vec2(2.,5.)*x-vec2(.8,1.3), da);
    d0 = min(d0,da/10.);
    dear(vec2(2.,5.)*x+vec2(.8,-1.3), da);
    d0 = min(d0,da/10.);
    
    // Nose
    dcircle(6.*x-vec2(0.,-.5),da);
    d0 = max(d0,-da/6.);
    dcircle(24.*x-vec2(-1.5,-2.),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(1.5,-2.),da);
    d0 = min(d0,da/24.);
    
    // Eyes
    dcircle(16.*x-vec2(-3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(16.*x-vec2(3.5,2.5),da);
    d0 = max(d0,-da/16.);
    dcircle(24.*x-vec2(-5.,3.5),da);
    d0 = min(d0,da/24.);
    dcircle(24.*x-vec2(5.,3.5),da);
    d0 = min(d0,da/24.);
    
    d = max(d, -d0);
}

// Distance to kewlers logo in hexagon
void dkewlers(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    x *= 1.2;
    
    dbox(x-vec2(0.,-.3),vec2(.6,.1),d0);
    dbox(x-vec2(-.5,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    dbox(x-vec2(-.5+1./3.,.25),vec2(.1,.5),da);
    d0 = min(d0,da);
    dbox(x-vec2(-.5+2./3.,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    dbox(x-vec2(.5,-.0),vec2(.1,.25),da);
    d0 = min(d0,da);
    
    d = max(d, -d0);
}

// Distance to farbrausch logo in hexagon
void dfarbrausch(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    x += vec2(.1,0.);
    x *= 1.2;
    
    dlinesegment(x,vec2(-.65,.05),vec2(-.5,.05),d0);
    dlinesegment(x,vec2(-.5,.05),vec2(-.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.2,-.49),vec2(-.0,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.0,-.49),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(-.07,0.),vec2(-.27,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.2,-.49),vec2(-.07,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.4,-.49),vec2(.2,-.49),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.13,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.33,0.),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.6,-.15),vec2(.51,-.33),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.53,0.),vec2(.6,-.15),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.53,.0),da);
    d0 = min(d0, da);
    dlinesegment(x,vec2(.7,0.),vec2(.68,-.04),da);
    d0 = min(d0, da);
    dpolygon(5.*(x+vec2(.3,.65)),6.,da);
    d0 = min(d0, da/5.);
    dpolygon(5.*(x+vec2(-.5,.65)),6.,da);
    d0 = min(d0, da/5.);
    
    stroke(d0,.035, d0);
    d = max(d, -d0);
}

// Distance to haujobb logo in hexagon
void dhaujobb(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    mat2 m;
	rot(.3,m);
    x = 1.1*x*m;
    x.x *= 1.1;
        
    x += vec2(-.05,.2);
    
    // Left leg
    dbox(x+.35*c.xx,vec2(.1,.05),d0);
    dbox(x+vec2(.3,.25),vec2(.05,.15),da);
    d0 = min(d0,da);
    dbox(x+vec2(.2,.15),vec2(.1,.05),da);
    d0 = min(d0,da);
    dbox(x+vec2(.15,.05),vec2(.05,.15),da);
    d0 = min(d0,da);
    
    // Right leg
    dbox(x-vec2(.65,.35),vec2(.05,.15),da);
    d0 = min(d0,da);

    // Torso
    rot(.2, m);
    dbox(m*(x-vec2(.25,.15)),vec2(.45,.05),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(-.15,.35)),vec2(.45,.05),da);
    d0 = min(d0,da);
    rot(pi/8.,m);
    dbox(m*(x-vec2(.0,.25)),vec2(.1,.15),da);
    d0 = min(d0,da);
    
    // Penis
    dbox(m*(x-vec2(.1,-.0)),vec2(.025,.1),da);
    d0 = min(d0,da);
    
    // Left hand
    rot(.3,m);
    dbox(m*(x-vec2(.235,.535)),vec2(.035,.15),da);
    d0 = min(d0,da);
    dbox(m*(x-vec2(.225,.7)),vec2(.075,.025),da);
    d0 = min(d0,da);
    
    // Right hand
    rot(-.2,m);
    dbox(m*(x+vec2(.585,-.2)),vec2(.0375,.1),da);
    d0 = min(d0,da);
    
    // Head
    dcircle(6.*(x-vec2(-.15,.58)),da);
    d0 = min(d0,da/6.);
    
    d0 -= .05*(abs(x.x)+abs(x.y)-.2);
    d = max(d,-d0);
}

// Distance to mercury logo in hexagon
void dmercury(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da;

    x += .1*c.yx;

    // Upper part
    dbox(x-.35*c.yx,vec2(.4,.35), da);
    d = max(d, -da);
    dbox(x-.7*c.yx, vec2(.2,.2), da);
    d = min(d,da);
    dbox(x-.25*c.yx,vec2(.2,.05),da);
    d = min(d,da);
    
    // Lower part
    dbox(x+.2*c.yx,vec2(.1,.4),da);
    d = max(d, -da);
    dbox(x+.2*c.yx, vec2(.4,.1),da);
    d = max(d, -da);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));
    col = length(col)/sqrt(1.)*c.xyy;
    float d, index = (iTime-mod(iTime, 4.))/4.;
    index = mod(index, 5.);
    if(index == 0.)
	    dmercury(2.*uv, d);
    else if(index == 1.)
        dhaujobb(2.*uv, d);
    else if(index == 2.)
        dfarbrausch(2.*uv, d);
    else if(index == 3.)
        dkewlers(2.*uv, d);
    else if(index == 4.)
        dspacepigs(2.*uv, d);
    col *= step(d,0.);
    stroke(d+.02,.01,d);
    col = mix(col, c.xxy, step(d,0.));
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
