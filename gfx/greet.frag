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

// Setup camera
void camerasetup(in vec2 uv, out vec3 ro, out vec3 dir)
{
    vec3 right = c.xyy, up = c.yxy, target = c.yyy;
    ro = c.yyx;
    dir = normalize(target + uv.x * right + uv.y * up - ro);
}

// Compute distance to stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0) - s;
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
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

// 3D Effect on text in intro (210 logo)
void texteffect(in vec3 x, out vec2 sdf)
{
    sdf = vec2(x.z,7.);
    vec2 uv = x.xy;
    
    //uv -= mod(iTime, 1.5)*c.xy;
    //uv += mod(iTime,1.5)*c.yx;
    uv /= 2.+mod(iTime, 4.);
    
    float time_index = (iTime-mod(iTime, 4.))/4.;
    lfnoise(.5*iTime*c.xx, time_index);
    time_index *= pi/3.;
    uv = mat2(cos(time_index),-sin(time_index), sin(time_index), cos(time_index))*uv;
    
    
    for(float i=0.; i<6.; i+=1.)
    {
        float d, d0, res = 4.*pow(3., i), bres = res/3.;
        vec2 ind;
        dhexagonpattern(res*uv, d0, ind);
        d = -d0;
        stroke(d, 0.1, d);
        
        vec2 cind = ind/res;
        
        float big_hexagons, big_border;
        vec2 big_ind;
        dhexagonpattern(bres*cind, big_hexagons, big_ind);
        stroke(big_hexagons, 0.2, big_hexagons);
        big_border = big_hexagons;
        stroke(big_border, .1/bres, big_border);
        
        vec2 dt;
        lfnoise(res*15.0*cind+2.0, dt.x);
        lfnoise(res*15.0*cind+3.0, dt.y);
        dt *= 2.;
        float dm, dm2;
        lfnoise(50.0*cind, dm);
        dm = 0.5+0.5*dm;
        lfnoise(6.5*cind-dt-2.0*iTime*c.xx, dm2);
        // change sign here for different effect
        dm2 = clamp(0.5-0.5*dm2,0.,1.);
        
        d = mix(d,-d, dm);
        
        zextrude(x.z, big_border, dm2-.5*i/7., d);
        //stroke(d, .01, d);
        d -= .1;
        sdf.x = min(sdf.x, d);
    }
    
    // Add guard objects for debugging
    
    float dr = .05;
    vec3 y = mod(x,dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    sdf.x = min(sdf.x, guard);

}

// Perform raymarching for actual object
void marchscene(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)
{
    flag = false;
    for(int ia=0; ia<max(N,0); ++ia)
	{
        x = ro + d*dir;
        texteffect(x,s);
        if(s.x < eps)
        {
            flag = true;
            break;
        }
        d += s.x;
	}
}

void calcnormal(in vec3 x, in float eps, out vec3 n)
{
    vec2 s, sp;
    texteffect(x, s);
    texteffect(x+eps*c.xyy, sp);
    n.x = sp.x-s.x;
    texteffect(x+eps*c.yxy, sp);
    n.y = sp.x-s.x;
    texteffect(x+eps*c.yyx, sp);
    n.z = sp.x-s.x;
    n = normalize(n);
}

// 3D rotational matrix
void rot3(in vec3 p, out mat3 m)
{
    m = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// Initial intro
void background2(in vec2 uv, out vec3 col)
{
    uv /= 2.+mod(iTime, 4.);
    float time_index = (iTime-mod(iTime, 4.))/4., rt;
    rt = time_index;
    
    lfnoise(.5*iTime*c.xx, time_index);
    time_index *= pi/3.;
    uv = mat2(cos(time_index),-sin(time_index), sin(time_index), cos(time_index))*uv;
    
    mat3 m;
    rot3(.1*iTime*vec3(1.1515,1.3321,1.5123) + rt, m);
    
    vec3 dark_green = abs(vec3(0.0,0.78,.52)),
        light_green = abs(vec3(.9, 0.98, 0.28)), 
        gray = .5*length(light_green)*c.xxx/sqrt(3.);
    
    col = c.yyy;
    for(float i=0.; i<7.; i+=1.)
    {
        rot3(.1*iTime*vec3(1.1515,1.3321,1.5123) + 3.*rt + i, m);
        
        float d, d0, res = 4.*pow(3., i), bres = res/3.;
        vec2 ind;
        dhexagonpattern(res*uv, d0, ind);
        d = -d0;
        stroke(d, 0.1, d);
        vec2 cind = ind/res;
        
        float big_hexagons, big_border;
        vec2 big_ind;
        dhexagonpattern(bres*cind, big_hexagons, big_ind);
        stroke(big_hexagons, 0.2, big_hexagons);
        big_border = big_hexagons;
        stroke(big_border, .1/bres, big_border);
        
        vec2 dt;
        lfnoise(res*15.0*cind+2.0, dt.x);
        lfnoise(res*15.0*cind+3.0, dt.y);
        dt *= 2.;
        float dm, dm2;
        lfnoise(50.0*cind, dm);
        dm = 0.5+0.5*dm;
        lfnoise(26.5*cind-dt-2.0*iTime*c.xx, dm2);
        dm2 = clamp(0.7+0.5*dm2, 0., 1.);
        
        light_green = mix(light_green, abs(m*m*m*m*light_green), dm2);
        dark_green = mix(c.yyy, abs(m*m*m*dark_green), dm2);
        
        gray = .5*length(light_green)*c.xxx/sqrt(3.);
        col += mix(light_green, dark_green, step(0.,big_border));
        //col = mix(col, light_green, step(big_border,0.));
        
        col = smoothstep(0.,12., iTime)*clamp(col*step(0.,d),0.,1.);
    }
    
    col += .2*dark_green;
}

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

// Distance to a triangle
void dtriangle(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float d)
{
    vec2 d1 = c.xz*(p1-p0).yx, d2 = c.xz*(p2-p1).yx, d3 = c.xz*(p0-p2).yx;
    d = -min(
        dot(p0-x,d1)/length(d1),
        min(
            dot(p1-x,d2)/length(d2),
            dot(p2-x,d3)/length(d3)
        )
    );
}

// Distance to schnappsgirls logo in hexagon
void dschnappsgirls(in vec2 x, out float d)
{
    dpolygon(.5*x,6.0,d);
    float da, d0;
    
    // Dress
    dtriangle(x, vec2(-.1,-.3), vec2(.5,-.3), vec2(.2, .6), d0);
    dlinesegment(x, vec2(-.1,.325), vec2(.5,.325), da);
    stroke(da,.06,da);
    d0 = max(d0,-da);
    
    // Head
    dcircle(7.*(x-vec2(.2,.5)), da);
    d0 = max(d0, -da+.5);
    d0 = min(d0, da/7.);
    
    // Legs
    dlinesegment(x, vec2(.125,-.3), vec2(.125,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.275,-.3), vec2(.275,-.6), da);
    stroke(da, .06, da);
    d0 = min(d0, da);
    
    // Shoulders
    dlinesegment(x, vec2(0.05,.25), vec2(.35,.25), da);
    stroke(da, .085, da);
    d0 = min(d0, da);
    
    // Arms
    dlinesegment(x, vec2(.385,.25), vec2(.5, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(.017,.25), vec2(-.1, -.1), da);
    stroke(da, .055, da);
    d0 = min(d0, da);
    
    // Glass
    dtriangle(x, vec2(-.6,.3), vec2(-.4,.1), vec2(-.2,.3), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dlinesegment(x, vec2(-.4,.15), vec2(-.4,-.1), da);
    stroke(da, .0125, da);
    d0 = min(d0, da);
    dtriangle(x, vec2(-.5,-.15), vec2(-.3,-.15), vec2(-.4,-.1), da);
    d0 = min(d0, da);
    
    // Liquid
    dtriangle(x, vec2(-.55,.25), vec2(-.4,.1), vec2(-.25,.25), da);
    d0 = min(d0, da);
    
    // Salad
    dlinesegment(x, vec2(-.4,.1), vec2(-.2,.5), da);
    stroke(da, .01, da);
    d0 = min(d0, da);
    dcircle(24.*(x-vec2(-.3,.3)), da);
    d0 = min(d0, da/24.);
    dcircle(24.*(x-vec2(-.25,.4)), da);
    d0 = min(d0, da/24.);
    
    d = max(d, -d0);
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
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), s = c.xy;

	vec3 ro, x, dir;
    
    float d = 0.;
    bool hit = false;
    
    vec3 col = c.yyy;
                
	camerasetup(uv, ro, dir);
    d = (.5-ro.z)/dir.z;
    marchscene(ro, dir, 500, 1.0e-4, x, s, d, hit);
    
    if(hit)
        background2(x.xy, col);
    else
        background2((ro-ro.z/dir.z*dir).xy, col);

    uv*=24.;
    uv /= 2.+mod(iTime, 4.);
    float time_index = (iTime-mod(iTime, 4.))/4., rt;
    rt = time_index;
    
    lfnoise(.5*iTime*c.xx, time_index);
    time_index *= pi/3.;
    time_index -= pi/2.-pi/6.*2.;
    uv = mat2(cos(time_index),-sin(time_index), sin(time_index), cos(time_index))*uv;
    
    float index = (iTime-mod(iTime, 4.))/4.;
    index = mod(index, 6.);
    if(index == 0.)
	    dmercury(uv, d);
    else if(index == 1.)
        dhaujobb(uv, d);
    else if(index == 2.)
        dfarbrausch(uv, d);
    else if(index == 3.)
        dkewlers(uv, d);
    else if(index == 4.)
        dspacepigs(uv, d);
    else if(index == 5.)
        dschnappsgirls(uv,d);
        
    col = mix(col, mix(col, c.xxx, .8), step(d,0.));
    
//     stroke(d+.02,.01,d);
//     col = mix(col, c.yyy, step(d,0.));
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
