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
    vec3 right = c.xyy, up = c.yxy, target = c.yyy+.05*vec3(cos(iTime), sin(iTime), 0.);
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
        dm2 = clamp(0.5+0.5*dm2,0.,1.);
        
        d = mix(d,-d, dm);
        
        zextrude(x.z, big_border, dm2-.5*i/7., d);
        //stroke(d, .01, d);
        d -= .1;
        sdf.x = min(sdf.x, d);
    }
    
    // Add guard objects for debugging
    
    float dr = .03;
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
void rot(in vec3 p, out mat3 m)
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
    rot(.1*iTime*vec3(1.1515,1.3321,1.5123) + rt, m);
    
    vec3 dark_green = abs(vec3(0.0,0.78,.52)),
        light_green = abs(vec3(.9, 0.98, 0.28)), 
        gray = .5*length(light_green)*c.xxx/sqrt(3.);
    
    col = c.yyy;
    for(float i=0.; i<7.; i+=1.)
    {
        rot(.1*iTime*vec3(1.1515,1.3321,1.5123) + 3.*rt + i, m);
        
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
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
