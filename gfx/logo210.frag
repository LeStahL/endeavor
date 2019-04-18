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

// Setup camera
void camerasetup(in vec2 uv, out vec3 ro, out vec3 dir)
{
    vec3 right = c.xyy, up = c.yxy, target = c.yyy;
    ro = c.yyx+.3*vec3(cos(iTime), sin(iTime), 0.)*(1.-smoothstep(11., 13., iTime));
    dir = normalize(target + uv.x * right + uv.y * up - ro);
}

// Compute distance to stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0) - s;
}

// Distance to circle segment
void dcirclesegment(in vec2 x, in float p0, in float p1, out float d)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.0*pi < philo.x && p+2.0*pi > philo.y) || (p-2.0*pi < philo.x && p-2.0*pi > philo.y))
        d = abs(length(x)-1.0);
    else d = min(
        length(x-vec2(cos(p0), sin(p0))),
        length(x-vec2(cos(p1), sin(p1)))
        );
}

// Distance to circle
void dcircle(in vec2 x, out float d)
{
    d = abs(length(x)-1.0);
}

// Distance to line segment
void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Distance to 210 logo
void dlogo210(in vec2 x, out float d)
{
    float d2;
    dcircle(x+c.zy, d);
    dlinesegment(x, c.yz, c.yx, d2);
    d = min(d, d2);
    dcirclesegment(x+c.xy, -.5*pi, .5*pi, d2);
    d = min(d, d2);
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

void box(in vec3 x, in vec3 b, out float d)
{
    d = length(max(abs(x)-b,0.0));
}

// graph traversal for 210 logo effect
void textpre(in vec3 x, out vec2 sdf)
{
    float blend = smoothstep(2.0, 6.0, iTime)*(1.0-smoothstep(6.0,12.0,iTime));
    //blend *= step(-x.x-2.*smoothstep(2.,8.,iTime),-1.);
    box(x, vec3(1., .5, .01+blend), sdf.x);
}

// Perform raymarching for bounding object
void marchbounds(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)
{
    flag = false;
    for(int ia=0; ia<max(N,0); ++ia)
	{
        x = ro + d*dir;
        textpre(x,s);
        if(s.x < eps)
        {
            flag = true;
            break;
        }
        d += s.x;
	}
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
    // Start with z=0 plane
    sdf = vec2(x.z, 7.0);
    vec2 ind;
    float hex;
    dhexagonpattern(48.0*x.xy, hex, ind);
    
    // compute hexagon indices in cartesian coordinates
    vec2 cind = ind/48.0;
    
    // build up team210 logo (t < 12.)
    float inner_logo, logo_border; 
    dlogo210(3.5*cind, inner_logo);
    stroke(inner_logo, 0.35, inner_logo);

    float blend = smoothstep(2.0, 6.0, iTime)*(1.0-smoothstep(6.0,12.0,iTime));
    if(inner_logo < 0.0 && blend >= 1.0e-3)
    {
        float noise;
        lfnoise(24.0*cind.xy-iTime, noise);
        zextrude(x.z,
                 1.5*x.z-inner_logo, 
                 .5*(0.5+0.5*noise)*blend*step(-cind.x-2.*smoothstep(2.,8.,iTime),-1.),
                 sdf.x);
        stroke(sdf.x, 0.05*blend, sdf.x);
        sdf.y = 7.0;
    }
    stroke(sdf.x,0.1,sdf.x);
    
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

// Initial intro
void background2(in vec2 uv, out vec3 col)
{
    col = c.yyy;
    
    // hexagonal grid
    float d, d0;
    vec2 ind;
    dhexagonpattern(48.0*uv, d0, ind);
    d = -d0;
    stroke(d, 0.1, d);
    vec2 cind = ind/48.0;
    
    // build up team210 logo (t < 12.)
    float inner_logo, logo_border; 
    dlogo210(3.5*cind, inner_logo);
    stroke(inner_logo, 0.35, inner_logo);
    stroke(inner_logo, 0.08, logo_border);
    
    // blend back to structure (t < 16., t > 12.)
    float blend = clamp(.25*(iTime-12.), 0., 1.);
    inner_logo = mix(inner_logo, d0, blend);
    logo_border = mix(logo_border, d0, blend);

    // make background change the color with time
    vec2 dt;
    lfnoise(15.0*cind+2.0, dt.x);
    lfnoise(15.0*cind+3.0, dt.y);
    dt *= 2.;
    float dm, dm2;
    lfnoise(50.0*cind, dm);
    dm = 0.5+0.5*dm;
    lfnoise(6.5*cind-dt-2.0*iTime*c.xx, dm2);
    dm2 = 0.5+0.5*dm2;
    
    // Colors
    vec3 orange = vec3(1.,0.27,0.);
    orange = mix(c.yyy, orange, dm2);
    vec3 gray = .5*length(orange)*c.xxx/sqrt(3.);
  
    col = mix(mix(orange,c.xxx,step(-1.,-cind.x-2.*smoothstep(2.,8.,iTime)+.024)), gray, step(-1.,-cind.x-2.*smoothstep(2.,8.,iTime)));
    col = mix(col, gray, step(0.,inner_logo));
    col = mix(col, c.yyy, step(logo_border,0.));
    
    // blend to black at the end
    col = mix(col, c.yyy, clamp(iTime-27., 0., 1.));
    
    // blend in at the beginning
    col = smoothstep(0.,12., iTime)*clamp(col*step(0.,d),0.,1.);
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
    marchbounds(ro, dir, 150, 2.0e-4, x, s, d, hit);

    if(hit) hit = false;
    else d = -ro.z/dir.z;
    marchscene(ro, dir, 500, 2.0e-4, x, s, d, hit);
    
    if(hit)
    {
        vec3 n;
        calcnormal(x, 2.0e-4, n);

        float rs = 1.9;
        vec3 l = x+1.*c.yyx,
        	re = normalize(reflect(-l,n));
        float rev = abs(dot(re,dir)), ln = abs(dot(l,n));
		background2(x.xy, col);
    }
    else
        background2((ro-ro.z/dir.z*dir).xy, col);

    col = clamp(col, 0., 1.);
    col = mix(c.yyy, col, smoothstep(0.,.5,iTime)*(1.-smoothstep(14.5,15.,iTime)));
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
