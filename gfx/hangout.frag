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

void rot(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// Distance to regular voronoi
void dvoronoi(in vec2 x, out float d, out vec2 ind)
{
    vec2 y = floor(x);
   	float ret = 1.;
    
    //find closest control point. ("In which cell am I?")
    vec2 pf=c.xx, p;
    float df=100., oo;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            vec2 pa;
            rand(p, pa.x);
            rand(p+.1,pa.y);
            p += pa;
            
            d = length(x-p);
            
            df = min(d,df);
            pf = mix(pf,p,step(d,df));
        }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            vec2 pa;
            rand(p, pa.x);
            rand(p+.1,pa.y);
            p += pa;

            vec2 o = p - pf;
            oo = length(o);
            if(oo < 1.e-4) continue;
            d = abs(.5*oo-dot(x-pf, o)/oo);
            ret = min(ret, d);
        }
    
    d = ret;
    ind = pf;
}

// 2D box
void dbox(in vec2 x, in vec2 b, out float d)
{
	vec2 da = abs(x)-b;
	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);
}

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

// Add distance functions with materials
void add(in vec2 sda, in vec2 sdb, out vec2 dst)
{
    dst = mix(sda, sdb, step(sdb.x, sda.x));
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

vec2 ind;
void scene(in vec3 x, out vec2 sdf)
{
    x.z *= 3.;
    
    sdf.y = 1.;
    mfnoise(x.xy, 1., 400., .45, sdf.x);
    stroke(sdf.x, .5, sdf.x);
    sdf.x = x.z+.1*sdf.x;
    
    float R = .07+.1*sdf.x, dis;
    lfnoise(.5*x.y*c.xx, dis);
    vec2 sdb;
    //if(abs(x.x) < 1.)
    {
        vec2 ya = abs(vec2(x.x,x.z)-.4*dis*c.xy)-.015*c.yx + .005*c.xy,
            yi = ya;
        
        ya = vec2(mod(ya.x, mix(2.,.2,smoothstep(0.,3.,x.y)))-.5*mix(2.,.2,smoothstep(0.,3.,x.y)), ya.y);
        yi -= ya;
        yi /= mix(2.,.2,smoothstep(0.,3.,x.y));
        float da;
        if(yi.x < 4.)
        {
            zextrude(x.y, -length(ya)+R, 1.e4, sdb.x);
            stroke(sdb.x, .003, da);
        } else da = 1.e-3;
        sdb.y = 2.;

        float phi = atan(ya.y, ya.x);
        dhexagonpattern(vec2(56.,12.)*vec2(x.y, phi), dis, ind);
        stroke(dis, .2, dis);
        stroke(sdb.x, .003*step(dis,0.), sdb.x);
        sdf.x = max(sdf.x,-da);

        // Add Station
        R = .9;
        vec2 sdc;
        zextrude(x.y-3., -length(x.xz)+R, 1.5, sdc.x);
        stroke(sdc.x, .05, sdc.x);
        sdc.x = mix(sdc.x, 1.e-3, step(-dis,0.));

        // Add normal walls
        float dd = x.y - 2.7, dis1;
        stroke(dd, .02, dd);
        vec2 ind0;
        dhexagonpattern(vec2(56.,22.)*x.xz, dis1, ind);
        stroke(dis1, .2, dis1);
        dd = mix(dd, 1.e-3, step(-dis1,0.));
        sdc.x = min(sdc.x, dd);


        // Remove exterior
        zextrude(x.y-3., -length(x.xz)+1.02*R, 1.05*.5, dd);
        stroke(dd, .03, dd);
        sdc.x = max(sdc.x, dd);

        // Remove interior
        zextrude(x.y-3.,-length(x.xz)+.95*R, 1.*.5, dd);
        stroke(dd, .03, dd);
        sdc.x = max(sdc.x, -dd);

        // Add Floor
        float fl;
        dbox(x.xy-3.*c.yx, vec2(R,.6*.5), fl);
        zextrude(x.z,-fl,.1, fl);
        sdc.x = min(sdc.x, fl);
        
        sdb.x = max(sdb.x, -dd);
        add(sdf, sdb, sdf);
        sdc.y = 2.;
        add(sdf, sdc, sdf);
        
        // Add guard objects for debugging
        float dr = .1;
        vec3 y = mod(x,dr)-.5*dr;
        float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
        guard = abs(guard)+dr*.1;
        sdf.x = min(sdf.x, guard);
    }
}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 5.e-3;
    vec2 s, na;
    
    scene(x,s);
    scene(x+dx*c.xyy, na);
    n.x = na.x;
    scene(x+dx*c.yxy, na);
    n.y = na.x;
    scene(x+dx*c.yyx, na);
    n.z = na.x;
    n = normalize(n-s.x);
}

void planet_texture(in vec2 x, out vec3 col)
{
    vec3 light_orange = vec3(1.00,0.69,0.05),
        orange = vec3(0.95,0.45,0.01),
        dark_orange = vec3(0.98,0.73,0.01);
    
    //rock like appearance
    float d;
    mfnoise(x, 1.,400.,.65,d);
	col = mix(vec3(0.19,0.02,0.00), vec3(0.91,0.45,0.02), .5+.5*d);
    
    // big structures
    stroke(d,.01, d);
    col = mix(mix(.8*vec3(0.99,0.49,0.02),c.yyy,d*clamp(.2-.5*x.y/12.,0.,1.)), col, smoothstep(-.05,.05,d));
    col = mix(col, vec3(0.15,0.05,0.00), clamp(.2-.5*x.y/12.,0.,1.));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0),
        s;
    vec3 col = c.yyy, 
        o = .5*c.yyx + mix(.3*iTime*c.yxy,.3*6.*c.yxy,smoothstep(0.,6.,clamp(iTime,0.,6.))) 
        	+ smoothstep(6.,12.,clamp(iTime, 6., 12.))*c.yyx*.15, 
        t = vec3(uv,0.) + mix(.3*iTime*c.yxy,.3*6.*c.yxy,smoothstep(0.,6.,clamp(iTime,0.,6.))) 
        	+ smoothstep(6.,12.,clamp(iTime, 6., 12.))*c.yyx*.15	
        	+ c.yxy,
        dir = normalize(t-o),
        x;
    float d = 0.;//(.14-o.z)/dir.z;
    
    int N = 450, i;
    for(i=0; i<N; ++i)
    {
        x = o + d * dir;
        scene(x,s);
        if(s.x < 1.e-4) break;
        d += s.x;
    }
    
    if(i<N)
    {
        vec3 n,
            l = normalize(x+c.yyx);
        normal(x,n);
        vec3 c1 = vec3(0.99,0.64,0.02);
        if(s.y == 1.)
        {
            planet_texture(x.xy, c1);
            col = .3*c1
                + (.3*c1)*abs(dot(l,n))
                + (1.3*c1+.1*c.xyy)*pow(abs(dot(reflect(-l,n),dir)),3.);
        }
        else if(s.y == 2.)
        {
            planet_texture(x.xy, c1);
            col = .3*c1
                + (.3*c1)*abs(dot(l,n))
                + (1.3*c1+.1*c.xyy)*pow(abs(dot(reflect(-l,n),dir)),3.);
            col = mix(col,.4*length(col)*c.xyy,.7);
            mat3 RR;
            rot(.01*ind.x*c.xxy, RR);
            col = abs(RR * col);
        }
    }
    /*
    vec3 ddd;
    rand(uv-iTime*c.xx, ddd.x);
    rand(uv-iTime*c.xx, ddd.y);
    rand(uv-iTime*c.xx, ddd.z);
    col -=.1*ddd;
    */
    col = clamp(col, 0.,1.);
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
