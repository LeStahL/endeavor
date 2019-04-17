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

// Scene constants
const float R = .5;

// Hash function
void rand(in vec2 x, out float num)
{
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

// 2D box
void dbox(in vec2 p, in vec2 b, out float dst)
{
  	vec2 d = abs(p) - b;
  	dst = length(max(d,0.0)) + min(max(d.x,d.y),0.0); 
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

void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)
{
    vec3 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Compute distance to regular polygon
void dpolygon(in vec2 x, in float N, in float R, out float d)
{
    d = 2.0*pi/N;
    float t = mod(acos(x.x/length(x)), d)-0.5*d;
    d = R-length(x)*cos(t)/cos(0.5*d);
}

// compute distance to regular star
void dstar(in vec2 x, in float N, in vec2 R, out float ds)
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
    ds = dot(x-p1,ff)/length(ff);
}

void rot(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

void dcapcylinder(in vec3 x, in float R, in float h, out float dst)
{
    dst = length(x.xy)-R;
    stroke(dst, .05*R, dst);
    zextrude(x.z, -dst, h, dst);
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

void add(in vec2 sda, in vec2 sdb, out vec2 sdf)
{
    sdf = mix(sda, sdb, step(sdb.x, sda.x));
}

void scene(in vec3 x, out vec2 sdf)
{
    float ra, tm = clamp((iTime-10.)/5., 0.,1.);
    lfnoise(iTime*c.xx, ra);
    vec2 dxy = mix(.04,.2,.5+.5*ra)*vec2(cos(x.z-4.*iTime),sin(x.z-4.*iTime));
    x.xy -= mix(c.yy, dxy, tm);
        mat3 r;
    rot(vec3(0.,0.,ra), r);
    x = mix(x, r*x, tm);
    
    // Rings
    vec3 y = vec3(x.xy, mod(x.z,.5)-.25), 
        ind = (x-y)/.5;
    dcapcylinder(y, 1.1*R, .05*R, sdf.x);
    sdf.y = 1.;
    
    // Vertical stripes
    float phi = atan(x.y, x.x),
        p0 = mod(phi, pi/4.)-pi/8.,
        pp0 = mod(phi-pi/8., pi/4.)-pi/8.;
    vec2 pa = 1.*R*vec2(cos(phi-pp0), sin(phi-pp0)),
        sda;
    sda.x = length(x-vec3(pa,x.z))-.04*R;
    sda.y = 1.;
    add(sdf, sda, sdf);
    pa = 1.1*R*vec2(cos(phi-p0), sin(phi-p0));
    sda.x = length(x.xy-pa)-.03*R;
    sda.y = 1.;
    add(sdf, sda, sdf);
    // Random stripes
    float p1 = p0 - pi/4.;// (sin(ind.z))*pi/4.;
    vec2 pb = 1.1*R*vec2(cos(phi-p1), sin(phi-p1));
    dlinesegment3(y, vec3(pa,.0), vec3(pb,.5), sda.x);
    p1 = p0+pi/4.;
    pb = 1.1*R*vec2(cos(phi-p1), sin(phi-p1));
    float dl;
    dlinesegment3(y, vec3(pa,.0), vec3(pb,-.5), dl);
    sda.x = min(sda.x,dl);
    stroke(sda.x, .025, sda.x);
    add(sdf, sda, sdf);
    sda.x = length(y-vec3(pa,0.))-.05;
    add(sdf, sda, sdf);
    
    // Mountains / Floor
    dbox(x.xz, vec2(.6*R, 1.e5*R), sda.x);
    zextrude(x.y+.9*R, -sda.x, .3*R, sda.x);
    sda.y = 2.;
    add(sdf, sda, sdf);
    
    dbox(vec2(abs(x.x)-.6*R,x.z), vec2(.05*R, 1.e5*R), sda.x);
    zextrude(x.y+.7*R, -sda.x, .001*R, sda.x);
    sda.y = 2.;
    add(sdf, sda, sdf);
    
    // Honeycomb
    float dh;
    vec2 pd = vec2(atan(x.y,x.x), x.z),
        hind;
    dhexagonpattern(2.*pi*vec2(3.,6.)*pd,dh, hind);
    stroke(dh, .3, dh);
    zextrude(length(x.xy)-1.1*R, -dh, .01, sda.x);
    sda.y = 1.;
    add(sdf, sda, sdf);
    
    // Add rails
    dbox(vec2(abs(x.x)-.3*R,x.z), vec2(.02*R, 1.e5*R), sda.x);
    zextrude(x.y+.7*R, -sda.x, .005*R, sda.x);
    sda.y = 1.;
    add(sdf, sda, sdf);
    dbox(y.xz, vec2(.5*R, .1*R), sda.x);
    zextrude(x.y+.73*R, -sda.x, .03*R, sda.x);
    sda.y = 3.;
    add(sdf, sda, sdf);
                
    // Add guard objects for debugging
    float dr = .1;
    y = mod(x,dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    sdf.x = min(sdf.x, guard);

}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 2.e-3;
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

void scene2(in vec3 x, out vec2 sdf)
{
}

void floor_pattern(in vec2 x, out vec3 col)
{
    float w = .02;
    vec2 y = mod(x, w)-.5*w;
    float d;
    dlinesegment(y, -.25*w*c.xx, .25*w*c.xx, d);
    stroke(d,.15*w, d);
    float dx;
    mfnoise(x, 10.,1000., .45, dx);
    col = mix(.4*c.xyy, .5*c.xyy, step(0.,d));
    col = mix(col, .4*c.xyy, .5+.5*dx);
}

// Distance to circle
void dcircle(in vec2 x, in float R, out float d)
{
    d = length(x)-R;
}

void colorize(in vec3 old, in vec2 uv, out vec3 col)
{
    uv.y += .3;
    
    // Main skull
    float d, da, db;
    dcircle(uv-.25*c.yx, .25, d);
    dcircle(uv-.03*c.yx, .1, da);
    dbox(uv-.03*c.yx,.1*c.xx, db);
    da = mix(da,db,.2);
    stroke(da, .03, da);
    d = min(d,da);
    
    // Remove eyes
    dcircle(vec2(abs(uv.x)-.1,1.6*uv.y-.3*abs(uv.x)-.3),.1, da);
    d = max(d,-da);
   	dcircle(vec2(abs(uv.x)-.1,uv.y-.3/1.6-.02),.025, da);
    d = min(d,da);
    dcircle(vec2(abs(uv.x)-.1,uv.y-.3/1.6-.02),.015, da);
    d = max(d,-da);
   	
    // Add hair
	float R = .25, phi;
    for(int i=0; i<10; ++i)
    {
        float ra;
        rand(float(i)*c.xx, ra);
        phi = .2 + 1.*ra;
        
        vec2 p = .25*c.yx+R*vec2(cos(phi), sin(phi));
        dlinesegment(vec2(abs(uv.x),uv.y), p, p+.2*c.yx, da);
        lfnoise(40.*uv.y*c.xx-4.*iTime-12.*ra,db);
        lfnoise(20.*uv.y*c.xx-7.*iTime-12.*ra-.5*db,db);
        stroke(da,.005+.01*db, da);
        d = min(d, da);
    }
    
    // Stripes
    vec3 col1 = mix(c.xxx, vec3(0.70,0.13,0.20), mod(round(14.*(uv.x+.25)),2.));
    
    // Blue box
    float d_, da_;
    dbox(uv-vec2(-.25*a,.5), .25*c.xx*vec2(4.*a,1.), d_);
    col1 = mix(col1, vec3(0.24,0.23,0.43), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d_));
    
    // Stars
    vec2 x = mod(uv,a/12.)-a/24.;
    dstar(x, 5., vec2(.02,.045), da_);
    if(d_ < 0.)
	    col1 = mix(col1, c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, -da_));
    x = mod(uv-a/24.,a/12.)-a/24.;
    dstar(x, 5., vec2(.02,.045), da_);
    if(d_ < 0.)
	    col1 = mix(col1, c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, -da_));
    
    dbox(uv-.27*c.yx, vec2(2.,.03), da_);
    col1 = mix(col1, c.xxx, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, da_));
    
    col1 = .5*length(col1)*c.xxx;
    
    col = mix(old, mix(old,col1,.9), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
    stroke(d-.01,.011,d);
    col = mix(col, c.yyy, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
}

// void colorize(in vec3 old, in vec2 uv, out vec3 col)
// {
//     uv.y += .3;
//     
//     // Main skull
//     float d, da, db;
//     dcircle(uv-.25*c.yx, .25, d);
//     dcircle(uv-.03*c.yx, .1, da);
//     dbox(uv-.03*c.yx,.1*c.xx, db);
//     da = mix(da,db,.2);
//     stroke(da, .03, da);
//     d = min(d,da);
//     
//     // Remove eyes
//     dcircle(vec2(abs(uv.x)-.1,1.6*uv.y-.3*abs(uv.x)-.3),.1, da);
//     d = max(d,-da);
//    	dcircle(vec2(abs(uv.x)-.1,uv.y-.3/1.6-.02),.025, da);
//     d = min(d,da);
//     dcircle(vec2(abs(uv.x)-.1,uv.y-.3/1.6-.02),.015, da);
//     d = max(d,-da);
//    	
//     // Add hair
// 	float R = .25, phi;
//     for(int i=0; i<10; ++i)
//     {
//         float ra;
//         rand(float(i)*c.xx, ra);
//         phi = .2 + 1.*ra;
//         
//         vec2 p = .25*c.yx+R*vec2(cos(phi), sin(phi));
//         dlinesegment(vec2(abs(uv.x),uv.y), p, p+.2*c.yx, da);
//         lfnoise(40.*uv.y*c.xx-4.*iTime-12.*ra,db);
//         lfnoise(20.*uv.y*c.xx-7.*iTime-12.*ra-.5*db,db);
//         stroke(da,.005+.01*db, da);
//         d = min(d, da);
//     }
//     
//     col = mix(old, mix(old,.1*c.xxx,.9), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
//     stroke(d-.01,.011,d);
//     col = mix(col, c.yyy, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d));
// }

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float tm = clamp((iTime-10.)/5., 0.,1.);
    float ra;
    mat3 r;
    lfnoise(iTime*c.xx, ra);
    rot(vec3(0.,0.,12.*ra), r);
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    vec2 dxy = mix(.04,.2,.5+.5*ra)*vec2(cos(uv.y-4.*iTime),sin(uv.y-4.*iTime));
    uv.xy -= mix(c.yy, dxy, tm);
    
    float N = 5.*(.5+.5*sin(1.e-3*iTime)), dd = 2.*pi/N*sin(.1*iTime);
    vec2 k = vec2(cos(dd),sin(dd)),
        tt = c.zx, e;
    mat2 Ra = mat2(k.x,k.y,-k.y,k.x);
    
    for(float i = 2.; i > .05; i = i*(.85+.1*sin(.3*iTime)))
    {
        vec2 uuv = mix(Ra*uv, floor(25.*Ra*uv)/25., clamp((iTime-15.)/5., 0.,1.));
        float dd;
        dstar(uuv, N, i*vec2(1.+.5*cos(3.4221*iTime),1.+.5*sin(2.153*iTime)), dd);
        e = vec2(dd,i);
        tt = mix(tt,e,step(-3./iResolution.y,e.x));
    }
    e = vec2(.025-length(uv),pi); 
    uv = mix(uv, floor(100.*uv)/100., clamp((iTime-20.)/5., 0., 1.));
    
    vec3 col = mix(.3*c.xxx*cos(2.*length(uv)),.8 + .5*cos(1.+uv.xyx+tt.y*1.5e1+iTime+vec3(0.,2.,4.)),tm), 
        o = c.yyx+(.04*iTime*iTime-1.e1)*c.yyx, dir = normalize(vec3(uv,1.));
    mat3 RR;
                        rot(vec3(1.1,1.4,1.6)*iTime ,RR);
                        col = abs(RR*col);
    col = clamp(col, 0., 1.);
    float d;
    vec2 s;
    vec3 x;
    
    fragColor = vec4(col, 1.);
    //if(length(uv) < 1.5*R)
    {
        // Analytical distance to infinite cylinder
        d =.6*R/length(dir.xy); 
		float dout = 1.15*R/length(dir.xy);

        x = o + d * dir;

        // Outer cylindrical structure
        int N = 300, i;
        for(i=0; i<N; ++i)
        {
            x = o + d * dir;
            scene(x, s);
            if(d > dout) 
                return;
            if(s.x < 1.e-3) break;
            d += s.x;
        }
        if(i < N)
        {
            // Calc normal
            vec3 n;
            normal(x, n);

            // Colorize
            vec3 l = normalize(x+c.yyx*sin(x.z));
            if(s.y == 1.)
            {
                col = .1*c.yyx
                    + .3*c.yyx*abs(dot(l, n))
                    + .9*c.xyx*pow(abs(dot(reflect(-l,n), dir)), 2.);
                
            }
            else if(s.y == 2.)
            {
                floor_pattern(x.xz, col);
                
                dir = n;
                o = x;
                d = 2.e-4;
                for(i=0; i<500; ++i)
                {
                    x = o + d * dir;
                    scene(x, s);
                    if(d > dout) 
                    {
    					col = mix(col, .2*c.xxx, tanh(.1*d));
                        
                        vec3 cc;
                        colorize(col, uv, cc);

                        col = mix(col, cc, smoothstep(1.,2.,iTime)*(1.-smoothstep(8.,9.,iTime)));
                        
                        vec3 bw = length(col)*c.xxx;
                        col = mix(bw, col, tm);
                        fragColor = vec4(col,1.0);
                        mat3 RR;
                        rot(vec3(1.1,1.4,1.6)*iTime ,RR);
                        col = abs(RR*col);
                        
                        return;
                    }
                    if(s.x < 1.e-4) break;
                    d += s.x;
                }
                
                if(i<N)
                {
                    normal(x, n);
                    vec3 c1 = .1*c.xyx
                        + .3*c.xyx*abs(dot(l, n))
                        + .9*c.xyx*pow(abs(dot(reflect(-l,n), dir)), 2.);
                    col = mix(col, c1, .3);
                    
                }
                
            }
            else if(s.y == 3.)
            {
                col = .1*c.yxy
                    + .2*c.yxy*abs(dot(l, n))
                    + .8*c.yxy*pow(abs(dot(reflect(-l,n), dir)), 2.);
                float da, db;
                lfnoise(4.*x.x*c.xx, db);
                mfnoise(x.z*c.xx-.03*db*c.xy, 1.e2, 1.e4, .45, da);
                col = mix(col, .4*col, .5+.5*da);
                mfnoise(x.y*c.xx-.03*db*c.xy, 1.e2, 1.e4, .45, da);
                col = mix(col, .4*col, .5+.5*da);
            }
        }
    }
    
    //float daa = 4.*abs(s.x-.4)-1.3;
    //col = mix(col, c.xyx, clamp(daa, 0., 1.));
	rot(vec3(1.1,1.4,1.6)*iTime ,RR);
	col = abs(RR*col);
    vec3 bw = .8*length(col*col)*c.xxx;
    col = mix(bw, col, tm);
    
    col = mix(col, c.yyy, tanh(2.e-1*(x.z-o.z)));
    
    vec3 cc;
    colorize(col, uv, cc);
    
    col = mix(col, cc, smoothstep(1.,2.,iTime)*(1.-smoothstep(8.,9.,iTime)));
    
    // Output to screen
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
