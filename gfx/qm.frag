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

// Hash function
void rand(in vec2 x, out float num)
{
    num = fract(sin(dot(abs(x) ,vec2(12.9898,78.233)))*43758.5453);
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
void dcircle(in vec2 x, in float r, out float d)
{
    d = abs(length(x)-r);
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

// Distance to circle segment
void circlesegment(in vec2 x, in float r, in float p0, in float p1, out float d)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))
    {
        d = abs(length(x)-r);
        return;
    }
    d = min(
        length(x-r*vec2(cos(p0), sin(p0))),
        length(x-r*vec2(cos(p1), sin(p1)))
        );
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

// Extrusion
void zextrude(in float z, in float d2d, in float h, out float d)
{
    vec2 w = vec2(-d2d, abs(z)-0.5*h);
    d = length(max(w,0.0));
}

// Distance to qm
void dqm(in vec2 x, out float d)
{
    // Q
    x.x += .1;
    float s = .25, d0;
    dcircle(x+.35*c.xy, s, d);
    dlinesegment(x+.35*c.xy, s*vec2(1./3., -1./3.), s*vec2(1.,-1.), d0);
    d = min(d,d0);
    
    // Connection
    dlinesegment(x-.075*c.xy, s*vec2(-.7,-1.), s*vec2(.7,-1.), d0);
    d = min(d,d0);
    
    // M
    vec2 y = x-.5*c.xy;
    circlesegment(y, s, 0., pi, d0);
    d = min(d,d0);
    dlinesegment(y, c.yy, s*vec2(0.,1.), d0);
    d = min(d,d0);
    y.x = -abs(y.x);
    dlinesegment(y, s*vec2(-1.,-1.), s*vec2(-1.,0.), d0);
    d = min(d,d0);
}

// Distance to regular voronoi
void dvoronoi(in vec2 x, out float d, out vec2 ind)
{
    vec2 y = floor(x);
       float ret = 1.;
    
    //find closest control point. ("In which cell am I?")
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    
    //compute voronoi distance: minimum distance to any edge
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    ind = pf;
}

void scene(in vec3 x, out float d)
{
        vec2 dis = c.yy, ind;
        
        if(iTime > 6.) // Distort
        {
            vec2 da;
            
            dvoronoi(12.*x.xy-iTime*c.xy, dis.x, ind);
            dvoronoi(12.*x.xy-iTime*c.yx, dis.y, ind);
            
            dvoronoi(23.*x.xy-iTime*c.xy, da.x, ind);
            dvoronoi(23.*x.xy-iTime*c.yx, da.x, ind);
            dis += .5*da;
            
            dis *= 3.;
        }

    x.xy += mix(c.yy,.1*dis,clamp((iTime-6.)/2.,0.,1.));
    
    dqm(x.xy,d);
    stroke(d, .12, d);

    d = -d;
    zextrude(x.z+1.1*clamp((iTime-10.)/2.,0.,1.),d,.3*clamp((iTime-10.)/2.,0.,1.),d);
    
    // Add guard objects for debugging
    float dr = .2;
    vec3 y = mod(x,dr)-.5*dr;
    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));
    guard = abs(guard)+dr*.1;
    d = min(d, guard);
}

void normal(in vec3 x, out vec3 n)
{
    const float dx = 1.e-4;
    float s;
    scene(x,s);
    scene(x+dx*c.xyy, n.x);
    scene(x+dx*c.yxy, n.y);
    scene(x+dx*c.yyx, n.z);
    n = normalize(n-s);
}

// Distance to regular star
void dstar(in vec2 x, in float N, in vec2 R, out float dst)
{
    float d = pi/N,
        p0 = atan(x.y,x.x),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 aa = mix(R,R.yx,i),
        p1 = aa.x*c.xy,
        ff = aa.y*vec2(cos(d),sin(d))-p1;
       ff = ff.yx*c.zx;
    dst = dot(x-p1,ff)/length(ff);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), ind;
    vec3 col = c.yyy, maincol = vec3(0.,0.27,1.);
    vec3 ca = length(maincol)/sqrt(3.)*c.xxx;
    uv.y -= .07;
    float d4, dborder, dborder2, dc;
    if(iTime > 0.) // Add grid after 2 seconds
    {
        vec2 y = mod(uv, .02)-.01,
            a = abs(y)-.0005;
        vec3 newcol = mix(mix(c.yyy, .1*maincol, sin(pi/2.*length(uv))), .5*maincol, step(a.x,0.));
        newcol = mix(newcol, .5*maincol, step(a.y,0.));
        y = mod(uv-.01,.08)-.04;
        a = abs(y)-.001;
        newcol = mix(newcol, maincol, step(a.x,0.));
        newcol = mix(newcol,maincol, step(a.y,0.));
        col = mix(col, newcol, clamp(iTime/2.,0.,1.)-clamp((iTime-14.)/2.,0.,1.));
    }
    if(iTime > 12.) // Add background scales
    {
        float db, dc;
        dhexagonpattern(32.*uv,db,ind);
        //dvoronoi(16.*uv,db, ind);
        db = mix(1.,db,clamp((iTime-12.)/2.,0.,1.));
        stroke(db,.1,db);
        rand(ind, dc);
        dc = mix(0., dc, clamp((iTime-12.)/2.,0.,1.));
            ca = mix(mix(ca,clamp(1.4*length(vec3(1.,dc*0.47,dc*.1))/sqrt(3.)*c.xxx*.5,0.,1.),clamp((iTime-12.)/2.,0.,1.)), c.yyy, step(db,0.));
            ca = mix(ca,mix(ca, length(vec3(.3,.01, .01))/sqrt(3.)*.1*c.xxx, clamp(abs((ind.y)/16.+.07),0.,1.)), clamp((iTime-12.)/2.,0.,1.));
        col = mix(col, mix(col, ca, step(-d4,0.)),clamp((iTime-12.)/2.,0.,1.));
    }
    vec2 dis = c.yy;
    if(iTime > 2.) // Add qm after another 2
    {
        
        
        if(iTime > 6.) // Distort
        {
            vec2 da;
            
            dvoronoi(12.*uv-iTime*c.xy, dis.x, ind);
            dvoronoi(12.*uv-iTime*c.yx, dis.y, ind);
            
            dvoronoi(23.*uv-iTime*c.xy, da.x, ind);
            dvoronoi(23.*uv-iTime*c.yx, da.x, ind);
            dis += .2*da;
            
            dis *= 2.;
        }
        dqm(uv+mix(c.yy,.1*dis,clamp((iTime-6.)/2.,0.,1.)), d4);
        stroke(d4, mix(.01,.12,clamp((iTime-2.)/2.,0.,1.)), d4);
        ca = maincol;
        
        if(iTime > 8.) // Mix color to voronoi profile
        {
            float db,db0;
            vec2 ind, ind0;
            //dvoronoi(16.*uv+dis,db, ind);
            dhexagonpattern(16.*(uv+.05*dis),db,ind);
            db = -db;
            stroke(db,.1,db);
            dhexagonpattern(32.*(uv+.05*dis),db0,ind0);
            db0 = -db0;
            stroke(db0,.1,db0);
            db = min(db,db0);
            ind += .1*ind0;
            dhexagonpattern(64.*(uv+.05*dis),db0,ind0);
            db0 = -db0;
            stroke(db0,.1,db0);
            db = min(db,db0);
            ind += .1*ind0;
            db = mix(1.,db,clamp((iTime-8.)/2.,0.,1.));
            rand(ind, dc);
            dc = mix(0., dc, clamp((iTime-8.)/2.,0.,1.));
            ca = mix(mix(ca,clamp(1.4*vec3(.1,dc*0.47,dc*1.),0.,1.),clamp((iTime-8.)/2.,0.,1.)), clamp(0.*vec3(0.27,.05,.48),0.,1.), step(db,0.));
            ca = mix(ca,mix(ca, vec3(.01,.01, .3), clamp(abs((ind.y-.07)/10.),0.,1.)), clamp((iTime-8.)/2.,0.,1.));
        }
        col = mix(col, ca, step(d4,0.));
    }
    if(iTime > 4.) // Add border
    {
        stroke(d4,mix(.01,.035,clamp((iTime-4.)/2.,0.,1.)),dborder);
        col = mix(col, clamp(.0*vec3(.17,0.17,.68),0.,1.), step(dborder, 0.));
        stroke(dborder,mix(.001,.0035,clamp((iTime-4.)/2.,0.,1.)),dborder2);
        col = mix(col, clamp(1.5*vec3(.0,0.51,.91),0.,1.), step(dborder2, 0.));
    }
    
    col = clamp(col, 0., 1.);

    if(iTime > 10.) // Yes, add that 3d already!
    {
        float s = 0., depth;
        vec3 o = c.yyx-1.*c.yxy + .5*mix(0.,.1, clamp((iTime-10.)/2.,0.,1.))*vec3(cos(iTime), sin(iTime),0.);
        vec3 dir = normalize(c.yyz+uv.x*c.xyy+uv.y*c.yxy-o), x,
            l = normalize(c.yxx), n;
           int i;
        depth = (-1.-o.z)/dir.z;
        for(i=0; i<750; ++i)
        {
            x = o+depth * dir;
            scene(x, s);
            if(s<1.e-4) break;
            depth += s;
            if(x.z<-1.3)break;
        }
        if(x.z > -1.29 && x.z < -1.01)
        {
            normal(x,n);
            col = .4*mix(c.yyy,clamp(1.3*maincol, 0., 1.), clamp((x.z+1.15)/.15,0.,1.)) 
                + 1.5*maincol*abs(dot(l,n))
                + maincol*pow(abs(dot(reflect(-l,n),dir)), 3.);
            col = clamp(col, 0., 1.);
        }
    }
    
    if(iTime > 14.) 
    {
        float d;
        vec2 y = mod(uv, .2)-.1, yi = floor((uv-y)/.2);
        dqm((yi-.1)*.2,d);
        stroke(d,.15,d);
        if(d<0.)
        {
            vec2 da;
            lfnoise(16.*yi*.2-2.e0*iTime*c.xy, da.x);
            lfnoise(16.*yi*.2-2.3e0*iTime*c.yx, da.y);
            da = .5+.5*da;
            
            float dg;
            vec2 ri;
            rand(yi+.1, ri.x);
            rand(yi+.2, ri.y);
            if(da.x < .5 && da.y <.4)
            {
                mat2 m;
                rot(ri.y-iTime,m);
                dstar(m*y-.05*da, 5.+floor(4.*ri.y), .2*ri.x*vec2(.01, .04), dg);
                col = mix(col, mix(col,c.xxx,.5), step(-dg,0.));
                stroke(dg, .0025, dg);
                col = mix(col, c.xxx, step(dg,0.));
            }
        }
    }
    
    vec3 rd;
    rand(uv-iTime,rd.x);
    rand(uv*2.-iTime,rd.y);
    rand(uv*3.-iTime,rd.z);
    col += .2*(-c.xxx+rd);
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
