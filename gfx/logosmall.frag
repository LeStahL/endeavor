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

const float r = .5;

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

// Stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0)-s;
}

// Add distance functions with materials
void add(in vec2 sda, in vec2 sdb, out vec2 dst)
{
    dst = mix(sda, sdb, step(sdb.x, sda.x));
}

void scene_bounds(in vec3 x, out vec2 s)
{
    vec2 sdmars = c.yx; // Mars has material 1
    sdmars.x = length(x)-r;
    
    vec2 sdphobos = 2.*c.yx; // Phobos has material 2
    float phi = .04*iTime, 
        theta = .4*iTime, 
        rtraj = r*1.2, rmoon = .07*r;
    sdphobos.x = length(x-rtraj*vec3(cos(theta)*cos(phi),
                                     cos(theta)*sin(phi),
                                     sin(theta)))-rmoon;
    add(sdmars, sdphobos, s);
    
    vec2 sddeimos = 3.*c.yx; // Deimos has material 3
    theta = (theta + 13.1)*1.5;
    phi = (phi - 6.22)*1.5;
    sddeimos.x = length(x-rtraj*vec3(cos(theta)*cos(phi),
                                     cos(theta)*sin(phi),
                                     sin(theta)))-.5*rmoon;
    add(s, sddeimos, s);
}

void planet(in vec3 x, out float d)
{
    vec2 texcoord = vec2(sqrt(2./(1./r-x.z)),sqrt(2./(1./r-x.z)))*x.xy;
    mfnoise(12.*texcoord-.2*iTime*c.xy, .9,250.,.45,d);
    stroke(d,.1, d);
    
    d = length(x)-r-mix(0.,.0003,smoothstep(-.05,.05,d))-.002*d;
}

void planet_normal(in vec3 x, out vec3 n)
{
    float s, dx=1.e-4;
    planet(x,s);
    planet(x+dx*c.xyy,n.x);
    planet(x+dx*c.yxy,n.y);
    planet(x+dx*c.yyx,n.z);
    n = normalize(n-s);
}

void phobos(in vec3 x, out float d)
{
    float phi = .04*iTime, 
        theta = .4*iTime, 
        rtraj = r*1.2, rmoon = .07*r;
    x -= rtraj*vec3(cos(theta)*cos(phi),
                                     cos(theta)*sin(phi),
                                     sin(theta));
    vec2 texcoord = vec2(sqrt(2./(1./rmoon-x.z)),sqrt(2./(1./rmoon-x.z)))*x.xy;
    mfnoise(12.*texcoord-.01*iTime*c.xy, 17.4,550.,.45,d);
    stroke(d,.01, d);
    
    d = length(x)-rmoon-.01*d;
}

void phobos_normal(in vec3 x, out vec3 n)
{
    float s, dx=1.e-4;
    phobos(x,s);
    phobos(x+dx*c.xyy,n.x);
    phobos(x+dx*c.yxy,n.y);
    phobos(x+dx*c.yyx,n.z);
    n = normalize(n-s);
}

void deimos(in vec3 x, out float d)
{
    float phi = .04*iTime, 
        theta = .4*iTime, 
        rtraj = r*1.2, rmoon = .07*r*.5;
    theta = (theta + 13.1)*1.5;
    phi = (phi - 6.22)*1.5;
    x -= rtraj*vec3(cos(theta)*cos(phi),
                                     cos(theta)*sin(phi),
                                     sin(theta));
    vec2 texcoord = vec2(sqrt(2./(1./rmoon-x.z)),sqrt(2./(1./rmoon-x.z)))*x.xy;
    mfnoise(12.*texcoord-.01*iTime*c.xy, 17.9,550.,.45,d);
    stroke(d,.01, d);
    
    d = length(x)-rmoon-.01*d;
}

void deimos_normal(in vec3 x, out vec3 n)
{
    float s, dx=1.e-4;
    deimos(x,s);
    deimos(x+dx*c.xyy,n.x);
    deimos(x+dx*c.yxy,n.y);
    deimos(x+dx*c.yyx,n.z);
    n = normalize(n-s);
}

void planet_texture(in vec2 x, out vec3 col)
{
    vec3 light_orange = vec3(1.00,0.69,0.05),
        orange = vec3(0.95,0.45,0.01),
        dark_orange = vec3(0.98,0.73,0.01);
    
    //rock like appearance
    float d;
    mfnoise(-50.+x-.2*iTime*c.xy, 1.,250.,.65,d);
	col = mix(vec3(0.19,0.02,0.00), vec3(0.91,0.45,0.02), .5+.5*d);
    
    // big structures
    mfnoise(x-.2*iTime*c.xy, .9,250.,.45,d);
    stroke(d,.04, d);
    col = mix(mix(.8*vec3(0.99,0.49,0.02),c.yyy,clamp(.2-.5*x.y/12.,0.,1.)), col, smoothstep(-.05,.05,d));
    
    col = mix(col, vec3(0.15,0.05,0.00), clamp(.2-.5*x.y/12.,0.,1.));
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

void rot(in vec3 p, out mat3 rot)
{
    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))
        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))
        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
     a = iResolution.x/iResolution.y;
    float dv, dvmin, vsize = 20.;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0),
        vind, vindmin;
    vec3 col = c.yyy;
    
    // Stars
    dvoronoi(vsize*uv, dv, vind);
    vec2 x = uv-vind/vsize;
    float starsize;
    rand(vind, starsize);
    starsize /= vsize*80.;
    vec3 starcol = vec3(1.00,1.00,0.99), nebulacol = mix(c.xxy,c.yxy,length(uv));
    vec3 starseed;
    rand(vind, starseed.x);
    rand(vind+.1, starseed.y);
    rand(vind+.2, starseed.z);
    mat3 RR;
    rot(.4*(-.5+starseed), RR);
    starcol = abs(RR * starcol);
    col = mix(col, starcol, smoothstep(6.*starseed.x*1.5/iResolution.y, -6.*starseed.x*1.5/iResolution.y, length(x)-starsize));
    col = mix(col, starcol, smoothstep(1.5/iResolution.y, -1.5/iResolution.y, length(x)-starsize));
    
    // Nebula
    float lfdensity, hfdensity;
    lfnoise(uv, lfdensity);
    mfnoise(uv, 6.e0,1.e4, .6, hfdensity);
    float density = .5*(lfdensity+hfdensity);
    col = mix(col, mix(nebulacol,c.xyy, 1.-density), density);
    
    fragColor = vec4(col,1.0);
    //if(length(uv)-r*1.154<0.)
    {	
        vec3 o = c.yyx, t = vec3(uv,0.), dir = normalize(t-o), x, n;
        float d;
        vec2 s;
        
        // trace planet bounding sphere. Nobody will notice 
        // the wrong surface at the border of the projection
        for(int i=0; i<250; ++i)
        {
            x = o + d * dir;
            scene_bounds(x, s);
            if(s.x<1.e-4)break;
            if(d>1.5)
            {
                col += .15*vec3(0.97,0.58,0.00)*(1.-smoothstep(r*1.1,r*1.24, length(uv)));
                fragColor = vec4(col, 1.);
                return;
            }
            d += s.x;
        }

        vec3 light = x+normalize(vec3(0.,1.,1.));
        if(s.y == 1.)
	        planet_normal(x,n);
        else if(s.y == 2.)
            phobos_normal(x,n);
        else if(s.y == 3.)
            deimos_normal(x,n);
        float dln = max(1., dot(light,n)),
            drv = max(0., dot(reflect(light,n),dir));
        if(s.y == 1.)
        {
            vec2 texcoord = vec2(sqrt(2./(1./r-x.z)),sqrt(2./(1./r-x.z)))*x.xy;
            planet_texture(12.*texcoord, col);
            col = -vec3(0.0,0.13,0.10)
            + .4*col*dln
            + 1.3*col*pow(drv,1.);
            
            // atmosphere glow
            //d = abs(d-.49)+.5;
            
            
        }
        else if(s.y == 2.)
            col = .05*vec3(0.66,0.54,0.46)*dln+.2*vec3(0.97,0.80,0.67)*pow(drv,1.);
        else if(s.y == 3.)
            col = .05*vec3(0.66,0.54,0.46)*dln+.2*vec3(0.36,0.36,0.36)*pow(drv,1.);
    }
    
    //else
    {
        // Stars
        
        // Glow effect
    }
    
    col = clamp(col, 0., 1.);
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
