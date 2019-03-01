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
    num = fract(sin(dot(x-15. ,vec2(12.9898,78.233)))*43758.5453);
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), ind;
    vec3 col = c.yyy;
    
    // Background grid
    float d = 1., da;
    dhexagonpattern(25.*uv, d, ind);
    d/=5.;
    stroke(d,.005,d);
    d = mix(d, 1., step(0., min(length(ind/25.)-.34, abs(ind.y/25.)-.1)));
    col = mix(col, .5*c.xxx, step(d,0.));
    d = mix(d, -1., step(0., min(length(ind/25.)-.34, abs(ind.y/25.)-.1)));
    col = mix(col, .1*c.xxx, step(-d+.05, 0.));
    
    // Rotating stuff
    mat2 m, mt;
    rot(iTime, mt);
    rot(pi/6.,m);
    
    // stripes with hexagons
    float w = .1, l, d0, d1;
    vec2 y = mod(uv-1.5*c.xy*iTime, vec2(w,.6*w))-.5*vec2(w,.6*w), yi = (uv-1.5*c.xy*iTime-y)/vec2(w,.6*w);
    rand(yi, l);
    rand(yi+1., d0);
    dlinesegment(y, -.4*l*w*c.xy, .4*w*l*c.xy, d);
    stroke(d, .01, d);
    dpolygon(mt*y/w*2.,6.,d1); 
    if(d0 < .5 && abs(yi.y) < 2.+1.3*sin(23.*yi.x-3.*iTime))
    {
	    col = mix(col, mix(vec3(0.78,.5*d0,.52),c.yyy,.2), step(mix(d,d1,step(d0,.25)),0.));
        stroke(d, .004, d);
        stroke(d1*w/2., .004, d1);
		col = mix(col, vec3(0.78,.5*d0,.52), step(mix(d,d1,step(d0,.25)),0.));
    }
    y = mod(.3*w+uv-1.5*c.xy*iTime, vec2(w,.6*w))-.5*vec2(w,.6*w);
        yi = (.5*w+uv-1.5*c.xy*iTime-y)/vec2(w,.6*w);
    rand(yi+5.*c.yx, l);
    rand(yi+15.*c.yx, d0);
    dlinesegment(y, -.4*l*w*c.xy, .4*w*l*c.xy, d);
    stroke(d, .01, d);
    dpolygon(mt*y/w*2.,6.,d1); 
    if(d0 < .3 && abs(yi.y) < 2.+1.5*sin(18.*yi.x-3.*iTime))
    {
	    col = mix(col, mix(vec3(0.78,.5*d0,.78),c.yyy,.2), step(mix(d,d1,step(d0,.25)),0.));
        stroke(d, .004, d);
        stroke(d1*w/2., .004, d1);
		col = mix(col, vec3(0.78,.5*d0,.78), step(mix(d,d1,step(d0,.25)),0.));
    }
    
    uv = 2.* mt * uv;
    dpolygon(.8*uv,6.,d);
    col = mix(col, mix(vec3(0.78,0.,.78),col,.7),step(d,0.));
    stroke(d,.01,d);
    col = mix(col, vec3(0.78,0.,.78),step(d,0.));
    
    // White grid
    dhexagonpattern(5.*uv, d, ind);
    d/=5.;
    stroke(d,.01,d);
    d = mix(d, 1., step(0., length(uv)-.37));
    dpolygon(1.02*m*(uv-.215*vec2(a,1.)*c.xz), 6., da);
    d = mix(d, 1., step(da,0.));
    col = mix(col, c.xxx, step(d,0.));
    
    // Planets
    float t = mod(iTime, 1.);
    d = length(uv-mix(c.yy,vec2(.35,0.),t))-mix(.1, .15, t);
    col = mix(col, mix(c.yyy,vec3(0.52,0.52,.52),.5), step(d,0.));
    stroke(d,.00475,d);
    col = mix(col, vec3(0.52,0.52,.52), step(d,0.));
    d = length(uv-mix(vec2(.35,0.), vec2(.2,-.3),t))-mix(.15, .07,t);
    col = mix(col, mix(c.yyy,vec3(0.52,0.52,.52),.5), step(d,0.));
    stroke(d,.00475,d);
    col = mix(col, vec3(0.52,0.52,.52), step(d,0.));
    d = length(uv-mix(vec2(.2,-.3), c.yy, t))-mix(.07, .1, t);
    col = mix(col, mix(c.yyy,vec3(0.52,0.52,.52),.5), step(d,0.));
    stroke(d,.00475,d);
    col = mix(col, vec3(0.52,0.52,.52), step(d,0.));
    
	// Hexagons
    dpolygon(3.5*m*(uv+.35*c.xy), 6., d);
    col = mix(col, mix(vec3(0.78,0.,.52),c.yyy,.2), step(d,0.));
    stroke(d/3.5,.00675, d);
    col = mix(col, vec3(0.98,0.,.72), step(d,0.));
    dpolygon(3.5*m*(uv)+vec2(1.075,.65), 6., d);
    col = mix(col, mix(vec3(0.78,0.,.52),c.yyy,.2), step(d,0.));
    stroke(d/3.5,.00675, d);
    col = mix(col, vec3(0.98,0.,.72), step(d,0.));
    dpolygon(3.5*m*(uv)-vec2(1.075,.65), 6., d);
    col = mix(col, mix(vec3(0.78,0.,.52),c.yyy,.2), step(d,0.));
    stroke(d/3.5,.00675, d);
    col = mix(col, vec3(0.98,0.,.72), step(d,0.));
    dpolygon(3.5*m*(uv)+vec2(.0,-1.25), 6., d);
    col = mix(col, mix(vec3(0.78,0.,.52),c.yyy,.2), step(d,0.));
    stroke(d/3.5,.00675, d);
    col = mix(col, vec3(0.98,0.,.72), step(d,0.));
    
    col = clamp(col, 0., 1.);
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
