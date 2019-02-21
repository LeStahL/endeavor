/* Generated with shader-compressor by NR4/Team210. */
#ifndef LOAD_H
#define LOAD_H
const char * load_frag =
"/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
"* Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
"*\n"
"* This program is free software: you can redistribute it and/or modify\n"
"* it under the terms of the GNU General Public License as published by\n"
"* the Free Software Foundation, either version 3 of the License, or\n"
"* (at your option) any later version.\n"
"*\n"
"* This program is distributed in the hope that it will be useful,\n"
"* but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"* GNU General Public License for more details.\n"
"*\n"
"* You should have received a copy of the GNU General Public License\n"
"* along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
"*/\n"
"\n"
"#version 130\n"
"\n"
"uniform float iTime, iProgress;\n"
"uniform vec2 iResolution;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"// Hash function\n"
"void rand(in vec2 x, out float num)\n"
"{\n"
"    num = fract(sin(dot(x-1. ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"// Arbitrary-frequency 2D noise\n"
"void lfnoise(in vec2 t, out float num)\n"
"{\n"
"    vec2 i = floor(t);\n"
"    t = fract(t);\n"
"    //t = ((6.*t-15.)*t+10.)*t*t*t;  // TODO: add this for slower perlin noise\n"
"    t = smoothstep(c.yy, c.xx, t); // TODO: add this for faster value noise\n"
"    vec2 v1, v2;\n"
"    rand(i, v1.x);\n"
"    rand(i+c.xy, v1.y);\n"
"    rand(i+c.yx, v2.x);\n"
"    rand(i+c.xx, v2.y);\n"
"    v1 = c.zz+2.*mix(v1, v2, t.y);\n"
"    num = mix(v1.x, v1.y, t.x);\n"
"}\n"
"\n"
"// Multi-frequency 2D noise\n"
"void mfnoise(in vec2 x, in float fmin, in float fmax, in float alpha, out float num)\n"
"{\n"
"    num = 0.;\n"
"    float a = 1., nf = 0., buf;\n"
"    for(float f = fmin; f<fmax; f = f*2.)\n"
"    {\n"
"        lfnoise(f*x, buf);\n"
"        num += a*buf;\n"
"        a *= alpha;\n"
"        nf += 1.;\n"
"    }\n"
"    num *= (1.-alpha)/(1.-pow(alpha, nf));\n"
"}\n"
"\n"
"// Compute distance to regular polygon\n"
"void dpolygon(in vec2 x, in float N, out float d)\n"
"{\n"
"    d = 2.0*pi/N;\n"
"    float t = mod(acos(x.x/length(x)), d)-0.5*d;\n"
"    d = 0.5-length(x)*cos(t)/cos(0.5*d);\n"
"}\n"
"\n"
"// Distance to hexagon pattern\n"
"void dhexagonpattern(in vec2 p, out float d, out vec2 ind) \n"
"{\n"
"    vec2 q = vec2( p.x*1.2, p.y + p.x*0.6 );\n"
"    \n"
"    vec2 pi = floor(q);\n"
"    vec2 pf = fract(q);\n"
"\n"
"    float v = mod(pi.x + pi.y, 3.0);\n"
"\n"
"    float ca = step(1.,v);\n"
"    float cb = step(2.,v);\n"
"    vec2  ma = step(pf.xy,pf.yx);\n"
"    \n"
"    d = dot( ma, 1.0-pf.yx + ca*(pf.x+pf.y-1.0) + cb*(pf.yx-2.0*pf.xy) );\n"
"    ind = pi + ca - cb*ma;\n"
"    ind = vec2(ind.x/1.2, ind.y);\n"
"    ind = vec2(ind.x, ind.y-ind.x*.6);\n"
"}\n"
"\n"
"// Compute distance to stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0) - s;\n"
"}\n"
"\n"
"// Apply distance function to colors\n"
"void render(in vec3 c1, in vec3 c2, in float d, out vec3 col)\n"
"{\n"
"    col = mix(c1, c2, step(0.0,d));\n"
"}\n"
"\n"
"// Image function\n"
"void mainImage(out vec4 fragColor, in vec2 fragCoord)\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    \n"
"    float d, d0, d1, rn, p, ln;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), ind;\n"
"    vec3 col, col1, col2;\n"
"    \n"
"    uv *= 2.0;\n"
"    \n"
"    // Pattern\n"
"    dhexagonpattern((30.0)*uv, d, ind);\n"
"    stroke(d, 0.2, d);\n"
"    d = min(d, d0);\n"
"    \n"
"    // Make colors\n"
"    rand(ind, rn);\n"
"    col1 = mix((.85+.15*rn)*vec3(1.0, 0., 0.), vec3(1.0, 1.0, .2), clamp(rn, 0.0, 1.0));\n"
"    col = col1;\n"
"    \n"
"    // Big hexagon\n"
"    dpolygon(2.0*ind/66.0, 6.0, d0);\n"
"    stroke(d0, 0.2, d0);\n"
"    d0 = -d0;\n"
"    \n"
"    // Fill hexagon with colorful small hexagons that scale with progress\n"
"    mfnoise(ind, .5,100., .35, ln);\n"
"    ln = .5+.5*ln;\n"
"    p = atan(ind.y, ind.x);\n"
"    col2 = .2*col;//length(col)/length(c.xxx)*c.xxx;\n"
"    col = mix(c.yyy, col2, ln);\n"
"    col = mix(col, clamp(3.*col1,0.,1.), step(p, -pi+2.*pi*iProgress));\n"
"    render(c.yyy, col, d, col);\n"
"    render(c.yyy, col, d0, col);\n"
"    \n"
"    // Big hexagon border\n"
"    stroke(d0+.1, 0.05, d0);\n"
"    render(col, mix(c.yyy, mix(col1, col2, iProgress), .75+.25*ln), -d0, col);\n"
"    stroke(d, 0.1, d);\n"
"    render(col, c.yyy, d, col);\n"
"    \n"
"    // Add lightning\n"
"    if(d0 > .1 && length(ind/66.) > .4)\n"
"    {\n"
"        float p0 = mod(p, pi/3.)-pi/6., pi = (p-p0)*3./pi;\n"
"        float d3, off;\n"
"        rand(pi*c.xx, off);\n"
"        lfnoise(.5*length(ind)*c.xx+50.*off-5.*iTime, d3);\n"
"        d3 = p0-.15*d3;\n"
"        stroke(d3, .1, d3);\n"
"        render(col, .7*col1, -d3, col);\n"
"        stroke(d3, .05, d3);\n"
"        render(col, .2*col1, -d3, col);\n"
"    }\n"
"    render(col, c.yyy, d, col);\n"
"    \n"
"    // Small hexagons\n"
"    for(float i=0.; i<30.; i+=1.)\n"
"    {\n"
"        float size;\n"
"        vec2 velocity;\n"
"        rand(i*c.xx, size);\n"
"        size = .1+.5*size;\n"
"        rand((i+1.3)*c.xy, velocity.x);\n"
"        rand((i+2.1)*c.yx, velocity.y);\n"
"        velocity = -1.+2.*velocity;\n"
"        \n"
"        vec2 location = -vec2(a,1.)+vec2(a,2.)*mod(velocity*.5*iTime, vec2(a,1.));\n"
"        \n"
"        dpolygon((uv-location)/size, 6.0, d0);\n"
"        stroke(d0, 0.01, d);\n"
"        d = -d;\n"
"        render(col, mix(col, mix(c.xyy, c.xxy, .2+.5*velocity.x), .2), d0, col);\n"
"        render(col, mix(c.xyy, c.xxy, clamp(.3+.5*velocity.x,0.,1.)), d, col);\n"
"    }\n"
"    \n"
"    // Output to screen\n"
"    fragColor = vec4(col,1.0);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
