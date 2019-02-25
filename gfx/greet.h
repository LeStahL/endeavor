/* Generated with shader-compressor by NR4/Team210. */
#ifndef GREET_H
#define GREET_H
const char * greet_frag =
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
"uniform float iTime;\n"
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
"// Setup camera\n"
"void camerasetup(in vec2 uv, out vec3 ro, out vec3 dir)\n"
"{\n"
"    vec3 right = c.xyy, up = c.yxy, target = c.yyy+.05*vec3(cos(iTime), sin(iTime), 0.);\n"
"    ro = c.yyx;\n"
"    dir = normalize(target + uv.x * right + uv.y * up - ro);\n"
"}\n"
"\n"
"// Compute distance to stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0) - s;\n"
"}\n"
"\n"
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
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
"// 3D Effect on text in intro (210 logo)\n"
"void texteffect(in vec3 x, out vec2 sdf)\n"
"{\n"
"    sdf = vec2(x.z,7.);\n"
"    vec2 uv = x.xy;\n"
"    \n"
"    //uv -= mod(iTime, 1.5)*c.xy;\n"
"    //uv += mod(iTime,1.5)*c.yx;\n"
"    uv /= 2.+mod(iTime, 4.);\n"
"    \n"
"    float time_index = (iTime-mod(iTime, 4.))/4.;\n"
"    lfnoise(.5*iTime*c.xx, time_index);\n"
"    time_index *= pi/3.;\n"
"    uv = mat2(cos(time_index),-sin(time_index), sin(time_index), cos(time_index))*uv;\n"
"    \n"
"    \n"
"    for(float i=0.; i<6.; i+=1.)\n"
"    {\n"
"        float d, d0, res = 4.*pow(3., i), bres = res/3.;\n"
"        vec2 ind;\n"
"        dhexagonpattern(res*uv, d0, ind);\n"
"        d = -d0;\n"
"        stroke(d, 0.1, d);\n"
"        \n"
"        vec2 cind = ind/res;\n"
"        \n"
"        float big_hexagons, big_border;\n"
"        vec2 big_ind;\n"
"        dhexagonpattern(bres*cind, big_hexagons, big_ind);\n"
"        stroke(big_hexagons, 0.2, big_hexagons);\n"
"        big_border = big_hexagons;\n"
"        stroke(big_border, .1/bres, big_border);\n"
"        \n"
"        vec2 dt;\n"
"        lfnoise(res*15.0*cind+2.0, dt.x);\n"
"        lfnoise(res*15.0*cind+3.0, dt.y);\n"
"        dt *= 2.;\n"
"        float dm, dm2;\n"
"        lfnoise(50.0*cind, dm);\n"
"        dm = 0.5+0.5*dm;\n"
"        lfnoise(6.5*cind-dt-2.0*iTime*c.xx, dm2);\n"
"        // change sign here for different effect\n"
"        dm2 = clamp(0.5-0.5*dm2,0.,1.);\n"
"        \n"
"        d = mix(d,-d, dm);\n"
"        \n"
"        zextrude(x.z, big_border, dm2-.5*i/7., d);\n"
"        //stroke(d, .01, d);\n"
"        d -= .1;\n"
"        sdf.x = min(sdf.x, d);\n"
"    }\n"
"    \n"
"    // Add guard objects for debugging\n"
"    \n"
"    float dr = .05;\n"
"    vec3 y = mod(x,dr)-.5*dr;\n"
"    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));\n"
"    guard = abs(guard)+dr*.1;\n"
"    sdf.x = min(sdf.x, guard);\n"
"\n"
"}\n"
"\n"
"// Perform raymarching for actual object\n"
"void marchscene(in vec3 ro, in vec3 dir, in int N, in float eps, out vec3 x, out vec2 s, out float d, out bool flag)\n"
"{\n"
"    flag = false;\n"
"    for(int ia=0; ia<max(N,0); ++ia)\n"
"	{\n"
"        x = ro + d*dir;\n"
"        texteffect(x,s);\n"
"        if(s.x < eps)\n"
"        {\n"
"            flag = true;\n"
"            break;\n"
"        }\n"
"        d += s.x;\n"
"	}\n"
"}\n"
"\n"
"void calcnormal(in vec3 x, in float eps, out vec3 n)\n"
"{\n"
"    vec2 s, sp;\n"
"    texteffect(x, s);\n"
"    texteffect(x+eps*c.xyy, sp);\n"
"    n.x = sp.x-s.x;\n"
"    texteffect(x+eps*c.yxy, sp);\n"
"    n.y = sp.x-s.x;\n"
"    texteffect(x+eps*c.yyx, sp);\n"
"    n.z = sp.x-s.x;\n"
"    n = normalize(n);\n"
"}\n"
"\n"
"// 3D rotational matrix\n"
"void rot(in vec3 p, out mat3 m)\n"
"{\n"
"    m = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"// Initial intro\n"
"void background2(in vec2 uv, out vec3 col)\n"
"{\n"
"    uv /= 2.+mod(iTime, 4.);\n"
"    float time_index = (iTime-mod(iTime, 4.))/4., rt;\n"
"    rt = time_index;\n"
"    \n"
"    lfnoise(.5*iTime*c.xx, time_index);\n"
"    time_index *= pi/3.;\n"
"    uv = mat2(cos(time_index),-sin(time_index), sin(time_index), cos(time_index))*uv;\n"
"    \n"
"    mat3 m;\n"
"    rot(.1*iTime*vec3(1.1515,1.3321,1.5123) + rt, m);\n"
"    \n"
"    vec3 dark_green = abs(vec3(0.0,0.78,.52)),\n"
"        light_green = abs(vec3(.9, 0.98, 0.28)), \n"
"        gray = .5*length(light_green)*c.xxx/sqrt(3.);\n"
"    \n"
"    col = c.yyy;\n"
"    for(float i=0.; i<7.; i+=1.)\n"
"    {\n"
"        rot(.1*iTime*vec3(1.1515,1.3321,1.5123) + 3.*rt + i, m);\n"
"        \n"
"        float d, d0, res = 4.*pow(3., i), bres = res/3.;\n"
"        vec2 ind;\n"
"        dhexagonpattern(res*uv, d0, ind);\n"
"        d = -d0;\n"
"        stroke(d, 0.1, d);\n"
"        vec2 cind = ind/res;\n"
"        \n"
"        float big_hexagons, big_border;\n"
"        vec2 big_ind;\n"
"        dhexagonpattern(bres*cind, big_hexagons, big_ind);\n"
"        stroke(big_hexagons, 0.2, big_hexagons);\n"
"        big_border = big_hexagons;\n"
"        stroke(big_border, .1/bres, big_border);\n"
"        \n"
"        vec2 dt;\n"
"        lfnoise(res*15.0*cind+2.0, dt.x);\n"
"        lfnoise(res*15.0*cind+3.0, dt.y);\n"
"        dt *= 2.;\n"
"        float dm, dm2;\n"
"        lfnoise(50.0*cind, dm);\n"
"        dm = 0.5+0.5*dm;\n"
"        lfnoise(26.5*cind-dt-2.0*iTime*c.xx, dm2);\n"
"        dm2 = clamp(0.7+0.5*dm2, 0., 1.);\n"
"        \n"
"        light_green = mix(light_green, abs(m*m*m*m*light_green), dm2);\n"
"        dark_green = mix(c.yyy, abs(m*m*m*dark_green), dm2);\n"
"        \n"
"        gray = .5*length(light_green)*c.xxx/sqrt(3.);\n"
"        col += mix(light_green, dark_green, step(0.,big_border));\n"
"        //col = mix(col, light_green, step(big_border,0.));\n"
"        \n"
"        col = smoothstep(0.,12., iTime)*clamp(col*step(0.,d),0.,1.);\n"
"    }\n"
"    \n"
"    col += .2*dark_green;\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), s = c.xy;\n"
"\n"
"	vec3 ro, x, dir;\n"
"    \n"
"    float d = 0.;\n"
"    bool hit = false;\n"
"    \n"
"    vec3 col = c.yyy;\n"
"                \n"
"	camerasetup(uv, ro, dir);\n"
"    d = (.5-ro.z)/dir.z;\n"
"    marchscene(ro, dir, 500, 1.0e-4, x, s, d, hit);\n"
"    \n"
"    if(hit)\n"
"        background2(x.xy, col);\n"
"    else\n"
"        background2((ro-ro.z/dir.z*dir).xy, col);\n"
"    \n"
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
