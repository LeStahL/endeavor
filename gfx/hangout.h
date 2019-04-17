/* Generated with shader-compressor by NR4/Team210. */
#ifndef HANGOUT_H
#define HANGOUT_H
const char * hangout_frag =
"/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" *\n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
" */\n"
" \n"
"#version 130\n"
"\n"
"uniform float iTime;\n"
"uniform vec2 iResolution;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"// Hash function\n"
"void rand(in vec2 x, out float num)\n"
"{\n"
"    x += 400.;\n"
"    num = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
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
"void rot(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"// Distance to regular voronoi\n"
"void dvoronoi(in vec2 x, out float d, out vec2 ind)\n"
"{\n"
"    vec2 y = floor(x);\n"
"   	float ret = 1.;\n"
"    \n"
"    //find closest control point. (\"In which cell am I?\")\n"
"    vec2 pf=c.xx, p;\n"
"    float df=100., oo;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            vec2 pa;\n"
"            rand(p, pa.x);\n"
"            rand(p+.1,pa.y);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            df = min(d,df);\n"
"            pf = mix(pf,p,step(d,df));\n"
"        }\n"
"    \n"
"    //compute voronoi distance: minimum distance to any edge\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            vec2 pa;\n"
"            rand(p, pa.x);\n"
"            rand(p+.1,pa.y);\n"
"            p += pa;\n"
"\n"
"            vec2 o = p - pf;\n"
"            oo = length(o);\n"
"            if(oo < 1.e-4) continue;\n"
"            d = abs(.5*oo-dot(x-pf, o)/oo);\n"
"            ret = min(ret, d);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    ind = pf;\n"
"}\n"
"\n"
"// 2D box\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"	vec2 da = abs(x)-b;\n"
"	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\n"
"// Add distance functions with materials\n"
"void add(in vec2 sda, in vec2 sdb, out vec2 dst)\n"
"{\n"
"    dst = mix(sda, sdb, step(sdb.x, sda.x));\n"
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
"vec2 ind;\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    x.z *= 3.;\n"
"    \n"
"    sdf.y = 1.;\n"
"    mfnoise(x.xy, 1., 400., .45, sdf.x);\n"
"    stroke(sdf.x, .5, sdf.x);\n"
"    sdf.x = x.z+.1*sdf.x;\n"
"    \n"
"    float R = .07+.1*sdf.x, dis;\n"
"    lfnoise(.5*x.y*c.xx, dis);\n"
"    vec2 sdb;\n"
"    //if(abs(x.x) < 1.)\n"
"    {\n"
"        vec2 ya = abs(vec2(x.x,x.z)-.4*dis*c.xy)-.015*c.yx + .005*c.xy,\n"
"            yi = ya;\n"
"        \n"
"        ya = vec2(mod(ya.x, mix(2.,.2,smoothstep(0.,3.,x.y)))-.5*mix(2.,.2,smoothstep(0.,3.,x.y)), ya.y);\n"
"        yi -= ya;\n"
"        yi /= mix(2.,.2,smoothstep(0.,3.,x.y));\n"
"        float da;\n"
"        if(yi.x < 4.)\n"
"        {\n"
"            zextrude(x.y, -length(ya)+R, 1.e4, sdb.x);\n"
"            stroke(sdb.x, .003, da);\n"
"        } else da = 1.e-3;\n"
"        sdb.y = 2.;\n"
"\n"
"        float phi = atan(ya.y, ya.x);\n"
"        dhexagonpattern(vec2(56.,12.)*vec2(x.y, phi), dis, ind);\n"
"        stroke(dis, .2, dis);\n"
"        stroke(sdb.x, .003*step(dis,0.), sdb.x);\n"
"        sdf.x = max(sdf.x,-da);\n"
"\n"
"        // Add Station\n"
"        R = .9;\n"
"        vec2 sdc;\n"
"        zextrude(x.y-3., -length(x.xz)+R, 1.5, sdc.x);\n"
"        stroke(sdc.x, .05, sdc.x);\n"
"        sdc.x = mix(sdc.x, 1.e-3, step(-dis,0.));\n"
"\n"
"        // Add normal walls\n"
"        float dd = x.y - 2.7, dis1;\n"
"        stroke(dd, .02, dd);\n"
"        vec2 ind0;\n"
"        dhexagonpattern(vec2(56.,22.)*x.xz, dis1, ind);\n"
"        stroke(dis1, .2, dis1);\n"
"        dd = mix(dd, 1.e-3, step(-dis1,0.));\n"
"        sdc.x = min(sdc.x, dd);\n"
"\n"
"\n"
"        // Remove exterior\n"
"        zextrude(x.y-3., -length(x.xz)+1.02*R, 1.05*.5, dd);\n"
"        stroke(dd, .03, dd);\n"
"        sdc.x = max(sdc.x, dd);\n"
"\n"
"        // Remove interior\n"
"        zextrude(x.y-3.,-length(x.xz)+.95*R, 1.*.5, dd);\n"
"        stroke(dd, .03, dd);\n"
"        sdc.x = max(sdc.x, -dd);\n"
"\n"
"        // Add Floor\n"
"        float fl;\n"
"        dbox(x.xy-3.*c.yx, vec2(R,.6*.5), fl);\n"
"        zextrude(x.z,-fl,.1, fl);\n"
"        sdc.x = min(sdc.x, fl);\n"
"        \n"
"        sdb.x = max(sdb.x, -dd);\n"
"        add(sdf, sdb, sdf);\n"
"        sdc.y = 2.;\n"
"        add(sdf, sdc, sdf);\n"
"        \n"
"        // Add guard objects for debugging\n"
"        float dr = .1;\n"
"        vec3 y = mod(x,dr)-.5*dr;\n"
"        float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));\n"
"        guard = abs(guard)+dr*.1;\n"
"        sdf.x = min(sdf.x, guard);\n"
"    }\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 5.e-3;\n"
"    vec2 s, na;\n"
"    \n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, na);\n"
"    n.x = na.x;\n"
"    scene(x+dx*c.yxy, na);\n"
"    n.y = na.x;\n"
"    scene(x+dx*c.yyx, na);\n"
"    n.z = na.x;\n"
"    n = normalize(n-s.x);\n"
"}\n"
"\n"
"void planet_texture(in vec2 x, out vec3 col)\n"
"{\n"
"    vec3 light_orange = vec3(1.00,0.69,0.05),\n"
"        orange = vec3(0.95,0.45,0.01),\n"
"        dark_orange = vec3(0.98,0.73,0.01);\n"
"    \n"
"    //rock like appearance\n"
"    float d;\n"
"    mfnoise(x, 1.,400.,.65,d);\n"
"	col = mix(vec3(0.19,0.02,0.00), vec3(0.91,0.45,0.02), .5+.5*d);\n"
"    \n"
"    // big structures\n"
"    stroke(d,.01, d);\n"
"    col = mix(mix(.8*vec3(0.99,0.49,0.02),c.yyy,d*clamp(.2-.5*x.y/12.,0.,1.)), col, smoothstep(-.05,.05,d));\n"
"    col = mix(col, vec3(0.15,0.05,0.00), clamp(.2-.5*x.y/12.,0.,1.));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0),\n"
"        s;\n"
"    vec3 col = c.yyy, \n"
"        o = .5*c.yyx + mix(.3*iTime*c.yxy,.3*6.*c.yxy,smoothstep(0.,6.,clamp(iTime,0.,6.))) \n"
"        	+ smoothstep(6.,12.,clamp(iTime, 6., 12.))*c.yyx*.15, \n"
"        t = vec3(uv,0.) + mix(.3*iTime*c.yxy,.3*6.*c.yxy,smoothstep(0.,6.,clamp(iTime,0.,6.))) \n"
"        	+ smoothstep(6.,12.,clamp(iTime, 6., 12.))*c.yyx*.15	\n"
"        	+ c.yxy,\n"
"        dir = normalize(t-o),\n"
"        x;\n"
"    float d = 0.;//(.14-o.z)/dir.z;\n"
"    \n"
"    int N = 450, i;\n"
"    for(i=0; i<N; ++i)\n"
"    {\n"
"        x = o + d * dir;\n"
"        scene(x,s);\n"
"        if(s.x < 1.e-4) break;\n"
"        d += s.x;\n"
"    }\n"
"    \n"
"    if(i<N)\n"
"    {\n"
"        vec3 n,\n"
"            l = normalize(x+c.yyx);\n"
"        normal(x,n);\n"
"        vec3 c1 = vec3(0.99,0.64,0.02);\n"
"        if(s.y == 1.)\n"
"        {\n"
"            planet_texture(x.xy, c1);\n"
"            col = .3*c1\n"
"                + (.3*c1)*abs(dot(l,n))\n"
"                + (1.3*c1+.1*c.xyy)*pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"        }\n"
"        else if(s.y == 2.)\n"
"        {\n"
"            planet_texture(x.xy, c1);\n"
"            col = .3*c1\n"
"                + (.3*c1)*abs(dot(l,n))\n"
"                + (1.3*c1+.1*c.xyy)*pow(abs(dot(reflect(-l,n),dir)),3.);\n"
"            col = mix(col,.4*length(col)*c.xyy,.7);\n"
"            mat3 RR;\n"
"            rot(.01*ind.x*c.xxy, RR);\n"
"            col = abs(RR * col);\n"
"        }\n"
"    }\n"
"    /*\n"
"    vec3 ddd;\n"
"    rand(uv-iTime*c.xx, ddd.x);\n"
"    rand(uv-iTime*c.xx, ddd.y);\n"
"    rand(uv-iTime*c.xx, ddd.z);\n"
"    col -=.1*ddd;\n"
"    */\n"
"    col = clamp(col, 0.,1.);\n"
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
