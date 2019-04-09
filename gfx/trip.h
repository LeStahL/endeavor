/* Generated with shader-compressor by NR4/Team210. */
#ifndef TRIP_H
#define TRIP_H
const char * trip_frag =
"/*\n"
" * Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
" * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>\n"
" * \n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" * \n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" * \n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n"
" */\n"
"\n"
"#version 130\n"
"\n"
"uniform float iTime, iProgress;\n"
"uniform vec2 iResolution;\n"
"\n"
"// Global constants\n"
"const float pi = acos(-1.);\n"
"const vec3 c = vec3(1.0, 0.0, -1.0);\n"
"float a = 1.0;\n"
"\n"
"// Scene constants\n"
"const float R = .5;\n"
"\n"
"// Hash function\n"
"void rand(in vec2 x, out float num)\n"
"{\n"
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
"// 2D box\n"
"void dbox(in vec2 p, in vec2 b, out float dst)\n"
"{\n"
"  	vec2 d = abs(p) - b;\n"
"  	dst = length(max(d,0.0)) + min(max(d.x,d.y),0.0); \n"
"}\n"
"\n"
"// Stroke\n"
"void stroke(in float d0, in float s, out float d)\n"
"{\n"
"    d = abs(d0)-s;\n"
"}\n"
"\n"
"// Distance to line segment\n"
"void dlinesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"void dlinesegment3(in vec3 x, in vec3 p1, in vec3 p2, out float d)\n"
"{\n"
"    vec3 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
"}\n"
"\n"
"// Compute distance to regular polygon\n"
"void dpolygon(in vec2 x, in float N, in float R, out float d)\n"
"{\n"
"    d = 2.0*pi/N;\n"
"    float t = mod(acos(x.x/length(x)), d)-0.5*d;\n"
"    d = R-length(x)*cos(t)/cos(0.5*d);\n"
"}\n"
"\n"
"// compute distance to regular star\n"
"void dstar(in vec2 x, in float N, in vec2 R, out float ds)\n"
"{\n"
"    float d = pi/N,\n"
"        p0 = acos(x.x/length(x)),\n"
"        p = mod(p0, d),\n"
"        i = mod(round((p-p0)/d),2.);\n"
"    x = length(x)*vec2(cos(p),sin(p));\n"
"    vec2 a = mix(R,R.yx,i),\n"
"    	p1 = a.x*c.xy,\n"
"        ff = a.y*vec2(cos(d),sin(d))-p1;\n"
"   	ff = ff.yx*c.zx;\n"
"    ds = dot(x-p1,ff)/length(ff);\n"
"}\n"
"\n"
"void rot(in vec3 p, out mat3 rot)\n"
"{\n"
"    rot = mat3(c.xyyy, cos(p.x), sin(p.x), 0., -sin(p.x), cos(p.x))\n"
"        *mat3(cos(p.y), 0., -sin(p.y), c.yxy, sin(p.y), 0., cos(p.y))\n"
"        *mat3(cos(p.z), -sin(p.z), 0., sin(p.z), cos(p.z), c.yyyx);\n"
"}\n"
"\n"
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\n"
"void dcapcylinder(in vec3 x, in float R, in float h, out float dst)\n"
"{\n"
"    dst = length(x.xy)-R;\n"
"    stroke(dst, .05*R, dst);\n"
"    zextrude(x.z, -dst, h, dst);\n"
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
"void add(in vec2 sda, in vec2 sdb, out vec2 sdf)\n"
"{\n"
"    sdf = mix(sda, sdb, step(sdb.x, sda.x));\n"
"}\n"
"\n"
"void scene(in vec3 x, out vec2 sdf)\n"
"{\n"
"    float ra, tm = clamp((iTime-10.)/5., 0.,1.);\n"
"    lfnoise(iTime*c.xx, ra);\n"
"    vec2 dxy = mix(.04,.2,.5+.5*ra)*vec2(cos(x.z-4.*iTime),sin(x.z-4.*iTime));\n"
"    x.xy -= mix(c.yy, dxy, tm);\n"
"        mat3 r;\n"
"    rot(vec3(0.,0.,ra), r);\n"
"    x = mix(x, r*x, tm);\n"
"    \n"
"    // Rings\n"
"    vec3 y = vec3(x.xy, mod(x.z,.5)-.25), \n"
"        ind = (x-y)/.5;\n"
"    dcapcylinder(y, 1.1*R, .05*R, sdf.x);\n"
"    sdf.y = 1.;\n"
"    \n"
"    // Vertical stripes\n"
"    float phi = atan(x.y, x.x),\n"
"        p0 = mod(phi, pi/4.)-pi/8.,\n"
"        pp0 = mod(phi-pi/8., pi/4.)-pi/8.;\n"
"    vec2 pa = 1.*R*vec2(cos(phi-pp0), sin(phi-pp0)),\n"
"        sda;\n"
"    sda.x = length(x-vec3(pa,x.z))-.04*R;\n"
"    sda.y = 1.;\n"
"    add(sdf, sda, sdf);\n"
"    pa = 1.1*R*vec2(cos(phi-p0), sin(phi-p0));\n"
"    sda.x = length(x.xy-pa)-.03*R;\n"
"    sda.y = 1.;\n"
"    add(sdf, sda, sdf);\n"
"    // Random stripes\n"
"    float p1 = p0 - pi/4.;// (sin(ind.z))*pi/4.;\n"
"    vec2 pb = 1.1*R*vec2(cos(phi-p1), sin(phi-p1));\n"
"    dlinesegment3(y, vec3(pa,.0), vec3(pb,.5), sda.x);\n"
"    p1 = p0+pi/4.;\n"
"    pb = 1.1*R*vec2(cos(phi-p1), sin(phi-p1));\n"
"    float dl;\n"
"    dlinesegment3(y, vec3(pa,.0), vec3(pb,-.5), dl);\n"
"    sda.x = min(sda.x,dl);\n"
"    stroke(sda.x, .025, sda.x);\n"
"    add(sdf, sda, sdf);\n"
"    sda.x = length(y-vec3(pa,0.))-.05;\n"
"    add(sdf, sda, sdf);\n"
"    \n"
"    // Mountains / Floor\n"
"    dbox(x.xz, vec2(.6*R, 1.e5*R), sda.x);\n"
"    zextrude(x.y+.9*R, -sda.x, .3*R, sda.x);\n"
"    sda.y = 2.;\n"
"    add(sdf, sda, sdf);\n"
"    \n"
"    dbox(vec2(abs(x.x)-.6*R,x.z), vec2(.05*R, 1.e5*R), sda.x);\n"
"    zextrude(x.y+.7*R, -sda.x, .001*R, sda.x);\n"
"    sda.y = 2.;\n"
"    add(sdf, sda, sdf);\n"
"    \n"
"    // Honeycomb\n"
"    float dh;\n"
"    vec2 pd = vec2(atan(x.y,x.x), x.z),\n"
"        hind;\n"
"    dhexagonpattern(2.*pi*vec2(3.,6.)*pd,dh, hind);\n"
"    stroke(dh, .3, dh);\n"
"    zextrude(length(x.xy)-1.1*R, -dh, .01, sda.x);\n"
"    sda.y = 1.;\n"
"    add(sdf, sda, sdf);\n"
"    \n"
"    // Add rails\n"
"    dbox(vec2(abs(x.x)-.3*R,x.z), vec2(.02*R, 1.e5*R), sda.x);\n"
"    zextrude(x.y+.7*R, -sda.x, .005*R, sda.x);\n"
"    sda.y = 1.;\n"
"    add(sdf, sda, sdf);\n"
"    dbox(y.xz, vec2(.5*R, .1*R), sda.x);\n"
"    zextrude(x.y+.73*R, -sda.x, .03*R, sda.x);\n"
"    sda.y = 3.;\n"
"    add(sdf, sda, sdf);\n"
"                \n"
"    // Add guard objects for debugging\n"
"    float dr = .1;\n"
"    y = mod(x,dr)-.5*dr;\n"
"    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));\n"
"    guard = abs(guard)+dr*.1;\n"
"    sdf.x = min(sdf.x, guard);\n"
"\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 2.e-3;\n"
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
"void floor_pattern(in vec2 x, out vec3 col)\n"
"{\n"
"    float w = .02;\n"
"    vec2 y = mod(x, w)-.5*w;\n"
"    float d;\n"
"    dlinesegment(y, -.25*w*c.xx, .25*w*c.xx, d);\n"
"    stroke(d,.15*w, d);\n"
"    float dx;\n"
"    mfnoise(x, 10.,1000., .45, dx);\n"
"    col = mix(.4*c.xyy, .5*c.xyy, step(0.,d));\n"
"    col = mix(col, .4*c.xyy, .5+.5*dx);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    float tm = clamp((iTime-10.)/5., 0.,1.);\n"
"    float ra;\n"
"    mat3 r;\n"
"    lfnoise(iTime*c.xx, ra);\n"
"    rot(vec3(0.,0.,12.*ra), r);\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec2 dxy = mix(.04,.2,.5+.5*ra)*vec2(cos(uv.y-4.*iTime),sin(uv.y-4.*iTime));\n"
"    uv.xy -= mix(c.yy, dxy, tm);\n"
"    \n"
"    float N = 5.*(.5+.5*sin(1.e-3*iTime)), dd = 2.*pi/N*sin(.1*iTime);\n"
"    vec2 k = vec2(cos(dd),sin(dd)),\n"
"        tt = c.zx, e;\n"
"    mat2 Ra = mat2(k.x,k.y,-k.y,k.x);\n"
"    \n"
"    for(float i = 2.; i > .05; i = i*(.85+.1*sin(.3*iTime)))\n"
"    {\n"
"        vec2 uuv = mix(Ra*uv, floor(25.*Ra*uv)/25., clamp((iTime-15.)/5., 0.,1.));\n"
"        float dd;\n"
"        dstar(uuv, N, i*vec2(1.+.5*cos(3.4221*iTime),1.+.5*sin(2.153*iTime)), dd);\n"
"        e = vec2(dd,i);\n"
"        tt = mix(tt,e,step(-3./iResolution.y,e.x));\n"
"    }\n"
"    e = vec2(.025-length(uv),pi); \n"
"    uv = mix(uv, floor(100.*uv)/100., clamp((iTime-20.)/5., 0., 1.));\n"
"    \n"
"    vec3 col = mix(.3*c.xxx*cos(2.*length(uv)),.8 + .5*cos(1.+uv.xyx+tt.y*1.5e1+iTime+vec3(0.,2.,4.)),tm), \n"
"        o = c.yyx+(.04*iTime*iTime-1.e1)*c.yyx, dir = normalize(vec3(uv,1.));\n"
"    mat3 RR;\n"
"                        rot(vec3(1.1,1.4,1.6)*iTime ,RR);\n"
"                        col = abs(RR*col);\n"
"    col = clamp(col, 0., 1.);\n"
"    float d;\n"
"    vec2 s;\n"
"    vec3 x;\n"
"    \n"
"    fragColor = vec4(col, 1.);\n"
"    //if(length(uv) < 1.5*R)\n"
"    {\n"
"        // Analytical distance to infinite cylinder\n"
"        d =.6*R/length(dir.xy); \n"
"		float dout = 1.15*R/length(dir.xy);\n"
"\n"
"        x = o + d * dir;\n"
"\n"
"        // Outer cylindrical structure\n"
"        int N = 300, i;\n"
"        for(i=0; i<N; ++i)\n"
"        {\n"
"            x = o + d * dir;\n"
"            scene(x, s);\n"
"            if(d > dout) \n"
"                return;\n"
"            if(s.x < 1.e-3) break;\n"
"            d += s.x;\n"
"        }\n"
"        if(i < N)\n"
"        {\n"
"            // Calc normal\n"
"            vec3 n;\n"
"            normal(x, n);\n"
"\n"
"            // Colorize\n"
"            vec3 l = normalize(x+c.yyx*sin(x.z));\n"
"            if(s.y == 1.)\n"
"            {\n"
"                col = .1*c.yyx\n"
"                    + .3*c.yyx*abs(dot(l, n))\n"
"                    + .9*c.xyx*pow(abs(dot(reflect(-l,n), dir)), 2.);\n"
"                \n"
"            }\n"
"            else if(s.y == 2.)\n"
"            {\n"
"                floor_pattern(x.xz, col);\n"
"                \n"
"                dir = n;\n"
"                o = x;\n"
"                d = 2.e-4;\n"
"                for(i=0; i<500; ++i)\n"
"                {\n"
"                    x = o + d * dir;\n"
"                    scene(x, s);\n"
"                    if(d > dout) \n"
"                    {\n"
"    					col = mix(col, .2*c.xxx, tanh(.1*d));\n"
"                        \n"
"                        vec3 bw = length(col)*c.xxx;\n"
"                        col = mix(bw, col, tm);\n"
"                        fragColor = vec4(col,1.0);\n"
"                        mat3 RR;\n"
"                        rot(vec3(1.1,1.4,1.6)*iTime ,RR);\n"
"                        col = abs(RR*col);\n"
"                        return;\n"
"                    }\n"
"                    if(s.x < 1.e-4) break;\n"
"                    d += s.x;\n"
"                }\n"
"                \n"
"                if(i<N)\n"
"                {\n"
"                    normal(x, n);\n"
"                    vec3 c1 = .1*c.xyx\n"
"                        + .3*c.xyx*abs(dot(l, n))\n"
"                        + .9*c.xyx*pow(abs(dot(reflect(-l,n), dir)), 2.);\n"
"                    col = mix(col, c1, .3);\n"
"                    \n"
"                }\n"
"                \n"
"            }\n"
"            else if(s.y == 3.)\n"
"            {\n"
"                col = .1*c.yxy\n"
"                    + .2*c.yxy*abs(dot(l, n))\n"
"                    + .8*c.yxy*pow(abs(dot(reflect(-l,n), dir)), 2.);\n"
"                float da, db;\n"
"                lfnoise(4.*x.x*c.xx, db);\n"
"                mfnoise(x.z*c.xx-.03*db*c.xy, 1.e2, 1.e4, .45, da);\n"
"                col = mix(col, .4*col, .5+.5*da);\n"
"                mfnoise(x.y*c.xx-.03*db*c.xy, 1.e2, 1.e4, .45, da);\n"
"                col = mix(col, .4*col, .5+.5*da);\n"
"            }\n"
"        }\n"
"    }\n"
"    \n"
"    //float daa = 4.*abs(s.x-.4)-1.3;\n"
"    //col = mix(col, c.xyx, clamp(daa, 0., 1.));\n"
"	rot(vec3(1.1,1.4,1.6)*iTime ,RR);\n"
"	col = abs(RR*col);\n"
"    vec3 bw = .8*length(col*col)*c.xxx;\n"
"    col = mix(bw, col, tm);\n"
"    \n"
"    col = mix(col, c.yyy, tanh(2.e-1*(x.z-o.z)));\n"
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