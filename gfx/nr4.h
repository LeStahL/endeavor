/* Generated with shader-compressor by NR4/Team210. */
#ifndef NR4_H
#define NR4_H
const char * nr4_frag =
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
"    num = fract(sin(dot(abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
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
"// 2D box\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"	vec2 da = abs(x)-b;\n"
"	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"// Compute distance to regular polygon\n"
"void dpolygon(in vec2 x, in float N, out float d)\n"
"{\n"
"    d = 2.0*pi/N;\n"
"    float t = mod(acos(x.x/length(x)), d)-0.5*d;\n"
"    d = -0.5+length(x)*cos(t)/cos(0.5*d);\n"
"}\n"
"\n"
"// Distance to circle\n"
"void dcircle(in vec2 x, out float d)\n"
"{\n"
"    d = length(x)-1.0;\n"
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
"// 2D rotational matrix\n"
"void rot(in float phi, out mat2 m)\n"
"{\n"
"    vec2 cs = vec2(cos(phi), sin(phi));\n"
"    m = mat2(cs.x, -cs.y, cs.y, cs.x);\n"
"}\n"
"\n"
"// Distance to pig ear\n"
"void dear(in vec2 x, out float d)\n"
"{\n"
"    d = abs(2.*x.y)\n"
"        -.95+smoothstep(0.,.5,clamp(abs(x.x),0.,1.))\n"
"        -.5*min(-abs(x.x),.01);\n"
"}\n"
"\n"
"// Distance to a triangle\n"
"void dtriangle(in vec2 x, in vec2 p0, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 d1 = c.xz*(p1-p0).yx, d2 = c.xz*(p2-p1).yx, d3 = c.xz*(p0-p2).yx;\n"
"    d = -min(\n"
"        dot(p0-x,d1)/length(d1),\n"
"        min(\n"
"            dot(p1-x,d2)/length(d2),\n"
"            dot(p2-x,d3)/length(d3)\n"
"        )\n"
"    );\n"
"}\n"
"\n"
"// Distance to circle segment\n"
"void circlesegment(in vec2 x, in float r, in float p0, in float p1, out float d)\n"
"{\n"
"    float p = atan(x.y, x.x);\n"
"    vec2 philo = vec2(max(p0, p1), min(p0, p1));\n"
"    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))\n"
"    {\n"
"        d = abs(length(x)-r);\n"
"        return;\n"
"    }\n"
"    d = min(\n"
"        length(x-r*vec2(cos(p0), sin(p0))),\n"
"        length(x-r*vec2(cos(p1), sin(p1)))\n"
"        );\n"
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
"// Extrusion\n"
"void zextrude(in float z, in float d2d, in float h, out float d)\n"
"{\n"
"    vec2 w = vec2(-d2d, abs(z)-0.5*h);\n"
"    d = length(max(w,0.0));\n"
"}\n"
"\n"
"// Distance to nr4\n"
"void dnr4(in vec2 x, out float d)\n"
"{\n"
"    circlesegment(x+.586*c.xy, .18, -3.*pi/8., 11.*pi/8., d);\n"
"    float d0;\n"
"    circlesegment(x-.01*c.xy, .18, pi/8., pi, d0);\n"
"    d = min(d,d0);\n"
"    dlinesegment(x-.01*c.xy, .18*vec2(-1.,0.), .18*vec2(-1.,-1.), d0);\n"
"    d = min(d,d0);\n"
"    circlesegment(x-.575*c.xy, .18, -11.*pi/8., 0., d0);\n"
"    d = min(d,d0);\n"
"    dlinesegment(x-.575*c.xy, .18*vec2(1.,0.), .18*vec2(1.,-1.), d0);\n"
"    d = min(d,d0);\n"
"}\n"
"\n"
"// Distance to regular voronoi\n"
"void dvoronoi(in vec2 x, out float d, out vec2 ind)\n"
"{\n"
"    vec2 y = floor(x);\n"
"   	float ret = 1.;\n"
"    \n"
"    //find closest control point. (\"In which cell am I?\")\n"
"    vec2 pf=c.yy, p;\n"
"    float df=10.;\n"
"    \n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            d = length(x-p);\n"
"            \n"
"            if(d < df)\n"
"            {\n"
"                df = d;\n"
"                pf = p;\n"
"            }\n"
"        }\n"
"    \n"
"    //compute voronoi distance: minimum distance to any edge\n"
"    for(int i=-1; i<=1; i+=1)\n"
"        for(int j=-1; j<=1; j+=1)\n"
"        {\n"
"            p = y + vec2(float(i), float(j));\n"
"            float pa;\n"
"            rand(p, pa);\n"
"            p += pa;\n"
"            \n"
"            vec2 o = p - pf;\n"
"            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);\n"
"            ret = min(ret, d);\n"
"        }\n"
"    \n"
"    d = ret;\n"
"    ind = pf;\n"
"}\n"
"\n"
"void scene(in vec3 x, out float d)\n"
"{\n"
"        vec2 dis = c.yy;\n"
"        \n"
"        if(iTime > 6.) // Distort\n"
"        {\n"
"            lfnoise(2.*x.xy-iTime*c.xy, dis.x);\n"
"            lfnoise(2.*x.xy-iTime*c.yx, dis.y);\n"
"            \n"
"            vec2 da;\n"
"            lfnoise(26.*x.xy-iTime*c.xy, da.x);\n"
"            lfnoise(26.*x.xy-iTime*c.yx, da.y);\n"
"            dis += .5*da;\n"
"        }\n"
"    dnr4(x.xy+mix(c.yy,.1*dis,clamp((iTime-6.)/2.,0.,1.)),d);\n"
"    stroke(d, .12, d);\n"
"\n"
"    d = -d;\n"
"	zextrude(x.z+1.1*clamp((iTime-10.)/2.,0.,1.),d,.3*clamp((iTime-10.)/2.,0.,1.),d);\n"
"    \n"
"    // Add guard objects for debugging\n"
"    float dr = .2;\n"
"    vec3 y = mod(x,dr)-.5*dr;\n"
"    float guard = -length(max(abs(y)-vec3(.5*dr*c.xx, .6),0.));\n"
"    guard = abs(guard)+dr*.1;\n"
"    d = min(d, guard);\n"
"}\n"
"\n"
"void normal(in vec3 x, out vec3 n)\n"
"{\n"
"    const float dx = 1.e-4;\n"
"    float s;\n"
"    scene(x,s);\n"
"    scene(x+dx*c.xyy, n.x);\n"
"    scene(x+dx*c.yxy, n.y);\n"
"    scene(x+dx*c.yyx, n.z);\n"
"    n = normalize(n-s);\n"
"}\n"
"\n"
"// Distance to regular star\n"
"void dstar(in vec2 x, in float N, in vec2 R, out float dst)\n"
"{\n"
"    float d = pi/N,\n"
"        p0 = atan(x.y,x.x),\n"
"        p = mod(p0, d),\n"
"        i = mod(round((p-p0)/d),2.);\n"
"    x = length(x)*vec2(cos(p),sin(p));\n"
"    vec2 aa = mix(R,R.yx,i),\n"
"    	p1 = aa.x*c.xy,\n"
"        ff = aa.y*vec2(cos(d),sin(d))-p1;\n"
"   	ff = ff.yx*c.zx;\n"
"    dst = dot(x-p1,ff)/length(ff);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0), ind;\n"
"    vec3 col = c.yyy;\n"
"    vec3 ca = length(vec3(1.,0.27,0.))/sqrt(3.)*c.xxx;\n"
"    \n"
"    float d4, dborder, dborder2, dc;\n"
"    if(iTime > 0.) // Add grid after 2 seconds\n"
"    {\n"
"        vec2 y = mod(uv, .02)-.01,\n"
"            a = abs(y)-.0005;\n"
"        vec3 newcol = mix(mix(c.yyy, .1*vec3(1.,0.27,0.), sin(pi/2.*length(uv))), .5*vec3(1.,0.27,0.), step(a.x,0.));\n"
"        newcol = mix(newcol, .5*vec3(1.,0.27,0.), step(a.y,0.));\n"
"        y = mod(uv-.01,.08)-.04;\n"
"        a = abs(y)-.001;\n"
"        newcol = mix(newcol, vec3(1.,0.27,0.), step(a.x,0.));\n"
"        newcol = mix(newcol, vec3(1.,0.27,0.), step(a.y,0.));\n"
"        col = mix(col, newcol, clamp(iTime/2.,0.,1.)-clamp((iTime-14.)/2.,0.,1.));\n"
"    }\n"
"    if(iTime > 12.) // Add background voronoi\n"
"    {\n"
"        float db, dc;\n"
"        vec2 ind;\n"
"		dvoronoi(16.*uv,db, ind);\n"
"		db = mix(1.,db,clamp((iTime-12.)/2.,0.,1.));\n"
"        stroke(db,.1,db);\n"
"        rand(ind, dc);\n"
"        dc = mix(0., dc, clamp((iTime-12.)/2.,0.,1.));\n"
"            ca = mix(mix(ca,clamp(1.4*length(vec3(1.,dc*0.47,dc*.1))/sqrt(3.)*c.xxx*.5,0.,1.),clamp((iTime-12.)/2.,0.,1.)), c.yyy, step(db,0.));\n"
"            ca = mix(ca,mix(ca, length(vec3(.3,.01, .01))/sqrt(3.)*.1*c.xxx, clamp(abs(ind.y/8.),0.,1.)), clamp((iTime-12.)/2.,0.,1.));\n"
"        col = mix(col, mix(col, ca, step(-d4,0.)),clamp((iTime-12.)/2.,0.,1.));\n"
"    }\n"
"    vec2 dis = c.yy;\n"
"    if(iTime > 2.) // Add nr4 after another 2\n"
"    {\n"
"        \n"
"        \n"
"        if(iTime > 6.) // Distort\n"
"        {\n"
"            lfnoise(2.*uv-iTime*c.xy, dis.x);\n"
"            lfnoise(2.*uv-iTime*c.yx, dis.y);\n"
"            \n"
"            vec2 da;\n"
"            lfnoise(26.*uv-iTime*c.xy, da.x);\n"
"            lfnoise(26.*uv-iTime*c.yx, da.y);\n"
"            dis += .5*da;\n"
"        }\n"
"        dnr4(uv+mix(c.yy,.1*dis,clamp((iTime-6.)/2.,0.,1.)), d4);\n"
"        stroke(d4, mix(.01,.12,clamp((iTime-2.)/2.,0.,1.)), d4);\n"
"        ca = vec3(1.,0.27,0.);\n"
"        \n"
"        if(iTime > 8.) // Mix color to voronoi profile\n"
"        {\n"
"            float db;\n"
"            vec2 ind;\n"
"            dvoronoi(16.*uv+dis,db, ind);\n"
"			db = mix(1.,db,clamp((iTime-8.)/2.,0.,1.));\n"
"            stroke(db,.1,db);\n"
"            rand(ind, dc);\n"
"            dc = mix(0., dc, clamp((iTime-8.)/2.,0.,1.));\n"
"            ca = mix(mix(ca,clamp(1.4*vec3(1.,dc*0.47,dc*.1),0.,1.),clamp((iTime-8.)/2.,0.,1.)), clamp(0.*vec3(.48,0.27,.1),0.,1.), step(db,0.));\n"
"            ca = mix(ca,mix(ca, vec3(.3,.01, .01), clamp(abs(ind.y/4.),0.,1.)), clamp((iTime-8.)/2.,0.,1.));\n"
"        }\n"
"        col = mix(col, ca, step(d4,0.));\n"
"    }\n"
"    if(iTime > 4.) // Add border\n"
"    {\n"
"        stroke(d4,mix(.01,.035,clamp((iTime-4.)/2.,0.,1.)),dborder);\n"
"        col = mix(col, clamp(.0*vec3(.68,0.17,.17),0.,1.), step(dborder, 0.));\n"
"    	stroke(dborder,mix(.001,.0035,clamp((iTime-4.)/2.,0.,1.)),dborder2);\n"
"        col = mix(col, clamp(vec3(.91,0.51,.0),0.,1.), step(dborder2, 0.));\n"
"    }\n"
"    \n"
"    col = clamp(col, 0., 1.);\n"
"\n"
"    if(iTime > 10.) // Yes, add that 3d already!\n"
"    {\n"
"        float s = 0., depth;\n"
"        vec3 o = c.yyx-1.*c.yxy + .5*mix(0.,.1, clamp((iTime-10.)/2.,0.,1.))*vec3(cos(iTime), sin(iTime),0.);\n"
"		vec3 dir = normalize(c.yyz+uv.x*c.xyy+uv.y*c.yxy-o), x,\n"
"	        l = normalize(c.yxx), n;\n"
"       	int i;\n"
"        depth = (-1.-o.z)/dir.z;\n"
"        for(i=0; i<350; ++i)\n"
"        {\n"
"            x = o+depth * dir;\n"
"            scene(x, s);\n"
"            if(s<1.e-4) break;\n"
"            depth += s;\n"
"            if(x.z<-1.3)break;\n"
"        }\n"
"        if(x.z > -1.29 && x.z < -1.01)\n"
"        {\n"
"            normal(x,n);\n"
"            col = .6*mix(c.yyy,c.xyy, clamp((x.z+1.15)/.15,0.,1.))*vec3(1.,0.27,0.) \n"
"                + .1*c.xxy*abs(dot(l,n))\n"
"                + vec3(1.,0.27,0.)*pow(abs(dot(reflect(-l,n),dir)), 3.);\n"
"            col = clamp(col, 0., 1.);\n"
"        }\n"
"    }\n"
"    \n"
"    if(iTime > 14.) \n"
"    {\n"
"        float d;\n"
"        vec2 y = mod(uv, .2)-.1, yi = floor((uv-y)/.2);\n"
"        dnr4((yi-.1)*.2,d);\n"
"        stroke(d,.15,d);\n"
"        if(d<0.)\n"
"        {\n"
"            vec2 da;\n"
"            lfnoise(16.*yi*.2-5.e0*iTime*c.xy, da.x);\n"
"            lfnoise(16.*yi*.2-6.3e0*iTime*c.yx, da.y);\n"
"            da = .5+.5*da;\n"
"            \n"
"            float dg;\n"
"            vec2 ri;\n"
"            rand(yi+.1, ri.x);\n"
"            rand(yi+.2, ri.y);\n"
"            if(da.x < .5 && da.y <.7)\n"
"            {\n"
"                mat2 m;\n"
"                rot(ri.y-iTime,m);\n"
"                dstar(m*(y-.05*da+.1*ri), 5.+floor(4.*ri.y), .05*ri.x*vec2(.01, .04), dg);\n"
"                col = mix(col, mix(col,c.xxy,.5), step(-dg,0.));\n"
"                stroke(dg, .0025, dg);\n"
"                col = mix(col, c.xxy, step(dg,0.));\n"
"            }\n"
"            lfnoise(11.*yi*.2-5.e0*iTime*c.xy, da.x);\n"
"            lfnoise(5.*yi*.2-6.3e0*iTime*c.yx, da.y);\n"
"            da = .5+.5*da;\n"
"            rand(yi+.3, ri.x);\n"
"            rand(yi+.4, ri.y);\n"
"            if(da.x < .5 && da.y <.7)\n"
"            {\n"
"                mat2 m;\n"
"                rot(ri.y-iTime,m);\n"
"                dstar(m*(y-.05*da+.1*ri), 5.+floor(4.*ri.y), .05*ri.x*vec2(.01, .04), dg);\n"
"                col = mix(col, mix(col,c.xxy,.5), step(-dg,0.));\n"
"                stroke(dg, .0025, dg);\n"
"                col = mix(col, c.xxy, step(dg,0.));\n"
"            }\n"
"        }\n"
"    }\n"
"    \n"
"    vec3 rd;\n"
"    rand(uv-iTime,rd.x);\n"
"    rand(uv*2.-iTime,rd.y);\n"
"    rand(uv*3.-iTime,rd.z);\n"
"    col += .2*(-c.xxx+rd);\n"
"    \n"
"    fragColor = vec4(col,1.0);\n"
"}\n"
"\n"
"\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
