/* Generated with shader-compressor by NR4/Team210. */
#ifndef FOURTWENTY_H
#define FOURTWENTY_H
const char * fourtwenty_frag =
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
"// compute distance to regular star\n"
"void dstar(in vec2 x, in float N, in vec2 R, out float dst)\n"
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
"    dst = dot(x-p1,ff)/length(ff);\n"
"}\n"
"\n"
"// 2D box\n"
"void dbox(in vec2 x, in vec2 b, out float d)\n"
"{\n"
"	vec2 da = abs(x)-b;\n"
"	d = length(max(da,c.yy)) + min(max(da.x,da.y),0.0);\n"
"}\n"
"\n"
"// Distance to circle\n"
"void dcircle(in vec2 x, in float R, out float d)\n"
"{\n"
"    d = length(x)-R;\n"
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
"void colorize(in vec3 old, in float ind, in vec2 uv, out vec3 col)\n"
"{\n"
"    col = old;\n"
"    float d=1., da, R = .5, dd = 1.;\n"
"    vec2 p0 = -.3*c.yx, dis;\n"
"    \n"
"    lfnoise(2.*uv.x*c.xx, dis.x);\n"
"    lfnoise(2.*uv.y*c.xx, dis.y);\n"
"    \n"
"    // Draw Leaves\n"
"    for(float phi=-pi/2.+pi/5.; phi<3.*pi/2.; phi += pi/5.)\n"
"    {\n"
"        vec2 p = R*vec2(cos(phi), sin(phi));\n"
"        dlinesegment(uv-.1*dis,p0,p,da);\n"
"        d = min(d,da);\n"
"        \n"
"        vec2 dad = p-p0;\n"
"        float dpd = dot(uv-p0, dad)/dot(dad,dad),\n"
"        	arg = 2.*pi*dpd*36.*length(dad),\n"
"            q = .3;\n"
"        stroke(da, (.08*length(dad)+.01-.005*sin(arg)/(1.+q*q-2.*q*cos(arg)))*smoothstep(0.,.5,dpd)*(1.-smoothstep(.5,1.,dpd)), da);\n"
"        dd = min(dd, da);\n"
"    }\n"
"    \n"
"    float scale;\n"
"    rand(ind*c.xx, scale);\n"
"    \n"
"    col = mix(col, (.5+scale)*vec3(0.24,0.54,0.01), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, dd)); \n"
"    stroke(dd,.002, dd);\n"
"    col = mix(col, (1.-scale)*vec3(0.16,0.40,0.05), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, dd)); \n"
"\n"
"    stroke(d,.002, d);\n"
"    col = mix(col, (1.-scale)*vec3(0.16,0.40,0.05), smoothstep(1.5/iResolution.y, -1.5/iResolution.y, d)); \n"
"    \n"
"    col = clamp(col, 0.,1.);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"    vec3 col = c.yyy;\n"
"    \n"
"    for(int i=0; i<80; ++i)\n"
"    {\n"
"        float phi;\n"
"        lfnoise(float(i)*c.xx-iTime, phi);\n"
"        mat2 RR = mat2(cos(phi),-sin(phi),sin(phi), cos(phi));\n"
"        float scale;\n"
"        rand(float(i)*c.xx+.1, scale);\n"
"        vec2 delta;\n"
"        rand(float(i)*c.xx+.3, delta.x);\n"
"        delta.x *= a;\n"
"        rand(float(i)*c.xx+.5, delta.y);\n"
"        delta -= .5*vec2(a,1.);\n"
"        colorize(col, float(i), (1.+3.*scale)*RR*(uv-delta), col);\n"
"    }\n"
"    col = clamp(col, 0.,1.);\n"
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
