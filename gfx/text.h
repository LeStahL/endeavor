/* Generated with shader-compressor by NR4/Team210. */
#ifndef TEXT_H
#define TEXT_H
const char * text_frag =
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
"uniform float iFontWidth, iExecutableSize, iTime;\n"
"uniform vec2 iResolution;\n"
"uniform sampler2D iChannel0, iFont;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"// Global constants\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"const float pi = acos(-1.);\n"
"float a; // Aspect ratio\n"
"\n"
"// Read short value from texture at index off\n"
"void rshort(in float off, out float val)\n"
"{\n"
"    // Parity of offset determines which byte is required.\n"
"    float hilo = mod(off, 2.);\n"
"    // Find the pixel offset your data is in (2 unsigned shorts per pixel).\n"
"    off *= .5;\n"
"    // - Determine texture coordinates.\n"
"    //     offset = i*iFontWidth+j for (i,j) in [0,iFontWidth]^2\n"
"    //     floor(offset/iFontWidth) = floor((i*iFontwidth+j)/iFontwidth)\n"
"    //                              = floor(i)+floor(j/iFontWidth) = i\n"
"    //     mod(offset, iFontWidth) = mod(i*iFontWidth + j, iFontWidth) = j\n"
"    // - For texture coordinates (i,j) has to be rescaled to [0,1].\n"
"    // - Also we need to add an extra small offset to the texture coordinate\n"
"    //   in order to always \"hit\" the right pixel. Pixel width is\n"
"    //     1./iFontWidth.\n"
"    //   Half of it is in the center of the pixel.\n"
"    vec2 ind = (vec2(mod(off, iFontWidth), floor(off/iFontWidth))+.05)/iFontWidth;\n"
"    // Get 4 bytes of data from the texture\n"
"    vec4 block = texture(iFont, ind);\n"
"    // Select the appropriate word\n"
"    vec2 data = mix(block.rg, block.ba, hilo);\n"
"    // Convert bytes to unsigned short. The lower bytes operate on 255,\n"
"    // the higher bytes operate on 65280, which is the maximum range \n"
"    // of 65535 minus the lower 255.\n"
"    val = round(dot(vec2(255., 65280.), data));\n"
"}\n"
"\n"
"// Read float value from texture at index off\n"
"void rfloat(in float off, out float val)\n"
"{\n"
"    // Convert the bytes to unsigned short as first step.\n"
"    float d;\n"
"    rshort(off, d);\n"
"    \n"
"    // Convert bytes to IEEE 754 float16. That is\n"
"    // 1 sign bit, 5 bit exponent, 11 bit mantissa.\n"
"    // Also it has a weird conversion rule that is not evident at all.\n"
"    float sign = floor(d/32768.),\n"
"        exponent = floor(d/1024.-sign*32.),\n"
"        significand = d-sign*32768.-exponent*1024.;\n"
"\n"
"    // Return full float16\n"
"    if(exponent == 0.)\n"
"        val = mix(1., -1., sign) * 5.960464477539063e-08 * significand;\n"
"    else val = mix(1., -1., sign) * (1. + significand * 9.765625e-4) * pow(2.,exponent-15.);\n"
"}\n"
"\n"
"\n"
"// 2D box\n"
"void box(in vec2 x, in vec2 b, out float dst)\n"
"{\n"
"    vec2 d = abs(x) - b;\n"
"    dst = length(max(d,c.yy)) + min(max(d.x,d.y),0.);\n"
"}\n"
"\n"
"// Distance to circle\n"
"void circle(in vec2 x, out float d)\n"
"{\n"
"    d = abs(length(x)-1.0);\n"
"}\n"
"\n"
"// Distance to line segment\n"
"void lineseg(in vec2 x, in vec2 p1, in vec2 p2, out float d)\n"
"{\n"
"    vec2 da = p2-p1;\n"
"    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));\n"
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
"// Get glyph data from texture\n"
"void dglyph(in vec2 x, in float ordinal, in float size, out float dst)\n"
"{\n"
"    float dis;\n"
"    box(x, 2.*size*c.xx, dis);\n"
"    if(dis > 0.)\n"
"    {\n"
"        dst = dis+.5*size;\n"
"        return;\n"
"    }\n"
"\n"
"    // Find glyph offset in glyph index\n"
"    float nglyphs, offset = 0;\n"
"    rfloat(1., nglyphs);\n"
"        \n"
"    for(float i=0.; i<nglyphs; i+=1.)\n"
"    {\n"
"        float ord;\n"
"        rfloat(2.+2.*i, ord);\n"
"        ord = floor(ord);\n"
"        \n"
"        if(ord == ordinal)\n"
"        {\n"
"            rfloat(2.+2.*i+1., offset);\n"
"            offset = floor(offset);\n"
"            break;\n"
"        }\n"
"    }\n"
"    \n"
"    if(offset == 0.) \n"
"    {\n"
"        dst = 1.;\n"
"        return;\n"
"    }\n"
"    \n"
"    // Get distance from glyph data\n"
"    float d = 1.;\n"
"    \n"
"    // Lines\n"
"    float nlines;\n"
"    rfloat(offset, nlines);\n"
"    nlines = floor(nlines);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<nlines; i+=1.)\n"
"    {\n"
"        float x1;\n"
"        rfloat(offset, x1);\n"
"        offset += 1.;\n"
"        float y1;\n"
"        rfloat(offset, y1);\n"
"        offset += 1.;\n"
"        float x2;\n"
"        rfloat(offset, x2);\n"
"        offset += 1.;\n"
"        float y2;\n"
"        rfloat(offset, y2);\n"
"        offset += 1.;\n"
"        float da;\n"
"        lineseg(x, size*vec2(x1,y1), size*vec2(x2, y2), da);\n"
"        d = min(d,da);\n"
"    }\n"
"    \n"
"    // Circles\n"
"    float ncircles;\n"
"    rfloat(offset, ncircles);\n"
"    ncircles = floor(ncircles);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<max(ncircles,0); i+=1.)\n"
"    {\n"
"        float xc;\n"
"        rfloat(offset, xc);\n"
"        offset += 1.;\n"
"        float yc;\n"
"        rfloat(offset, yc);\n"
"        offset += 1.;\n"
"        float r;\n"
"        rfloat(offset, r);\n"
"        offset += 1.;\n"
"        float da;\n"
"        circle( size*r*(x-size*vec2(xc, yc)),da);\n"
"        d = min(d, da/size/r);\n"
"    }\n"
"    \n"
"    // Circle segments\n"
"    float nsegments;\n"
"    rfloat(offset, nsegments);\n"
"    nsegments = floor(nsegments);\n"
"    offset += 1.;\n"
"    for(float i=0.; i<max(nsegments,0); i+=1.)\n"
"    {\n"
"        float xc;\n"
"        rfloat(offset, xc);\n"
"        offset += 1.;\n"
"        float yc;\n"
"        rfloat(offset, yc);\n"
"        offset += 1.;\n"
"        float r;\n"
"        rfloat(offset, r);\n"
"        offset += 1.;\n"
"        float phi0;\n"
"        rfloat(offset, phi0);\n"
"        offset += 1.;\n"
"        float phi1;\n"
"        rfloat(offset, phi1);\n"
"        offset += 1.;\n"
"        float da;\n"
"        circlesegment(x-size*vec2(xc,yc), size*r, phi0, phi1, da);\n"
"        d = min(d, da);\n"
"    }\n"
"    \n"
"    if(nlines+ncircles+nsegments == 0.)\n"
"        dst = dis;\n"
"    else dst = d;\n"
"}\n"
"\n"
"// Get distance to string from database\n"
"void dstring(in vec2 x, in float ordinal, in float size, out float dst)\n"
"{\n"
"    // Get string database offset\n"
"    float stroff0;\n"
"    rfloat(0., stroff0);\n"
"    stroff0 = floor(stroff0);\n"
"    \n"
"    // Return 1 if wrong ordinal is supplied\n"
"    float nstrings;\n"
"    rfloat(stroff0, nstrings);\n"
"    nstrings = floor(nstrings);\n"
"    if(ordinal >= nstrings)\n"
"    {\n"
"        dst = 1.;\n"
"        return;\n"
"    }\n"
"    \n"
"    // Get offset and length of string from string database index\n"
"    float stroff;\n"
"    rfloat(stroff0+1.+2.*ordinal, stroff);\n"
"    stroff = floor(stroff);\n"
"    float len;\n"
"    rfloat(stroff0+2.+2.*ordinal, len);\n"
"    len = floor(len);\n"
"    \n"
"    // Draw glyphs\n"
"    vec2 dx = mod(x-size, 2.*size)-size, \n"
"        ind = ceil((x-dx+size)/2./size);\n"
"    \n"
"    // Bounding box\n"
"    float bound;\n"
"    box(x-size*(len-3.)*c.xy, vec2(size*len, 1.*size), bound);\n"
"    if(bound > 0.)\n"
"    {\n"
"        dst = bound+.5*size;\n"
"        return;\n"
"    }\n"
"    \n"
"    float da;\n"
"    rfloat(stroff+ind.x, da);\n"
"    da = floor(da);\n"
"    dglyph(dx, da, .7*size, dst);\n"
"}\n"
"\n"
"// distance to a floating point number string\n"
"// for debugging stuff while shader is loaded\n"
"void dfloat(in vec2 x, in float num, in float size, out float dst)\n"
"{\n"
"    float d = 1., index = 0.;\n"
"    \n"
"    // Determine sign and output it if present\n"
"    float sign = sign(num), exp = 0.;\n"
"    if(sign<0.)\n"
"    {\n"
"        float da;\n"
"        dglyph(x, 45., .7*size, da);\n"
"        d = min(d, da);\n"
"        index += 1.;\n"
"        num *= -1.;\n"
"    }\n"
"    \n"
"    // The first power of ten that floors num to anything not zero is the exponent\n"
"    for(exp = -15.; exp < max(15., -32.+sign); exp += 1.)\n"
"        if(floor(num*pow(10.,exp)) != 0.)\n"
"            break;\n"
"    exp *= -1.;\n"
"    // Determine the significand and output it\n"
"    for(float i = exp; i >= max(exp-5.,-33); i -= 1.)\n"
"    {\n"
"        float po = pow(10.,i);\n"
"        float ca = floor(num/po);\n"
"        num -= ca*po;\n"
"        \n"
"        float da;\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, da);\n"
"        d = min(d, da);\n"
"        index += 1.;\n"
"        if(i == exp) // decimal point\n"
"        {\n"
"            dglyph(x-2.*index*size*c.xy, 46., .7*size, da);\n"
"            d = min(d, da);\n"
"            index += 1.;\n"
"        }\n"
"    }\n"
"    \n"
"    // Output the exponent\n"
"    float db;\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 101., .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    if(exp < 0.) // Sign\n"
"    {\n"
"        dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 45., .7*size,db);\n"
"        d = min(d, db);\n"
"        index += 1.;\n"
"        exp *= -1.;\n"
"    }\n"
"    float ca = floor(exp/10.);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    ca = floor(exp-10.*ca);\n"
"    dglyph(x+.7*size*c.xy-2.*index*size*c.xy, 48.+ca, .7*size, db);\n"
"    d = min(d, db);\n"
"    index += 1.;\n"
"    \n"
"    dst = d;\n"
"}\n"
"\n"
"// Add scene contents\n"
"void add(in vec4 src1, in vec4 src2, out vec4 dst)\n"
"{\n"
"    dst = mix(src1, src2, step(src2.x, src1.x));\n"
"}\n"
"\n"
"void blend(in vec4 src, in float tlo, in float thi, out vec4 dst)\n"
"{\n"
"    vec4 added;\n"
"    add(dst, src, added);\n"
"    dst = mix(src, added, smoothstep(tlo-.5,tlo+.5,iTime)*(1.-smoothstep(thi-.5,thi+.5,iTime)));\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    a = iResolution.x/iResolution.y;\n"
"    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);\n"
"\n"
"    float d;\n"
"    dstring(uv, 3., .1, d); // Team210 present\n"
"    vec4 sdf;\n"
"    sdf.gba = mix(texture(iChannel0, uv).xyz, c.xyy, step(0.,d));\n"
"    sdf.x = d;\n"
"    blend(sdf, 5., 20., sdf);\n"
"\n"
"    fragColor = vec4(sdf.gba, 1.);\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
