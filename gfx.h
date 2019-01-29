/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef GFX_H_
# define GFX_H_

const char *gfx_frag =
 "#version 130\n"
 "float iScale,iNBeats;"
 "uniform float iTime;"
 "uniform vec2 iResolution;"
 "uniform sampler2D iFont;"
 "uniform float iFontWidth;"
 "uniform sampler2D iSequence;"
 "uniform float iSequenceWidth,iExecutableSize;"
 "const vec3 c=vec3(1.,0.,-1.);"
 "const float pi=acos(-1.);"
 "float a;\n"
 "#define FSAA 2\n"
 "float rshort(float off)"
 "{"
   "float hilo=mod(off,2.);"
   "off*=.5;"
   "vec2 ind=(vec2(mod(off,iFontWidth),floor(off/iFontWidth))+.05)/iFontWidth;"
   "vec4 block=texture(iFont,ind);"
   "vec2 data=mix(block.xy,block.zw,hilo);"
   "return round(dot(vec2(255.,65280.),data));"
 "}"
 "float rfloat(float off)"
 "{"
   "float d=rshort(off),sign=floor(d/32768.),exponent=floor(d/1024.-sign*32.),significand=d-sign*32768.-exponent*1024.;"
   "if(exponent==0.)"
     "return mix(1.,-1.,sign)*5.96046e-08*significand;"
   "return mix(1.,-1.,sign)*(1.+significand*.000976563)*pow(2.,exponent-15.);"
 "}"
 "float rshorts(float off)"
 "{"
   "float hilo=mod(off,2.);"
   "off*=.5;"
   "vec2 ind=(vec2(mod(off,iSequenceWidth),floor(off/iSequenceWidth))+.05)/iSequenceWidth;"
   "vec4 block=texture(iSequence,ind);"
   "vec2 data=mix(block.xy,block.zw,hilo);"
   "return round(dot(vec2(255.,65280.),data));"
 "}"
 "float rfloats(float off)"
 "{"
   "float d=rshorts(off),sign=floor(d/32768.),exponent=floor(d/1024.-sign*32.),significand=d-sign*32768.-exponent*1024.;"
   "if(exponent==0.)"
     "return mix(1.,-1.,sign)*5.96046e-08*significand;"
   "return mix(1.,-1.,sign)*(1.+significand*.000976563)*pow(2.,exponent-15.);"
 "}\n"
 "#define NTRK 2\n"
 "#define NMOD 5\n"
 "#define NPTN 2\n"
 "#define NNOT 32\n"
 "int trk_sep(int index)"
 "{"
   "return int(rfloats(index));"
 "}"
 "int trk_syn(int index)"
 "{"
   "return int(rfloats(index+1+1*NTRK));"
 "}"
 "float trk_norm(int index)"
 "{"
   "return rfloats(index+1+2*NTRK);"
 "}"
 "float trk_rel(int index)"
 "{"
   "return rfloats(index+1+3*NTRK);"
 "}"
 "float mod_on(int index)"
 "{"
   "return rfloats(index+1+4*NTRK);"
 "}"
 "float mod_off(int index)"
 "{"
   "return rfloats(index+1+4*NTRK+1*NMOD);"
 "}"
 "int mod_ptn(int index)"
 "{"
   "return int(rfloats(index+1+4*NTRK+2*NMOD));"
 "}"
 "float mod_transp(int index)"
 "{"
   "return rfloats(index+1+4*NTRK+3*NMOD);"
 "}"
 "int ptn_sep(int index)"
 "{"
   "return int(rfloats(index+1+4*NTRK+4*NMOD));"
 "}"
 "float note_on(int index)"
 "{"
   "return rfloats(index+2+4*NTRK+4*NMOD+NPTN);"
 "}"
 "float note_off(int index)"
 "{"
   "return rfloats(index+2+4*NTRK+4*NMOD+NPTN+1*NNOT);"
 "}"
 "float note_pitch(int index)"
 "{"
   "return rfloats(index+2+4*NTRK+4*NMOD+NPTN+2*NNOT);"
 "}"
 "float note_vel(int index)"
 "{"
   "return rfloats(index+2+4*NTRK+4*NMOD+NPTN+3*NNOT);"
 "}"
 "const float BPM=35.,BPS=BPM/60.,SPB=60./BPM;"
 "float scale()"
 "{"
   "float max_mod_off=4.;"
   "int drum_index=24;"
   "float drum_synths=10.,r=0.,d=0.,BT=mod(BPS*iTime*44100.,max_mod_off);"
   "if(BT>max_mod_off)"
     "return r;"
   "float time=SPB*BT,r_sidechain=1.,Bon=0.,Boff=0.;"
   "for(int trk=0;trk<NTRK;trk++)"
     "{"
       "if(trk==drum_index)"
         "continue;"
       "int tsep=trk_sep(trk),tlen=trk_sep(trk+1)-tsep,_modU=tlen-1;"
       "for(int i=0;i<tlen-1;i++)"
         "if(BT<mod_on(tsep+i))"
           "{"
             "_modU=i;"
             "break;"
           "}"
       "int _modL=tlen-1;"
       "for(int i=0;i<tlen-1;i++)"
         "if(BT<mod_off(tsep+i)+trk_rel(trk))"
           "{"
             "_modL=i;"
             "break;"
           "}"
       "for(int _mod=_modL;_mod<=_modU;_mod++)"
         "{"
           "float B=BT-mod_on(tsep+_mod);"
           "int ptn=mod_ptn(tsep+_mod),psep=ptn_sep(ptn),plen=ptn_sep(ptn+1)-psep,_noteU=plen-1;"
           "for(int i=0;i<plen-1;i++)"
             "if(B<note_on(psep+i+1))"
               "{"
                 "_noteU=i;"
                 "break;"
               "}"
           "int _noteL=plen-1;"
           "for(int i=0;i<plen-1;i++)"
             "if(B<=note_off(psep+i)+trk_rel(trk))"
               "{"
                 "_noteL=i;"
                 "break;"
               "}"
           "for(int _note=_noteL;_note<=_noteU;_note++)"
             "{"
               "Bon=note_on(psep+_note);"
               "Boff=note_off(psep+_note);"
               "int Bdrum=int(note_pitch(psep+_note));"
               "d=max(d,clamp(smoothstep(Bon,mix(Bon,3.*Boff-2.*Bon,.5),time)*(1.-smoothstep(mix(Bon,3.*Boff-2.*Bon,.5),3.*Boff-2.*Bon,time)),0.,1.));"
             "}"
           "return d;"
         "}"
     "}"
   "return 0.;"
 "}"
 "float rand(vec2 x)"
 "{"
   "return fract(sin(dot(x-1.,vec2(12.9898,78.233)))*43758.5);"
 "}"
 "vec3 taylorInvSqrt(vec3 r)"
 "{"
   "return 1.79284-.853735*r;"
 "}"
 "vec3 permute(vec3 x)"
 "{"
   "return mod((x*34.+1.)*x,289.);"
 "}"
 "float snoise(vec2 P)"
 "{"
   "const vec2 C=vec2(.211325,.366025);"
   "vec2 i=floor(P+dot(P,C.yy)),x0=P-i+dot(i,C.xx),i1;"
   "i1.x=step(x0.y,x0.x);"
   "i1.y=1.-i1.x;"
   "vec4 x12=x0.xyxy+vec4(C.xx,C.xx*2.-1.);"
   "x12.xy-=i1;"
   "i=mod(i,289.);"
   "vec3 p=permute(permute(i.y+vec3(0.,i1.y,1.))+i.x+vec3(0.,i1.x,1.)),m=max(.5-vec3(dot(x0,x0),dot(x12.xy,x12.xy),dot(x12.zw,x12.zw)),0.);"
   "m=m*m;"
   "m=m*m;"
   "vec3 x=fract(p*(1./41.))*2.-1.,gy=abs(x)-.5,ox=floor(x+.5),gx=x-ox;"
   "m*=taylorInvSqrt(gx*gx+gy*gy);"
   "vec3 g;"
   "g.x=gx.x*x0.x+gy.x*x0.y;"
   "g.yz=gx.yz*x12.xz+gy.yz*x12.yw;"
   "return-1.+2.*(130.*dot(m,g));"
 "}"
 "float mfsnoise(vec2 x,float f0,float f1,float phi)"
 "{"
   "float sum=0.,a=1.2;"
   "for(float f=f0;f<f1;f=f*2.)"
     "sum=a*snoise(f*x)+sum,a=a*phi;"
   "return sum;"
 "}"
 "mat3 rot(vec3 p)"
 "{"
   "return mat3(c.xyyy,cos(p.x),sin(p.x),0.,-sin(p.x),cos(p.x))*mat3(cos(p.y),0.,-sin(p.y),c.yxy,sin(p.y),0.,cos(p.y))*mat3(cos(p.z),-sin(p.z),0.,sin(p.z),cos(p.z),c.yyyx);"
 "}"
 "vec2 add(vec2 sda,vec2 sdb)"
 "{"
   "return mix(sda,sdb,step(sdb.x,sda.x));"
 "}"
 "vec2 sub(vec2 sda,vec2 sdb)"
 "{"
   "return mix(-sda,sdb,step(sda.x,sdb.x));"
 "}"
 "float lineseg(vec2 x,vec2 p1,vec2 p2)"
 "{"
   "vec2 d=p2-p1;"
   "return length(x-mix(p1,p2,clamp(dot(x-p1,d)/dot(d,d),0.,1.)));"
 "}"
 "float lineseg(vec3 x,vec3 p1,vec3 p2)"
 "{"
   "vec3 d=p2-p1;"
   "return length(x-mix(p1,p2,clamp(dot(x-p1,d)/dot(d,d),0.,1.)));"
 "}"
 "float circle(vec2 x,float r)"
 "{"
   "return abs(length(x)-r);"
 "}"
 "float circlesegment(vec2 x,float r,float p0,float p1)"
 "{"
   "float p=atan(x.y,x.x);"
   "vec2 philo=vec2(max(p0,p1),min(p0,p1));"
   "if(p<philo.x&&p>philo.y||p+2.*pi<philo.x&&p+2.*pi>philo.y||p-2.*pi<philo.x&&p-2.*pi>philo.y)"
     "return abs(length(x)-r);"
   "return min(length(x-r*vec2(cos(p0),sin(p0))),length(x-r*vec2(cos(p1),sin(p1))));"
 "}"
 "float dpoly_min(vec2 x,float N,float R)"
 "{"
   "float d=2.*pi/N,t=mod(acos(x.x/length(x)),d)-.5*d;"
   "return R-length(x)*cos(t)/cos(.5*d);"
 "}"
 "float box(vec2 x,vec2 b)"
 "{"
   "vec2 d=abs(x)-b;"
   "return length(max(d,c.yy))+min(max(d.x,d.y),0.);"
 "}"
 "float dglyph(vec2 x,float ordinal,float size)"
 "{"
   "float dis=box(x,2.*size*c.xx);"
   "if(dis>0.)"
     "return dis+.5*size;"
   "float nglyphs=rfloat(1.),offset=0;"
   "for(float i=0.;i<nglyphs;i+=1.)"
     "{"
       "float ord=floor(rfloat(2.+2.*i));"
       "if(ord==ordinal)"
         "{"
           "offset=floor(rfloat(2.+2.*i+1.));"
           "break;"
         "}"
     "}"
   "if(offset==0.)"
     "return 1.;"
   "float d=1.,nlines=floor(rfloat(offset));"
   "offset+=1.;"
   "for(float i=0.;i<nlines;i+=1.)"
     "{"
       "float x1=rfloat(offset);"
       "offset+=1.;"
       "float y1=rfloat(offset);"
       "offset+=1.;"
       "float x2=rfloat(offset);"
       "offset+=1.;"
       "float y2=rfloat(offset);"
       "offset+=1.;"
       "d=min(d,lineseg(x,size*vec2(x1,y1),size*vec2(x2,y2)));"
     "}"
   "float ncircles=floor(rfloat(offset));"
   "offset+=1.;"
   "for(float i=0.;i<ncircles;i+=1.)"
     "{"
       "float xc=rfloat(offset);"
       "offset+=1.;"
       "float yc=rfloat(offset);"
       "offset+=1.;"
       "float r=rfloat(offset);"
       "offset+=1.;"
       "d=min(d,circle(x-size*vec2(xc,yc),size*r));"
     "}"
   "float nsegments=floor(rfloat(offset));"
   "offset+=1.;"
   "for(float i=0.;i<nsegments;i+=1.)"
     "{"
       "float xc=rfloat(offset);"
       "offset+=1.;"
       "float yc=rfloat(offset);"
       "offset+=1.;"
       "float r=rfloat(offset);"
       "offset+=1.;"
       "float phi0=rfloat(offset);"
       "offset+=1.;"
       "float phi1=rfloat(offset);"
       "offset+=1.;"
       "d=min(d,circlesegment(x-size*vec2(xc,yc),size*r,phi0,phi1));"
     "}"
   "if(nlines+ncircles+nsegments==0.)"
     "return dis;"
   "return d;"
 "}"
 "float dstring(vec2 x,float ordinal,float size)"
 "{"
   "float stroff0=floor(rfloat(0.)),nstrings=floor(rfloat(stroff0));"
   "if(ordinal>=nstrings)"
     "return 1.;"
   "float stroff=floor(rfloat(stroff0+1.+2.*ordinal)),len=floor(rfloat(stroff0+2.+2.*ordinal));"
   "vec2 dx=mod(x-size,2.*size)-size,ind=ceil((x-dx+size)/2./size);"
   "float bound=box(x-size*(len-3.)*c.xy,vec2(size*len,size));"
   "if(bound>0.)"
     "return bound+.5*size;"
   "return dglyph(dx,floor(rfloat(stroff+ind.x)),.7*size);"
 "}"
 "float dfloat(vec2 x,float num,float size)"
 "{"
   "float d=1.,index=0.,sign=sign(num),exp=0.;"
   "if(sign<0.)"
     "d=min(d,dglyph(x,45.,.7*size)),index+=1.,num*=-1.;"
   "for(exp=-15.;exp<15.;exp+=1.)"
     "if(floor(num*pow(10.,exp))!=0.)"
       "break;"
   "exp*=-1.;"
   "for(float i=exp;i>=exp-5.;i-=1.)"
     "{"
       "float po=pow(10.,i),ca=floor(num/po);"
       "num-=ca*po;"
       "d=min(d,dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size));"
       "index+=1.;"
       "if(i==exp)"
         "d=min(d,dglyph(x-2.*index*size*c.xy,46.,.7*size)),index+=1.;"
     "}"
   "d=min(d,dglyph(x+.7*size*c.xy-2.*index*size*c.xy,101.,.7*size));"
   "index+=1.;"
   "if(exp<0.)"
     "d=min(d,dglyph(x+.7*size*c.xy-2.*index*size*c.xy,45.,.7*size)),index+=1.,exp*=-1.;"
   "float ca=floor(exp/10.);"
   "d=min(d,dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size));"
   "index+=1.;"
   "ca=floor(exp-10.*ca);"
   "d=min(d,dglyph(x+.7*size*c.xy-2.*index*size*c.xy,48.+ca,.7*size));"
   "index+=1.;"
   "return d;"
 "}"
 "float logo(vec2 x,float r)"
 "{"
   "return min(min(circle(x+r*c.zy,r),lineseg(x,r*c.yz,r*c.yx)),circlesegment(x+r*c.xy,r,-.5*pi,.5*pi));"
 "}"
 "float stroke(float d,float w)"
 "{"
   "return abs(d)-w;"
 "}"
 "vec2 ind;"
 "float hexagon(vec2 p)"
 "{"
   "vec2 q=vec2(p.x*1.2,p.y+p.x*.6),pi=floor(q),pf=fract(q);"
   "float v=mod(pi.x+pi.y,3.),ca=step(1.,v),cb=step(2.,v);"
   "vec2 ma=step(pf.xy,pf.yx);"
   "ind=pi+ca-cb*ma;"
   "return dot(ma,1.-pf.yx+ca*(pf.x+pf.y-1.)+cb*(pf.yx-2.*pf.xy));"
 "}"
 "float zextrude(float z,float d2d,float h)"
 "{"
   "vec2 w=vec2(-d2d,abs(z)-.5*h);"
   "return length(max(w,0.));"
 "}"
 "float box(vec3 x,vec3 b)"
 "{"
   "return length(max(abs(x)-b,0.));"
 "}"
 "vec2 inset(vec3 x)"
 "{"
   "float rs=1.9;"
   "return vec2(-.11+min(x.y+.4,abs(length(x)-rs)),2.);"
 "}"
 "vec2 scene(vec3 x)"
 "{"
   "vec2 sdf=vec2(x.y+.4,1.);"
   "float rs=1.9;"
   "sdf=add(sdf,vec2(abs(length(x)-2.*rs),0.));"
   "vec2 pt=vec2(atan(x.x,x.y+.4),-acos(x.z/length(x+.4*c.yxy)));"
   "float d=stroke(zextrude(length(x)-rs,-stroke(hexagon(vec2(5.,10.)*pt),.1),.1),.05);"
   "sdf=add(sdf,vec2(d,3.));"
   "if(rand(ind)<.5)"
     "{"
       "float d=stroke(zextrude(length(x)-rs,stroke(hexagon(vec2(5.,10.)*pt),.1),.1),.01);"
       "sdf=add(sdf,vec2(d,2.));"
     "}"
   "d=stroke(zextrude(x.y+.4,-stroke(length(x.xz)-rs,.2),.1),.05);"
   "sdf=add(sdf,vec2(d,4.));"
   "float dr=.1;"
   "vec3 y=mod(x,dr)-.5*dr;"
   "float guard=-length(max(abs(y)-vec3(.5*dr*c.xx,.6),0.));"
   "guard=abs(guard)+dr*.1;"
   "sdf.x=min(sdf.x,guard);"
   "return sdf;"
 "}"
 "vec2 greetings(vec3 x)"
 "{"
   "vec2 sdf=c.xy;"
   "return sdf;"
 "}"
 "vec2 textpre(vec3 x)"
 "{"
   "vec2 sdf=vec2(x.z,7.);"
   "float structure=stroke(logo(x.xy+.3*c.xy,.6),.25),blend=smoothstep(2.,6.,iTime)*(1.-smoothstep(6.,12.,iTime));"
   "if(structure<0.&&blend>=.001)"
     "{"
       "float blend=smoothstep(2.,6.,iTime)*(1.-smoothstep(6.,12.,iTime));"
       "sdf=vec2(stroke(zextrude(x.z,1.5*x.z-stroke(logo(x.xy+.3*c.xy,.6),.25),blend*clamp(1.-exp(-(x.x-34.)-8.*iTime),0.,.5)),.05*blend),7.);"
     "}"
   "sdf.x=abs(sdf.x)-.3;"
   "return sdf;"
 "}"
 "vec2 texteffect(vec3 x)"
 "{"
   "vec2 sdf=vec2(x.z,7.);"
   "float hex=hexagon(18.*x.xy);"
   "vec2 cind=ind/18.;"
   "cind=vec2(cind.x/1.2,cind.y);"
   "cind=vec2(cind.x,cind.y-cind.x*.6);"
   "float structure=stroke(logo(cind+.3*c.xy,.6),.25),blend=smoothstep(2.,6.,iTime)*(1.-smoothstep(6.,12.,iTime));"
   "if(structure<0.&&blend>=.001)"
     "{"
       "float blend=smoothstep(2.,6.,iTime)*(1.-smoothstep(6.,12.,iTime));"
       "sdf=vec2(stroke(zextrude(x.z,2.*x.z-stroke(logo(cind.xy+.3*c.xy,.6),.25),(.6+.15*snoise(4.*cind.xy-iTime))*blend*clamp(1.-exp(-(ind.x-34.)-8.*iTime),0.,1.)),.05*blend),7.);"
     "}"
   "float dr=.03;"
   "vec3 y=mod(x,dr)-.5*dr;"
   "float guard=-length(max(abs(y)-vec3(.5*dr*c.xx,.6),0.));"
   "guard=abs(guard)+dr*.1;"
   "sdf.x=min(sdf.x,guard);"
   "return sdf;"
 "}"
 "vec2 texteffect2(vec3 x)"
 "{"
   "vec2 sdf=vec2(x.z,7.);"
   "return sdf;"
 "}"
 "vec3 post1(vec2 uv,vec3 col)"
 "{"
   "if(uv.y<.8)"
     "return col+=vec3(0.,.05,.1)*sin(uv.y*1050.+5.*iTime),col=clamp(col,c.yyy,c.xxx),col;"
   "vec3 blu=vec3(.2,.68,1.);"
   "float px=1.5/iResolution.y,dt0=logo(uv-2.*vec2(-.45*a,.45),.04);"
   "dt0=stroke(dt0,.01);"
   "col=mix(col,mix(col,blu,.5),smoothstep(px,-px,dt0));"
   "dt0=stroke(dt0,.0025);"
   "col=mix(col,blu,smoothstep(px,-px,dt0));"
   "dt0=stroke(lineseg(uv-2.*vec2(-.45*a,.45)-.2*c.xy,c.yy,c.xy),.05);"
   "col=mix(col,mix(col,blu,.5),smoothstep(px,-px,dt0));"
   "float dt1=stroke(dt0,.0025);"
   "col=mix(col,blu,smoothstep(px,-px,dt1));"
   "float dta=dstring(uv-2.*vec2(-.45*a,.45)-.3*c.xy,1.,.025);"
   "dta=min(dta,dfloat(uv-2.*vec2(-.45*a,.45)-.7*c.xy,iTime,.025));"
   "dta=stroke(dta,.0025);"
   "col=mix(col,clamp(blu,0.,1.),smoothstep(px,-px,dta));"
   "dt0=stroke(lineseg(uv-2.*vec2(-.45*a,.45)-1.4*c.xy,c.yy,c.xy),.05);"
   "col=mix(col,mix(col,blu,.5),smoothstep(px,-px,dt0));"
   "dt1=stroke(dt0,.0025);"
   "col=mix(col,blu,smoothstep(px,-px,dt1));"
   "dta=dstring(uv-2.*vec2(-.45*a,.45)-1.5*c.xy,2.,.025);"
   "dta=min(dta,dfloat(uv-2.*vec2(-.45*a,.45)-1.7*c.xy,iExecutableSize,.025));"
   "dta=stroke(dta,.0025);"
   "col=mix(col,clamp(blu,0.,1.),smoothstep(px,-px,dta));"
   "col+=vec3(0.,.05,.1)*sin(uv.y*1050.+5.*iTime);"
   "col=clamp(col,c.yyy,c.xxx);"
   "return col;"
 "}\n"
 "#define raymarch(scene,xc,ro,d,dir,s,N,eps,flag)flag=false;for(int ia=0;ia<N;++ia){xc=ro+d*dir;s=scene(xc);if(s.x<eps){flag=true;break;}d+=s.x;}\n"
 "#define calcnormal(scene,n,eps,xc){float ss=scene(xc).x;n=normalize(vec3(scene(xc+eps*c.xyy).x-ss,scene(xc+eps*c.yxy).x-ss,scene(xc+eps*c.yyx).x-ss));}\n"
 "#define camerasetup(camera,ro,r,u,t,uv,dir){camera(ro,r,u,t);t+=uv.x*r+uv.y*u;dir=normalize(t-ro);}\n"
 "#define post(color,uv){color=post1(uv,color);}\n"
 "void camera1(out vec3 ro,out vec3 r,out vec3 u,out vec3 t)"
 "{"
   "ro=c.yyx-c.yyx+.1*(iTime-28.)*c.yyx,r=c.xyy,u=c.yxy,t=c.yyy-c.yyx+.1*(iTime-28.)*c.yyx;"
 "}"
 "void camera0(out vec3 ro,out vec3 r,out vec3 u,out vec3 t)"
 "{"
   "float blend=0.;"
   "ro=c.yyx-.5*c.yxy*blend;"
   "r=c.xyy;"
   "u=c.yxy+.5*c.yyx*blend;"
   "t=.1*c.yyx;"
 "}"
 "vec3 stdcolor(vec2 x)"
 "{"
   "return.5+.5*cos(iTime+x.xyx+vec3(0,2,4));"
 "}"
 "float star(vec2 x,float r0)"
 "{"
   "return 1.-smoothstep(.5*r0,r0,length(x));"
 "}"
 "vec3 background(vec2 x)"
 "{"
   "vec3 col=mix(vec3(1.,.56,.44),vec3(1.,1.,.87),abs(x.y));"
   "float d=length(vec2(x.x,abs(x.y+.15))-.3*c.yx)-.15;"
   "col=mix(col,c.xxx,smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d));"
   "float da=.5*snoise(2.*vec2(x.x,abs(x.y+.15))-.4*iTime),dx=da+mfsnoise(vec2(x.x-.2*iTime,abs(x.y+.15)),1.,1000.,.55);"
   "col=mix(col,.8*vec3(1.,.7,.57),clamp(2.5+dx,0.,1.));"
   "col=mix(col,.9*vec3(1.,.7,.57),clamp(2.3+dx,0.,1.));"
   "col=mix(col,vec3(1.,.7,.57),clamp(2.1+dx,0.,1.));"
   "da=.5*snoise(2.*vec2(x.x,abs(x.y+.15))-.4*iTime-15.);"
   "dx=da+mfsnoise(vec2(x.x-.1*iTime-15.,abs(x.y+.15)),1.,1000.,.55);"
   "col=mix(col,.8*vec3(1.,1.,.87),clamp(2.5+dx,0.,1.));"
   "col=mix(col,.9*vec3(1.,1.,.87),clamp(2.3+dx,0.,1.));"
   "col=mix(col,vec3(1.,1.,.87),clamp(1.6+dx,0.,1.));"
   "col=mix(col,c.xxx,.9*smoothstep(1.5/iResolution.y,-1.5/iResolution.y,d));"
   "return col;"
 "}"
 "vec3 background2(vec2 uv)"
 "{"
   "float d=stroke(-hexagon(18.*uv),.1);"
   "vec2 cind=ind/18.;"
   "cind=vec2(cind.x/1.2,cind.y);"
   "cind=vec2(cind.x,cind.y-cind.x*.6);"
   "float structure=exp(-(ind.x-34.)-8.*iTime)+stroke(logo(cind+.3*c.xy,.6),.25);"
   "structure=mix(structure,hexagon(18.*uv),clamp(.25*(iTime-12.),0.,1.));"
   "vec2 dind=cind;"
   "float endeavor=dstring(cind+2.*(-6.+1.2*iTime-16.8)*c.xy,0.,.8);"
   "endeavor=stroke(endeavor,.2);"
   "structure=mix(structure,endeavor,clamp(.25*(iTime-14.),0.,1.));"
   "d=mix(d,stroke(-hexagon(8.*uv),.1),clamp(.5*(iTime-24.),0.,1.));"
   "structure=mix(structure,hexagon(8.*uv),clamp(.5*(iTime-24.),0.,1.));"
   "dind=ind/8.;"
   "dind=vec2(dind.x/1.2,dind.y);"
   "dind=vec2(dind.x,dind.y-dind.x*.6);"
   "cind=mix(cind,dind,clamp(.5*(iTime-24.),0.,1.));"
   "vec2 dt=vec2(snoise(cind+2.),snoise(cind+3.));"
   "float m=1.5+.5*snoise(10.*cind)+mix(-2.,clamp(.5+.5*snoise(.05*cind-dt-iTime*c.xx),0.,1.),clamp(.125*(iTime-1.),0.,1.));"
   "vec3 c1=mix(c.yyy,c.yyy,m)*smoothstep(-1.5/iResolution.y,1.5/iResolution.y,d);"
   "c1=mix(c1,mix(c.yyy,vec3(1.,.27,0.),m),smoothstep(-1.5/iResolution.y,1.5/iResolution.y,stroke(structure,.05)))*smoothstep(-1.5/iResolution.y,1.5/iResolution.y,d);"
   "c1=clamp(c1,0.,1.);"
   "if(structure>0.)"
     "c1=mix(.3*length(c1)*c.xxx/sqrt(3.),c1,clamp(.5*(iTime-24.),0.,1.));"
   "c1=mix(c1,c.yyy,clamp(iTime-27.,0.,1.));"
   "return clamp(c1,0.,1.);"
 "}"
 "vec3 color(float rev,float ln,float mat,vec2 uv,vec3 x)"
 "{"
   "if(mat==7.)"
     "return background2(x.xy);"
   "if(mat==6.)"
     "return clamp(.7*c.xxx+.7*c.xxy*ln+c.xxx*abs(pow(rev,8.)),0.,1.);"
   "if(mat==2.)"
     "{"
       "vec3 col=.1*c.xyy+.3*c.xyy*abs(ln)+.8*c.xxy*abs(pow(rev,8.));"
       "vec2 pt=vec2(atan(x.x,x.y+.4),-acos(x.z/length(x+.4*c.yxy)));"
       "float d=stroke(-hexagon(6.*vec2(5.,10.)*pt),.1),m=rand(ind);"
       "col=mix(col,mix(col,vec3(1.,.27,0.),m),smoothstep(-1.5/iResolution.y,1.5/iResolution.y,d));"
       "return col;"
     "}"
   "return.1*c.xyy+.3*c.xyy*abs(ln)+.8*c.xxy*abs(pow(rev,8.));"
 "}"
 "void mainImage(out vec4 fragColor,in vec2 fragCoord)"
 "{"
   "iScale=scale();"
   "a=iResolution.x/iResolution.y;"
   "vec3 ro,r,u,t,x,dir;"
   "vec2 s,uv;"
   "float d=0.;"
   "bool hit;"
   "vec3 col=c.yyy;"
   "\n#if FSAA!=1\n"
   "for(int jii=0;jii<FSAA;++jii)"
     "for(int jij=0;jij<FSAA;++jij)"
       "{"
         "vec2 o=vec2(float(jii),float(jij))/float(FSAA)-.5;"
         "uv=(-iResolution.xy+2.*(fragCoord+o))/iResolution.y;"
         "\n#else\n"
         "uv=(-iResolution.xy+2.*fragCoord)/iResolution.y;"
         "\n#endif\n"
         "if(iTime<1000.)"
           "col+=step(uv.y,iScale)*c.xxy;"
         "else"
           " if(iTime<28.)"
             "{"
               "vec3 c1=c.yyy;"
               "camerasetup(camera0,ro,r,u,t,uv,dir);"
               "raymarch(textpre,x,ro,d,dir,s,100,2e-05,hit);"
               "if(hit)"
                 "hit=false;"
               "else"
                 " d=-ro.z/dir.z;"
               "raymarch(texteffect,x,ro,d,dir,s,200,2e-05,hit);"
               "if(hit)"
                 "{"
                   "vec3 n;"
                   "calcnormal(scene,n,.0002,x);"
                   "float rs=1.9;"
                   "vec3 l=x+c.yyx,re=normalize(reflect(-l,n)),v=normalize(x-ro);"
                   "float rev=dot(re,v),ln=dot(l,n);"
                   "c1=color(rev,ln,s.y,uv,x);"
                 "}"
               "else"
                 " c1=background2((ro-ro.z/dir.z*dir).xy);"
               "col+=c1;"
             "}"
           "else"
             " if(iTime<10000.)"
               "{"
                 "vec3 c1=c.yyy;"
                 "camerasetup(camera1,ro,r,u,t,uv,dir);"
                 "d=0.;"
                 "raymarch(inset,x,ro,d,dir,s,40,.0001,hit);"
                 "if(hit)"
                   "hit=false;"
                 "raymarch(scene,x,ro,d,dir,s,300,.0001,hit);"
                 "if(hit)"
                   "{"
                     "vec3 n;"
                     "calcnormal(scene,n,.0002,x);"
                     "float rs=1.9;"
                     "vec3 l=-1.*c.yxy+1.5*c.yyx,re=normalize(reflect(-l,n)),v=normalize(x-ro);"
                     "float rev=dot(re,v),ln=dot(l,n);"
                     "c1=color(rev,ln,s.y,uv,x);"
                     "if(s.y==1.)"
                       "{"
                         "for(float k=.7;k>=.7;k-=.1)"
                           "{"
                             "dir=reflect(dir,n);"
                             "d=.0002;"
                             "ro=x;"
                             "raymarch(inset,x,ro,d,dir,s,20,.0001,hit);"
                             "raymarch(scene,x,ro,d,dir,s,300,.0001,hit);"
                             "if(hit)"
                               "{"
                                 "vec3 n2;"
                                 "calcnormal(scene,n2,.0002,x);"
                                 "re=normalize(reflect(-l,n));"
                                 "v=normalize(x-ro);"
                                 "rev=abs(dot(re,v));"
                                 "ln=abs(dot(l,n));"
                                 "c1=mix(c1,color(rev,ln,s.y,uv,x),k);"
                               "}"
                             "else"
                               " c1=mix(c1,background(uv),k);"
                           "}"
                       "}"
                   "}"
                 "else"
                   " c1=background(uv);"
                 "for(float k=0.;k<8.;k+=1.)"
                   "{"
                     "vec2 dx=.15*vec2(-1.+2.*rand(c.xx+k),-1.+2.*rand(c.xx+k+1.));"
                     "vec3 cx=c.xxx-.2*vec3(rand(c.xx+k+2.),rand(c.xx+k+3.),rand(c.xx+k+4.));"
                     "float sx=.05+.05*rand(c.xx+k+5.),da=dpoly_min(uv-.15*c.yx+dx,6.,sx);"
                     "c1=mix(c1,mix(c1,cx,.5),smoothstep(-1.5/iResolution.y,1.5/iResolution.y,da));"
                   "}"
                 "col+=mix(mix(c.yyy,c1,clamp(iTime-29.,0.,1.)),c.yyy,clamp(iTime-44.,0.,1.));"
               "}"
         "\n#if FSAA!=1\n"
       "}"
   "col/=float(FSAA*FSAA);"
   "\n#else\n"
   "\n#endif\n"
   "post(col,uv);"
   "fragColor=vec4(col,1.);"
 "}"
 "void main()"
 "{"
   "mainImage(gl_FragColor,gl_FragCoord.xy);"
 "}";

#endif // GFX_H_
