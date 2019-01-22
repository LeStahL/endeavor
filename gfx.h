/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef GFX_H_
# define GFX_H_

const char *gfx_frag =
 "const vec3 c=vec3(1.,0.,-1.);"
 "const float pi=acos(-1.);"
 "float a;\n"
 "#define FSAA 2\n"
 "vec3 col=c.yyy;"
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
   "return length(x)-r;"
 "}"
 "float circlesegment(vec2 x,float r,float p0,float p1)"
 "{"
   "float p=atan(x.y,x.x);"
   "p=clamp(p,p0,p1);"
   "return length(x-r*vec2(cos(p),sin(p)));"
 "}"
 "float dglyph(int ascii)"
 "{"
   "return 1.;"
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
 "float dpoly_min(vec2 x,float N,float R)"
 "{"
   "float d=2.*pi/N,t=mod(acos(x.x/length(x)),d)-.5*d;"
   "return R-length(x)*cos(t)/cos(.5*d);"
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
 "}\n"
 "#define raymarch(scene,xc,ro,d,dir,s,N,eps,flag)flag=false;for(int i=0;i<N;++i){xc=ro+d*dir;s=scene(xc);if(s.x<eps){flag=true;break;}d+=s.x;}\n"
 "#define calcnormal(scene,n,eps,xc){float ss=scene(xc).x;n=normalize(vec3(scene(xc+eps*c.xyy).xc-ss,scene(xc+eps*c.yxy).xc-ss,scene(xc+eps*c.yyx).xc-ss));}\n"
 "#define camerasetup(camera,ro,r,u,t,uv,dir){camera(ro,r,u,t);t+=uv.x*r+uv.y*u;dir=normalize(t-ro);}\n"
 "#define post(color,uv){col=mix(clamp(col,c.yyy,c.xxx),c.xxx,smoothstep(1.5/iResolution.y,-1.5/iResolution.y,stroke(logo(uv-2.*vec2(-.45*a,.45),.04),.01)));col+=vec3(0.,0.05,0.1)*sin(uv.y*1050.+5.*iTime);}\n"
 "void camera1(out vec3 ro,out vec3 r,out vec3 u,out vec3 t)"
 "{"
   "ro=c.yyx-c.yyx+.1*(iTime-28.)*c.yyx,r=c.xyy,u=c.yxy,t=c.yyy-c.yyx+.1*(iTime-28.)*c.yyx;"
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
   "d=mix(d,stroke(-hexagon(38.*uv),.1),clamp(.5*(iTime-16.),0.,1.));"
   "structure=mix(structure,hexagon(38.*uv),clamp(.5*(iTime-16.),0.,1.));"
   "vec2 dind=ind/38.;"
   "dind=vec2(dind.x/1.2,dind.y);"
   "dind=vec2(dind.x,dind.y-dind.x*.6);"
   "cind=mix(cind,dind,clamp(.5*(iTime-16.),0.,1.));"
   "float dp=6.*pi/8.,endeavor=circlesegment(cind*c.zx-vec2(1.2,.5),.3,-dp,dp);"
   "endeavor=min(endeavor,lineseg(cind,vec2(-1.1,.5),vec2(-1.3,.5)));"
   "endeavor=min(endeavor,circlesegment(cind.yx-vec2(-.4,.5).yx,.3,-dp,dp));"
   "endeavor=min(endeavor,circlesegment(cind-vec2(.4,.5),.3,-dp,dp));"
   "endeavor=min(endeavor,lineseg(cind,vec2(.3,.45),vec2(.3,.55)));"
   "endeavor=min(endeavor,circlesegment(cind*c.zx-vec2(-1.3,.5),.3,-dp,dp));"
   "endeavor=min(endeavor,lineseg(cind,vec2(1.2,.5),vec2(1.4,.5)));"
   "endeavor=min(endeavor,circlesegment(cind.yx-vec2(-1.2,-.5).yx,.3,-dp,dp));"
   "endeavor=min(endeavor,lineseg(cind,vec2(-1.1,-.5),vec2(-1.3,-.5)));"
   "endeavor=min(endeavor,circlesegment((cind.yx-vec2(-.3,-.5).yx)*c.zx,.3,-dp,dp));"
   "endeavor=min(endeavor,circle(cind-vec2(.6,-.5),.3));"
   "endeavor=min(endeavor,circlesegment(cind*c.zx-vec2(-1.5,-.5),.3,0.,.66*dp));"
   "endeavor=min(endeavor,lineseg(cind,vec2(1.2,-.5),vec2(1.2,-.8)));"
   "endeavor=stroke(endeavor,.13);"
   "structure=mix(structure,endeavor,clamp(.25*(iTime-14.),0.,1.));"
   "d=mix(d,stroke(-hexagon(8.*uv),.1),clamp(.5*(iTime-24.),0.,1.));"
   "structure=mix(structure,hexagon(8.*uv),clamp(.5*(iTime-24.),0.,1.));"
   "dind=ind/8.;"
   "dind=vec2(dind.x/1.2,dind.y);"
   "dind=vec2(dind.x,dind.y-dind.x*.6);"
   "cind=mix(cind,dind,clamp(.5*(iTime-24.),0.,1.));"
   "vec2 dt=vec2(snoise(cind+2.),snoise(cind+3.));"
   "float m=1.5+.5*snoise(10.*cind)+mix(-2.,clamp(.5+.5*snoise(.05*cind-dt-iTime*c.xx),0.,1.),clamp(.125*(iTime-1.),0.,1.));"
   "vec3 c1=mix(c.yyy,.8*c.yyy,m)*smoothstep(-1.5/iResolution.y,1.5/iResolution.y,d);"
   "c1=mix(c1,mix(c.yyy,vec3(1.,.27,0.),m),smoothstep(-1.5/iResolution.y,1.5/iResolution.y,stroke(structure,.05)))*smoothstep(-1.5/iResolution.y,1.5/iResolution.y,d);"
   "c1=clamp(c1,0.,1.);"
   "if(structure>0.)"
     "c1=mix(.7*length(c1)*c.xxx/sqrt(3.),c1,clamp(.5*(iTime-24.),0.,1.));"
   "c1=mix(c1,c.yyy,clamp(iTime-27.,0.,1.));"
   "return clamp(c1,0.,1.);"
 "}"
 "vec3 color(float rev,float ln,float mat,vec2 uv,vec3 x)"
 "{"
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
   "a=iResolution.x/iResolution.y;"
   "vec3 ro,r,u,t,x,dir;"
   "vec2 s,uv;"
   "float d=0.;"
   "bool hit;"
   "\n#if FSAA!=1\n"
   "for(int i=0;i<FSAA;++i)"
     "for(int j=0;j<FSAA;++j)"
       "{"
         "vec2 o=vec2(float(i),float(j))/float(FSAA)-.5;"
         "uv=(-iResolution.xy+2.*(fragCoord+o))/iResolution.y;"
         "\n#else\n"
         "uv=(-iResolution.xy+2.*fragCoord)/iResolution.y;"
         "\n#endif\n"
         "if(iTime<28.)"
           "col+=background2(uv);"
         "else"
           " if(iTime<10000.)"
             "{"
               "vec3 c1=c.yyy;"
               "camerasetup(camera1,ro,r,u,t,uv,dir);"
               "raymarch(inset,x,ro,d,dir,s,40,.0001,hit);"
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
                             "calcnormal(scene,n,.0002,x),l=-1.*c.yxy+1.5*c.yyx,re=normalize(reflect(-l,n)),v=normalize(x-ro),rev=abs(dot(re,v)),ln=abs(dot(l,n)),c1=mix(c1,color(rev,ln,s.y,uv,x),k);"
                           "else"
                             " c1=mix(c1,background(uv),k);"
                         "}"
                     "}"
                 "}"
               "else"
                 " c1=background(uv);"
               "for(float i=0.;i<8.;i+=1.)"
                 "{"
                   "vec2 dx=.15*vec2(-1.+2.*rand(c.xx+i),-1.+2.*rand(c.xx+i+1.));"
                   "vec3 cx=c.xxx-.2*vec3(rand(c.xx+i+2.),rand(c.xx+i+3.),rand(c.xx+i+4.));"
                   "float sx=.05+.05*rand(c.xx+i+5.),da=dpoly_min(uv-.15*c.yx+dx,6.,sx);"
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
 "}";

#endif // GFX_H_
