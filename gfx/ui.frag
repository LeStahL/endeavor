// Global constants
const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);
float a; // Aspect ratio

// 2D box
void box(in vec2 x, in vec2 b, out float dst)
{
    vec2 d = abs(x) - b;
    dst = length(max(d,c.yy)) + min(max(d.x,d.y),0.);
}

// 2D rhomboid
void rhomboid(in vec2 x, in vec2 b, in float tilt, out float dst)
{
    x.x -= tilt/2./b.y*x.y;
    box(x,b,dst);
}


// Distance to circle segment
void circlesegment(in vec2 x, in float r, in float p0, in float p1, out float d)
{
    float p = atan(x.y, x.x);
    vec2 philo = vec2(max(p0, p1), min(p0, p1));
    if((p < philo.x && p > philo.y) || (p+2.*pi < philo.x && p+2.*pi > philo.y) || (p-2.*pi < philo.x && p-2.*pi > philo.y))
    {
        d = abs(length(x)-r);
        return;
    }
    d = min(
        length(x-r*vec2(cos(p0), sin(p0))),
        length(x-r*vec2(cos(p1), sin(p1)))
        );
}

// Distance to line segment
void linesegment(in vec2 x, in vec2 p1, in vec2 p2, out float d)
{
    vec2 da = p2-p1;
    d = length(x-mix(p1, p2, clamp(dot(x-p1, da)/dot(da,da),0.,1.)));
}

// Compute distance to stroke
void stroke(in float d0, in float s, out float d)
{
    d = abs(d0) - s;
}

// Add scene contents
void add(in vec4 src1, in vec4 src2, out vec4 dst)
{
    dst = mix(src1, src2, smoothstep(0., 1.5/iResolution.y, -src2.x));
}

// UI Window Control
void window(in vec2 x, in vec2 size, in vec3 bg, in float title_index, out vec4 col)
{
    size.x *= .5;
    col = vec4(1., bg);
    
    const float cellsize = .015, bordersize = .0025;
    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),
        bordercolor = vec3(1.00,0.71,0.02);
    vec4 c2 = vec4(1., titlecolor);
    
    // Window background
    box(x+.5*size*c.yx,size*vec2(1.,.5),c2.x);
    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/size.y), .5);
    add(col, c2, col);
    
    // Title bar
    c2.gba = titlecolor;
    rhomboid(x+.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    rhomboid(x, vec2(.65*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    rhomboid(x-.8*size.x*c.xy, vec2(.1*size.x,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    
    // Border of title bar
    c2 = vec4(1., bordercolor);
    stroke(col.x,bordersize,c2.x);
    add(col,c2,col);
    
    // Window Border
    linesegment(x, -.9*size.x*c.xy, -size.x*c.xy, c2.x);
    float d;
    linesegment(x, -size.x*c.xy, -size, d);
    c2.x = min(c2.x, d);
    linesegment(x, -size, size*c.xz, d);
    c2.x = min(c2.x, d);
    linesegment(x, size*c.xz, size*c.xy, d);
    c2.x = min(c2.x, d);
    linesegment(x, .9*size.x*c.xy, size.x*c.xy, d);
    c2.x = min(c2.x, d);
    stroke(c2.x,.5*bordersize,c2.x);
    add(col, c2, col);
}

void progressbar(in vec2 x, in float width, in float progress, out vec4 col)
{
    const float cellsize = .015, bordersize = .0025;
    vec3 titlecolor = mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),.5-.5*x.y/cellsize),
        bordercolor = vec3(1.00,0.71,0.02), bg = c.yyy;
    vec4 c2 = vec4(1., titlecolor);
    
    // Window background
    box(x+.5*width*c.yx,width*c.xy,c2.x);
    c2.gba = mix(bg, mix(vec3(0.82,0.00,0.09),vec3(0.45,0.00,0.06),-x.y/cellsize), .5);
    add(col, c2, col);
    
    // Bar background
    c2.gba = titlecolor;
    rhomboid(x, vec2(.5*width,cellsize), cellsize, c2.x);
   	add(col, c2, col);
    
    // Border
    c2.gba = bordercolor;
    stroke(c2.x,.5*bordersize,c2.x);
    add(col, c2, col);
    
    // Progress
    float wc = width/cellsize;
    x.x -= .5*x.y;
    vec2 y = vec2(mod(x.x, 1.2*cellsize)-.6*cellsize, x.y),
        index = (x-y)/.6/cellsize;
    if(abs(index.x) < .8*wc && -index.x > .8*wc*(1.-2.*progress))
    {
        box(y, vec2(.5*cellsize, .8*cellsize), c2.x);
        add(col, c2, col);
    }
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	a = iResolution.x/iResolution.y;
    vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
    
    vec3 old = texture(iChannel0, uv).rgb;
    
    vec4 col, c1;
    window(uv, vec2(.4,.2), old, 0., col);
    
    progressbar(uv+.05*c.yx, .3, .2+.2*sin(iTime), c1);
    add(col, c1, col);
    progressbar(uv+.1*c.yx, .3, .5+.5*sin(iTime), c1);
    add(col, c1, col);
	progressbar(uv+.15*c.yx, .3, .9+.1*sin(iTime), c1);
    add(col, c1, col);

    fragColor = vec4(col.gba, 1.);
}
