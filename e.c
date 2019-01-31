/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
 * Copyright (C) 2018 Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// define this to get some extra bytes
//#define SUPER_SMALL

// WIN32 code
#ifdef _MSC_VER

const char *demoname = "Endeavour/Team210";
unsigned int muted = 0.;

int _fltused = 0;

#define ABS(x) ((x)<0?(-x):(x))
#define sign(x) ((x)<0?-1.:1.)

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#include "commctrl.h"
#include <mmsystem.h>
#include <Mmreg.h>

#include <GL/gl.h>
#include "glext.h"

// Standard library and CRT rewrite for saving executable size
void *memset(void *ptr, int value, size_t num)
{
    for(int i=num-1; i>=0; i--)
        ((unsigned char *)ptr)[i] = value;
    return ptr;
}

size_t strlen(const char *str)
{
    int len = 0;
    while(str[len] != '\0') ++len;
    return len;
}

void *malloc( unsigned int size )
{
    return GlobalAlloc(GMEM_ZEROINIT, size) ;
}
#else 

#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>

#include <GL/gl.h>
#include <GL/glx.h>
#include "glext.h"

#include <time.h>
#include <sys/time.h>

#include <dlfcn.h>

// #include <stdlib.h>
#include <string.h>
#endif

// TODO: remove
#include <stdio.h>

#include <math.h>

// fonts
#include "font/font.h"
#include "sequence.h"

// OpenGL extensions
PFNGLGETPROGRAMIVPROC glGetProgramiv;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLCREATESHADERPROC glCreateShader;
PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLUNIFORM2FPROC glUniform2f;
PFNGLUNIFORM1FPROC glUniform1f;
PFNGLGENFRAMEBUFFERSPROC glGenFramebuffers;
PFNGLBINDFRAMEBUFFERPROC glBindFramebuffer;
PFNGLFRAMEBUFFERTEXTURE2DPROC glFramebufferTexture2D;
PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC glNamedRenderbufferStorageEXT;
PFNGLUNIFORM1IPROC glUniform1i;

#ifdef _MSC_VER
PFNGLACTIVETEXTUREPROC glActiveTexture;
#endif

// TODO: remove below
void debug(int shader_handle)
{
    printf("debugging shader.\n");
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len = 4618;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar CompileLog[4618];
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
    }
    else 
        printf("shader compilation successful.\n");
//     Sleep(20000); TODO ADDN
}
// TODO: remove above

//TODO: use full hd again
// Graphics shader globals
// int w = 1366, h = 768,
// int w = 800, h = 450,
int w = 1920, h = 1080,
    gfx_handle, gfx_program, 
    time_location, resolution_location, 
    font_texture_location, font_width_location,
    sfx_program, sfx_blockoffset_location, 
    sfx_samplerate_location, sfx_volumelocation, 
    sfx_texs_location,
    sfx_sequence_texture_location, sfx_sequence_width_location,
    gfx_sequence_texture_location, gfx_sequence_width_location,
    gfx_executable_size_location;
    
// Demo globals
double t_start = 0., 
    t_now = 0., 
    t_end = 180.; // TODO: set to sensible end
unsigned int font_texture_handle, sequence_texture_handle;
float executable_size = 0.;

// Music shader globals
int sample_rate = 44100, channels = 2;
double duration1 = 312.*.43; //3 min running time
float *smusic1;
int music1_size;
#define texs 1024
int block_size = texs*texs;
unsigned int paused = 0;

#ifdef _MSC_VER
HWAVEOUT hWaveOut;
#endif
double t_paused;

// Pure opengl drawing code, essentially cross-platform
void draw()
{
    // End demo when it is over
    if(t_now-t_start > t_end)
#ifdef _MSC_VER
        ExitProcess(0);
#else
        exit(0)
#endif
    
    glUniform1i(font_texture_location, 0);
    glUniform1f(font_width_location, font_texture_size);
    glUniform1i(gfx_sequence_texture_location, 1);
    glUniform1f(gfx_sequence_width_location, sequence_texture_size);
    glUniform1f(gfx_executable_size_location, executable_size);
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, font_texture_handle);
    
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, sequence_texture_handle);
    
    glUniform1f(time_location, t_now-t_start);
    glUniform2f(resolution_location, w, h);
    
    glBegin(GL_QUADS);
    
    glVertex3f(-1,-1,0);
    glVertex3f(-1,1,0);
    glVertex3f(1,1,0);
    glVertex3f(1,-1,0);
    
    glEnd();
    
    glFlush();
}

#ifdef _MSC_VER
LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_KEYDOWN:
            switch(wParam)
            {
                case VK_ESCAPE:
                    ExitProcess(0);
                    break;
#ifndef SUPER_SMALL
                case VK_SPACE:
                    // pause/unpaused render timer
                    paused = !paused;
                    if(paused)
                    {
                        SetTimer(hwnd, 1, USER_TIMER_MAXIMUM, NULL);
                        if(!muted)
                            waveOutPause(hWaveOut);
                        
                        t_paused = t_now;
                    }
                    else
                    {
                        SetTimer(hwnd, 1, 1000./60., NULL);
                        if(!muted)
                            waveOutRestart(hWaveOut);
                        
                        SYSTEMTIME st_now;
                        GetSystemTime(&st_now);            
                        t_now = (float)st_now.wMinute*60.+(float)st_now.wSecond+(float)st_now.wMilliseconds/1000.;
                        t_start += t_now - t_paused;
                    }
                    break;
#endif
            }
            break;
        
        case WM_RBUTTONDOWN:
            ExitProcess(0);
            break;
            
        case WM_TIMER:
            HDC hdc = GetDC(hwnd);
                 
            SYSTEMTIME st_now;
            GetSystemTime(&st_now);            
            t_now = (float)st_now.wMinute*60.+(float)st_now.wSecond+(float)st_now.wMilliseconds/1000.;
            
            draw();
            
            SwapBuffers(hdc);
            
            break;
            
        default:
            break;
            
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
#endif

#ifdef _MSC_VER
LRESULT CALLBACK DialogProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch(uMsg)
    {
        case WM_COMMAND:
            UINT id =  LOWORD(wParam);
            HWND hSender = (HWND)lParam;
            switch(id)
            {
                case 5:
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    if(index == 0)
                    {
                        w = 1920;
                        h = 1080;
                    }
                    else if(index == 1)
                    {
                        w = 960;
                        h = 540;
                    }
                    break;
                case 6:
                    muted = !muted;
                    if(muted)
                        SendMessage(hSender, BM_SETCHECK, BST_CHECKED, 0);
                    else
                        SendMessage(hSender, BM_SETCHECK, BST_UNCHECKED, 0);
                    break;
                case 7:
                    DestroyWindow(hwnd);
                    PostQuitMessage(0);
                    break;
            }
            break;
            
        case WM_CLOSE:
            ExitProcess(0);
            break;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}
#endif

#ifdef _MSC_VER
int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    //TODO: remove
    AllocConsole();
    freopen("CONIN$", "r", stdin);
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    
    // Show executable size on top
#define MAX_PATH 1024
    char szFileName[MAX_PATH];
    GetModuleFileName(NULL, szFileName, MAX_PATH);
    printf("%s\n", szFileName);
    
    FILE *ex = fopen(szFileName, "rb");
    fseek(ex, 0, SEEK_END);
    executable_size = (float)ftell(ex)/1024.;
    
    // Display settings selector
    WNDCLASS wca = { 0 };
    wca.lpfnWndProc   = DialogProc;
    wca.hInstance     = hInstance;
    wca.lpszClassName = L"Settings";
    RegisterClass(&wca);
    HWND lwnd = CreateWindowEx(
        0,                              // Optional window styles.
        L"Settings",                     // Window class
        demoname,    // Window text
        WS_OVERLAPPEDWINDOW,            // Window style

        // Size and position
        200, 200, 300, 300,

        NULL,       // Parent window    
        NULL,       // Menu
        hInstance,  // Instance handle
        NULL        // Additional application data
        );
    
    // Add "Resolution: " text
    HWND hResolutionText = CreateWindow(WC_STATIC, "Resolution: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,15,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add resolution Combo box
    HWND hResolutionComboBox = CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 10, 150, 80, lwnd, (HMENU)5, hInstance,
     NULL);
    
    // Add items to resolution combo box and select full HD
    const char *fullhd = "1920 x 1080",
        *halfhd = "960 x 540";
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fullhd)); 
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (halfhd));
    SendMessage(hResolutionComboBox, CB_SETCURSEL, 0, 0);
    
    // Add mute checkbox
    HWND hMuteCheckbox = CreateWindow(WC_BUTTON, TEXT("Mute"),
                     WS_VISIBLE | WS_CHILD | BS_CHECKBOX,
                     10, 40, 100, 20,        
                     lwnd, (HMENU) 6, hInstance, NULL);
    
    // Add start button
    HWND hwndButton = CreateWindow( 
    WC_BUTTON,  // Predefined class; Unicode assumed 
    "Abfahrt!",      // Button text 
    WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,  // Styles 
    185,         // x position 
    165,         // y position 
    90,        // Button width
    90,        // Button height
    lwnd,     // Parent window
    (HMENU)7,       // No menu.
    hInstance, 
    NULL);      // Pointer not needed.
    
    // Show the selector
    ShowWindow(lwnd, TRUE);
    UpdateWindow(lwnd);
    
    MSG msg;
    while(GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg); 
    }
    
    printf("Rendering Demo with:\nSound ");
    if(muted)printf("muted");
    else printf("playing");
    printf("\nResolution: %d * %d\n", w, h);

    // Display demo window
    CHAR WindowClass[]  = "Team210 Demo Window";
    
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
    wc.lpfnWndProc = &WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO); 
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = NULL;
    wc.lpszClassName = WindowClass;
    wc.hIconSm = NULL;
    
    RegisterClassEx(&wc);
    
    // Get full screen information
    HMONITOR hmon = MonitorFromWindow(0, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi = { sizeof(mi) };
    GetMonitorInfo(hmon, &mi);
    
    // Create the window.
    HWND hwnd = CreateWindowEx(
        0,                                                          // Optional window styles.
        WindowClass,                                                // Window class
        ":: NR4^QM/Team210 :: GO - MAKE A DEMO ::",                                 // Window text
        WS_POPUP | WS_VISIBLE,                                      // Window style
        mi.rcMonitor.left,
        mi.rcMonitor.top,
        mi.rcMonitor.right - mi.rcMonitor.left,
        mi.rcMonitor.bottom - mi.rcMonitor.top,                     // Size and position
        
        NULL,                                                       // Parent window    
        NULL,                                                       // Menu
        hInstance,                                                  // Instance handle
        0                                                           // Additional application data
    );
    
    // Show it
    ShowWindow(hwnd, TRUE);
    UpdateWindow(hwnd);
    
    // Create OpenGL context
    PIXELFORMATDESCRIPTOR pfd =
    {
        sizeof(PIXELFORMATDESCRIPTOR),
        1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    //Flags
        PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
        32,                   // Colordepth of the framebuffer.
        0, 0, 0, 0, 0, 0,
        0,
        0,
        0,
        0, 0, 0, 0,
        24,                   // Number of bits for the depthbuffer
        8,                    // Number of bits for the stencilbuffer
        0,                    // Number of Aux buffers in the framebuffer.
        PFD_MAIN_PLANE,
        0,
        0, 0, 0
    };
    
    HDC hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    HGLRC glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // Draw black screen with loading bar
    glClearColor(0.,0.,0.,1.);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(.2, .2, .2);
    glBegin(GL_QUADS);            
    glVertex2f(-.5f, -.05f);
    glVertex2f(-.5f, .05f);
    glVertex2f(.5f, .05f);
    glVertex2f(.5f, -0.05f);
    glEnd();
    
    SwapBuffers(hdc);
    
    // OpenGL extensions
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) wglGetProcAddress("glGetShaderiv");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) wglGetProcAddress("glGetShaderInfoLog");
    glCreateShader = (PFNGLCREATESHADERPROC) wglGetProcAddress("glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) wglGetProcAddress("glCreateProgram");
    glShaderSource = (PFNGLSHADERSOURCEPROC) wglGetProcAddress("glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) wglGetProcAddress("glCompileShader");
    glAttachShader = (PFNGLATTACHSHADERPROC) wglGetProcAddress("glAttachShader");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) wglGetProcAddress("glLinkProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC) wglGetProcAddress("glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) wglGetProcAddress("glGetUniformLocation");
    glUniform2f = (PFNGLUNIFORM2FPROC) wglGetProcAddress("glUniform2f");
    glUniform1f = (PFNGLUNIFORM1FPROC) wglGetProcAddress("glUniform1f");
    glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) wglGetProcAddress("glGenFramebuffers");
    glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) wglGetProcAddress("glBindFramebuffer");
    glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) wglGetProcAddress("glFramebufferTexture2D");
    glNamedRenderbufferStorageEXT = (PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC) wglGetProcAddress("glNamedRenderbufferStorage");
    glActiveTexture = (PFNGLACTIVETEXTUREPROC) wglGetProcAddress("glActiveTexture");
    glUniform1i = (PFNGLUNIFORM1IPROC) wglGetProcAddress("glUniform1i");
    
    // Set render timer to 60 fps
    UINT_PTR t = SetTimer(hwnd, 1, 1000./60., NULL);
#else
int main(int argc, char **args)
{
    XInitThreads();
    
    Display                 *display;
    Window                  root;
    GLint                   att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    XVisualInfo             *vi;
    Colormap                cmap;
    XSetWindowAttributes    swa;
    Window                  win;
    GLXContext              glc;
    XWindowAttributes       gwa;
    XEvent                  xevent;

    display = XOpenDisplay(NULL);
    root = DefaultRootWindow(display);
    vi = glXChooseVisual(display, 0, att);
    cmap = XCreateColormap(display, root, vi->visual, AllocNone);
    swa.colormap = cmap;
    
    swa.event_mask = ExposureMask | KeyPressMask;
    win = XCreateWindow(display, root, 0, 0, w, h, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
 
    Atom atoms[2] = { XInternAtom(display, "_NET_WM_STATE_FULLSCREEN", True), None };
    XChangeProperty(
        display, 
        win, 
        XInternAtom(display, "_NET_WM_STATE", True),
                    XA_ATOM, 32, PropModeReplace,(unsigned char*) atoms, 1
    );
    XMapWindow(display, win);
    glc = glXCreateContext(display, vi, NULL, GL_TRUE);
    glXMakeCurrent(display, win, glc);

    // OpenGL extensions
    void *gl = (void*)dlopen("libGL.so", RTLD_LAZY | RTLD_GLOBAL);
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) dlsym(gl, "glGetProgramiv");
    glGetShaderiv = (PFNGLGETSHADERIVPROC) dlsym(gl, "glGetShaderiv");
    glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) dlsym(gl, "glGetShaderInfoLog");
    glCreateShader = (PFNGLCREATESHADERPROC) dlsym(gl, "glCreateShader");
    glCreateProgram = (PFNGLCREATEPROGRAMPROC) dlsym(gl, "glCreateProgram");
    glShaderSource = (PFNGLSHADERSOURCEPROC) dlsym(gl, "glShaderSource");
    glCompileShader = (PFNGLCOMPILESHADERPROC) dlsym(gl, "glCompileShader");
    glAttachShader = (PFNGLATTACHSHADERPROC) dlsym(gl, "glAttachShader");
    glLinkProgram = (PFNGLLINKPROGRAMPROC) dlsym(gl, "glLinkProgram");
    glUseProgram = (PFNGLUSEPROGRAMPROC) dlsym(gl, "glUseProgram");
    glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) dlsym(gl, "glGetUniformLocation");
    glUniform2f = (PFNGLUNIFORM2FPROC) dlsym(gl, "glUniform2f");
    glUniform1f = (PFNGLUNIFORM1FPROC) dlsym(gl, "glUniform1f");
    glGenFramebuffers = (PFNGLGENFRAMEBUFFERSPROC) dlsym(gl, "glGenFramebuffers");
    glBindFramebuffer = (PFNGLBINDFRAMEBUFFERPROC) dlsym(gl, "glBindFramebuffer");
    glFramebufferTexture2D = (PFNGLFRAMEBUFFERTEXTURE2DPROC) dlsym(gl, "glFramebufferTexture2D");
    glNamedRenderbufferStorageEXT = (PFNGLNAMEDRENDERBUFFERSTORAGEEXTPROC) dlsym(gl, "glNamedRenderbufferStorage");
    glUniform1i = (PFNGLUNIFORM1IPROC) dlsym(gl, "glUniform1i");
    
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t_start = (double)tv.tv_sec+(double)tv.tv_usec/1.e6;
    
#endif
    
     // Load sfx shader
#undef VAR_IBLOCKOFFSET
#undef VAR_ISAMPLERATE
#undef VAR_IVOLUME
#include "sfx.h"
#ifndef VAR_IVOLUME
    #define VAR_IVOLUME "iVolume"
#endif
#ifndef VAR_ISAMPLERATE
    #define VAR_ISAMPLERATE "iSampleRate"
#endif
#ifndef VAR_IBLOCKOFFSET
    #define VAR_IBLOCKOFFSET "iBlockOffset"
#endif
#ifndef VAR_ITEXSIZE
    #define VAR_ITEXSIZE "iTexSize"
#endif
#ifndef VAR_ISEQUENCE
    #define VAR_ISEQUENCE "iSequence"
#endif
#ifndef VAR_ISEQUENCEWIDTH
    #define VAR_ISEQUENCEWIDTH "iSequenceWidth"
#endif
    int sfx_size = strlen(sfx_frag),
        sfx_handle = glCreateShader(GL_FRAGMENT_SHADER);
    sfx_program = glCreateProgram();
    glShaderSource(sfx_handle, 1, (GLchar **)&sfx_frag, &sfx_size);
    glCompileShader(sfx_handle);
    debug(sfx_handle);
    printf("sfx_debugged");
    glAttachShader(sfx_program, sfx_handle);
    glLinkProgram(sfx_program);
    glUseProgram(sfx_program);
    sfx_samplerate_location = glGetUniformLocation(sfx_program, VAR_ISAMPLERATE);
    sfx_blockoffset_location = glGetUniformLocation(sfx_program, VAR_IBLOCKOFFSET);
    sfx_volumelocation = glGetUniformLocation(sfx_program, VAR_IVOLUME);
    sfx_texs_location = glGetUniformLocation(sfx_program, VAR_ITEXSIZE);
    sfx_sequence_texture_location = glGetUniformLocation(sfx_program, VAR_ISEQUENCE);
    sfx_sequence_width_location = glGetUniformLocation(sfx_program, VAR_ISEQUENCEWIDTH);
    
    // Initialize sequence texture
    printf("sequence texture width is: %d\n", sequence_texture_size); // TODO: remove
    glGenTextures(1, &sequence_texture_handle);
    glBindTexture(GL_TEXTURE_2D, sequence_texture_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, sequence_texture_size, sequence_texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, sequence_texture);
    
    int nblocks1 = sample_rate*duration1/block_size+1;
    music1_size = nblocks1*block_size; 
    smusic1 = (float*)malloc(4*music1_size);
    
    unsigned int snd_framebuffer;
    glGenFramebuffers(1, &snd_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
    glPixelStorei(GL_PACK_ALIGNMENT,  4);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
    
    unsigned int snd_texture;
    glGenTextures(1, &snd_texture);
    glBindTexture(GL_TEXTURE_2D, snd_texture);
    glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, texs, texs, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, snd_texture, 0);
    
    // Render sfx
    printf("nblocks: %d\n", nblocks1);
    for(int i=0; i<nblocks1; ++i)
    {
        double tstart = (double)(i*block_size);
        
        glViewport(0,0,texs,texs);
        
        glUniform1f(sfx_volumelocation, 1.);
        glUniform1f(sfx_samplerate_location, (float)sample_rate);
        glUniform1f(sfx_blockoffset_location, (float)tstart);
        glUniform1i(sfx_texs_location, texs);
        glUniform1i(sfx_sequence_texture_location, 0);
        glUniform1f(sfx_sequence_width_location, sequence_texture_size);
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, sequence_texture_handle);
        
//         glActiveTexture(GL_TEXTURE1);
        
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();

        glFlush();

        glReadPixels(0, 0, texs, texs, GL_RGBA, GL_UNSIGNED_BYTE, smusic1+i*block_size);
                
        printf("Block: %d/%d\n", i, nblocks1);
    }
    glFlush();
    
    unsigned short *buf = (unsigned short*)smusic1;
    short *dest = (short*)smusic1;
    for(int j=0; j<2*nblocks1*block_size; ++j)
    {
        dest[j] = (buf[j]-(1<<15));
    }
    
    // Reset everything for rendering gfx again
    glViewport(0, 0, w, h);
    glUseProgram(0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    
    glClearColor(0.,0.,0.,1.);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(.2, .2, .2);
    glBegin(GL_QUADS);            
    glVertex2f(-.5f, -.05f);
    glVertex2f(-.5f, .05f);
    glVertex2f(.5f, .05f);
    glVertex2f(.5f, -0.05f);
    glEnd();
    
    glColor3f(1., 1., 1.);
    glBegin(GL_QUADS);            
    glVertex2f(-.5f, -.05f);
    glVertex2f(-.5f, .05f);
    glVertex2f(0.f, .05f);
    glVertex2f(0.f, -0.05f);
    glEnd();
    
    SwapBuffers(hdc);
    
    // Load gfx shader
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#undef VAR_IFONT
#undef VAR_IFONTWIDTH
#undef VAR_ISEQUENCE
#undef VAR_ISEQUENCEWIDTH
#include "gfx.h"
#ifndef VAR_IRESOLUTION
    #define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ITIME
    #define VAR_ITIME "iTime"
#endif
#ifndef VAR_IFONT
    #define VAR_IFONT "iFont"
#endif
#ifndef VAR_IFONTWIDTH
    #define VAR_IFONTWIDTH "iFontWidth"
#endif
#ifndef VAR_ISEQUENCE
    #define VAR_ISEQUENCE "iSequence"
#endif
#ifndef VAR_ISEQUENCEWIDTH
    #define VAR_ISEQUENCEWIDTH "iSequenceWidth"
#endif
#ifndef VAR_IEXECUTABLESIZE
    #define VAR_IEXECUTABLESIZE "iExecutableSize"
#endif
    int gfx_size = strlen(gfx_frag),
        gfx_handle = glCreateShader(GL_FRAGMENT_SHADER);
    gfx_program = glCreateProgram();
    glShaderSource(gfx_handle, 1, (GLchar **)&gfx_frag, &gfx_size);
    glCompileShader(gfx_handle);
    debug(gfx_handle);
    glAttachShader(gfx_program, gfx_handle);
    glLinkProgram(gfx_program);
    glUseProgram(gfx_program);
    time_location =  glGetUniformLocation(gfx_program, VAR_ITIME);
    resolution_location = glGetUniformLocation(gfx_program, VAR_IRESOLUTION);
    font_texture_location = glGetUniformLocation(gfx_program, VAR_IFONT);
    font_width_location = glGetUniformLocation(gfx_program, VAR_IFONTWIDTH);
    gfx_sequence_texture_location = glGetUniformLocation(gfx_program, VAR_ISEQUENCE);
    gfx_sequence_width_location = glGetUniformLocation(gfx_program, VAR_ISEQUENCEWIDTH);
    gfx_executable_size_location = glGetUniformLocation(gfx_program, VAR_IEXECUTABLESIZE);
    
    glUseProgram(gfx_program);
    glViewport(0, 0, w, h);
    
    glClearColor(0.,0.,0.,1.);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(1., 1., 1.);
    glBegin(GL_QUADS);            
    glVertex2f(-.5f, -.05f);
    glVertex2f(-.5f, .05f);
    glVertex2f(0.5f, .05f);
    glVertex2f(0.5f, -0.05f);
    glEnd();
    
    SwapBuffers(hdc);
    
    // Initialize font texture
    printf("font texture width is: %d\n", font_texture_size); // TODO: remove
    glGenTextures(1, &font_texture_handle);
    glBindTexture(GL_TEXTURE_2D, font_texture_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, font_texture_size, font_texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, font_texture);
    
    // Play sound
#ifdef _MSC_VER
    hWaveOut = 0;
    int n_bits_per_sample = 16;
	WAVEFORMATEX wfx = { WAVE_FORMAT_PCM, channels, sample_rate, sample_rate*channels*n_bits_per_sample/8, channels*n_bits_per_sample/8, n_bits_per_sample, 0 };
	waveOutOpen(&hWaveOut, WAVE_MAPPER, &wfx, 0, 0, CALLBACK_NULL);
	
	WAVEHDR header = { smusic1, 4*music1_size, 0, 0, 0, 0, 0, 0 };
	waveOutPrepareHeader(hWaveOut, &header, sizeof(WAVEHDR));
	if(!muted)
        waveOutWrite(hWaveOut, &header, sizeof(WAVEHDR));
#endif
    
    // Main loop
#ifdef _MSC_VER
    // Get start time for relative time sync
    SYSTEMTIME st_start;
    GetSystemTime(&st_start);
    t_start = (float)st_start.wMinute*60.+(float)st_start.wSecond+(float)st_start.wMilliseconds/1000.;
    
//     MSG msg;
    while(GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg); 
    }
    
    return msg.wParam;
}
#else
    int x_file_descriptor = ConnectionNumber(display);
    fd_set x_file_descriptor_set;
    
    // Main loop
    while(1)
    {
        // Exit demo if any key is pressed.
        while(XPending(display))
        {
            XNextEvent(display, &xevent);
            if(xevent.type == KeyPress) 
                exit(0);
        }
        
        FD_ZERO(&x_file_descriptor_set);
        FD_SET(x_file_descriptor, &x_file_descriptor_set);
        
        struct timeval t;
        t.tv_usec = 1.e6/60.;
        t.tv_sec = 0;
        
        int num_ready_fds = select(x_file_descriptor + 1, &x_file_descriptor_set, NULL, NULL, &t);
        if (num_ready_fds == 0)    
        {            
            struct timeval tv_now;
            gettimeofday(&tv_now, NULL);
            t_now = (double)tv_now.tv_sec+(double)tv_now.tv_usec/1.e6;

            draw();
            
            glXSwapBuffers(display, win);
        }
    }
    return 0;
}
#endif
