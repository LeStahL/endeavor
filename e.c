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

// TODO: remove
#include <stdio.h>

#include <math.h>

#include "sequence.h"

// OpenGL extensions
PFNGLGETPROGRAMIVPROC glGetProgramiv;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
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
PFNGLACTIVETEXTUREPROC glActiveTexture;


// TODO: remove below
void debug(int shader_handle)
{
    printf("    Debugging shader with handle %d.\n", shader_handle);
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("    FAILED.\n");
        GLint len;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("    Log length: %d\n", len);
        GLchar *CompileLog = (GLchar*)malloc(len*sizeof(GLchar));
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("    Error messages:\n%s\n", CompileLog);
        free(CompileLog);
    }
    else 
        printf("    Shader compilation successful.\n");
}

void debugp(int program)
{
    printf("    Debugging program.\n");
    int compile_status = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("    FAILED.\n");
        GLint len;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
        printf("    Log length: %d\n", len);
        GLchar *CompileLog = (GLchar*)malloc(len*sizeof(GLchar));
        glGetProgramInfoLog(program, len, NULL, CompileLog);
        printf("    Error messages:\n%s\n", CompileLog);
        free(CompileLog);
    }
    else 
        printf("    Program linking successful.\n");
}
// TODO: remove above

// Graphics shader globals
int w = 1920, h = 1080,
    gfx_handle, gfx_program, 
    time_location, resolution_location, 
    font_texture_location, font_width_location,
    sfx_program, sfx_handle, sfx_blockoffset_location, 
    sfx_samplerate_location, sfx_volumelocation, 
    sfx_texs_location,
    sfx_sequence_texture_location, sfx_sequence_width_location,
    gfx_sequence_texture_location, gfx_sequence_width_location,
    gfx_executable_size_location,
    load_program, load_resolution_location, load_time_location,
    load_progress_location,
    post_program, post_resolution_location, post_fsaa_location,
    post_channel0_location,
    fsaa = 25, txaa = 1,
    gfx_fsaa_location, gfx_txaa_location,
    logo210_time_location, logo210_resolution_location,
    logo210_program, logo210_handle,
    greet_time_location, greet_resolution_location,
    greet_program, greet_handle,
    text_time_location, text_resolution_location,
    text_font_width_location, text_executable_size_location,
    text_channel0_location, text_font_location,
    text_program, text_handle,
    logo_small_handle, logo_small_program,
    logo_small_time_location, logo_small_resolution_location,
    nr4_handle,nr4_program,
    nr4_time_location, nr4_resolution_location;
    
float font_texture_width;
    
// Demo globals
double t_start = 0., 
    t_now = 0., 
    t_end = 180.; // TODO: set to sensible end
unsigned int font_texture_handle, sequence_texture_handle;
float executable_size = 0.;
unsigned int loading = 1, music_loading = 0;
int music_block = 0;
unsigned int snd_framebuffer;
unsigned int scene_override = 0, override_index = 0;

// Music shader globals
int sample_rate = 44100, channels = 2;
double duration1 = 312.*.43; //3 min running time
float *smusic1;
int music1_size;
float texs = 256;
int block_size = 256*256, 
    nblocks1;
unsigned int paused = 0;
float progress = 0.;

HWAVEOUT hWaveOut;
WAVEHDR silence_header = {0, 0, 0, 0, 0, 0, 0, 0 };
double t_paused;
HANDLE load_music_thread;
DWORD load_music_thread_id;
HANDLE load_gfx_thread;
DWORD load_gfx_thread_id;

GLuint first_pass_framebuffer = 0, first_pass_texture;

HDC hdc;
HGLRC glrc;
GLenum error;

float t_load_end = 0.;

DWORD WINAPI LoadMusicThread( LPVOID lpParam)
{
#undef VAR_IBLOCKOFFSET
#undef VAR_ISAMPLERATE
#undef VAR_IVOLUME
#undef VAR_ITEXSIZE
#undef VAR_ISEQUENCE
#undef VAR_ISEQUENCEWIDTH
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
    printf("---> SFX shader:\n");
    debug(sfx_handle);
    glAttachShader(sfx_program, sfx_handle);
    glLinkProgram(sfx_program);
    printf("---> SFX program:\n");
    debugp(sfx_program);
    glUseProgram(sfx_program);
    sfx_samplerate_location = glGetUniformLocation(sfx_program, VAR_ISAMPLERATE);
    sfx_blockoffset_location = glGetUniformLocation(sfx_program, VAR_IBLOCKOFFSET);
    sfx_volumelocation = glGetUniformLocation(sfx_program, VAR_IVOLUME);
    sfx_texs_location = glGetUniformLocation(sfx_program, VAR_ITEXSIZE);
    sfx_sequence_texture_location = glGetUniformLocation(sfx_program, VAR_ISEQUENCE);
    sfx_sequence_width_location = glGetUniformLocation(sfx_program, VAR_ISEQUENCEWIDTH);
    printf("++++ SFX shader created.\n");
    
    music_loading = 1;
    progress += .1; //TODO: add better value here as soon as the real time is known
    
    return 0;
}

DWORD WINAPI LoadLogo210Thread( LPVOID lpParam)
{
    // Load gfx shader
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx/logo210.h"
#ifndef VAR_IRESOLUTION
    #define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ITIME
    #define VAR_ITIME "iTime"
#endif
    int logo210_size = strlen(logo210_frag),
        logo210_handle = glCreateShader(GL_FRAGMENT_SHADER);
    logo210_program = glCreateProgram();
    glShaderSource(logo210_handle, 1, (GLchar **)&logo210_frag, &logo210_size);
    glCompileShader(logo210_handle);
    printf("---> Logo 210 shader:\n");
    debug(logo210_handle);
    glAttachShader(logo210_program, logo210_handle);
    glLinkProgram(logo210_program);
    printf("---> Logo 210 program:\n");
    debugp(logo210_program);
    glUseProgram(logo210_program);
    logo210_time_location =  glGetUniformLocation(logo210_program, VAR_ITIME);
    logo210_resolution_location = glGetUniformLocation(logo210_program, VAR_IRESOLUTION);
    printf("++++ Logo 210 shader created.\n");
    
    progress += .1;
    
    return 0;
}

DWORD WINAPI LoadGreetThread( LPVOID lpParam)
{
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx/greet.h"
#ifndef VAR_IRESOLUTION
    #define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ITIME
    #define VAR_ITIME "iTime"
#endif
    int greet_size = strlen(greet_frag),
        greet_handle = glCreateShader(GL_FRAGMENT_SHADER);
    greet_program = glCreateProgram();
    glShaderSource(greet_handle, 1, (GLchar **)&greet_frag, &greet_size);
    glCompileShader(greet_handle);
    printf("---> Greetings shader:\n");
    debug(greet_handle);
    glAttachShader(greet_program, greet_handle);
    glLinkProgram(greet_program);
    printf("---> Greetings program:\n");
    debugp(greet_program);
    glUseProgram(greet_program);
    greet_time_location =  glGetUniformLocation(greet_program, VAR_ITIME);
    greet_resolution_location = glGetUniformLocation(greet_program, VAR_IRESOLUTION);
    printf("++++ Greetings shader created.\n");
    
    progress += .1;
    
    return 0;
}

#include "font/font.h"
DWORD WINAPI LoadTextThread(LPVOID lpParam)
{
#undef VAR_IFONTWIDTH
#undef VAR_IEXECUTABLESIZE
#undef VAR_IRESOLUTION
#undef VAR_ICHANNEL0
#undef VAR_IFONT
#undef VAR_ITIME
#include "gfx/text.h"
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
#ifndef VAR_IEXECUTABLESIZE
    #define VAR_IEXECUTABLESIZE "iExecutableSize"
#endif
#ifndef VAR_ICHANNEL0
    #define VAR_ICHANNEL0 "iChannel0"
#endif
    
    int text_size = strlen(text_frag);
    text_handle = glCreateShader(GL_FRAGMENT_SHADER);
    text_program = glCreateProgram();
    glShaderSource(text_handle, 1, (GLchar **)&text_frag, &text_size);
    glCompileShader(text_handle);
    printf("---> Text shader:\n");
    debug(text_handle);
    glAttachShader(text_program, text_handle);
    glLinkProgram(text_program);
    printf("---> Text program:\n");
    debugp(text_program);
    glUseProgram(text_program);
    text_time_location =  glGetUniformLocation(text_program, VAR_ITIME);
    text_resolution_location = glGetUniformLocation(text_program, VAR_IRESOLUTION);
    text_font_location= glGetUniformLocation(text_program, VAR_IFONT);
    text_font_width_location = glGetUniformLocation(text_program, VAR_IFONTWIDTH);
    text_executable_size_location = glGetUniformLocation(text_program, VAR_IEXECUTABLESIZE);
    text_channel0_location = glGetUniformLocation(text_program, VAR_ICHANNEL0);
    printf("++++ Text shader created.\n");
    
    // Initialize font texture
    printf("font texture width is: %d\n", font_texture_size); // TODO: remove
    glGenTextures(1, &font_texture_handle);
    glBindTexture(GL_TEXTURE_2D, font_texture_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, font_texture_size, font_texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, font_texture);
    
    font_texture_width = font_texture_size;
    
    progress += .1;
    
    return 0;
}

DWORD WINAPI LoadLogoSmallThread( LPVOID lpParam)
{
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx/logosmall.h"
#ifndef VAR_IRESOLUTION
    #define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ITIME
    #define VAR_ITIME "iTime"
#endif
    int logo_small_size = strlen(logosmall_frag),
        logo_small_handle = glCreateShader(GL_FRAGMENT_SHADER);
    logo_small_program = glCreateProgram();
    glShaderSource(logo_small_handle, 1, (GLchar **)&logosmall_frag, &logo_small_size);
    glCompileShader(logo_small_handle);
    printf("---> Endeavour logo shader:\n");
    debug(logo_small_handle);
    glAttachShader(logo_small_program, logo_small_handle);
    glLinkProgram(logo_small_program);
    printf("---> Endeavour logo program:\n");
    debugp(logo_small_program);
    glUseProgram(logo_small_program);
    logo_small_time_location =  glGetUniformLocation(logo_small_program, VAR_ITIME);
    logo_small_resolution_location = glGetUniformLocation(logo_small_program, VAR_IRESOLUTION);
    printf("++++ Endeavour logo shader created.\n");
    
    progress += .1;
    
    return 0;
}

DWORD WINAPI LoadNr4Thread( LPVOID lpParam)
{
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx/nr4.h"
#ifndef VAR_IRESOLUTION
    #define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ITIME
    #define VAR_ITIME "iTime"
#endif
    int nr4_size = strlen(nr4_frag),
        nr4_handle = glCreateShader(GL_FRAGMENT_SHADER);
    nr4_program = glCreateProgram();
    glShaderSource(nr4_handle, 1, (GLchar **)&nr4_frag, &nr4_size);
    glCompileShader(nr4_handle);
    printf("---> Endeavour NR4 shader:\n");
    debug(nr4_handle);
    glAttachShader(nr4_program, nr4_handle);
    glLinkProgram(nr4_program);
    printf("---> Endeavour NR4 program:\n");
    debugp(nr4_program);
    glUseProgram(nr4_program);
    nr4_time_location =  glGetUniformLocation(nr4_program, VAR_ITIME);
    nr4_resolution_location = glGetUniformLocation(nr4_program, VAR_IRESOLUTION);
    printf("++++ Endeavour NR4 shader created.\n");
    
    progress += .1;
    
    return 0;
}

void quad()
{
    glBegin(GL_QUADS);
    glVertex3f(-1,-1,0);
    glVertex3f(-1,1,0);
    glVertex3f(1,1,0);
    glVertex3f(1,-1,0);
    glEnd();
    glFlush();
}

// Pure opengl drawing code, essentially cross-platform
void draw()
{
        
    if(progress == 1.)
        t_load_end = t_now;
    
    if(t_now-t_load_end > 5.)
    {
        loading = 0;
        glUseProgram(0.);
    }
    
    // Render first pass
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    glViewport(0,0,w,h);
    glClear(GL_COLOR_BUFFER_BIT);
    
    if(loading)
    {
        glUseProgram(load_program);
        glUniform1f(load_progress_location, progress);
        glUniform2f(load_resolution_location, w, h);
        glUniform1f(load_time_location, t_now-t_start);
    }
    else //TODO: arrange scenes in the right order
    {
        float t = t_now-t_load_end-5.;
        if(scene_override)
        {
            if(override_index == 1)
            {
                glUseProgram(logo210_program);
                glUniform1f(logo210_time_location, t);
                glUniform2f(logo210_resolution_location, w, h);
            }
            else if(override_index == 2)
            {
                glUseProgram(logo_small_program);
                glUniform1f(logo_small_time_location, t);
                glUniform2f(logo_small_resolution_location, w, h);
            }
            else if(override_index == 3)
            {
                glUseProgram(nr4_program);
                glUniform1f(nr4_time_location, t);
                glUniform2f(nr4_resolution_location, w, h);
            }
            else if(override_index == 4)
            {
                glUseProgram(greet_program);
                glUniform1f(greet_time_location, t);
                glUniform2f(greet_resolution_location, w, h);
            }
        }
        else
        {
            if(t < 15.)
            {
                glUseProgram(logo210_program);
                glUniform1f(logo210_time_location, t);
                glUniform2f(logo210_resolution_location, w, h);
            }
            else if(t<25.)
            {
                glUseProgram(logo_small_program);
                glUniform1f(logo_small_time_location, t);
                glUniform2f(logo_small_resolution_location, w, h);
            }
            else if(t<50.) //TODO: move to end
            {
                glUseProgram(nr4_program);
                glUniform1f(nr4_time_location, t-25.);
                glUniform2f(nr4_resolution_location, w, h);
            }
            else if(t < 5000.)
            {
                glUseProgram(greet_program);
                glUniform1f(greet_time_location, t);
                glUniform2f(greet_resolution_location, w, h);
            }
            else ExitProcess(0);
        }
        quad();
        
        // Render text to buffer
        glUseProgram(text_program);
        glUniform2f(text_resolution_location, w, h);
        glUniform1f(text_font_width_location, font_texture_size);
        glUniform1f(text_executable_size_location, executable_size);
        glUniform1f(text_time_location, t);
        glUniform1i(text_channel0_location, 0);
        glUniform1i(text_font_location, 1);
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, first_pass_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, font_texture_handle);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, font_texture_size, font_texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    }
   
    quad();
    
    // Render second pass (Post processing) to screen
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0,0,w,h);
    
    glUseProgram(post_program);
    glUniform2f(post_resolution_location, w, h);
    glUniform1f(post_fsaa_location, fsaa);
    glUniform1i(post_channel0_location, 0);
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, first_pass_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    
    quad();
    
    glBindTexture(GL_TEXTURE_2D, 0);
    
    
    if(music_loading)
    {
        glBindFramebuffer(GL_FRAMEBUFFER, snd_framebuffer);
        glUseProgram(sfx_program);
        
        if(music_block >= nblocks1) 
        {
            if(!muted)
                music_loading = 0;
        }
        else
        {
            if(!muted)
            {
                printf("Rendering SFX block %d/%d -> %le\n", music_block, nblocks1, .5*(float)music_block/(float)nblocks1);
                double tstart = (double)(music_block*block_size);
                
                glViewport(0,0,texs,texs);
                
                glUniform1f(sfx_volumelocation, 1.);
                glUniform1f(sfx_samplerate_location, (float)sample_rate);
                glUniform1f(sfx_blockoffset_location, (float)tstart);
                glUniform1i(sfx_texs_location, texs);
                glUniform1i(sfx_sequence_texture_location, 0);
                glUniform1f(sfx_sequence_width_location, sequence_texture_size);
                
                glActiveTexture(GL_TEXTURE0);
                glBindTexture(GL_TEXTURE_2D, sequence_texture_handle);
                
                quad();

                glReadPixels(0, 0, texs, texs, GL_RGBA, GL_UNSIGNED_BYTE, smusic1+music_block*block_size);
                glFlush();
                
                unsigned short *buf = (unsigned short*)smusic1;
                short *dest = (short*)smusic1;
                for(int j=2*music_block*block_size; j<2*(music_block+1)*block_size; ++j)
                    dest[j] = (buf[j]-(1<<15));
            }
            
            progress += .5*1./(float)nblocks1;
            
            ++music_block;
        }
    }
}

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
                case VK_SPACE:
                    // pause/unpaused render timer
                    paused = !paused;
                    if(paused)
                        waveOutPause(hWaveOut);
                    else
                        waveOutRestart(hWaveOut);
                    break;
            }
            break;
        
        case WM_RBUTTONDOWN:
            ExitProcess(0);
            break;
            
        default:
            break;
            
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

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
                {
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
                    else if(index == 2)
                    {
                        w = 1024;
                        h = 768;
                    }
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
                case 8: // Full screen Antialiasing
                {
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    fsaa = (index + 1)*(index + 1);
                }
                    break;
                case 9: // Texture buffer size
                {
                    int index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    texs = 128;
                    for(int i=0; i<index; ++i)
                        texs *= 2;
                    block_size = texs*texs;
                }
                    break;
                case 10:
                {
                    override_index = SendMessage(hSender, CB_GETCURSEL, 0, 0);
                    scene_override = override_index > 0;
                } 
                break;
            }
            break;
            
        case WM_CLOSE:
            ExitProcess(0);
            break;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

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
     100, 10, 175, 80, lwnd, (HMENU)5, hInstance,
     NULL);
    
    // Add items to resolution combo box and select full HD
    const char *fullhd = "1920*1080",
        *halfhd = "960*540",
        *normal = "1024*768";
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fullhd)); 
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (halfhd));
    SendMessage(hResolutionComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (normal));
    SendMessage(hResolutionComboBox, CB_SETCURSEL, 0, 0);
    
    // Add mute checkbox
    HWND hMuteCheckbox = CreateWindow(WC_BUTTON, TEXT("Mute"),
                     WS_VISIBLE | WS_CHILD | BS_CHECKBOX,
                     10, 40, 100, 20,        
                     lwnd, (HMENU) 6, hInstance, NULL);
    
    // Add "Antialiasing: " text
    HWND hAntialiasingText = CreateWindow(WC_STATIC, "FSAA: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,65,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add Fullscreen Antialiasing combo box
    HWND hFSAAComboBox= CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 60, 175, 280, lwnd, (HMENU)8, hInstance,
     NULL);
    
    // Populate with entries
    const char *fsaa1= "None",
        *fsaa4 = "4*FSAA",
        *fsaa9 = "9*FSAA",
        *fsaa16 = "16*FSAA",
        *fsaa25 = "25*FSAA";
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa1)); 
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa4));
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa9)); 
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa16));
    SendMessage(hFSAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (fsaa25));
    SendMessage(hFSAAComboBox, CB_SETCURSEL, 4, 0);
    
    // Add "SFX Buffer: " text
    HWND hTXAAText = CreateWindow(WC_STATIC, "SFX Buffer: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,95,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add SFX buffer size combo box
    HWND hTXAAComboBox= CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 90, 175, 280, lwnd, (HMENU)9, hInstance,
     NULL);
    
    // Populate with entries
    const char *buf128= "128^2 px",
        *buf256 = "256^2 px",
        *buf512 = "512^2 px",
        *buf1024 = "1024^2 px";
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf128)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf256));
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf512)); 
    SendMessage(hTXAAComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (buf1024));
    SendMessage(hTXAAComboBox, CB_SETCURSEL, 1, 0);

    // Add "Antialiasing: " text
    HWND hSceneText = CreateWindow(WC_STATIC, "Scene: ", WS_VISIBLE | WS_CHILD | SS_LEFT, 10,125,100,100, lwnd, NULL, hInstance, NULL);
    
    // Add scene selector
    HWND hSceneComboBox = CreateWindow(WC_COMBOBOX, TEXT(""), 
     CBS_DROPDOWN | CBS_HASSTRINGS | WS_CHILD | WS_OVERLAPPED | WS_VISIBLE,
     100, 120, 175, 280, lwnd, (HMENU)10, hInstance,
     NULL);
    
    // Populate with entries
    const char *all_scenes = "All scenes",
        *logo210_scene= "Team210 Logo",
        *logoendeavor_scene = "Endeavour Logo",
        *nr4_scene = "NR4 Graffiti build-up",
        *greet_scene = "Greetings";
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (all_scenes)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (logo210_scene)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (logoendeavor_scene));
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (nr4_scene)); 
    SendMessage(hSceneComboBox,(UINT) CB_ADDSTRING,(WPARAM) 0,(LPARAM) (greet_scene));
    SendMessage(hSceneComboBox, CB_SETCURSEL, 0, 0);
    
    // Add start button
    HWND hwndButton = CreateWindow(WC_BUTTON,"Abfahrt!",WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,185,165,90,90,lwnd,(HMENU)7,hInstance,NULL);
    
    // Show the selector
    ShowWindow(lwnd, TRUE);
    UpdateWindow(lwnd);
    
    MSG msg = { 0 };
    while(GetMessage(&msg, NULL, 0, 0) > 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg); 
    }
    
    printf("Rendering Demo with:\nSound ");
    if(muted)printf("muted");
    else printf("playing");
    printf("\nResolution: %d * %d\n", w, h);
    printf("FSAA: %d*\n", fsaa);

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
    
    hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // OpenGL extensions
    glGetProgramiv = (PFNGLGETPROGRAMIVPROC) wglGetProcAddress("glGetProgramiv");
    glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) wglGetProcAddress("glGetProgramInfoLog");
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
    
    // Load loading bar shader
    printf("++++ Creating Loading bar.\n");
#undef VAR_ITIME
#undef VAR_IPROGRESS
#undef VAR_IRESOLUTION
#include "gfx/load.h"
#ifndef VAR_ITIME
#define VAR_ITIME "iTime"
#endif
#ifndef VAR_IPROGRESS
#define VAR_IPROGRESS "iProgress"
#endif
#ifndef VAR_IRESOLUTION
#define VAR_IRESOLUTION "iResolution"
#endif
    int load_size = strlen(load_frag),
        load_handle = glCreateShader(GL_FRAGMENT_SHADER);
    load_program = glCreateProgram();
    glShaderSource(load_handle, 1, (GLchar **)&load_frag, &load_size);
    glCompileShader(load_handle);
    printf("---> Load shader:\n");
    debug(load_handle);
    glAttachShader(load_program, load_handle);
    glLinkProgram(load_program);
    printf("---> Load Program:\n");
    debugp(load_program);
    glUseProgram(load_program);
    load_progress_location = glGetUniformLocation(load_program, VAR_IPROGRESS);
    load_time_location = glGetUniformLocation(load_program, VAR_ITIME);
    load_resolution_location = glGetUniformLocation(load_program, VAR_IRESOLUTION);
    printf("++++ Loading bar created.\n");
    
    // Load post processing shader
    printf("++++ Creating Post Shader.\n");
#undef VAR_IFSAA
#undef VAR_IRESOLUTION
#undef VAR_ICHANNEL0
#include "gfx/post.h"
#ifndef VAR_IFSAA
#define VAR_IFSAA "iFSAA"
#endif
#ifndef VAR_IRESOLUTION
#define VAR_IRESOLUTION "iResolution"
#endif
#ifndef VAR_ICHANNEL0
#define VAR_ICHANNEL0 "iChannel0"
#endif
    int post_size = strlen(post_frag),
        post_handle = glCreateShader(GL_FRAGMENT_SHADER);
    post_program = glCreateProgram();
    glShaderSource(post_handle, 1, (GLchar **)&post_frag, &post_size);
    glCompileShader(post_handle);
    printf("---> Post shader:\n");
    debug(post_handle);
    glAttachShader(post_program, post_handle);
    glLinkProgram(post_program);
    printf("---> Post Program:\n");
    debugp(post_program);
    glUseProgram(post_program);
    post_channel0_location = glGetUniformLocation(post_program, VAR_ICHANNEL0);
    post_fsaa_location = glGetUniformLocation(post_program, VAR_IFSAA);
    post_resolution_location = glGetUniformLocation(post_program, VAR_IRESOLUTION);
    printf("++++ Post shader created.\n");
    
    // Create framebuffer for rendering first pass to
    glGenFramebuffers(1, &first_pass_framebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, first_pass_framebuffer);
    glGenTextures(1, &first_pass_texture);
    glBindTexture(GL_TEXTURE_2D, first_pass_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, first_pass_texture, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    
    // Generate music framebuffer
    // Initialize sequence texture
    printf("sequence texture width is: %d\n", sequence_texture_size); // TODO: remove
    glGenTextures(1, &sequence_texture_handle);
    glBindTexture(GL_TEXTURE_2D, sequence_texture_handle);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, sequence_texture_size, sequence_texture_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, sequence_texture);
    
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
    
    // Music allocs
    nblocks1 = sample_rate*duration1/block_size+1;
    music1_size = nblocks1*block_size; 
    smusic1 = (float*)malloc(4*music1_size);
    short *dest = (short*)smusic1;
    for(int i=0; i<2*music1_size; ++i)
        dest[i] = 0;
    
    draw();
    SwapBuffers(hdc);

    // Load music shader
    LoadMusicThread(0);
    draw();
    SwapBuffers(hdc);
    
    // Load Logo 210 shader
    LoadLogo210Thread(0);
    draw();
    SwapBuffers(hdc);
    
    LoadGreetThread(0);
    draw();
    SwapBuffers(hdc);
    
    LoadTextThread(0);
    draw();
    SwapBuffers(hdc);
    
    LoadLogoSmallThread(0);
    draw();
    SwapBuffers(hdc);
    
    LoadNr4Thread(0);
    draw();
    SwapBuffers(hdc);
    
    hWaveOut = 0;
    int n_bits_per_sample = 16;
    WAVEFORMATEX wfx = { WAVE_FORMAT_PCM, channels, sample_rate, sample_rate*channels*n_bits_per_sample/8, channels*n_bits_per_sample/8, n_bits_per_sample, 0 };
    waveOutOpen(&hWaveOut, WAVE_MAPPER, &wfx, 0, 0, CALLBACK_NULL);
    
    WAVEHDR header = { smusic1, 4*music1_size, 0, 0, 0, 0, 0, 0 };
    waveOutPrepareHeader(hWaveOut, &header, sizeof(WAVEHDR));
    waveOutWrite(hWaveOut, &header, sizeof(WAVEHDR));
    
    // Main loop
    t_start = 0.;
    while(1)
    {
        while ( PeekMessageA( &msg, NULL, 0, 0, PM_REMOVE ) ) 
        {
            if ( msg.message == WM_QUIT ) {
                return 0;
            }
            TranslateMessage( &msg );
            DispatchMessageA( &msg );
        }
        
        static MMTIME MMTime = { TIME_SAMPLES, 0};
        waveOutGetPosition(hWaveOut, &MMTime, sizeof(MMTIME));
        t_now = ((double)MMTime.u.sample)/( 44100.0);
        
        draw();
        SwapBuffers(hdc);

    }
    return msg.wParam;
}

