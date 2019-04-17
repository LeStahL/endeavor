# Puffin Cure by Team210 - 64k Demo at Under Construction 2k18
# Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

project(endeavor)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CTRL_ALT_TEST_MINIFIER OFF)

if(WIN32)
    find_program(KKRUNCHY_K7 kkrunchy_k7.exe HINTS ${PROJECT_BINARY_DIR})
    
    if(${CTRL_ALT_TEST_MINIFIER})
        find_program(MINIFIER shader_minifier.exe HINTS ${PROJECT_BINARY_DIR})
        add_custom_target(minification
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMAND ${MINIFIER} "${PROJECT_SOURCE_DIR}/gfx.frag" -o "${PROJECT_SOURCE_DIR}/gfx.h"  --preserve-all-globals --preserve-externals
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/sfx.h" "${PROJECT_SOURCE_DIR}/sfx.frag" --preserve-externals
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/load.h" "${PROJECT_SOURCE_DIR}/gfx/load.frag" --preserve-all-globals --preserve-externals
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/post.h" "${PROJECT_SOURCE_DIR}/gfx/post.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/logo210.h" "${PROJECT_SOURCE_DIR}/gfx/logo210.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/greet.h" "${PROJECT_SOURCE_DIR}/gfx/greet.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/text.h" "${PROJECT_SOURCE_DIR}/gfx/text.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/logosmall.h" "${PROJECT_SOURCE_DIR}/gfx/logosmall.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/nr4.h" "${PROJECT_SOURCE_DIR}/gfx/nr4.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/qm.h" "${PROJECT_SOURCE_DIR}/gfx/qm.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/trip.h" "${PROJECT_SOURCE_DIR}/gfx/trip.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/surface.h" "${PROJECT_SOURCE_DIR}/gfx/surface.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/hangout.h" "${PROJECT_SOURCE_DIR}/gfx/hangout.frag"
            COMMAND ${MINIFIER} -o "${PROJECT_SOURCE_DIR}/gfx/fourtwenty.h" "${PROJECT_SOURCE_DIR}/gfx/fourtwenty.frag"
        )
    else()
        add_custom_target(minification
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/minification/
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/sfx.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/load.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/post.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/logo210.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/greet.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/text.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/logosmall.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/nr4.h\""
#             COMMAND del "\"${PROJECT_SOURCE_DIR}/gfx/qm.h\""
            COMMAND python shader-compressor.py "${PROJECT_SOURCE_DIR}/gfx.frag" -o "${PROJECT_SOURCE_DIR}/gfx.h"  --no-minification 
            COMMAND python shader-compressor.py  -o "${PROJECT_SOURCE_DIR}/sfx.h" "${PROJECT_SOURCE_DIR}/sfx.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/load.h" "${PROJECT_SOURCE_DIR}/gfx/load.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/post.h" "${PROJECT_SOURCE_DIR}/gfx/post.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/logo210.h" "${PROJECT_SOURCE_DIR}/gfx/logo210.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/greet.h" "${PROJECT_SOURCE_DIR}/gfx/greet.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/text.h" "${PROJECT_SOURCE_DIR}/gfx/text.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/logosmall.h" "${PROJECT_SOURCE_DIR}/gfx/logosmall.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/nr4.h" "${PROJECT_SOURCE_DIR}/gfx/nr4.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/qm.h" "${PROJECT_SOURCE_DIR}/gfx/qm.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/trip.h" "${PROJECT_SOURCE_DIR}/gfx/trip.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/surface.h" "${PROJECT_SOURCE_DIR}/gfx/surface.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/hangout.h" "${PROJECT_SOURCE_DIR}/gfx/hangout.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/fourtwenty.h" "${PROJECT_SOURCE_DIR}/gfx/fourtwenty.frag" --no-minification
        )
    endif()
    
    add_custom_target(texture_packing
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/font/
#         COMMAND del "\"${PROJECT_SOURCE_DIR}/font/font.h\""
        COMMAND python pack.py
    )
    
    add_definitions(/Gy /O1 /fp:fast /GR- /GS- /MT)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:demo /SUBSYSTEM:Windows  /OPT:REF  /OPT:ICF=10 /VERBOSE /NODEFAULTLIB") 
    add_executable(endeavor e.c)
    
    # for debug output: User32.lib ucrt.lib msvcrt.lib
    target_link_libraries(endeavor opengl32.lib Winmm.lib User32.lib ucrt.lib msvcrt.lib)
    
    add_dependencies(endeavor minification)
    add_dependencies(endeavor texture_packing)
    add_custom_command(TARGET endeavor COMMAND ${KKRUNCHY_K7} --best ${CMAKE_CURRENT_BINARY_DIR}/Release/endeavor.exe --out ${CMAKE_CURRENT_BINARY_DIR}/Release/endeavorc.exe POST_BUILD)
endif()

