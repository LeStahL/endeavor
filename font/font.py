# Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
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

import numpy

# Specification of data:
# cicle: [x, y, r]
# circlesegment: [x, y, r, plow, phigh], plow/phigh counter-clockwise
# line: [x1, y1, x2, y2 ]

# Glyph data
def glyph(char):
    pi = numpy.pi;
    circles = []
    segments = []
    lines = []
    if char == 'a':
        segments = [ [ 0., 0., 1., -3.*pi/8., 11.*pi/8.] ]
        lines = [ [ -1., 0., 1., 0. ] ]
    elif char == 'b':
        segments = [ [ 0., 0., 1., -7.*pi/8., 7.*pi/8.] ]
        lines =  [ [ -1./3., -1./3., -1./3., 1./3. ], [ -1./3., 0., .5, 0. ] ]
    elif char == 'c':
        segments = [ [ 0., 0., 1., pi/8., 15.*pi/8. ] ]
    elif char == 'd':
        segments = [ [ 0., 0., 1., -7.*pi/8., 7.*pi/8.] ]
        lines = [ [ -1./3,-1./3., -1./3., 1./3.] ]
    elif char == 'e':
        segments = [ [ 0., 0., 1., pi/8., 15.*pi/8. ] ]
        lines = [ [ -1., 0., 0., 0. ] ]
    elif char == 'f':
        segments = [ [ 0., 0., 1., pi/8., pi ] ]
        lines = [ [ -1., 0., 0., 0. ], [ -1., 0., -1., -1. ] ]
    elif char == 'g':
        segments = [ [ 0., 0., 1., pi/4., 2.*pi ] ]
        lines = [ [ 0., 0., 1., 0. ] ]
    elif char == 'h':
        segments = [ [ 0., 0., 1., 5.*pi/8., 11.*pi/8. ], [ 0., 0., 1., -3.*pi/8., 3.*pi/8. ] ]
        lines = [ [ -1./3., 0., 1./3., 0. ] ]
    elif char == 'i':
        segments = [ [ 0., 0., 1., pi/8., 7.*pi/8. ], [ 0., 0., 1., 9.*pi/8., 15.*pi/8. ] ]
        lines = [ [ 0., 1., 0., -1. ] ]
    elif char == 'j':
        #segments = [ [ 0., 0., 1., pi/8., 7.*pi/8. ], [ 0., 0., 1., 9.*pi/8., 3.*pi/2. ] ]
        #lines = [ [ 0., -1., 0., 1./3. ] ]
        segments = [ [ 0., 0., 1., 9.*pi/8., 2.*pi ] ]
        lines = [ [ 1., 0., 1., 1. ] ]
    elif char == 'k':
        segments = [ [ 0., 0., 1., pi/8., 15.*pi/8. ] ]
        lines = [ [ -1., -1., -1., 1. ] ]
    elif char == 'l':
        segments = [ [ 0., 0., 1., pi, 15.*pi/8. ] ]
        lines = [ [ -1., 0., -1., 1. ] ]
    elif char == 'm':
        segments = [ [ 0., 0., 1., -3.*pi/8., 11.*pi/8. ] ]
        lines = [ [ 0., 0., 0., 1. ] ]
    elif char == 'n':
        segments = [ [ 0., 0., 1., -3.*pi/8., 11.*pi/8. ] ]
    elif char == 'o':
        circles = [ [ 0., 0., 1. ] ]
    elif char == 'p':
        segments = [ [ 0., 0., 1., 0., 7.*pi/8. ] ]
        lines = [ [ -1./3., -1., -1./3., 1./3. ],  [ -1./3., 0., 1., 0. ] ]
    elif char == 'q':
        circles = [ [ 0., 0., 1. ] ]
        lines = [ [1./3., -1./3., 1., -1.] ]
    elif char == 'r':
        segments = [ [ 0., 0., 1., pi/8., pi ] ]
        lines = [ [ -1., 0., -1., -1. ] ]
    elif char == 's':
        segments = [ [ 0., 0., 1., pi/8., 7.*pi/8. ], [ 0., 0., 1., 9.*pi/8., 15.*pi/8. ] ]
        lines = [ [ -0.9238795325112868, 0.3826834323650899, 0.9238795325112868, -0.3826834323650899 ] ]
    elif char == 't':
        segments = [ [ 0., 0., 1., pi, 15.*pi/8. ] ]
        lines = [ [ -1., 0., -1., 1. ], [ -1., 0., 0., 0. ] ]
    elif char == 'u':
        segments = [ [ 0., 0., 1., -pi, 0. ] ]
        lines = [ [ -1., 0., -1., 1. ], [ 1., 0., 1., 1. ] ]
    elif char == 'v':
        segments = [ [ 0., 0., 1., -11.*pi/8., 3.*pi/8. ] ]
    elif char == 'w':
        segments = [ [ 0., 0., 1., -11.*pi/8., 3.*pi/8. ] ]
        lines = [ [ 0., -1., 0., 0. ] ]
    elif char == 'x':
        segments = [ [ -1., 0., 1., -pi/2., pi/2. ], [ 1., 0., 1., pi/2., 3.*pi/2. ] ]
    elif char == 'y':
        segments = [ [ -1., 0., 1., 0., pi/2. ], [ 1., 0., 1., pi/2., pi ] ]
        lines = [ [ 0., 0., 0., -1. ] ]
    elif char == 'z':
        segments = [ [ 0., 0., 1., pi/8., 7.*pi/8. ], [ 0., 0., 1., 9.*pi/8., 15.*pi/8. ] ]
        lines = [ [ -0.9238795325112868, -0.3826834323650899, 0.9238795325112868, 0.3826834323650899 ] ]
    elif char == '0':
        circles = [ [ 0., 0., 1. ] ]
        lines = [ [ -1./3, -1./3., 1./3., 1./3. ] ]
    elif char == '1':
        segments = [ [ 0., 0., 1., 0., 7.*pi/8. ] ]
        lines = [ [ 1., 0., 1., -1. ] ]
    elif char == '2':
        segments = [ [ 0., 0., 1., -pi/2., 7.*pi/8. ] ]
        lines = [ [ -.5, -1., 1., -1. ] ]
    elif char == '3':
        segments = [ [ 0., 0., 1., -7.*pi/8., 7.*pi/8. ] ]
        lines = [ [ 0., 0., 1./3., 0. ] ]
    elif char == '4':
        segments = [ [ 0., 0., 1., -11.*pi/8., 0. ] ]
        lines = [ [ 1., 0., 1., -1. ] ]
    elif char == '5':
        segments = [ [ 0., 0., 1., -7.*pi/8., 7.*pi/8. ] ]
        lines += [ [ -0.9238795325112868, 0.3826834323650899, -0.9238795325112868, 1. ], [ -0.9238795325112868, 1., 1., 1. ] ]
    elif char == '6':
        circles = [ [ 0., 0., 1. ] ]
        lines = [ [ 0., 1., 1., 1. ] ]
    elif char == '7':
        segments = [ [ 0., 0., 1., 0., 7.*pi/8. ] ]
        lines = [ [ 1., 0., 1., -1. ], [ 0., 0., 1./3., 0. ] ]
    elif char == '8':
        circles = [ [ 0., 0., 1. ] ]
        lines = [ [ -1./3., 0., 1./3., 0. ] ]
    elif char == '9':
        circles = [ [ 0., 0., 1. ] ]
        lines = [ [ 1., -1., 1., 0. ] ]
    elif char == '.':
        circles = [ [ -1., -1., 0. ] ]
    elif char == ':':
        circles = [ [ -1., -.5, 0. ], [-1., .5, 0.] ]
    elif char == ',':
        lines = [ [ -1., -1., -2./3., -2./3. ] ]
    elif char == '-':
        lines = [ [ -2./3., 0., 2./3., 0. ] ]
    elif char == '+':
        lines = [ [ -2./3., 0., 2./3., 0. ], [ 0., -2./3., 0., 2./3.] ]
    elif char == '?':
        segments = [ [ 0., 0., 1., 0., 7.*pi/8. ] ]
        lines = [ [ 1., 0., 1., -1./3. ] ]
        circles = [ [ 1., -1., 0. ] ]
    elif char == '!':
        lines = [ [ 0., 1., 0., -1./3. ] ]
        circles = [ [ 0., -1., 0. ] ]
    elif char == '/':
        lines = [ [ -1., -1., 1., 1. ] ]
        
    return [ circles, segments, lines ]
        
# Length of glyph data in texture (for building index with offsets)
def pack_length(char):
    gly = glyph(char)
    circles = gly[0]
    segments = gly[1]
    lines = gly[2]    
    return 3 + len(circles)*3 + len(segments)*5 + len(lines)*4
