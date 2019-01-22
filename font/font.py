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
# circlesegment: [x, y, r, plow, phigh]
# line: [x1, y1, x2, y2 ]

def glyph(char):
    pi = numpy.pi;
    circles = []
    segments = []
    lines = []
    if char == 'a':
        segments = [ [ 0., 0., 1., 13.*pi/8., 11.*pi/8.] ]
        lines = [ [ -.5, 0., .5, 0. ] ]
    return [ circles, segments, lines ]
        
