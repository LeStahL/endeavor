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

import struct
import font

# Pack everything as float. If executable size is a problem, this can be optimized slightly.
# Pack alignment:
# nglyphs
# for glyph in glyphs
#     ordinal
#     offset
# for glyph in glyphs
#     nlines
#     for line in lines
#         x1 y1 x2 y2
#     ncircles
#     for circle in circles
#         xc yc r
#     ncirclesegments
#     for segment in circlesegments
#         xc yc r phi0 phi1
# nstrings
# for string in nstrings
#     offset
# for string in nstrings
#     for char in string
#         char

# Read string database
strings = None
with open('strings.txt', 'rt') as f:
    strings = f.readlines()
    f.close()

# Font has only lowercase letters
for i in range(len(strings)):
    strings[i] = strings[i].lower()
    
# Get the list of unique contained ordinals
ordinals = list(set(''.join(strings)))

# Pack number of glyphs
texture = struct.pack('@e', len(ordinals));

# Pack the according glyph table
table = ""

# Pack the glyph data
data = ""



