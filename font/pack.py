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
nglyphs = len(ordinals)
print("Packing glyph data of: ", ordinals)

# Pack number of glyphs
fmt = '@e'
texture = struct.pack(fmt, float(len(ordinals)));

# Pack the according glyph table
pack_len = 1 + nglyphs
table = ""
for char in ordinals:
    # Pack ordinal
    texture += struct.pack(fmt, float(ord(char)))
    
    # Pack offset
    texture += struct.pack(fmt, float(pack_len))
    
    # Update offset
    pack_len += font.pack_length(char)

# Pack the glyph data
for char in ordinals:
    # Get glyph inlines
    glyph = font.glyph(char)
    
    # Pack number of lines
    lines = glyph[2]
    texture += struct.pack(fmt, float(len(lines)))
    
    # Pack lines
    for line in lines:
        for i in range(4):
            texture += struct.pack(fmt, float(line[i]))
            
    # Pack number of circles
    circles = glyph[0]
    texture += struct.pack(fmt, float(len(circles)))
    
    # Pack circles
    for circle in circles:
        for i in range(3):
            texture += struct.pack(fmt, float(circle[i]))
            
    # Pack number of circle segments
    segments = glyph[1]
    texture += struct.pack(fmt, float(len(segments)))
    
    # Pack segments
    for segment in segments:
        for i in range(5):
            texture += struct.pack(fmt, float(segment[i]))
print("Finished packing texture.")
    
#Generate C header text with texture data
# Write output to c header file or stdout
# Fill last 4-byte-block with zero
length = int(len(texture)) # in bytes
while ((length % 4) != 0):
    texture += bytes(10)
    length += 1
print("Packed font is "+str(length)+" bytes.")


