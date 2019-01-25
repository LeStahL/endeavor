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

import numpy as np

num = 0.337e-11

text = ""

# Sign
sign = np.sign(num)
if sign<0.:
    text += "-"
    num *= -1.;

# Exponent
exp = 0.
for exp in np.arange(-15., 15., 1.):
    if np.floor(num*np.power(10., exp)) != 0.:
        break
exp *= -1.

# Significand
for i in np.arange(exp, exp-7., -1.):
    ca = np.floor(num/pow(10., i))
    num = num - ca*pow(10., i)
    text += str(int(ca))
    if i==exp:
        text += '.'

# Pack exponent
text += 'e'
if exp < 0.:
    text += '-'
    exp *= -1.
ca = np.floor(exp/10.)
text += str(int(ca))
ca = np.floor(exp-ca*10.)
text += str(int(ca))

print(text)
