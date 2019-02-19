// Read short value from texture at index off
float rshort(float off)
{
    // Parity of offset determines which byte is required.
    float hilo = mod(off, 2.);
    // Find the pixel offset your data is in (2 unsigned shorts per pixel).
    off *= .5;
    // - Determine texture coordinates.
    //     offset = i*iFontWidth+j for (i,j) in [0,iFontWidth]^2
    //     floor(offset/iFontWidth) = floor((i*iFontwidth+j)/iFontwidth)
    //                              = floor(i)+floor(j/iFontWidth) = i
    //     mod(offset, iFontWidth) = mod(i*iFontWidth + j, iFontWidth) = j
    // - For texture coordinates (i,j) has to be rescaled to [0,1].
    // - Also we need to add an extra small offset to the texture coordinate
    //   in order to always "hit" the right pixel. Pixel width is
    //     1./iFontWidth.
    //   Half of it is in the center of the pixel.
    vec2 ind = (vec2(mod(off, iFontWidth), floor(off/iFontWidth))+.05)/iFontWidth;
    // Get 4 bytes of data from the texture
    vec4 block = texture(iFont, ind);
    // Select the appropriate word
    vec2 data = mix(block.rg, block.ba, hilo);
    // Convert bytes to unsigned short. The lower bytes operate on 255,
    // the higher bytes operate on 65280, which is the maximum range 
    // of 65535 minus the lower 255.
    return round(dot(vec2(255., 65280.), data));
}
