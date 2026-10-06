export const activeVoxels_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D tTriCount;
uniform sampler2D tVolumeData;

uniform float uIsoValue;
uniform vec3 uGridDim;
uniform vec3 uGridTexDim;
uniform vec2 uScale;

#include common

// cube corners (excluding origin)
const vec3 c1 = vec3(1., 0., 0.);
const vec3 c2 = vec3(1., 1., 0.);
const vec3 c3 = vec3(0., 1., 0.);
const vec3 c4 = vec3(0., 0., 1.);
const vec3 c5 = vec3(1., 0., 1.);
const vec3 c6 = vec3(1., 1., 1.);
const vec3 c7 = vec3(0., 1., 1.);

vec3 index3dFrom2d(vec2 coord) {
    vec2 gridTexPos = coord * uGridTexDim.xy;
    vec2 columnRow = ivec2Div(gridTexPos, uGridDim.xy);
    vec2 posXY = gridTexPos - columnRow * uGridDim.xy;
    float posZ = columnRow.y * intDiv(uGridTexDim.x, uGridDim.x) + columnRow.x;
    return vec3(posXY, posZ);
}

vec4 texture3dFrom2dNearest(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice = floor(pos.z * gridDim.z + 0.5); // round to nearest z-slice
    float column = intDiv(intMod(zSlice * gridDim.x, texDim.x), gridDim.x);
    float row = intDiv(zSlice * gridDim.x, texDim.x);
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / (texDim / uScale);
    return texture2D(tex, coord);
}

float voxelValue(vec3 pos) {
    pos = min(max(vec3(0.0), pos), uGridDim - vec3(1.0));
    vec4 v = texture3dFrom2dNearest(tVolumeData, pos / uGridDim, uGridDim, uGridTexDim.xy);
    #ifdef dValueChannel_red
        return v.r;
    #else
        return v.a;
    #endif
}

void main(void) {
    vec2 uv = gl_FragCoord.xy / uGridTexDim.xy;
    vec3 posXYZ = index3dFrom2d(uv);

    // get MC case as the sum of corners that are below the given iso level
    float c = step(voxelValue(posXYZ), uIsoValue)
        + 2. * step(voxelValue(posXYZ + c1), uIsoValue)
        + 4. * step(voxelValue(posXYZ + c2), uIsoValue)
        + 8. * step(voxelValue(posXYZ + c3), uIsoValue)
        + 16. * step(voxelValue(posXYZ + c4), uIsoValue)
        + 32. * step(voxelValue(posXYZ + c5), uIsoValue)
        + 64. * step(voxelValue(posXYZ + c6), uIsoValue)
        + 128. * step(voxelValue(posXYZ + c7), uIsoValue);
    c *= step(c, 254.);

    // handle out of bounds positions
    posXYZ += 1.0;
    posXYZ.xy += 1.0; // pixel padding (usually ok even if the texture has no padding)
    if (posXYZ.x >= uGridDim.x || posXYZ.y >= uGridDim.y || posXYZ.z >= uGridDim.z)
        c = 0.0;

    // get total triangles to generate for calculated MC case from triCount texture
    float totalTrianglesToGenerate = texture2D(tTriCount, vec2(intMod(c, 16.), floor(c / 16.)) / 16.).a;
    gl_FragColor = vec4(vec3(totalTrianglesToGenerate * 3.0), c / 255.0);
}
`;
