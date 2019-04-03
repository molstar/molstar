precision mediump float;
precision mediump sampler2D;

uniform sampler2D tTriCount;
uniform sampler2D tVolumeData;

uniform float uIsoValue;
uniform vec3 uGridDim;
uniform vec3 uGridTexDim;

varying vec2 vCoordinate;

vec3 index3dFrom2d(vec2 coord) {
    vec2 gridTexPos = coord * uGridTexDim.xy;
    vec2 columnRow = floor(gridTexPos / uGridDim.xy);
    vec2 posXY = gridTexPos - columnRow * uGridDim.xy;
    float posZ = columnRow.y * floor(uGridTexDim.x / uGridDim.x) + columnRow.x;
    vec3 posXYZ = vec3(posXY, posZ) / uGridDim;
    return posXYZ;
}

float intDiv(float a, float b) { return float(int(a) / int(b)); }
float intMod(float a, float b) { return a - b * float(int(a) / int(b)); }

vec4 texture3dFrom2dNearest(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice = floor(pos.z * gridDim.z + 0.5); // round to nearest z-slice
    float column = intMod(zSlice * gridDim.x, texDim.x) / gridDim.x;
    float row = floor(intDiv(zSlice * gridDim.x, texDim.x));
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / texDim;
    return texture2D(tex, coord);
}

vec4 voxel(vec3 pos) {
    return texture3dFrom2dNearest(tVolumeData, pos, uGridDim, uGridTexDim.xy);
}

void main(void) {
    vec3 posXYZ = index3dFrom2d(vCoordinate);

    // get MC case as the sum of corners that are below the given iso level
    float c = step(voxel(posXYZ).a, uIsoValue)
        + 2. * step(voxel(posXYZ + vec3(1., 0., 0.) / uGridDim).a, uIsoValue)
        + 4. * step(voxel(posXYZ + vec3(1., 1., 0.) / uGridDim).a, uIsoValue)
        + 8. * step(voxel(posXYZ + vec3(0., 1., 0.) / uGridDim).a, uIsoValue)
        + 16. * step(voxel(posXYZ + vec3(0., 0., 1.) / uGridDim).a, uIsoValue)
        + 32. * step(voxel(posXYZ + vec3(1., 0., 1.) / uGridDim).a, uIsoValue)
        + 64. * step(voxel(posXYZ + vec3(1., 1., 1.) / uGridDim).a, uIsoValue)
        + 128. * step(voxel(posXYZ + vec3(0., 1., 1.) / uGridDim).a, uIsoValue);
    c *= step(c, 254.);

    // get total triangles to generate for calculated MC case from triCount texture
    float totalTrianglesToGenerate = texture2D(tTriCount, vec2(intMod(c, 16.), floor(c / 16.)) / 16.).a;
    gl_FragColor = vec4(vec3(totalTrianglesToGenerate * 255.0 * 3.0), c);
}