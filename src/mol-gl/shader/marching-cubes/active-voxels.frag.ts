export default `
precision highp float;
precision highp sampler2D;

uniform sampler2D tTriCount;
uniform sampler2D tVolumeData;

uniform float uIsoValue;
uniform vec3 uGridDim;
uniform vec3 uGridTexDim;
uniform vec2 uScale;

// cube corners
const vec3 c0 = vec3(0., 0., 0.);
const vec3 c1 = vec3(1., 0., 0.);
const vec3 c2 = vec3(1., 1., 0.);
const vec3 c3 = vec3(0., 1., 0.);
const vec3 c4 = vec3(0., 0., 1.);
const vec3 c5 = vec3(1., 0., 1.);
const vec3 c6 = vec3(1., 1., 1.);
const vec3 c7 = vec3(0., 1., 1.);

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
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / (texDim / uScale);
    // return texture2D(tex, coord + 0.5 / texDim);
    return texture2D(tex, coord);
}

vec4 voxel(vec3 pos) {
    return texture3dFrom2dNearest(tVolumeData, pos, uGridDim, uGridTexDim.xy);
}

void main(void) {
    vec2 uv = gl_FragCoord.xy / uGridTexDim.xy;
    vec3 posXYZ = index3dFrom2d(uv);

    // get MC case as the sum of corners that are below the given iso level
    float c = step(voxel(posXYZ).a, uIsoValue)
        + 2. * step(voxel(posXYZ + c1 / uGridDim).a, uIsoValue)
        + 4. * step(voxel(posXYZ + c2 / uGridDim).a, uIsoValue)
        + 8. * step(voxel(posXYZ + c3 / uGridDim).a, uIsoValue)
        + 16. * step(voxel(posXYZ + c4 / uGridDim).a, uIsoValue)
        + 32. * step(voxel(posXYZ + c5 / uGridDim).a, uIsoValue)
        + 64. * step(voxel(posXYZ + c6 / uGridDim).a, uIsoValue)
        + 128. * step(voxel(posXYZ + c7 / uGridDim).a, uIsoValue);
    c *= step(c, 254.);

    // get total triangles to generate for calculated MC case from triCount texture
    float totalTrianglesToGenerate = texture2D(tTriCount, vec2(intMod(c, 16.), floor(c / 16.)) / 16.).a;
    gl_FragColor = vec4(vec3(floor(totalTrianglesToGenerate * 255.0 + 0.5) * 3.0), c);

    // gl_FragColor = vec4(255.0, 0.0, 0.0, voxel(posXYZ + c4 / uGridDim).a * 255.0);
    // gl_FragColor = vec4(255.0, 0.0, 0.0, voxel(posXYZ).a * 255.0);

    // vec2 uv = vCoordinate;
    // uv = gl_FragCoord.xy / uGridTexDim.xy;

    // if (uv.y < 0.91) discard;
    // gl_FragColor = vec4(vCoordinate * 255.0, 0.0, 255.0);
    // gl_FragColor = vec4(250.0, 0.0, 0.0, 255.0);
}
`;