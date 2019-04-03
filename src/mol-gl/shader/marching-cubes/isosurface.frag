precision highp float;
precision highp sampler2D;

uniform sampler2D tActiveVoxelsPyramid;
uniform sampler2D tActiveVoxelsBase;
uniform sampler2D tActiveVoxelsTotal;
uniform sampler2D tVolumeData;
uniform sampler2D tTriIndices;

uniform float uIsoValue;
uniform float uLevels;
uniform float uSize;

uniform vec3 uGridDim;
uniform vec3 uGridTexDim;
uniform mat4 uGridTransform;

// scale to volume data coord
uniform vec2 uScale;

varying vec2 vCoordinate;

// const vec3 c0  = vec3(0., 0., 0.);
// const vec3 c1  = vec3(1., 0., 0.);
// const vec3 c2  = vec3(1., 1., 0.);
// const vec3 c3  = vec3(0., 1., 0.);
// const vec3 c4  = vec3(0., 0., 1.);
// const vec3 c5  = vec3(1., 0., 1.);
// const vec3 c6  = vec3(1., 1., 1.);
// const vec3 c7  = vec3(0., 1., 1.);

const vec3 p0 = vec3(1., 0., 0.);
const vec3 p1 = vec3(1., 1., 0.);
const vec3 p2 = vec3(0., 1., 0.);
const vec3 p3 = vec3(0., 0., 1.);
const vec3 p4 = vec3(1., 0., 1.);
const vec3 p5 = vec3(1., 1., 1.);
const vec3 p6 = vec3(0., 1., 1.);

vec3 index3dFrom2d(vec2 coord) {
    vec2 gridTexPos = coord * uGridTexDim.xy;
    vec2 columnRow = floor(gridTexPos / uGridDim.xy);
    vec2 posXY = gridTexPos - columnRow * uGridDim.xy;
    float posZ = columnRow.y * floor(uGridTexDim.x / uGridDim.x) + columnRow.x;
    vec3 posXYZ = vec3(posXY, posZ);
    return posXYZ;
}

float intDiv(float a, float b) { return float(int(a) / int(b)); }
float intMod(float a, float b) { return a - b * float(int(a) / int(b)); }

vec4 texture3dFrom2dNearest(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice = floor(pos.z * gridDim.z + 0.5); // round to nearest z-slice
    float column = intMod(zSlice * gridDim.x, texDim.x) / gridDim.x;
    float row = floor(intDiv(zSlice * gridDim.x, texDim.x));
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / texDim;
    return texture2D(tex, coord + 0.5 / texDim);
}

vec4 voxel(vec3 pos) {
    return texture3dFrom2dNearest(tVolumeData, pos, uGridDim, uGridTexDim.xy);
}

void main(void) {
    // get 1D index
    float vI = dot(floor(uSize * vCoordinate), vec2(1.0, uSize));

    // ignore 1D indices outside of the grid
    if(vI >= texture2D(tActiveVoxelsTotal, vec2(0.5)).r) discard;

    float offset = uSize - 2.;
    float k = 1. / uSize;

    vec2 relativePosition = k * vec2(offset, 0.);
    vec4 partialSums = texture2D(tActiveVoxelsPyramid, relativePosition);
    float start = 0.;
    vec4 starts = vec4(0.);
    vec4 ends = vec4(0.);
    float diff = 2.;
    vec4 m = vec4(0.);
    vec2 position = vec2(0.);
    vec4 vI4 = vec4(vI);

    // traverse the different levels of the pyramid
    for(int i = 1; i < 12; i++) {
        if(float(i) >= uLevels) break;

        offset -= diff;
        diff *= 2.;
        relativePosition = position + k * vec2(offset, 0.);

        ends = partialSums.wzyx + vec4(start);
        starts = vec4(start, ends.xyz);
        m = vec4(greaterThanEqual(vI4, starts)) * vec4(lessThan(vI4, ends));
        relativePosition += m.y * vec2(k, 0.) + m.z * vec2(0., k) + m.w * vec2(k, k);

        start = dot(m, starts);
        position = 2. * (relativePosition - k * vec2(offset, 0.));
        partialSums = texture2D(tActiveVoxelsPyramid, relativePosition);
    }

    ends = partialSums.wzyx + vec4(start);
    starts = vec4(start, ends.xyz);
    m = vec4(greaterThanEqual(vI4, starts)) * vec4(lessThan(vI4, ends));
    position += m.y * vec2(k, 0.) + m.z * vec2(0., k) + m.w * vec2(k, k);

    vec2 coord2d = position * uScale;
    vec3 coord3d = floor(index3dFrom2d(coord2d) + 0.5);

    float edgeIndex = texture2D(tActiveVoxelsBase, position * uScale).a;

    // current vertex for the up to 15 MC cases
    float currentVertex = vI - dot(m, starts);

    // get index into triIndices table
    float mcIndex = 16. * edgeIndex + currentVertex;
    vec4 mcData = texture2D(tTriIndices, vec2(intMod(mcIndex, 64.), floor(mcIndex / 64.)) / 64.);
    mcIndex = floor(mcData.a * 255.0 + 0.5);

    // bit mask to avoid conditionals (see comment below) for getting MC case corner
    vec4 m0 = vec4(mcIndex);

    // get edge value masks
    vec4 m1 = vec4(equal(m0, vec4(0., 1., 2., 3.)));
    vec4 m2 = vec4(equal(m0, vec4(4., 5., 6., 7.)));
    vec4 m3 = vec4(equal(m0, vec4(8., 9., 10., 11.)));

    // apply bit masks
    vec3 b0 = coord3d +
                m1.y * p0 +
                m1.z * p1 +
                m1.w * p2 +
                m2.x * p3 +
                m2.y * p4 +
                m2.z * p5 +
                m2.w * p6 +
                m3.y * p0 +
                m3.z * p1 +
                m3.w * p2;
    vec3 b1 = coord3d +
                m1.x * p0 +
                m1.y * p1 +
                m1.z * p2 +
                m2.x * p4 +
                m2.y * p5 +
                m2.z * p6 +
                m2.w * p3 +
                m3.x * p3 +
                m3.y * p4 +
                m3.z * p5 +
                m3.w * p6;

    // vec3 b0 = coord3d;
    // vec3 b1 = coord3d;
    // if (mcIndex == 0.0) {
    //     b0 += c0; b1 += c1;
    // } else if (mcIndex == 1.0) {
    //     b0 += c1; b1 += c2;
    // } else if (mcIndex == 2.0) {
    //     b0 += c2; b1 += c3;
    // } else if (mcIndex == 3.0) {
    //     b0 += c3; b1 += c0;
    // } else if (mcIndex == 4.0) {
    //     b0 += c4; b1 += c5;
    // } else if (mcIndex == 5.0) {
    //     b0 += c5; b1 += c6;
    // } else if (mcIndex == 6.0) {
    //     b0 += c6; b1 += c7;
    // } else if (mcIndex == 7.0) {
    //     b0 += c7; b1 += c4;
    // } else if (mcIndex == 8.0) {
    //     b0 += c0; b1 += c4;
    // } else if (mcIndex == 9.0) {
    //     b0 += c1; b1 += c5;
    // } else if (mcIndex == 10.0) {
    //     b0 += c2; b1 += c6;
    // } else if (mcIndex == 11.0) {
    //     b0 += c3; b1 += c7;
    // }
    // b0 = floor(b0 + 0.5);
    // b1 = floor(b1 + 0.5);

    float v0 = voxel(b0 / uGridDim).a;
    float v1 = voxel(b1 / uGridDim).a;

    float t = (uIsoValue - v0) / (v0 - v1);
    gl_FragColor.xyz = b0 + t * (b0 - b1);
}