export default `
precision highp float;
precision highp sampler2D;

uniform sampler2D tActiveVoxelsPyramid;
uniform sampler2D tActiveVoxelsBase;
uniform sampler2D tVolumeData;
uniform sampler2D tTriIndices;

uniform float uIsoValue;
uniform float uLevels;
uniform float uSize;
uniform float uCount;

uniform vec3 uGridDim;
uniform vec3 uGridTexDim;
uniform mat4 uGridTransform;

// scale to volume data coord
uniform vec2 uScale;

// varying vec2 vCoordinate;

#include common

// cube corners
const vec3 c0 = vec3(0., 0., 0.);
const vec3 c1 = vec3(1., 0., 0.);
const vec3 c2 = vec3(1., 1., 0.);
const vec3 c3 = vec3(0., 1., 0.);
const vec3 c4 = vec3(0., 0., 1.);
const vec3 c5 = vec3(1., 0., 1.);
const vec3 c6 = vec3(1., 1., 1.);
const vec3 c7 = vec3(0., 1., 1.);

const float EPS = 0.00001;

vec3 index3dFrom2d(vec2 coord) {
    vec2 gridTexPos = coord * uGridTexDim.xy;
    vec2 columnRow = floor(gridTexPos / uGridDim.xy);
    vec2 posXY = gridTexPos - columnRow * uGridDim.xy;
    float posZ = columnRow.y * floor(uGridTexDim.x / uGridDim.x) + columnRow.x;
    vec3 posXYZ = vec3(posXY, posZ);
    return posXYZ;
}

vec4 texture3dFrom2dNearest(sampler2D tex, vec3 pos, vec3 gridDim, vec2 texDim) {
    float zSlice = floor(pos.z * gridDim.z + 0.5); // round to nearest z-slice
    float column = intMod(zSlice * gridDim.x, texDim.x) / gridDim.x;
    float row = floor(intDiv(zSlice * gridDim.x, texDim.x));
    vec2 coord = (vec2(column * gridDim.x, row * gridDim.y) + (pos.xy * gridDim.xy)) / (texDim / uScale);
    return texture2D(tex, coord + 0.5 / (texDim / uScale));
    // return texture2D(tex, coord);
}

vec4 voxel(vec3 pos) {
    return texture3dFrom2dNearest(tVolumeData, pos, uGridDim, uGridTexDim.xy);
}

void main(void) {
    // get 1D index
    float vI = dot(floor(uSize * (gl_FragCoord.xy / uSize)), vec2(1.0, uSize));

    // ignore 1D indices outside of the grid
    if(vI >= uCount) discard;

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

    vec2 coord2d = position / uScale;
    vec3 coord3d = floor(index3dFrom2d(coord2d) + 0.5);

    float edgeIndex = floor(texture2D(tActiveVoxelsBase, position).a + 0.5);

    // current vertex for the up to 15 MC cases
    float currentVertex = vI - dot(m, starts);

    // get index into triIndices table
    float mcIndex = 16. * edgeIndex + currentVertex;
    vec4 mcData = texture2D(tTriIndices, vec2(intMod(mcIndex, 64.), floor(mcIndex / 64.)) / 64.);

    // bit mask to avoid conditionals (see comment below) for getting MC case corner
    vec4 m0 = vec4(floor(mcData.a * 255.0 + 0.5));

    // get edge value masks
    vec4 m1 = vec4(equal(m0, vec4(0., 1., 2., 3.)));
    vec4 m2 = vec4(equal(m0, vec4(4., 5., 6., 7.)));
    vec4 m3 = vec4(equal(m0, vec4(8., 9., 10., 11.)));

    // apply bit masks
    vec3 b0 = coord3d +
                m1.y * c1 +
                m1.z * c2 +
                m1.w * c3 +
                m2.x * c4 +
                m2.y * c5 +
                m2.z * c6 +
                m2.w * c7 +
                m3.y * c1 +
                m3.z * c2 +
                m3.w * c3;
    vec3 b1 = coord3d +
                m1.x * c1 +
                m1.y * c2 +
                m1.z * c3 +
                m2.x * c5 +
                m2.y * c6 +
                m2.z * c7 +
                m2.w * c4 +
                m3.x * c4 +
                m3.y * c5 +
                m3.z * c6 +
                m3.w * c7;

    // the conditionals that are avoided by above bitmasks
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

    vec4 d0 = voxel(b0 / uGridDim);
    vec4 d1 = voxel(b1 / uGridDim);

    float v0 = d0.a;
    float v1 = d1.a;

    float t = (uIsoValue - v0) / (v0 - v1);
    // t = -0.5;
    gl_FragData[0].xyz = (uGridTransform * vec4(b0 + t * (b0 - b1), 1.0)).xyz;
    gl_FragData[0].w = decodeFloatRGB(d0.rgb); // group id

    // normals from gradients
    vec3 n0 = -normalize(vec3(
        voxel((b0 - c1) / uGridDim).a - voxel((b0 + c1) / uGridDim).a,
        voxel((b0 - c3) / uGridDim).a - voxel((b0 + c3) / uGridDim).a,
        voxel((b0 - c4) / uGridDim).a - voxel((b0 + c4) / uGridDim).a
    ));
    vec3 n1 = -normalize(vec3(
        voxel((b1 - c1) / uGridDim).a - voxel((b1 + c1) / uGridDim).a,
        voxel((b1 - c3) / uGridDim).a - voxel((b1 + c3) / uGridDim).a,
        voxel((b1 - c4) / uGridDim).a - voxel((b1 + c4) / uGridDim).a
    ));
    gl_FragData[1].xyz = -vec3(
        n0.x + t * (n0.x - n1.x),
        n0.y + t * (n0.y - n1.y),
        n0.z + t * (n0.z - n1.z)
    );

    mat3 normalMatrix = transpose3(inverse3(mat3(uGridTransform)));
    gl_FragData[1].xyz = normalMatrix * gl_FragData[1].xyz;
}
`;