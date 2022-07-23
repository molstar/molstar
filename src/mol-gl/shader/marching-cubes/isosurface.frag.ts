export const isosurface_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#if __VERSION__ == 100
    uniform sampler2D tActiveVoxelsPyramid;
#else
    precision highp isampler2D;
    uniform isampler2D tActiveVoxelsPyramid;
#endif

uniform sampler2D tActiveVoxelsBase;
uniform sampler2D tVolumeData;
uniform sampler2D tTriIndices;

uniform float uIsoValue;
uniform float uLevels;
uniform float uSize;
uniform float uCount;
uniform bool uInvert;

uniform vec3 uGridDim;
uniform vec3 uGridTexDim;
uniform mat4 uGridTransform;

// scale to volume data coord
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
    return texture2D(tex, coord + 0.5 / (texDim / uScale));
}

vec4 voxel(vec3 pos) {
    pos = min(max(vec3(0.0), pos), uGridDim - vec3(1.0));
    return texture3dFrom2dNearest(tVolumeData, pos / uGridDim, uGridDim, uGridTexDim.xy);
}

vec4 voxelPadded(vec3 pos) {
    pos = min(max(vec3(0.0), pos), uGridDim - vec3(vec2(2.0), 1.0)); // remove xy padding
    return texture3dFrom2dNearest(tVolumeData, pos / uGridDim, uGridDim, uGridTexDim.xy);
}

int idot2(const in ivec2 a, const in ivec2 b) {
    return a.x * b.x + a.y * b.y;
}

int idot4(const in ivec4 a, const in ivec4 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

#if __VERSION__ == 100
    int pyramidVoxel(vec2 pos) {
        return int(unpackRGBToInt(texture2D(tActiveVoxelsPyramid, pos / (vec2(1.0, 0.5) * uSize)).rgb));
    }
#else
    int pyramidVoxel(vec2 pos) {
        return texture2D(tActiveVoxelsPyramid, pos / (vec2(1.0, 0.5) * uSize)).r;
    }
#endif

vec4 baseVoxel(vec2 pos) {
    return texture2D(tActiveVoxelsBase, pos / uSize);
}

vec4 getGroup(const in vec3 p) {
    vec3 gridDim = uGridDim - vec3(1.0, 1.0, 0.0); // remove xy padding
    // note that we swap x and z because the texture is flipped around y
    #if defined(dAxisOrder_012)
        float group = p.z + p.y * gridDim.z + p.x * gridDim.z * gridDim.y; // 210
    #elif defined(dAxisOrder_021)
        float group = p.y + p.z * gridDim.y + p.x * gridDim.y * gridDim.z; // 120
    #elif defined(dAxisOrder_102)
        float group = p.z + p.x * gridDim.z + p.y * gridDim.z * gridDim.x; // 201
    #elif defined(dAxisOrder_120)
        float group = p.x + p.z * gridDim.x + p.y * gridDim.x * gridDim.z; // 021
    #elif defined(dAxisOrder_201)
        float group = p.y + p.x * gridDim.y + p.z * gridDim.y * gridDim.x; // 102
    #elif defined(dAxisOrder_210)
        float group = p.x + p.y * gridDim.x + p.z * gridDim.x * gridDim.y; // 012
    #endif
    return vec4(group > 16777215.5 ? vec3(1.0) : packIntToRGB(group), 1.0);
}

void main(void) {
    // get 1D index
    int vI = int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(uSize);

    // ignore 1D indices outside of the grid
    if(vI >= int(uCount)) discard;

    ivec2 offset = ivec2(int(uSize) - 2, 0);

    int start = 0;
    ivec4 starts = ivec4(0);
    ivec4 ends = ivec4(0);
    int diff = 2;
    ivec4 m = ivec4(0);
    ivec2 position = ivec2(0);
    ivec4 vI4 = ivec4(vI);

    ivec2 relativePosition = ivec2(0);
    int end = 0;
    ivec2 pos1 = ivec2(0);
    ivec2 pos2 = ivec2(0);
    ivec2 pos3 = ivec2(0);
    ivec2 pos4 = ivec2(0);
    ivec3 vI3 = ivec3(vI);
    ivec3 mask = ivec3(0);

    // traverse the different levels of the pyramid
    for(int i = 1; i < 14; i++) {
        if(float(i) >= uLevels) break;

        offset.x -= diff;
        diff *= 2;
        relativePosition = position + offset;

        end = start + pyramidVoxel(vec2(relativePosition));
        pos1 = ivec2(relativePosition);
        starts.x = start;
        ends.x = end;
        pos2 = ivec2(relativePosition + ivec2(1, 0));
        starts.y = ends.x;
        ends.y = ends.x + pyramidVoxel(vec2(pos2));
        pos3 = relativePosition + ivec2(0, 1);
        starts.z = ends.y;
        ends.z = ends.y + pyramidVoxel(vec2(pos3));
        pos4 = relativePosition + ivec2(1, 1);
        starts.w = ends.z;
        mask = ivec3(greaterThanEqual(vI3, starts.rgb)) * ivec3(lessThan(vI3, ends.rgb));
        m = ivec4(mask, 1 - int(any(bvec3(mask))));

        relativePosition = m.x * pos1 + m.y * pos2 + m.z * pos3 + m.w * pos4;
        start = idot4(m, starts);
        position = 2 * (relativePosition - offset);
    }

    end = start + int(baseVoxel(vec2(position)).r * 255.0);
    pos1 = position;
    starts.x = start;
    ends.x = end;
    pos2 = position + ivec2(1, 0);
    starts.y = ends.x;
    ends.y = ends.x + int(baseVoxel(vec2(pos2)).r * 255.0);
    pos3 = position + ivec2(0, 1);
    starts.z = ends.y;
    ends.z = ends.y + int(baseVoxel(vec2(pos3)).r * 255.0);
    pos4 = position + ivec2(1, 1);
    starts.w = ends.z;
    mask = ivec3(greaterThanEqual(vI3, starts.rgb)) * ivec3(lessThan(vI3, ends.rgb));
    m = ivec4(mask, 1 - int(any(bvec3(mask))));
    position = m.x * pos1 + m.y * pos2 + m.z * pos3 + m.w * pos4;

    vec2 coord2d = (vec2(position) / uSize) / uScale;
    vec3 coord3d = floor(index3dFrom2d(coord2d) + 0.5);

    float edgeIndex = floor(baseVoxel(vec2(position)).a * 255.0 + 0.5);

    // current vertex for the up to 15 MC cases
    int currentVertex = vI - idot4(m, starts);

    // ensure winding-order is the same for negative and positive iso-levels
    if (uInvert) {
        int v = imod(currentVertex + 1, 3);
        if (v == 1) currentVertex += 2;
        else if (v == 0) currentVertex -= 2;
    }

    // get index into triIndices table
    int mcIndex = 16 * int(edgeIndex) + currentVertex;
    vec4 mcData = texture2D(tTriIndices, vec2(imod(mcIndex, 64), mcIndex / 64) / 64.);

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
    //     b1 += c1;
    // } else if (mcIndex == 1.0) {
    //     b0 += c1; b1 += c2;
    // } else if (mcIndex == 2.0) {
    //     b0 += c2; b1 += c3;
    // } else if (mcIndex == 3.0) {
    //     b0 += c3;
    // } else if (mcIndex == 4.0) {
    //     b0 += c4; b1 += c5;
    // } else if (mcIndex == 5.0) {
    //     b0 += c5; b1 += c6;
    // } else if (mcIndex == 6.0) {
    //     b0 += c6; b1 += c7;
    // } else if (mcIndex == 7.0) {
    //     b0 += c7; b1 += c4;
    // } else if (mcIndex == 8.0) {
    //     b1 += c4;
    // } else if (mcIndex == 9.0) {
    //     b0 += c1; b1 += c5;
    // } else if (mcIndex == 10.0) {
    //     b0 += c2; b1 += c6;
    // } else if (mcIndex == 11.0) {
    //     b0 += c3; b1 += c7;
    // }
    // b0 = floor(b0 + 0.5);
    // b1 = floor(b1 + 0.5);

    vec4 d0 = voxel(b0);
    vec4 d1 = voxel(b1);

    float v0 = d0.a;
    float v1 = d1.a;

    float t = (uIsoValue - v0) / (v0 - v1);
    gl_FragData[0].xyz = (uGridTransform * vec4(b0 + t * (b0 - b1), 1.0)).xyz;

    // group id
    #if __VERSION__ == 100 || defined(dConstantGroup)
        // webgl1 does not support 'flat' interpolation (i.e. no interpolation)
        // ensure a constant group id per triangle as needed
        #ifdef dPackedGroup
            gl_FragData[1] = vec4(voxel(coord3d).rgb, 1.0);
        #else
            gl_FragData[1] = getGroup(coord3d);
        #endif
    #else
        #ifdef dPackedGroup
            gl_FragData[1] = vec4(t < 0.5 ? d0.rgb : d1.rgb, 1.0);
        #else
            gl_FragData[1] = getGroup(t < 0.5 ? b0 : b1);
        #endif
    #endif

    // normals from gradients
    vec3 n0 = -normalize(vec3(
        voxelPadded(b0 - c1).a - voxelPadded(b0 + c1).a,
        voxelPadded(b0 - c3).a - voxelPadded(b0 + c3).a,
        voxelPadded(b0 - c4).a - voxelPadded(b0 + c4).a
    ));
    vec3 n1 = -normalize(vec3(
        voxelPadded(b1 - c1).a - voxelPadded(b1 + c1).a,
        voxelPadded(b1 - c3).a - voxelPadded(b1 + c3).a,
        voxelPadded(b1 - c4).a - voxelPadded(b1 + c4).a
    ));
    gl_FragData[2].xyz = -vec3(
        n0.x + t * (n0.x - n1.x),
        n0.y + t * (n0.y - n1.y),
        n0.z + t * (n0.z - n1.z)
    );

    // ensure normal-direction is the same for negative and positive iso-levels
    if (uInvert) {
        gl_FragData[2].xyz *= -1.0;
    }

    // apply normal matrix
    gl_FragData[2].xyz *= transpose3(inverse3(mat3(uGridTransform)));
}
`;