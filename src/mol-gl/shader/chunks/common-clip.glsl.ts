/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export default `
#if dClipObjectCount != 0
    vec3 quaternionTransform(vec4 q, vec3 v) {
        vec3 t = 2.0 * cross(q.xyz, v);
        return v + q.w * t + cross(q.xyz, t);
    }

    vec4 computePlane(vec3 normal, vec3 inPoint) {
        return vec4(normalize(normal), -dot(normal, inPoint));
    }

    float planeSD(vec4 plane, vec3 center) {
        return -dot(plane.xyz, center - plane.xyz * -plane.w);
    }

    float sphereSD(vec3 position, vec4 rotation, vec3 size, vec3 center) {
        return (
            length(quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position) / size) - 1.0
        ) * min(min(size.x, size.y), size.z);
    }

    float cubeSD(vec3 position, vec4 rotation, vec3 size, vec3 center) {
        vec3 d = abs(quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position)) - size;
        return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
    }

    float cylinderSD(vec3 position, vec4 rotation, vec3 size, vec3 center) {
        vec3 t = quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position);

        vec2 d = abs(vec2(length(t.xz), t.y)) - size.xy;
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    }

    float infiniteConeSD(vec3 position, vec4 rotation, vec3 size, vec3 center) {
        vec3 t = quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position);

        float q = length(t.xy);
        return dot(size.xy, vec2(q, t.z));
    }

    float getSignedDistance(vec3 center, int type, vec3 position, vec4 rotation, vec3 scale) {
        if (type == 1) {
            vec3 normal = quaternionTransform(rotation, vec3(0.0, 1.0, 0.0));
            vec4 plane = computePlane(normal, position);
            return planeSD(plane, center);
        } else if (type == 2) {
            return sphereSD(position, rotation, scale * 0.5, center);
        } else if (type == 3) {
            return cubeSD(position, rotation, scale * 0.5, center);
        } else if (type == 4) {
            return cylinderSD(position, rotation, scale * 0.5, center);
        } else if (type == 5) {
            return infiniteConeSD(position, rotation, scale * 0.5, center);
        } else {
            return 0.1;
        }
    }

    #if __VERSION__ != 300
        // 8-bit
        int bitwiseAnd(int a, int b) {
            int d = 128;
            int result = 0;
            for (int i = 0; i < 8; ++i) {
                if (d <= 0) break;
                if (a >= d && b >= d) result += d;
                if (a >= d) a -= d;
                if (b >= d) b -= d;
                d /= 2;
            }
            return result;
        }

        bool hasBit(int mask, int bit) {
            return bitwiseAnd(mask, bit) == 0;
        }
    #else
        bool hasBit(int mask, int bit) {
            return (mask & bit) == 0;
        }
    #endif

    // flag is a bit-flag for clip-objects to ignore (note, object ids start at 1 not 0)
    bool clipTest(vec4 sphere, int flag) {
        for (int i = 0; i < dClipObjectCount; ++i) {
            if (flag == 0 || hasBit(flag, i + 1)) {
                // TODO take sphere radius into account?
                if (getSignedDistance(sphere.xyz, uClipObjectType[i], uClipObjectPosition[i], uClipObjectRotation[i], uClipObjectScale[i]) <= 0.0)
                    return true;
            }
        }
        return false;
    }
#endif
`;