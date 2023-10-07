/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export const common_clip = `
#if dClipObjectCount != 0
    vec3 quaternionTransform(const in vec4 q, const in vec3 v) {
        vec3 t = 2.0 * cross(q.xyz, v);
        return v + q.w * t + cross(q.xyz, t);
    }

    vec4 computePlane(const in vec3 normal, const in vec3 inPoint) {
        return vec4(normalize(normal), -dot(normal, inPoint));
    }

    float planeSD(const in vec4 plane, const in vec3 center) {
        return -dot(plane.xyz, center - plane.xyz * -plane.w);
    }

    float sphereSD(const in vec3 position, const in vec4 rotation, const in vec3 size, const in vec3 center) {
        return (
            length(quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position) / size) - 1.0
        ) * min(min(size.x, size.y), size.z);
    }

    float cubeSD(const in vec3 position, const in vec4 rotation, const in vec3 size, const in vec3 center) {
        vec3 d = abs(quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position)) - size;
        return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
    }

    float cylinderSD(const in vec3 position, const in vec4 rotation, const in vec3 size, const in vec3 center) {
        vec3 t = quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position);

        vec2 d = abs(vec2(length(t.xz), t.y)) - size.xy;
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    }

    float infiniteConeSD(const in vec3 position, const in vec4 rotation, const in vec3 size, const in vec3 center) {
        vec3 t = quaternionTransform(vec4(-rotation.x, -rotation.y, -rotation.z, rotation.w), center - position);

        float q = length(t.xy);
        return dot(size.xy, vec2(q, t.z));
    }

    float getSignedDistance(const in vec3 center, const in int type, const in vec3 position, const in vec4 rotation, const in vec3 scale) {
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

    #if __VERSION__ == 100
        // 8-bit
        int bitwiseAnd(in int a, in int b) {
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

        bool hasBit(const in int mask, const in int bit) {
            return bitwiseAnd(mask, bit) == 0;
        }
    #else
        bool hasBit(const in int mask, const in int bit) {
            return (mask & bit) == 0;
        }
    #endif

    bool clipTest(const in vec4 sphere) {
        // flag is a bit-flag for clip-objects to ignore (note, object ids start at 1 not 0)
        #if defined(dClipping)
            int flag = int(floor(vClipping * 255.0 + 0.5));
        #else
            int flag = 0;
        #endif

        #pragma unroll_loop_start
        for (int i = 0; i < dClipObjectCount; ++i) {
            if (flag == 0 || hasBit(flag, UNROLLED_LOOP_INDEX + 1)) {
                // TODO take sphere radius into account?
                bool test = getSignedDistance(sphere.xyz, uClipObjectType[i], uClipObjectPosition[i], uClipObjectRotation[i], uClipObjectScale[i]) <= 0.0;
                if ((!uClipObjectInvert[i] && test) || (uClipObjectInvert[i] && !test)) {
                    return true;
                }
            }
        }
        #pragma unroll_loop_end
        return false;
    }
#endif
`;