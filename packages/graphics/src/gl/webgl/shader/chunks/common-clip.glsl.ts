/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

export const common_clip = `
vec3 quaternionTransform(const in vec4 q, const in vec3 v) {
    vec3 t = 2.0 * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

// signed distance of \`t\` (in the shape's frame) to a shape of half size \`size\` (1 plane, 2 sphere, 3 cube, 4 cylinder,
// 5 infinite cone)
float localSD(const in int type, const in vec3 t, const in vec3 size) {
    if (type == 1) {
        return -t.y;
    } else if (type == 2) {
        return (length(t / size) - 1.0) * min(min(size.x, size.y), size.z);
    } else if (type == 3) {
        vec3 d = abs(t) - size;
        return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
    } else if (type == 4) {
        vec2 d = abs(vec2(length(t.xz), t.y)) - size.xy;
        return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
    } else if (type == 5) {
        return dot(size.xy, vec2(length(t.xy), t.z));
    }
    return 0.1;
}

// signed distance of \`center\` to a shape given by position, rotation, scale and the transform of the point (image trims)
float getSignedDistance(const in vec3 center, const in int type, const in vec3 position, const in vec4 rotation, const in vec3 scale, const in mat4 transform) {
    vec3 c = (transform * vec4(center, 1.0)).xyz;
    return localSD(type, quaternionTransform(vec4(-rotation.xyz, rotation.w), c - position), scale * 0.5);
}

// a clip object by its frame and shape (\`Clip.Objects\`): the type in \`shape.w\`, negative where inverted
int clipType(const in vec4 shape) {
    return int(abs(shape.w) + 0.5);
}

bool clipInverted(const in vec4 shape) {
    return shape.w < 0.0;
}

float clipObjectSD(const in vec3 center, const in mat4 frame, const in vec4 shape) {
    return localSD(clipType(shape), (frame * vec4(center, 1.0)).xyz, shape.xyz);
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

    // bit of the clip object \`index\` in the clipping groups flag
    int clipBit(const in int index) {
        int bit = 1;
        for (int i = 0; i < 16; ++i) {
            if (i == index) break;
            bit *= 2;
        }
        return bit;
    }
#else
    bool hasBit(const in int mask, const in int bit) {
        return (mask & bit) == 0;
    }

    // bit of the clip object \`index\` in the clipping groups flag
    int clipBit(const in int index) {
        return 1 << index;
    }
#endif

#if defined(dClipObjects)
    bool clipIgnored(const in int index) {
        #if defined(dClipping)
            int flag = int(floor(vClipping * 255.0 + 0.5));
            return flag != 0 && !hasBit(flag, clipBit(index));
        #else
            return false;
        #endif
    }

    bool clipTest(const in vec3 center) {
        // flag is a bit-flag for clip-objects to ignore
        #if defined(dClipping)
            int flag = int(floor(vClipping * 255.0 + 0.5));
        #else
            int flag = 0;
        #endif

        for (int i = 0; i < uClipObjectCount; ++i) {
            if (i != uSolidInteriorClip && (flag == 0 || hasBit(flag, clipBit(i)))) {
                bool test = clipObjectSD(center, uClipObjectFrame[i], uClipObjectShape[i]) <= 0.0;
                if (test != clipInverted(uClipObjectShape[i])) {
                    return true;
                }
            }
        }
        return false;
    }

    float clipDistance(const in vec3 center) {
        #if defined(dClipping)
            int flag = int(floor(vClipping * 255.0 + 0.5));
        #else
            int flag = 0;
        #endif

        float d = 1.0e10;
        for (int i = 0; i < uClipObjectCount; ++i) {
            if (i != uSolidInteriorClip && (flag == 0 || hasBit(flag, clipBit(i)))) {
                float sd = clipObjectSD(center, uClipObjectFrame[i], uClipObjectShape[i]);
                d = min(d, clipInverted(uClipObjectShape[i]) ? -sd : sd);
            }
        }
        return d;
    }

    vec3 clipNormal(const in vec3 center) {
        vec2 e = vec2(0.01, 0.0);
        return normalize(vec3(
            clipDistance(center + e.xyy) - clipDistance(center - e.xyy),
            clipDistance(center + e.yxy) - clipDistance(center - e.yxy),
            clipDistance(center + e.yyx) - clipDistance(center - e.yyx)
        ));
    }

    vec2 clipHalfLine(const in float f0, const in float f1) {
        float df = f1 - f0;
        if (abs(df) < 1.0e-12) return f0 <= 0.0 ? vec2(-1.0e6, 1.0e6) : vec2(1.0e6, -1.0e6);
        float sr = -f0 / df;
        return df < 0.0 ? vec2(sr, 1.0e6) : vec2(-1.0e6, sr);
    }

    vec2 clipSlab(const in float q0, const in float q1, const in float size) {
        vec2 a = clipHalfLine(q0 - size, q1 - size);
        vec2 b = clipHalfLine(-q0 - size, -q1 - size);
        return vec2(max(a.x, b.x), min(a.y, b.y));
    }

    vec2 clipQuadratic(const in float a, const in float b, const in float c) {
        if (a < 1.0e-12) return c <= 0.0 ? vec2(-1.0e6, 1.0e6) : vec2(1.0e6, -1.0e6);
        float disc = b * b - a * c;
        if (disc < 0.0) return vec2(1.0e6, -1.0e6);
        float sq = sqrt(disc);
        return vec2((-b - sq) / a, (-b + sq) / a);
    }

    vec2 clipObjectInterval(const in mat4 frame, const in vec4 shape, const in vec3 pFront, const in vec3 pBack) {
        vec3 qF = (frame * vec4(pFront, 1.0)).xyz;
        vec3 qB = (frame * vec4(pBack, 1.0)).xyz;
        int type = clipType(shape);
        if (type == 1) return clipHalfLine(-qF.y, -qB.y);
        vec3 size = shape.xyz;
        if (type == 2) {
            vec3 o = qF / size, d = (qB - qF) / size;
            return clipQuadratic(dot(d, d), dot(o, d), dot(o, o) - 1.0);
        } else if (type == 3) {
            vec2 x = clipSlab(qF.x, qB.x, size.x), y = clipSlab(qF.y, qB.y, size.y), z = clipSlab(qF.z, qB.z, size.z);
            return vec2(max(x.x, max(y.x, z.x)), min(x.y, min(y.y, z.y)));
        } else if (type == 4) {
            vec2 o = qF.xz, d = (qB - qF).xz;
            vec2 body = clipQuadratic(dot(d, d), dot(o, d), dot(o, o) - size.x * size.x);
            vec2 slab = clipSlab(qF.y, qB.y, size.y);
            return vec2(max(body.x, slab.x), min(body.y, slab.y));
        } else if (type == 5) {
            vec2 xyF = qF.xy, xyd = (qB - qF).xy;
            float hF = -size.y * qF.z, hd = -size.y * (qB.z - qF.z);
            vec2 nappe = clipHalfLine(-hF, -hF - hd);
            float a = size.x * size.x * dot(xyd, xyd) - hd * hd;
            float b = size.x * size.x * dot(xyF, xyd) - hF * hd;
            float c = size.x * size.x * dot(xyF, xyF) - hF * hF;
            if (a > 1.0e-12) {
                vec2 g = clipQuadratic(a, b, c);
                return vec2(max(g.x, nappe.x), min(g.y, nappe.y));
            }
            float disc = b * b - a * c;
            if (a < -1.0e-12 && disc >= 0.0) {
                float sq = sqrt(disc);
                float s1 = (-b + sq) / a, s2 = (-b - sq) / a;
                vec2 lo = vec2(max(-1.0e6, nappe.x), min(s1, nappe.y));
                vec2 hi = vec2(max(s2, nappe.x), min(1.0e6, nappe.y));
                return lo.x <= lo.y ? lo : hi;
            }
            vec2 lin = clipHalfLine(c, 2.0 * b + a + c);
            return vec2(max(lin.x, nappe.x), min(lin.y, nappe.y));
        }
        return vec2(1.0e6, -1.0e6);
    }

    float clipExit(const in vec3 pFront, const in vec3 pBack, const in float sStart) {
        #if defined(dClipping)
            int flag = int(floor(vClipping * 255.0 + 0.5));
        #else
            int flag = 0;
        #endif

        float s = sStart;
        for (int k = 0; k <= uClipObjectCount; ++k) {
            float sNext = s;
            for (int i = 0; i < uClipObjectCount; ++i) {
                if (i != uSolidInteriorClip && (flag == 0 || hasBit(flag, clipBit(i)))) {
                    vec2 iv = clipObjectInterval(uClipObjectFrame[i], uClipObjectShape[i], pFront, pBack);
                    if (clipInverted(uClipObjectShape[i])) {
                        if (iv.x > iv.y || s > iv.y) sNext = 1.0e6;
                        else if (s < iv.x) sNext = max(sNext, iv.x);
                    } else if (s > iv.x && s < iv.y) {
                        sNext = max(sNext, iv.y);
                    }
                }
            }
            if (sNext == s) return s;
            s = sNext;
            if (s >= 1.0) return -1.0;
        }
        return -1.0;
    }

    float clipCapExit(const in vec3 pFront, const in vec3 pBack, const in bool back) {
        float s = -1.0;
        for (int i = 0; i < uClipObjectCount; ++i) {
            if (i == uSolidInteriorClip) {
                vec2 iv = clipObjectInterval(uClipObjectFrame[i], uClipObjectShape[i], pFront, pBack);
                if (iv.x < iv.y) s = clipInverted(uClipObjectShape[i]) != back ? iv.x : iv.y;
            }
        }
        return s;
    }

    float clipCapDistance(const in vec3 center) {
        float d = 0.0;
        for (int i = 0; i < uClipObjectCount; ++i) {
            if (i == uSolidInteriorClip) {
                float sd = clipObjectSD(center, uClipObjectFrame[i], uClipObjectShape[i]);
                d = clipInverted(uClipObjectShape[i]) ? -sd : sd;
            }
        }
        return d;
    }

    vec3 clipCapNormal(const in vec3 center) {
        vec2 e = vec2(0.01, 0.0);
        return normalize(vec3(
            clipCapDistance(center + e.xyy) - clipCapDistance(center - e.xyy),
            clipCapDistance(center + e.yxy) - clipCapDistance(center - e.yxy),
            clipCapDistance(center + e.yyx) - clipCapDistance(center - e.yyx)
        ));
    }
#endif
`;
