export const clip_instance = `
#if defined(dClipVariant_instance) && dClipObjectCount != 0
    vec3 mCenter = (uModel * aTransform * vec4(uInvariantBoundingSphere.xyz, 1.0)).xyz;
    if (clipTest(mCenter)) {
        // move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
    }
#endif
`;