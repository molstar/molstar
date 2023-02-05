export const clip_instance = `
#if defined(dClipVariant_instance) && dClipObjectCount != 0
    vec4 mCenter = uModel * aTransform * vec4(uInvariantBoundingSphere.xyz, 1.0);
    if (clipTest(vec4(mCenter.xyz, uInvariantBoundingSphere.w)))
        // move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
#endif
`;