export default `
#if defined(dClipVariant_instance) && dClipObjectCount != 0
    int flag = 0;
    #if defined(dClipping)
        flag = int(floor(vClipping * 255.0 + 0.5));
    #endif

    vec4 mCenter = uModel * aTransform * vec4(uInvariantBoundingSphere.xyz, 1.0);
    if (clipTest(vec4(mCenter.xyz, uInvariantBoundingSphere.w), flag))
        // move out of [ -w, +w ] to 'discard' in vert shader
        gl_Position.z = 2.0 * gl_Position.w;
#endif
`;