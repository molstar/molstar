export default `
#if defined(dFlatShaded) && defined(enabledStandardDerivatives)
    vec3 fdx = dFdx(vViewPosition);
    vec3 fdy = dFdy(vViewPosition);
    vec3 normal = -normalize(cross(fdx, fdy));
#else
    vec3 normal = -normalize(vNormal);
    #ifdef dDoubleSided
        normal = normal * (float(frontFacing) * 2.0 - 1.0);
    #endif
#endif
`