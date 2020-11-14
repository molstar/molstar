export default `
#if defined(dRenderVariant_color)
    #if !defined(dRenderMode_volume) && !defined(dRenderMode_isosurface)
        uniform sampler2D tDepth;
        uniform vec2 uDrawingBufferSize;

        float getDepth(const in vec2 coords) {
            #ifdef dPackedDepth
                return unpackRGBAToDepth(texture2D(tDepth, coords));
            #else
                return texture2D(tDepth, coords).r;
            #endif
        }
    #endif
    uniform int uRenderWboit;
#endif

float calcDepth(const in vec3 pos) {
    vec2 clipZW = pos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}
`;