export const wboit_params = `
#if defined(dRenderVariant_colorWboit)
    #if !defined(dRenderMode_volume) && !defined(dRenderMode_isosurface)
        uniform sampler2D tDepth;
        uniform vec2 uDrawingBufferSize;

        float getDepth(const in vec2 coords) {
            // always packed due to merged depth from primitives and volumes
            return unpackRGBAToDepth(texture2D(tDepth, coords));
        }
    #endif
#endif

uniform bool uRenderWboit;

float calcDepth(const in vec3 pos) {
    vec2 clipZW = pos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}
`;