/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const wboit_write = `
#if defined(dRenderVariant_colorWboit)
    if (uRenderMask == MaskOpaque) {
        if (preFogAlpha < 1.0) {
            discard;
        }
    } else if (uRenderMask == MaskTransparent) {
        // opaque fragments are drawn in the opaque pass; direct volumes have none: drawn here also where opaque
        #ifdef dGeometryType_directVolume
            bool wboitOpaque = false;
        #else
            bool wboitOpaque = preFogAlpha == 1.0;
        #endif
        if (!wboitOpaque && fragmentDepth > getDepth(gl_FragCoord.xy / uDrawingBufferSize)) {
            #ifdef dTransparentBackfaces_off
                if (interior) discard;
            #endif
            float alpha = gl_FragColor.a;
            // depth is reversed: larger in front
            float wboitWeight = alpha * clamp(pow(fragmentDepth, 2.0), 0.01, 1.0);
            gl_FragColor = vec4(gl_FragColor.rgb * alpha * wboitWeight, alpha);
            // extra alpha is to handle pre-multiplied alpha
            #ifndef dGeometryType_directVolume
                gl_FragData[1] = vec4((uTransparentBackground ? alpha : 1.0) * alpha * wboitWeight);
            #else
                gl_FragData[1] = vec4(alpha * alpha * wboitWeight);
            #endif
        } else {
            discard;
        }
    }
#endif
`;
