/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const wboit_write = `
#if defined(dRenderVariant_colorWboit)
    if (!uRenderWboit) {
        if (preFogAlpha < 1.0) {
            discard;
        }
    } else if (uRenderWboit) {
        // the 'fragmentDepth > 0.99' check is to handle precision issues with packed depth
        if (preFogAlpha != 1.0 && !interior && (fragmentDepth < getDepth(gl_FragCoord.xy / uDrawingBufferSize) || fragmentDepth > 0.99)) {
            float alpha = gl_FragColor.a;
            float wboitWeight = alpha * clamp(pow(1.0 - fragmentDepth, 2.0), 0.01, 1.0);
            gl_FragColor = vec4(gl_FragColor.rgb * alpha * wboitWeight, alpha);
            // extra alpha is to handle pre-multiplied alpha
            #if !defined(dRenderMode_volume) && !defined(dRenderMode_isosurface)
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