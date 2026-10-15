/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Background at the far plane: skybox (cube texture), image, or horizontal and radial gradients. Fragment coordinates
 * in GL window space (bottom-left origin) as on WebGL; images and cube faces sampled in upload order.
 */
export const background_wgsl = /* wgsl */ `
struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) vPosition: vec4<f32>,
};

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.vPosition = vec4<f32>(in.aPosition, 1.0, 1.0);
    out.position = out.vPosition;
    return out;
}

fn screenSpaceDither(screenPos: vec2<f32>) -> vec3<f32> {
    var dither = vec3<f32>(dot(vec2<f32>(171.0, 231.0), screenPos));
    dither = fract(dither / vec3<f32>(103.0, 71.0, 97.0));
    return dither / 255.0;
}

fn saturateColor(c: vec3<f32>, amount: f32) -> vec3<f32> {
    // https://www.w3.org/TR/WCAG21/#dfn-relative-luminance
    let W = vec3<f32>(0.2125, 0.7154, 0.0721);
    let intensity = vec3<f32>(dot(c, W));
    return mix(intensity, c, 1.0 + amount);
}

fn lightenColor(c: vec3<f32>, amount: f32) -> vec3<f32> {
    return c + amount;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSize = uObject.uTexSize;
    let viewport = uObject.uViewport;
    let viewportAdjusted = uObject.uViewportAdjusted != 0u;
    // gl_FragCoord.xy of WebGL
    let fragCoord = vec2<f32>(in.position.x, texSize.y - in.position.y);

    #if defined(dVariant_skybox)
        let t = uObject.uViewDirectionProjectionInverse * in.vPosition;
        let color = textureSampleLevel(tSkybox, sLinearMipmap, uObject.uRotation * normalize(t.xyz / t.w), uObject.uBlur * 8.0);
        return vec4<f32>(lightenColor(saturateColor(color.rgb, uObject.uSaturation), uObject.uLightness), uObject.uOpacity);
    #elif defined(dVariant_image)
        var coords: vec2<f32>;
        if (viewportAdjusted) {
            coords = ((fragCoord - viewport.xy) * (texSize / viewport.zw) / uObject.uImageScale) + uObject.uImageOffset;
        } else {
            coords = (fragCoord / uObject.uImageScale) + uObject.uImageOffset;
        }
        let color = textureSampleLevel(tImage, sLinearMipmap, vec2<f32>(coords.x, 1.0 - coords.y), uObject.uBlur * 8.0);
        return vec4<f32>(lightenColor(saturateColor(color.rgb, uObject.uSaturation), uObject.uLightness), uObject.uOpacity);
    #elif defined(dVariant_horizontalGradient)
        var d: f32;
        if (viewportAdjusted) {
            d = ((fragCoord.y - viewport.y) * (texSize.y / viewport.w) / texSize.y) + 1.0 - (uObject.uGradientRatio * 2.0);
        } else {
            d = (fragCoord.y / texSize.y) + 1.0 - (uObject.uGradientRatio * 2.0);
        }
        let color = mix(uObject.uGradientColorB, uObject.uGradientColorA, clamp(d, 0.0, 1.0));
        return vec4<f32>(color + screenSpaceDither(fragCoord), 1.0);
    #elif defined(dVariant_radialGradient)
        var d: f32;
        if (viewportAdjusted) {
            d = distance(vec2<f32>(0.5), (fragCoord - viewport.xy) * (texSize / viewport.zw) / texSize) + uObject.uGradientRatio - 0.5;
        } else {
            d = distance(vec2<f32>(0.5), fragCoord / texSize) + uObject.uGradientRatio - 0.5;
        }
        let color = mix(uObject.uGradientColorB, uObject.uGradientColorA, 1.0 - clamp(d, 0.0, 1.0));
        return vec4<f32>(color + screenSpaceDither(fragCoord), 1.0);
    #else
        return vec4<f32>(0.0);
    #endif
}
`;
