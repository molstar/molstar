/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext } from './compat';

export type WebGLState = {
    currentProgramId: number
    currentMaterialId: number
    currentRenderItemId: number

    /**
     * specifies which WebGL capability to enable
     * - `gl.BLEND`: blending of the computed fragment color values
     * - `gl.CULL_FACE`: culling of polygons
     * - `gl.DEPTH_TEST`: depth comparisons and updates to the depth buffer
     * - `gl.DITHER`: dithering of color components before they get written to the color buffer
     * - `gl.POLYGON_OFFSET_FILL`: adding an offset to depth values of polygon's fragments
     * - `gl.SAMPLE_ALPHA_TO_COVERAGE`: computation of a temporary coverage value determined by the alpha value
     * - `gl.SAMPLE_COVERAGE`: ANDing the fragment's coverage with the temporary coverage value
     * - `gl.SCISSOR_TEST`: scissor test that discards fragments that are outside of the scissor rectangle
     * - `gl.STENCIL_TEST`: stencil testing and updates to the stencil buffer
     */
    enable: (cap: number) => void
    /**
     * specifies which WebGL capability to disable
     * - `gl.BLEND`: blending of the computed fragment color values
     * - `gl.CULL_FACE`: culling of polygons
     * - `gl.DEPTH_TEST`: depth comparisons and updates to the depth buffer
     * - `gl.DITHER`: dithering of color components before they get written to the color buffer
     * - `gl.POLYGON_OFFSET_FILL`: adding an offset to depth values of polygon's fragments
     * - `gl.SAMPLE_ALPHA_TO_COVERAGE`: computation of a temporary coverage value determined by the alpha value
     * - `gl.SAMPLE_COVERAGE`: ANDing the fragment's coverage with the temporary coverage value
     * - `gl.SCISSOR_TEST`: scissor test that discards fragments that are outside of the scissor rectangle
     * - `gl.STENCIL_TEST`: stencil testing and updates to the stencil buffer
     */
    disable: (cap: number) => void

    /** specifies whether polygons are front- or back-facing by setting a winding orientation */
    frontFace: (mode: number) => void
    /** specifies whether or not front- and/or back-facing polygons can be culled */
    cullFace: (mode: number) => void
    /** sets whether writing into the depth buffer is enabled or disabled */
    depthMask: (flag: boolean) => void
    /** sets which color components to enable or to disable */
    colorMask: (red: boolean, green: boolean, blue: boolean, alpha: boolean) => void
    /** specifies the color values used when clearing color buffers, used when calling `gl.clear`, clamped to [0, 1] */
    clearColor: (red: number, green: number, blue: number, alpha: number) => void

    /** defines which function is used for blending pixel arithmetic */
    blendFunc: (src: number, dst: number) => void
    /** defines which function is used for blending pixel arithmetic for RGB and alpha components separately */
    blendFuncSeparate: (srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number) => void

    /** set both the RGB blend equation and alpha blend equation to a single equation, determines how a new pixel is combined with an existing */
    blendEquation: (mode: number) => void
    /** set the RGB blend equation and alpha blend equation separately, determines how a new pixel is combined with an existing */
    blendEquationSeparate: (modeRGB: number, modeAlpha: number) => void

    reset: () => void
}

export function createState(gl: GLRenderingContext): WebGLState {
    let enabledCapabilities: { [k: number]: boolean } = {};

    let currentFrontFace = gl.getParameter(gl.FRONT_FACE);
    let currentCullFace = gl.getParameter(gl.CULL_FACE_MODE);
    let currentDepthMask = gl.getParameter(gl.DEPTH_WRITEMASK);
    let currentColorMask = gl.getParameter(gl.COLOR_WRITEMASK);
    let currentClearColor = gl.getParameter(gl.COLOR_CLEAR_VALUE);

    let currentBlendSrcRGB = gl.getParameter(gl.BLEND_SRC_RGB);
    let currentBlendDstRGB = gl.getParameter(gl.BLEND_DST_RGB);
    let currentBlendSrcAlpha = gl.getParameter(gl.BLEND_SRC_ALPHA);
    let currentBlendDstAlpha = gl.getParameter(gl.BLEND_DST_ALPHA);

    let currentBlendEqRGB = gl.getParameter(gl.BLEND_EQUATION_RGB);
    let currentBlendEqAlpha = gl.getParameter(gl.BLEND_EQUATION_ALPHA);

    return {
        currentProgramId: -1,
        currentMaterialId: -1,
        currentRenderItemId: -1,

        enable: (cap: number) => {
            if (enabledCapabilities[cap] !== true ) {
                gl.enable(cap);
                enabledCapabilities[cap] = true;
            }
        },
        disable: (cap: number) => {
            if (enabledCapabilities[cap] !== false) {
                gl.disable(cap);
                enabledCapabilities[cap] = false;
            }
        },

        frontFace: (mode: number) => {
            if (mode !== currentFrontFace) {
                gl.frontFace(mode);
                currentFrontFace = mode;
            }
        },
        cullFace: (mode: number) => {
            if (mode !== currentCullFace) {
                gl.cullFace(mode);
                currentCullFace = mode;
            }
        },
        depthMask: (flag: boolean) => {
            if (flag !== currentDepthMask) {
                gl.depthMask(flag);
                currentDepthMask = flag;
            }
        },
        colorMask: (red: boolean, green: boolean, blue: boolean, alpha: boolean) => {
            if (red !== currentColorMask[0] || green !== currentColorMask[1] || blue !== currentColorMask[2] || alpha !== currentColorMask[3]) {
                gl.colorMask(red, green, blue, alpha);
                currentColorMask[0] = red;
                currentColorMask[1] = green;
                currentColorMask[2] = blue;
                currentColorMask[3] = alpha;
            }
        },
        clearColor: (red: number, green: number, blue: number, alpha: number) => {
            if (red !== currentClearColor[0] || green !== currentClearColor[1] || blue !== currentClearColor[2] || alpha !== currentClearColor[3]) {
                gl.clearColor(red, green, blue, alpha);
                currentClearColor[0] = red;
                currentClearColor[1] = green;
                currentClearColor[2] = blue;
                currentClearColor[3] = alpha;
            }
        },

        blendFunc: (src: number, dst: number) => {
            if (src !== currentBlendSrcRGB || dst !== currentBlendDstRGB || src !== currentBlendSrcAlpha || dst !== currentBlendDstAlpha) {
                gl.blendFunc(src, dst);
                currentBlendSrcRGB = src;
                currentBlendDstRGB = dst;
                currentBlendSrcAlpha = src;
                currentBlendDstAlpha = dst;
            }
        },
        blendFuncSeparate: (srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number) => {
            if (srcRGB !== currentBlendSrcRGB || dstRGB !== currentBlendDstRGB || srcAlpha !== currentBlendSrcAlpha || dstAlpha !== currentBlendDstAlpha) {
                gl.blendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha);
                currentBlendSrcRGB = srcRGB;
                currentBlendDstRGB = dstRGB;
                currentBlendSrcAlpha = srcAlpha;
                currentBlendDstAlpha = dstAlpha;
            }
        },

        blendEquation: (mode: number) => {
            if (mode !== currentBlendEqRGB || mode !== currentBlendEqAlpha) {
                gl.blendEquation(mode);
                currentBlendEqRGB = mode;
                currentBlendEqAlpha = mode;
            }
        },
        blendEquationSeparate: (modeRGB: number, modeAlpha: number) => {
            if (modeRGB !== currentBlendEqRGB || modeAlpha !== currentBlendEqAlpha) {
                gl.blendEquationSeparate(modeRGB, modeAlpha);
                currentBlendEqRGB = modeRGB;
                currentBlendEqAlpha = modeAlpha;
            }
        },

        reset: () => {
            enabledCapabilities = {};

            currentFrontFace = gl.getParameter(gl.FRONT_FACE);
            currentCullFace = gl.getParameter(gl.CULL_FACE_MODE);
            currentDepthMask = gl.getParameter(gl.DEPTH_WRITEMASK);
            currentColorMask = gl.getParameter(gl.COLOR_WRITEMASK);
            currentClearColor = gl.getParameter(gl.COLOR_CLEAR_VALUE);

            currentBlendSrcRGB = gl.getParameter(gl.BLEND_SRC_RGB);
            currentBlendDstRGB = gl.getParameter(gl.BLEND_DST_RGB);
            currentBlendSrcAlpha = gl.getParameter(gl.BLEND_SRC_ALPHA);
            currentBlendDstAlpha = gl.getParameter(gl.BLEND_DST_ALPHA);

            currentBlendEqRGB = gl.getParameter(gl.BLEND_EQUATION_RGB);
            currentBlendEqAlpha = gl.getParameter(gl.BLEND_EQUATION_ALPHA);
        }
    };
}