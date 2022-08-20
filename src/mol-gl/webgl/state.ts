/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    /** specifies the depth value used when clearing depth buffer, used when calling `gl.clear` */
    clearDepth: (depth: number) => void
    /** sets the depth comparison function */
    depthFunc: (func: number) => void
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
    /** specifies the source and destination blending factors, clamped to [0, 1] */
    blendColor: (red: number, green: number, blue: number, alpha: number) => void

    enableVertexAttrib: (index: number) => void
    clearVertexAttribsState: () => void
    disableUnusedVertexAttribs: () => void

    viewport: (x: number, y: number, width: number, height: number) => void
    scissor: (x: number, y: number, width: number, height: number) => void

    reset: () => void
}

export function createState(gl: GLRenderingContext): WebGLState {
    let enabledCapabilities: Record<number, boolean> = {};

    let currentFrontFace = gl.getParameter(gl.FRONT_FACE);
    let currentCullFace = gl.getParameter(gl.CULL_FACE_MODE);
    let currentDepthMask = gl.getParameter(gl.DEPTH_WRITEMASK);
    let currentClearDepth = gl.getParameter(gl.DEPTH_CLEAR_VALUE);
    let currentDepthFunc = gl.getParameter(gl.DEPTH_FUNC);
    let currentColorMask = gl.getParameter(gl.COLOR_WRITEMASK);
    let currentClearColor = gl.getParameter(gl.COLOR_CLEAR_VALUE);

    let currentBlendSrcRGB = gl.getParameter(gl.BLEND_SRC_RGB);
    let currentBlendDstRGB = gl.getParameter(gl.BLEND_DST_RGB);
    let currentBlendSrcAlpha = gl.getParameter(gl.BLEND_SRC_ALPHA);
    let currentBlendDstAlpha = gl.getParameter(gl.BLEND_DST_ALPHA);
    let currentBlendColor = gl.getParameter(gl.BLEND_COLOR);

    let currentBlendEqRGB = gl.getParameter(gl.BLEND_EQUATION_RGB);
    let currentBlendEqAlpha = gl.getParameter(gl.BLEND_EQUATION_ALPHA);

    let maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
    const vertexAttribsState: number[] = [];

    let currentViewport: [number, number, number, number] = gl.getParameter(gl.VIEWPORT);
    let currentScissor: [number, number, number, number] = gl.getParameter(gl.SCISSOR_BOX);

    const clearVertexAttribsState = () => {
        for (let i = 0; i < maxVertexAttribs; ++i) {
            vertexAttribsState[i] = 0;
        }
    };
    clearVertexAttribsState();

    return {
        currentProgramId: -1,
        currentMaterialId: -1,
        currentRenderItemId: -1,

        enable: (cap: number) => {
            if (enabledCapabilities[cap] !== true) {
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
        clearDepth: (depth: number) => {
            if (depth !== currentClearDepth) {
                gl.clearDepth(depth);
                currentClearDepth = depth;
            }
        },
        depthFunc: (func: number) => {
            if (func !== currentDepthFunc) {
                gl.depthFunc(func);
                currentDepthFunc = func;
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
        blendColor: (red: number, green: number, blue: number, alpha: number) => {
            if (red !== currentBlendColor[0] || green !== currentBlendColor[1] || blue !== currentBlendColor[2] || alpha !== currentBlendColor[3]) {
                gl.blendColor(red, green, blue, alpha);
                currentBlendColor[0] = red;
                currentBlendColor[1] = green;
                currentBlendColor[2] = blue;
                currentBlendColor[3] = alpha;
            }
        },

        enableVertexAttrib: (index: number) => {
            gl.enableVertexAttribArray(index);
            vertexAttribsState[index] = 1;
        },
        clearVertexAttribsState,
        disableUnusedVertexAttribs: () => {
            for (let i = 0; i < maxVertexAttribs; ++i) {
                if (vertexAttribsState[i] === 0) gl.disableVertexAttribArray(i);
            }
        },

        viewport: (x: number, y: number, width: number, height: number) => {
            if (x !== currentViewport[0] || y !== currentViewport[1] || width !== currentViewport[2] || height !== currentViewport[3]) {
                gl.viewport(x, y, width, height);
                currentViewport[0] = x;
                currentViewport[1] = y;
                currentViewport[2] = width;
                currentViewport[3] = height;
            }
        },

        scissor: (x: number, y: number, width: number, height: number) => {
            if (x !== currentScissor[0] || y !== currentScissor[1] || width !== currentScissor[2] || height !== currentScissor[3]) {
                gl.scissor(x, y, width, height);
                currentScissor[0] = x;
                currentScissor[1] = y;
                currentScissor[2] = width;
                currentScissor[3] = height;
            }
        },

        reset: () => {
            enabledCapabilities = {};

            currentFrontFace = gl.getParameter(gl.FRONT_FACE);
            currentCullFace = gl.getParameter(gl.CULL_FACE_MODE);
            currentDepthMask = gl.getParameter(gl.DEPTH_WRITEMASK);
            currentClearDepth = gl.getParameter(gl.DEPTH_CLEAR_VALUE);
            currentDepthFunc = gl.getParameter(gl.DEPTH_FUNC);
            currentColorMask = gl.getParameter(gl.COLOR_WRITEMASK);
            currentClearColor = gl.getParameter(gl.COLOR_CLEAR_VALUE);

            currentBlendSrcRGB = gl.getParameter(gl.BLEND_SRC_RGB);
            currentBlendDstRGB = gl.getParameter(gl.BLEND_DST_RGB);
            currentBlendSrcAlpha = gl.getParameter(gl.BLEND_SRC_ALPHA);
            currentBlendDstAlpha = gl.getParameter(gl.BLEND_DST_ALPHA);
            currentBlendColor = gl.getParameter(gl.BLEND_COLOR);

            currentBlendEqRGB = gl.getParameter(gl.BLEND_EQUATION_RGB);
            currentBlendEqAlpha = gl.getParameter(gl.BLEND_EQUATION_ALPHA);

            maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
            vertexAttribsState.length = 0;
            for (let i = 0; i < maxVertexAttribs; ++i) {
                vertexAttribsState[i] = 0;
            }

            currentViewport = gl.getParameter(gl.VIEWPORT);
            currentScissor = gl.getParameter(gl.SCISSOR_BOX);
        }
    };
}