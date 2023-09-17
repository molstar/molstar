/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext } from './compat';
import { WebGLExtensions } from './extensions';

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
     * - `ext.CLIP_DISTANCE[0-7]`: clip distance 0 to 7 (with `ext` being `WEBGL_clip_cull_distance`)
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
     * - `ext.CLIP_DISTANCE[0-7]`: clip distance 0 to 7 (with `ext` being `WEBGL_clip_cull_distance`)
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

    /** sets the front and back function and reference value for stencil testing */
    stencilFunc: (func: number, ref: number, mask: number) => void
    /** sets the front and/or back function and reference value for stencil testing */
    stencilFuncSeparate: (face: number, func: number, ref: number, mask: number) => void
    /** controls enabling and disabling of both the front and back writing of individual bits in the stencil planes */
    stencilMask: (mask: number) => void
    /** controls enabling and disabling of both the front and back writing of individual bits in the stencil planes */
    stencilMaskSeparate: (face: number, mask: number) => void
    /** sets both the front and back-facing stencil test actions */
    stencilOp: (fail: number, zfail: number, zpass: number) => void
    /** sets the front and/or back-facing stencil test actions */
    stencilOpSeparate: (face: number, fail: number, zfail: number, zpass: number) => void

    enableVertexAttrib: (index: number) => void
    clearVertexAttribsState: () => void
    disableUnusedVertexAttribs: () => void

    viewport: (x: number, y: number, width: number, height: number) => void
    scissor: (x: number, y: number, width: number, height: number) => void

    /**
     * controls the clipping volume behavior
     * @param origin must be `ext.LOWER_LEFT` (default) or `ext.UPPER_LEFT`.
     * @param depth must be `ext.NEGATIVE_ONE_TO_ONE` (default) or `ext.ZERO_TO_ONE`.
     * with `ext` being `EXT_clip_control`
     */
    clipControl?: (origin: number, depth: number) => void

    reset: () => void
}

export function createState(gl: GLRenderingContext, e: WebGLExtensions): WebGLState {
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

    let currentStencilFunc = gl.getParameter(gl.STENCIL_FUNC);
    let currentStencilValueMask = gl.getParameter(gl.STENCIL_VALUE_MASK);
    let currentStencilRef = gl.getParameter(gl.STENCIL_REF);
    let currentStencilBackFunc = gl.getParameter(gl.STENCIL_BACK_FUNC);
    let currentStencilBackValueMask = gl.getParameter(gl.STENCIL_BACK_VALUE_MASK);
    let currentStencilBackRef = gl.getParameter(gl.STENCIL_BACK_REF);
    let currentStencilWriteMask = gl.getParameter(gl.STENCIL_WRITEMASK);
    let currentStencilBackWriteMask = gl.getParameter(gl.STENCIL_BACK_WRITEMASK);
    let currentStencilFail = gl.getParameter(gl.STENCIL_FAIL);
    let currentStencilPassDepthPass = gl.getParameter(gl.STENCIL_PASS_DEPTH_PASS);
    let currentStencilPassDepthFail = gl.getParameter(gl.STENCIL_PASS_DEPTH_FAIL);
    let currentStencilBackFail = gl.getParameter(gl.STENCIL_BACK_FAIL);
    let currentStencilBackPassDepthPass = gl.getParameter(gl.STENCIL_BACK_PASS_DEPTH_PASS);
    let currentStencilBackPassDepthFail = gl.getParameter(gl.STENCIL_BACK_PASS_DEPTH_FAIL);

    let maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
    const vertexAttribsState: number[] = [];

    let currentViewport: [number, number, number, number] = gl.getParameter(gl.VIEWPORT);
    let currentScissor: [number, number, number, number] = gl.getParameter(gl.SCISSOR_BOX);

    let currentClipOrigin = e.clipControl ? gl.getParameter(e.clipControl.CLIP_ORIGIN) : -1;
    let currentClipDepthMode = e.clipControl ? gl.getParameter(e.clipControl.CLIP_DEPTH_MODE) : -1;

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

        stencilFunc: (func: number, ref: number, mask: number) => {
            if (func !== currentStencilFunc || ref !== currentStencilRef || mask !== currentStencilValueMask || func !== currentStencilBackFunc || ref !== currentStencilBackRef || mask !== currentStencilBackValueMask) {
                gl.stencilFunc(func, ref, mask);
                currentStencilFunc = func;
                currentStencilRef = ref;
                currentStencilValueMask = mask;
                currentStencilBackFunc = func;
                currentStencilBackRef = ref;
                currentStencilBackValueMask = mask;
            }
        },
        stencilFuncSeparate: (face: number, func: number, ref: number, mask: number) => {
            if (face === gl.FRONT) {
                if (func !== currentStencilFunc || ref !== currentStencilRef || mask !== currentStencilValueMask) {
                    gl.stencilFuncSeparate(face, func, ref, mask);
                    currentStencilFunc = func;
                    currentStencilRef = ref;
                    currentStencilValueMask = mask;
                }
            } else if (face === gl.BACK) {
                if (func !== currentStencilBackFunc || ref !== currentStencilBackRef || mask !== currentStencilBackValueMask) {
                    gl.stencilFuncSeparate(face, func, ref, mask);
                    currentStencilBackFunc = func;
                    currentStencilBackRef = ref;
                    currentStencilBackValueMask = mask;
                }
            } else if (face === gl.FRONT_AND_BACK) {
                if (func !== currentStencilFunc || ref !== currentStencilRef || mask !== currentStencilValueMask || func !== currentStencilBackFunc || ref !== currentStencilBackRef || mask !== currentStencilBackValueMask) {
                    gl.stencilFuncSeparate(face, func, ref, mask);
                    currentStencilFunc = func;
                    currentStencilRef = ref;
                    currentStencilValueMask = mask;
                    currentStencilBackFunc = func;
                    currentStencilBackRef = ref;
                    currentStencilBackValueMask = mask;
                }
            }
        },
        stencilMask: (mask: number) => {
            if (mask !== currentStencilWriteMask || mask !== currentStencilBackWriteMask) {
                gl.stencilMask(mask);
                currentStencilWriteMask = mask;
                currentStencilBackWriteMask = mask;
            }
        },
        stencilMaskSeparate: (face: number, mask: number) => {
            if (face === gl.FRONT) {
                if (mask !== currentStencilWriteMask) {
                    gl.stencilMaskSeparate(face, mask);
                    currentStencilWriteMask = mask;
                }
            } else if (face === gl.BACK) {
                if (mask !== currentStencilBackWriteMask) {
                    gl.stencilMaskSeparate(face, mask);
                    currentStencilBackWriteMask = mask;
                }
            } else if (face === gl.FRONT_AND_BACK) {
                if (mask !== currentStencilWriteMask || mask !== currentStencilBackWriteMask) {
                    gl.stencilMaskSeparate(face, mask);
                    currentStencilWriteMask = mask;
                    currentStencilBackWriteMask = mask;
                }
            }
        },
        stencilOp: (fail: number, zfail: number, zpass: number) => {
            if (fail !== currentStencilFail || zfail !== currentStencilPassDepthFail || zpass !== currentStencilPassDepthPass || fail !== currentStencilBackFail || zfail !== currentStencilBackPassDepthFail || zpass !== currentStencilBackPassDepthPass) {
                gl.stencilOp(fail, zfail, zpass);
                currentStencilFail = fail;
                currentStencilPassDepthFail = zfail;
                currentStencilPassDepthPass = zpass;
                currentStencilBackFail = fail;
                currentStencilBackPassDepthFail = zfail;
                currentStencilBackPassDepthPass = zpass;
            }
        },
        stencilOpSeparate: (face: number, fail: number, zfail: number, zpass: number) => {
            if (face === gl.FRONT) {
                if (fail !== currentStencilFail || zfail !== currentStencilPassDepthFail || zpass !== currentStencilPassDepthPass) {
                    gl.stencilOpSeparate(face, fail, zfail, zpass);
                    currentStencilFail = fail;
                    currentStencilPassDepthFail = zfail;
                    currentStencilPassDepthPass = zpass;
                }
            } else if (face === gl.BACK) {
                if (fail !== currentStencilBackFail || zfail !== currentStencilBackPassDepthFail || zpass !== currentStencilBackPassDepthPass) {
                    gl.stencilOpSeparate(face, fail, zfail, zpass);
                    currentStencilBackFail = fail;
                    currentStencilBackPassDepthFail = zfail;
                    currentStencilBackPassDepthPass = zpass;
                }
            } else if (face === gl.FRONT_AND_BACK) {
                if (fail !== currentStencilFail || zfail !== currentStencilPassDepthFail || zpass !== currentStencilPassDepthPass || fail !== currentStencilBackFail || zfail !== currentStencilBackPassDepthFail || zpass !== currentStencilBackPassDepthPass) {
                    gl.stencilOpSeparate(face, fail, zfail, zpass);
                    currentStencilFail = fail;
                    currentStencilPassDepthFail = zfail;
                    currentStencilPassDepthPass = zpass;
                    currentStencilBackFail = fail;
                    currentStencilBackPassDepthFail = zfail;
                    currentStencilBackPassDepthPass = zpass;
                }
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

        clipControl: e.clipControl ? (origin: number, depth: number) => {
            if (origin !== currentClipOrigin || depth !== currentClipDepthMode) {
                e.clipControl!.clipControl(origin, depth);
                currentClipOrigin = origin;
                currentClipDepthMode = depth;
            }
        } : undefined,

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

            currentStencilFunc = gl.getParameter(gl.STENCIL_FUNC);
            currentStencilValueMask = gl.getParameter(gl.STENCIL_VALUE_MASK);
            currentStencilRef = gl.getParameter(gl.STENCIL_REF);
            currentStencilBackFunc = gl.getParameter(gl.STENCIL_BACK_FUNC);
            currentStencilBackValueMask = gl.getParameter(gl.STENCIL_BACK_VALUE_MASK);
            currentStencilBackRef = gl.getParameter(gl.STENCIL_BACK_REF);
            currentStencilWriteMask = gl.getParameter(gl.STENCIL_WRITEMASK);
            currentStencilBackWriteMask = gl.getParameter(gl.STENCIL_BACK_WRITEMASK);
            currentStencilFail = gl.getParameter(gl.STENCIL_FAIL);
            currentStencilPassDepthPass = gl.getParameter(gl.STENCIL_PASS_DEPTH_PASS);
            currentStencilPassDepthFail = gl.getParameter(gl.STENCIL_PASS_DEPTH_FAIL);
            currentStencilBackFail = gl.getParameter(gl.STENCIL_BACK_FAIL);
            currentStencilBackPassDepthPass = gl.getParameter(gl.STENCIL_BACK_PASS_DEPTH_PASS);
            currentStencilBackPassDepthFail = gl.getParameter(gl.STENCIL_BACK_PASS_DEPTH_FAIL);

            maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
            vertexAttribsState.length = 0;
            for (let i = 0; i < maxVertexAttribs; ++i) {
                vertexAttribsState[i] = 0;
            }

            currentViewport = gl.getParameter(gl.VIEWPORT);
            currentScissor = gl.getParameter(gl.SCISSOR_BOX);

            currentClipOrigin = e.clipControl ? gl.getParameter(e.clipControl.CLIP_ORIGIN) : -1;
            currentClipDepthMode = e.clipControl ? gl.getParameter(e.clipControl.CLIP_DEPTH_MODE) : -1;
        }
    };
}