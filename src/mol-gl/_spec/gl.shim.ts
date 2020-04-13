/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '../../mol-util/id-factory';

const c = {
    ACTIVE_ATTRIBUTE_MAX_LENGTH: 35722,
    ACTIVE_UNIFORM_MAX_LENGTH: 35719,
    INFO_LOG_LENGTH: 35716,
    NUM_COMPRESSED_TEXTURE_FORMATS: 34466,
    SHADER_COMPILER: 36346,
    SHADER_SOURCE_LENGTH: 35720,

    MAX_TEXTURE_MAX_ANISOTROPY_EXT: 0x84FF,
    MAX_TEXTURE_IMAGE_UNITS_NV: 0x8872
};

const gl = {
    ACTIVE_ATTRIBUTES: 35721,
    ACTIVE_TEXTURE: 34016,
    ACTIVE_UNIFORMS: 35718,
    ALIASED_LINE_WIDTH_RANGE: 33902,
    ALIASED_POINT_SIZE_RANGE: 33901,
    ALPHA: 6406,
    ALPHA_BITS: 3413,
    ALWAYS: 519,
    ARRAY_BUFFER: 34962,
    ARRAY_BUFFER_BINDING: 34964,
    ATTACHED_SHADERS: 35717,
    BACK: 1029,
    BLEND: 3042,
    BLEND_COLOR: 32773,
    BLEND_DST_ALPHA: 32970,
    BLEND_DST_RGB: 32968,
    BLEND_EQUATION: 32777,
    BLEND_EQUATION_ALPHA: 34877,
    BLEND_EQUATION_RGB: 32777,
    BLEND_SRC_ALPHA: 32971,
    BLEND_SRC_RGB: 32969,
    BLUE_BITS: 3412,
    BOOL: 35670,
    BOOL_VEC2: 35671,
    BOOL_VEC3: 35672,
    BOOL_VEC4: 35673,
    BROWSER_DEFAULT_WEBGL: 37444,
    BUFFER_SIZE: 34660,
    BUFFER_USAGE: 34661,
    BYTE: 5120,
    CCW: 2305,
    CLAMP_TO_EDGE: 33071,
    COLOR_ATTACHMENT0: 36064,
    COLOR_BUFFER_BIT: 16384,
    COLOR_CLEAR_VALUE: 3106,
    COLOR_WRITEMASK: 3107,
    COMPILE_STATUS: 35713,
    COMPRESSED_TEXTURE_FORMATS: 34467,
    CONSTANT_ALPHA: 32771,
    CONSTANT_COLOR: 32769,
    CONTEXT_LOST_WEBGL: 37442,
    CULL_FACE: 2884,
    CULL_FACE_MODE: 2885,
    CURRENT_PROGRAM: 35725,
    CURRENT_VERTEX_ATTRIB: 34342,
    CW: 2304,
    DECR: 7683,
    DECR_WRAP: 34056,
    DELETE_STATUS: 35712,
    DEPTH_ATTACHMENT: 36096,
    DEPTH_BITS: 3414,
    DEPTH_BUFFER_BIT: 256,
    DEPTH_CLEAR_VALUE: 2931,
    DEPTH_COMPONENT: 6402,
    DEPTH_COMPONENT16: 33189,
    DEPTH_FUNC: 2932,
    DEPTH_RANGE: 2928,
    DEPTH_STENCIL: 34041,
    DEPTH_STENCIL_ATTACHMENT: 33306,
    DEPTH_TEST: 2929,
    DEPTH_WRITEMASK: 2930,
    DITHER: 3024,
    DONT_CARE: 4352,
    DST_ALPHA: 772,
    DST_COLOR: 774,
    DYNAMIC_DRAW: 35048,
    ELEMENT_ARRAY_BUFFER: 34963,
    ELEMENT_ARRAY_BUFFER_BINDING: 34965,
    EQUAL: 514,
    FASTEST: 4353,
    FLOAT: 5126,
    FLOAT_MAT2: 35674,
    FLOAT_MAT3: 35675,
    FLOAT_MAT4: 35676,
    FLOAT_VEC2: 35664,
    FLOAT_VEC3: 35665,
    FLOAT_VEC4: 35666,
    FRAGMENT_SHADER: 35632,
    FRAMEBUFFER: 36160,
    FRAMEBUFFER_ATTACHMENT_OBJECT_NAME: 36049,
    FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE: 36048,
    FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE: 36051,
    FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL: 36050,
    FRAMEBUFFER_BINDING: 36006,
    FRAMEBUFFER_COMPLETE: 36053,
    FRAMEBUFFER_INCOMPLETE_ATTACHMENT: 36054,
    FRAMEBUFFER_INCOMPLETE_DIMENSIONS: 36057,
    FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: 36055,
    FRAMEBUFFER_UNSUPPORTED: 36061,
    FRONT: 1028,
    FRONT_AND_BACK: 1032,
    FRONT_FACE: 2886,
    FUNC_ADD: 32774,
    FUNC_REVERSE_SUBTRACT: 32779,
    FUNC_SUBTRACT: 32778,
    GENERATE_MIPMAP_HINT: 33170,
    GEQUAL: 518,
    GREATER: 516,
    GREEN_BITS: 3411,
    HIGH_FLOAT: 36338,
    HIGH_INT: 36341,
    IMPLEMENTATION_COLOR_READ_FORMAT: 0x8B9B,
    IMPLEMENTATION_COLOR_READ_TYPE: 0x8B9A,
    INCR: 7682,
    INCR_WRAP: 34055,
    INT: 5124,
    INT_VEC2: 35667,
    INT_VEC3: 35668,
    INT_VEC4: 35669,
    INVALID_ENUM: 1280,
    INVALID_FRAMEBUFFER_OPERATION: 1286,
    INVALID_OPERATION: 1282,
    INVALID_VALUE: 1281,
    INVERT: 5386,
    KEEP: 7680,
    LEQUAL: 515,
    LESS: 513,
    LINEAR: 9729,
    LINEAR_MIPMAP_LINEAR: 9987,
    LINEAR_MIPMAP_NEAREST: 9985,
    LINES: 1,
    LINE_LOOP: 2,
    LINE_STRIP: 3,
    LINE_WIDTH: 2849,
    LINK_STATUS: 35714,
    LOW_FLOAT: 36336,
    LOW_INT: 36339,
    LUMINANCE: 6409,
    LUMINANCE_ALPHA: 6410,
    MAX_COMBINED_TEXTURE_IMAGE_UNITS: 35661,
    MAX_CUBE_MAP_TEXTURE_SIZE: 34076,
    MAX_FRAGMENT_UNIFORM_VECTORS: 36349,
    MAX_RENDERBUFFER_SIZE: 34024,
    MAX_TEXTURE_IMAGE_UNITS: 34930,
    MAX_TEXTURE_SIZE: 3379,
    MAX_VARYING_VECTORS: 36348,
    MAX_VERTEX_ATTRIBS: 34921,
    MAX_VERTEX_TEXTURE_IMAGE_UNITS: 35660,
    MAX_VERTEX_UNIFORM_VECTORS: 36347,
    MAX_VIEWPORT_DIMS: 3386,
    MEDIUM_FLOAT: 36337,
    MEDIUM_INT: 36340,
    MIRRORED_REPEAT: 33648,
    NEAREST: 9728,
    NEAREST_MIPMAP_LINEAR: 9986,
    NEAREST_MIPMAP_NEAREST: 9984,
    NEVER: 512,
    NICEST: 4354,
    NONE: 0,
    NOTEQUAL: 517,
    NO_ERROR: 0,
    ONE: 1,
    ONE_MINUS_CONSTANT_ALPHA: 32772,
    ONE_MINUS_CONSTANT_COLOR: 32770,
    ONE_MINUS_DST_ALPHA: 773,
    ONE_MINUS_DST_COLOR: 775,
    ONE_MINUS_SRC_ALPHA: 771,
    ONE_MINUS_SRC_COLOR: 769,
    OUT_OF_MEMORY: 1285,
    PACK_ALIGNMENT: 3333,
    POINTS: 0,
    POLYGON_OFFSET_FACTOR: 32824,
    POLYGON_OFFSET_FILL: 32823,
    POLYGON_OFFSET_UNITS: 10752,
    RED_BITS: 3410,
    RENDERBUFFER: 36161,
    RENDERBUFFER_ALPHA_SIZE: 36179,
    RENDERBUFFER_BINDING: 36007,
    RENDERBUFFER_BLUE_SIZE: 36178,
    RENDERBUFFER_DEPTH_SIZE: 36180,
    RENDERBUFFER_GREEN_SIZE: 36177,
    RENDERBUFFER_HEIGHT: 36163,
    RENDERBUFFER_INTERNAL_FORMAT: 36164,
    RENDERBUFFER_RED_SIZE: 36176,
    RENDERBUFFER_STENCIL_SIZE: 36181,
    RENDERBUFFER_WIDTH: 36162,
    RENDERER: 7937,
    REPEAT: 10497,
    REPLACE: 7681,
    RGB: 6407,
    RGB5_A1: 32855,
    RGB565: 36194,
    RGBA: 6408,
    RGBA4: 32854,
    SAMPLER_2D: 35678,
    SAMPLER_CUBE: 35680,
    SAMPLES: 32937,
    SAMPLE_ALPHA_TO_COVERAGE: 32926,
    SAMPLE_BUFFERS: 32936,
    SAMPLE_COVERAGE: 32928,
    SAMPLE_COVERAGE_INVERT: 32939,
    SAMPLE_COVERAGE_VALUE: 32938,
    SCISSOR_BOX: 3088,
    SCISSOR_TEST: 3089,
    SHADER_TYPE: 35663,
    SHADING_LANGUAGE_VERSION: 35724,
    SHORT: 5122,
    SRC_ALPHA: 770,
    SRC_ALPHA_SATURATE: 776,
    SRC_COLOR: 768,
    STATIC_DRAW: 35044,
    STENCIL_ATTACHMENT: 36128,
    STENCIL_BACK_FAIL: 34817,
    STENCIL_BACK_FUNC: 34816,
    STENCIL_BACK_PASS_DEPTH_FAIL: 34818,
    STENCIL_BACK_PASS_DEPTH_PASS: 34819,
    STENCIL_BACK_REF: 36003,
    STENCIL_BACK_VALUE_MASK: 36004,
    STENCIL_BACK_WRITEMASK: 36005,
    STENCIL_BITS: 3415,
    STENCIL_BUFFER_BIT: 1024,
    STENCIL_CLEAR_VALUE: 2961,
    STENCIL_FAIL: 2964,
    STENCIL_FUNC: 2962,
    STENCIL_INDEX: 6401,
    STENCIL_INDEX8: 36168,
    STENCIL_PASS_DEPTH_FAIL: 2965,
    STENCIL_PASS_DEPTH_PASS: 2966,
    STENCIL_REF: 2967,
    STENCIL_TEST: 2960,
    STENCIL_VALUE_MASK: 2963,
    STENCIL_WRITEMASK: 2968,
    STREAM_DRAW: 35040,
    SUBPIXEL_BITS: 3408,
    TEXTURE: 5890,
    TEXTURE0: 33984,
    TEXTURE1: 33985,
    TEXTURE2: 33986,
    TEXTURE3: 33987,
    TEXTURE4: 33988,
    TEXTURE5: 33989,
    TEXTURE6: 33990,
    TEXTURE7: 33991,
    TEXTURE8: 33992,
    TEXTURE9: 33993,
    TEXTURE10: 33994,
    TEXTURE11: 33995,
    TEXTURE12: 33996,
    TEXTURE13: 33997,
    TEXTURE14: 33998,
    TEXTURE15: 33999,
    TEXTURE16: 34000,
    TEXTURE17: 34001,
    TEXTURE18: 34002,
    TEXTURE19: 34003,
    TEXTURE20: 34004,
    TEXTURE21: 34005,
    TEXTURE22: 34006,
    TEXTURE23: 34007,
    TEXTURE24: 34008,
    TEXTURE25: 34009,
    TEXTURE26: 34010,
    TEXTURE27: 34011,
    TEXTURE28: 34012,
    TEXTURE29: 34013,
    TEXTURE30: 34014,
    TEXTURE31: 34015,
    TEXTURE_2D: 3553,
    TEXTURE_BINDING_2D: 32873,
    TEXTURE_BINDING_CUBE_MAP: 34068,
    TEXTURE_CUBE_MAP: 34067,
    TEXTURE_CUBE_MAP_NEGATIVE_X: 34070,
    TEXTURE_CUBE_MAP_NEGATIVE_Y: 34072,
    TEXTURE_CUBE_MAP_NEGATIVE_Z: 34074,
    TEXTURE_CUBE_MAP_POSITIVE_X: 34069,
    TEXTURE_CUBE_MAP_POSITIVE_Y: 34071,
    TEXTURE_CUBE_MAP_POSITIVE_Z: 34073,
    TEXTURE_MAG_FILTER: 10240,
    TEXTURE_MIN_FILTER: 10241,
    TEXTURE_WRAP_S: 10242,
    TEXTURE_WRAP_T: 10243,
    TRIANGLES: 4,
    TRIANGLE_FAN: 6,
    TRIANGLE_STRIP: 5,
    UNPACK_ALIGNMENT: 3317,
    UNPACK_COLORSPACE_CONVERSION_WEBGL: 37443,
    UNPACK_FLIP_Y_WEBGL: 37440,
    UNPACK_PREMULTIPLY_ALPHA_WEBGL: 37441,
    UNSIGNED_BYTE: 5121,
    UNSIGNED_INT: 5125,
    UNSIGNED_SHORT: 5123,
    UNSIGNED_SHORT_4_4_4_4: 32819,
    UNSIGNED_SHORT_5_5_5_1: 32820,
    UNSIGNED_SHORT_5_6_5: 33635,
    VALIDATE_STATUS: 35715,
    VENDOR: 7936,
    VERSION: 7938,
    VERTEX_ATTRIB_ARRAY_BUFFER_BINDING: 34975,
    VERTEX_ATTRIB_ARRAY_ENABLED: 34338,
    VERTEX_ATTRIB_ARRAY_NORMALIZED: 34922,
    VERTEX_ATTRIB_ARRAY_POINTER: 34373,
    VERTEX_ATTRIB_ARRAY_SIZE: 34339,
    VERTEX_ATTRIB_ARRAY_STRIDE: 34340,
    VERTEX_ATTRIB_ARRAY_TYPE: 34341,
    VERTEX_SHADER: 35633,
    VIEWPORT: 2978,
    ZERO: 0
};
type gl = typeof gl

export function createGl(width: number, height: number, contextAttributes: WebGLContextAttributes): WebGLRenderingContext {
    const getNextId = idFactory();
    const items: { [k: number]: any } = {};
    const boundTextures: { [k: number]: any } = {};
    const viewport = { x: 0, y: 0, width, height };
    const depthRange = { zNear: 0, zFar: 1 };
    const colorClearValue = { r: 0, g: 0, b: 0, a: 0 };

    // let _activeFramebuffer: number
    // let _activeRenderbuffer: number

    return {
        ...gl,
        canvas: {
            width, height
        } as HTMLCanvasElement,
        getAttachedShaders: function(program: WebGLProgram) {
            return [] as WebGLShader[];
        },
        getBufferParameter: function(target: number, pname: number) {
            return 0;
        },
        getContextAttributes: function() { return contextAttributes; },
        getFramebufferAttachmentParameter: function() {},
        getProgramInfoLog: function() { return ''; },
        getShaderInfoLog: function() { return ''; },
        getRenderbufferParameter: function() {},
        getShaderPrecisionFormat: function(shadertype: number, precisiontype: number) {
            return {
                precision: 0,
                rangeMax: 0,
                rangeMin: 0
            };
        },
        getShaderSource: function(shader: WebGLShader | null) { return ''; },
        getTexParameter: function() {},
        getUniform: function() {},
        getVertexAttrib: function() {},
        getVertexAttribOffset: function(index: number, pname: number) { return 0; },
        hint: function(target: number, mode: number) {},
        isBuffer: function(buffer: WebGLBuffer | null) { return true; },
        isEnabled: function(cap: number) { return true; },
        isFramebuffer: function(framebuffer: WebGLFramebuffer | null) { return true; },
        isProgram: function(program: WebGLProgram | null) { return true; },
        isRenderbuffer: function(renderbuffer: WebGLRenderbuffer | null) { return true; },
        isShader: function(shader: WebGLShader | null) { return true; },
        isTexture: function(texture: WebGLTexture | null) { return true; },
        getExtension: function (extensionName: string): any {
            switch (extensionName) {
                case 'EXT_blend_minmax': return {
                    MAX_EXT: 0,
                    MIN_EXT: 0
                } as EXT_blend_minmax;
                case 'EXT_texture_filter_anisotropic': return {
                    MAX_TEXTURE_MAX_ANISOTROPY_EXT: 0,
                    TEXTURE_MAX_ANISOTROPY_EXT: 0
                } as EXT_texture_filter_anisotropic;
                case 'EXT_frag_depth': return {} as EXT_frag_depth;
                case 'EXT_shader_texture_lod': return {} as EXT_shader_texture_lod;
                case 'EXT_sRGB': return {
                    FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING_EXT: 0,
                    SRGB8_ALPHA8_EXT: 0,
                    SRGB_ALPHA_EXT: 0,
                    SRGB_EXT: 0
                } as EXT_sRGB;
                case 'OES_vertex_array_object': return {
                    VERTEX_ARRAY_BINDING_OES: 0,
                    bindVertexArrayOES: function(arrayObject: WebGLVertexArrayObjectOES) { },
                    createVertexArrayOES: function(): WebGLVertexArrayObjectOES { return {}; },
                    deleteVertexArrayOES: function(arrayObject: WebGLVertexArrayObjectOES) { },
                    isVertexArrayOES: function(value: any) { return true; }
                } as OES_vertex_array_object;
                case 'WEBGL_color_buffer_float': return {
                    FRAMEBUFFER_ATTACHMENT_COMPONENT_TYPE_EXT: 0,
                    RGB32F_EXT: 0,
                    RGBA32F_EXT: 0,
                    UNSIGNED_NORMALIZED_EXT: 0
                } as WEBGL_color_buffer_float;
                case 'WEBGL_compressed_texture_astc': return null;
                case 'WEBGL_compressed_texture_s3tc_srgb': return null;
                case 'WEBGL_debug_shaders': return {
                    getTranslatedShaderSource(shader: WebGLShader) { return ''; }
                } as WEBGL_debug_shaders;
                case 'WEBGL_draw_buffers': return null;
                case 'WEBGL_lose_context': return {
                    loseContext: function() { },
                    restoreContext: function() { },
                } as WEBGL_lose_context;
                case 'WEBGL_depth_texture': return {
                    UNSIGNED_INT_24_8_WEBGL: 0
                } as WEBGL_depth_texture;
                case 'WEBGL_debug_renderer_info': return {
                    UNMASKED_RENDERER_WEBGL: 0,
                    UNMASKED_VENDOR_WEBGL: 0
                } as WEBGL_debug_renderer_info;
                case 'WEBGL_compressed_texture_s3tc': return null;
                case 'OES_texture_half_float_linear': return {} as OES_texture_half_float_linear;
                case 'OES_texture_half_float': return {
                    HALF_FLOAT_OES: 0
                } as OES_texture_half_float;
                case 'OES_texture_float_linear': return {} as OES_texture_float_linear;
                case 'OES_texture_float': return {} as OES_texture_float;
                case 'OES_standard_derivatives': return {
                    FRAGMENT_SHADER_DERIVATIVE_HINT_OES: 0
                } as OES_standard_derivatives;
                case 'OES_element_index_uint': return {} as OES_element_index_uint;
                case 'ANGLE_instanced_arrays': return {
                    drawArraysInstancedANGLE: function(mode: number, first: number, count: number, primcount: number) {},
                    drawElementsInstancedANGLE: function(mode: number, count: number, type: number, offset: number, primcount: number) {},
                    vertexAttribDivisorANGLE: function(index: number, divisor: number) {},
                    VERTEX_ATTRIB_ARRAY_DIVISOR_ANGLE: 0
                } as ANGLE_instanced_arrays;
            }
            return null;
        },
        createBuffer: function () {
            const id = getNextId();
            items[id] = {
                which: 'buffer',
            };
            return id;
        },
        deleteBuffer: function () { },
        bindBuffer: function () { },
        bufferData: function () { },
        getParameter: function (pname: number) {
            switch (pname) {
                case c.MAX_TEXTURE_MAX_ANISOTROPY_EXT: return 16;
                case c.MAX_TEXTURE_IMAGE_UNITS_NV: return 16;

                case gl.ELEMENT_ARRAY_BUFFER_BINDING:
                case gl.ARRAY_BUFFER_BINDING:
                case gl.FRAMEBUFFER_BINDING:
                case gl.CURRENT_PROGRAM:
                case gl.RENDERBUFFER_BINDING:
                    return 0;
                    // return _activeFramebuffer
                    // return _activeRenderbuffer
                case gl.TEXTURE_BINDING_2D:
                case gl.TEXTURE_BINDING_CUBE_MAP:
                    return null;

                case gl.VERSION:
                    return '1.0.0';
                case gl.VENDOR:
                    return 'shim';
                case gl.RENDERER:
                    return 'shim-renderer';
                case gl.SHADING_LANGUAGE_VERSION:
                    return 'WebGL GLSL ES 1.0 shim';

                case gl.COMPRESSED_TEXTURE_FORMATS:
                    return new Uint32Array(0);

                // Int arrays
                case gl.MAX_VIEWPORT_DIMS:
                case gl.SCISSOR_BOX:
                    return new Int32Array([ 0, 0, 4096, 4096 ]);
                case gl.VIEWPORT:
                    const { x, y, width, height } = viewport;
                    return new Int32Array([ x, y, width, height ]);

                // Float arrays
                case gl.ALIASED_LINE_WIDTH_RANGE:
                    return new Float32Array([0, 1]);
                case gl.ALIASED_POINT_SIZE_RANGE:
                    return new Float32Array([0, 255]);
                case gl.DEPTH_RANGE:
                    return new Float32Array([ depthRange.zNear, depthRange.zFar ]);
                case gl.BLEND_COLOR:
                    return new Float32Array([0, 0, 0, 0]);
                case gl.COLOR_CLEAR_VALUE:
                    const { r, g, b, a } = colorClearValue;
                    return new Float32Array([ r, g, b, a ]);

                case gl.COLOR_WRITEMASK:
                    return 0;

                case gl.DEPTH_CLEAR_VALUE:
                case gl.LINE_WIDTH:
                case gl.POLYGON_OFFSET_FACTOR:
                case gl.POLYGON_OFFSET_UNITS:
                case gl.SAMPLE_COVERAGE_VALUE:
                    return 1;

                case gl.BLEND:
                case gl.CULL_FACE:
                case gl.DEPTH_TEST:
                case gl.DEPTH_WRITEMASK:
                case gl.DITHER:
                case gl.POLYGON_OFFSET_FILL:
                case gl.SAMPLE_COVERAGE_INVERT:
                case gl.SCISSOR_TEST:
                case gl.STENCIL_TEST:
                case gl.UNPACK_FLIP_Y_WEBGL:
                case gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL:
                    return false;

                case gl.MAX_TEXTURE_SIZE:
                case gl.MAX_CUBE_MAP_TEXTURE_SIZE:
                    return 16384;

                case gl.MAX_VERTEX_UNIFORM_VECTORS:
                case gl.MAX_FRAGMENT_UNIFORM_VECTORS:
                    return 4096;

                case gl.MAX_VARYING_VECTORS:
                case gl.MAX_TEXTURE_IMAGE_UNITS:
                case gl.MAX_COMBINED_TEXTURE_IMAGE_UNITS:
                case gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS:
                    return 32;

                case gl.MAX_RENDERBUFFER_SIZE:
                    return 2048;

                case gl.ACTIVE_TEXTURE:
                case gl.ALPHA_BITS:
                case gl.BLEND_DST_ALPHA:
                case gl.BLEND_DST_RGB:
                case gl.BLEND_EQUATION_ALPHA:
                case gl.BLEND_EQUATION_RGB:
                case gl.BLEND_SRC_ALPHA:
                case gl.BLEND_SRC_RGB:
                case gl.BLUE_BITS:
                case gl.CULL_FACE_MODE:
                case gl.DEPTH_BITS:
                case gl.DEPTH_FUNC:
                case gl.FRONT_FACE:
                case gl.GENERATE_MIPMAP_HINT:
                case gl.GREEN_BITS:
                case gl.MAX_VERTEX_ATTRIBS:
                case gl.PACK_ALIGNMENT:
                case gl.RED_BITS:
                case gl.SAMPLE_BUFFERS:
                case gl.SAMPLES:
                case gl.STENCIL_BACK_FAIL:
                case gl.STENCIL_BACK_FUNC:
                case gl.STENCIL_BACK_PASS_DEPTH_FAIL:
                case gl.STENCIL_BACK_PASS_DEPTH_PASS:
                case gl.STENCIL_BACK_REF:
                case gl.STENCIL_BACK_VALUE_MASK:
                case gl.STENCIL_BACK_WRITEMASK:
                case gl.STENCIL_BITS:
                case gl.STENCIL_CLEAR_VALUE:
                case gl.STENCIL_FAIL:
                case gl.STENCIL_FUNC:
                case gl.STENCIL_PASS_DEPTH_FAIL:
                case gl.STENCIL_PASS_DEPTH_PASS:
                case gl.STENCIL_REF:
                case gl.STENCIL_VALUE_MASK:
                case gl.STENCIL_WRITEMASK:
                case gl.SUBPIXEL_BITS:
                case gl.UNPACK_ALIGNMENT:
                case gl.UNPACK_COLORSPACE_CONVERSION_WEBGL:
                    return 0;

                default:
                    return 0;
            }
        },
        getSupportedExtensions: function () {
            return [
                'EXT_blend_minmax', 'EXT_texture_filter_anisotropic', 'EXT_frag_depth',
                'EXT_shader_texture_lod', 'EXT_sRGB', 'OES_vertex_array_object',
                'WEBGL_color_buffer_float', 'WEBGL_debug_shaders', 'WEBGL_lose_context',
                'WEBGL_depth_texture', 'WEBGL_debug_renderer_info', 'OES_texture_half_float_linear',
                'OES_texture_half_float', 'OES_texture_float_linear', 'OES_texture_float',
                'OES_standard_derivatives', 'OES_element_index_uint', 'ANGLE_instanced_arrays'
            ];
        },
        createShader: function (type: number) {
            const id = getNextId();
            items[id] = {
                which: 'shader',
                type: type,
            };
            return id as WebGLShader;
        },
        getShaderParameter: function (shader: WebGLShader, pname: number) {
            switch (pname) {
                case gl.SHADER_TYPE: return items[shader as number].type;
                case gl.COMPILE_STATUS: return true;
                default: throw `getShaderParameter ${pname}`;
            }
        },
        shaderSource: function () { },
        compileShader: function () { },
        createProgram: function () {
            const id = getNextId();
            items[id] = {
                which: 'program',
                shaders: [],
            };
            return id;
        },
        attachShader: function (program: WebGLProgram, shader: WebGLShader) {
            items[program as number].shaders.push(shader as number);
        },
        bindAttribLocation: function () { },
        linkProgram: function () { },
        getProgramParameter: function (program: number, pname: number) {
            switch (pname) {
                case gl.LINK_STATUS: return true;
                case gl.ACTIVE_UNIFORMS: return 4;
                case gl.ACTIVE_ATTRIBUTES: return 0;
                case gl.ACTIVE_UNIFORMS: return 0;
                case gl.DELETE_STATUS: return false;
                case gl.VALIDATE_STATUS: return true;
                case gl.ATTACHED_SHADERS: return 2;
                default: throw `getProgramParameter ${pname}`;
            }
        },
        deleteShader: function () { },
        deleteProgram: function () { },
        viewport: function (x: number, y: number, width: number, height: number) {
            viewport.x = x;
            viewport.y = y;
            viewport.width = width;
            viewport.height = height;
        },
        clearColor: function (red: number, green: number, blue: number, alpha: number) {
            colorClearValue.r = red;
            colorClearValue.g = green;
            colorClearValue.b = blue;
            colorClearValue.a = alpha;
        },
        clearDepth: function () { },
        depthFunc: function () { },
        enable: function () { },
        disable: function () { },
        frontFace: function () { },
        cullFace: function () { },
        activeTexture: function () { },
        createTexture: function () {
            const id = getNextId();
            items[id] = {
                which: 'texture',
            };
            return id;
        },
        deleteTexture: function () { },
        bindTexture: function (target: number, texture: WebGLTexture | null) {
            boundTextures[target as number] = texture;
        },
        texParameterf: function () { },
        texParameteri: function () { },
        pixelStorei: function () { },
        texImage2D: function () { },
        texSubImage2D: function () { },
        compressedTexImage2D: function () { },
        useProgram: function () { },
        getUniformLocation: function () {
            return 0;
        },
        getActiveUniform: function (program: WebGLProgram, index: number) {
            return {
                size: 1,
                type: gl.INT_VEC3,
                name: `__activeUniform${index}`,
            };
        },
        getActiveAttrib: function (program: WebGLProgram, index: number) {
            return {
                size: 1,
                type: gl.FLOAT,
                name: `__activeAttrib${index}`
            };
        },
        clear: function () { },
        uniform1f: function () { },
        uniform1fv: function () { },
        uniform1i: function () { },
        uniform1iv: function () { },
        uniform2f: function () { },
        uniform2fv: function () { },
        uniform2i: function () { },
        uniform2iv: function () { },
        uniform3f: function () { },
        uniform3fv: function () { },
        uniform3i: function () { },
        uniform3iv: function () { },
        uniform4f: function () { },
        uniform4fv: function () { },
        uniform4i: function () { },
        uniform4iv: function () { },
        uniformMatrix2fv: function () { },
        uniformMatrix3fv: function () { },
        uniformMatrix4fv: function () { },
        getAttribLocation: function () { return 1; },
        vertexAttribPointer: function () { },
        enableVertexAttribArray: function () { },
        disableVertexAttribArray: function () { },
        drawElements: function () { },
        drawArrays: function () { },
        depthMask: function () { },
        depthRange: function (zNear: number, zFar: number) {
            depthRange.zNear = zNear;
            depthRange.zFar = zFar;
        },
        bufferSubData: function () { },
        blendFunc: function () { },
        createFramebuffer: function () {
            const id = getNextId();
            items[id] = {
                which: 'framebuffer',
                shaders: [],
            };
            return id;
        },
        bindFramebuffer: function () { },
        framebufferTexture2D: function () { },
        checkFramebufferStatus: function () {
            return gl.FRAMEBUFFER_COMPLETE;
        },
        deleteFramebuffer: function () { },
        createRenderbuffer: function () {
            const id = getNextId();
            items[id] = {
                which: 'renderbuffer',
                shaders: [],
            };
            return id;
        },
        bindRenderbuffer: function () { },
        deleteRenderbuffer: function () { },
        renderbufferStorage: function () { },
        framebufferRenderbuffer: function () { },
        scissor: function () { },
        colorMask: function () { },
        lineWidth: function () { },
        vertexAttrib1f: function () { },
        vertexAttrib1fv: function () { },
        vertexAttrib2f: function () { },
        vertexAttrib2fv: function () { },
        vertexAttrib3f: function () { },
        vertexAttrib3fv: function () { },
        vertexAttrib4f: function () { },
        vertexAttrib4fv: function () { },
        validateProgram: function () { },
        generateMipmap: function () { },
        isContextLost: function () { return false; },
        drawingBufferWidth: 1024,
        drawingBufferHeight: 1024,
        blendColor: function () { },
        blendEquation: function () { },
        blendEquationSeparate: function () { },
        blendFuncSeparate: function () { },
        clearStencil: function () { },
        compressedTexSubImage2D: function () { },
        copyTexImage2D: function () { },
        copyTexSubImage2D: function () { },
        detachShader: function () { },
        finish: function () { },
        flush: function () { },
        getError: function () { return 0; },
        polygonOffset: function () { },
        readPixels: function () { },
        sampleCoverage: function () { },
        stencilFunc: function () { },
        stencilFuncSeparate: function () { },
        stencilMask: function () { },
        stencilMaskSeparate: function () { },
        stencilOp: function () { },
        stencilOpSeparate: function () { },
    };
}