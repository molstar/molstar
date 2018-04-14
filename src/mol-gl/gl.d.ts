// Type definitions for gl 4.0.4
// Project: headless-gl
// Definitions by: Ivan Perevezentsev (https://github.com/ip)

declare module 'gl' {

    /*~ Note that ES6 modules cannot directly export callable functions.
    *~ This file should be imported using the CommonJS-style:
    *~   import x = require('someLibrary');
    *~
    *~ Refer to the documentation to understand common
    *~ workarounds for this limitation of ES6 modules.
    */

    /*~ This declaration specifies that the function
    *~ is the exported object from the file
    */
    export = createWebglContext;

    function createWebglContext(
      width: number,
      height: number,
      options?: WebGLContextAttributes): WebGLRenderingContext;
  }

  interface WebGLContextAttributes {
    failIfMajorPerformanceCaveat?: boolean;
    alpha?: boolean;
    antialias?: boolean;
    depth?: boolean;
    premultipliedAlpha?: boolean;
    preserveDrawingBuffer?: boolean;
    stencil?: boolean;
  }

  interface WEBGL_compressed_texture_s3tc {
      readonly COMPRESSED_RGB_S3TC_DXT1_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT1_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT3_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT5_EXT: number;
  }

  declare var WEBGL_compressed_texture_s3tc: {
      prototype: WEBGL_compressed_texture_s3tc;
      new(): WEBGL_compressed_texture_s3tc;
      readonly COMPRESSED_RGB_S3TC_DXT1_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT1_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT3_EXT: number;
      readonly COMPRESSED_RGBA_S3TC_DXT5_EXT: number;
  };

  interface WEBGL_debug_renderer_info {
      readonly UNMASKED_RENDERER_WEBGL: number;
      readonly UNMASKED_VENDOR_WEBGL: number;
  }

  declare var WEBGL_debug_renderer_info: {
      prototype: WEBGL_debug_renderer_info;
      new(): WEBGL_debug_renderer_info;
      readonly UNMASKED_RENDERER_WEBGL: number;
      readonly UNMASKED_VENDOR_WEBGL: number;
  };

  interface WEBGL_depth_texture {
      readonly UNSIGNED_INT_24_8_WEBGL: number;
  }

  declare var WEBGL_depth_texture: {
      prototype: WEBGL_depth_texture;
      new(): WEBGL_depth_texture;
      readonly UNSIGNED_INT_24_8_WEBGL: number;
  };

  interface WebGLActiveInfo {
      readonly name: string;
      readonly size: number;
      readonly type: number;
  }

  declare var WebGLActiveInfo: {
      prototype: WebGLActiveInfo;
      new(): WebGLActiveInfo;
  };

  interface WebGLBuffer extends WebGLObject {
  }

  declare var WebGLBuffer: {
      prototype: WebGLBuffer;
      new(): WebGLBuffer;
  };

  // interface WebGLContextEvent extends Event {
  //     readonly statusMessage: string;
  // }

  interface WebGLFramebuffer extends WebGLObject {
  }

  declare var WebGLFramebuffer: {
      prototype: WebGLFramebuffer;
      new(): WebGLFramebuffer;
  };

  interface WebGLObject {
  }

  declare var WebGLObject: {
      prototype: WebGLObject;
      new(): WebGLObject;
  };

  interface WebGLProgram extends WebGLObject {
  }

  declare var WebGLProgram: {
      prototype: WebGLProgram;
      new(): WebGLProgram;
  };

  interface WebGLRenderbuffer extends WebGLObject {
  }

  declare var WebGLRenderbuffer: {
      prototype: WebGLRenderbuffer;
      new(): WebGLRenderbuffer;
  };

  interface WebGLRenderingContext {
      readonly drawingBufferHeight: number;
      readonly drawingBufferWidth: number;
      activeTexture(texture: number): void;
      attachShader(program: WebGLProgram | null, shader: WebGLShader | null): void;
      bindAttribLocation(program: WebGLProgram | null, index: number, name: string): void;
      bindBuffer(target: number, buffer: WebGLBuffer | null): void;
      bindFramebuffer(target: number, framebuffer: WebGLFramebuffer | null): void;
      bindRenderbuffer(target: number, renderbuffer: WebGLRenderbuffer | null): void;
      bindTexture(target: number, texture: WebGLTexture | null): void;
      blendColor(red: number, green: number, blue: number, alpha: number): void;
      blendEquation(mode: number): void;
      blendEquationSeparate(modeRGB: number, modeAlpha: number): void;
      blendFunc(sfactor: number, dfactor: number): void;
      blendFuncSeparate(srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number): void;
      bufferData(target: number, size: number | ArrayBufferView | ArrayBuffer, usage: number): void;
      bufferSubData(target: number, offset: number, data: ArrayBufferView | ArrayBuffer): void;
      checkFramebufferStatus(target: number): number;
      clear(mask: number): void;
      clearColor(red: number, green: number, blue: number, alpha: number): void;
      clearDepth(depth: number): void;
      clearStencil(s: number): void;
      colorMask(red: boolean, green: boolean, blue: boolean, alpha: boolean): void;
      compileShader(shader: WebGLShader | null): void;
      compressedTexImage2D(target: number, level: number, internalformat: number, width: number, height: number, border: number, data: ArrayBufferView): void;
      compressedTexSubImage2D(target: number, level: number, xoffset: number, yoffset: number, width: number, height: number, format: number, data: ArrayBufferView): void;
      copyTexImage2D(target: number, level: number, internalformat: number, x: number, y: number, width: number, height: number, border: number): void;
      copyTexSubImage2D(target: number, level: number, xoffset: number, yoffset: number, x: number, y: number, width: number, height: number): void;
      createBuffer(): WebGLBuffer | null;
      createFramebuffer(): WebGLFramebuffer | null;
      createProgram(): WebGLProgram | null;
      createRenderbuffer(): WebGLRenderbuffer | null;
      createShader(type: number): WebGLShader | null;
      createTexture(): WebGLTexture | null;
      cullFace(mode: number): void;
      deleteBuffer(buffer: WebGLBuffer | null): void;
      deleteFramebuffer(framebuffer: WebGLFramebuffer | null): void;
      deleteProgram(program: WebGLProgram | null): void;
      deleteRenderbuffer(renderbuffer: WebGLRenderbuffer | null): void;
      deleteShader(shader: WebGLShader | null): void;
      deleteTexture(texture: WebGLTexture | null): void;
      depthFunc(func: number): void;
      depthMask(flag: boolean): void;
      depthRange(zNear: number, zFar: number): void;
      detachShader(program: WebGLProgram | null, shader: WebGLShader | null): void;
      disable(cap: number): void;
      disableVertexAttribArray(index: number): void;
      drawArrays(mode: number, first: number, count: number): void;
      drawElements(mode: number, count: number, type: number, offset: number): void;
      enable(cap: number): void;
      enableVertexAttribArray(index: number): void;
      finish(): void;
      flush(): void;
      framebufferRenderbuffer(target: number, attachment: number, renderbuffertarget: number, renderbuffer: WebGLRenderbuffer | null): void;
      framebufferTexture2D(target: number, attachment: number, textarget: number, texture: WebGLTexture | null, level: number): void;
      frontFace(mode: number): void;
      generateMipmap(target: number): void;
      getActiveAttrib(program: WebGLProgram | null, index: number): WebGLActiveInfo | null;
      getActiveUniform(program: WebGLProgram | null, index: number): WebGLActiveInfo | null;
      getAttachedShaders(program: WebGLProgram | null): WebGLShader[] | null;
      getAttribLocation(program: WebGLProgram | null, name: string): number;
      getBufferParameter(target: number, pname: number): any;
      getContextAttributes(): WebGLContextAttributes;
      getError(): number;
      getExtension(name: string): any;
      getFramebufferAttachmentParameter(target: number, attachment: number, pname: number): any;
      getParameter(pname: number): any;
      getProgramInfoLog(program: WebGLProgram | null): string | null;
      getProgramParameter(program: WebGLProgram | null, pname: number): any;
      getRenderbufferParameter(target: number, pname: number): any;
      getShaderInfoLog(shader: WebGLShader | null): string | null;
      getShaderParameter(shader: WebGLShader | null, pname: number): any;
      getShaderPrecisionFormat(shadertype: number, precisiontype: number): WebGLShaderPrecisionFormat | null;
      getShaderSource(shader: WebGLShader | null): string | null;
      getSupportedExtensions(): string[] | null;
      getTexParameter(target: number, pname: number): any;
      getUniform(program: WebGLProgram | null, location: WebGLUniformLocation | null): any;
      getUniformLocation(program: WebGLProgram | null, name: string): WebGLUniformLocation | null;
      getVertexAttrib(index: number, pname: number): any;
      getVertexAttribOffset(index: number, pname: number): number;
      hint(target: number, mode: number): void;
      isBuffer(buffer: WebGLBuffer | null): boolean;
      isContextLost(): boolean;
      isEnabled(cap: number): boolean;
      isFramebuffer(framebuffer: WebGLFramebuffer | null): boolean;
      isProgram(program: WebGLProgram | null): boolean;
      isRenderbuffer(renderbuffer: WebGLRenderbuffer | null): boolean;
      isShader(shader: WebGLShader | null): boolean;
      isTexture(texture: WebGLTexture | null): boolean;
      lineWidth(width: number): void;
      linkProgram(program: WebGLProgram | null): void;
      pixelStorei(pname: number, param: number | boolean): void;
      polygonOffset(factor: number, units: number): void;
      readPixels(x: number, y: number, width: number, height: number, format: number, type: number, pixels: ArrayBufferView | null): void;
      renderbufferStorage(target: number, internalformat: number, width: number, height: number): void;
      sampleCoverage(value: number, invert: boolean): void;
      scissor(x: number, y: number, width: number, height: number): void;
      shaderSource(shader: WebGLShader | null, source: string): void;
      stencilFunc(func: number, ref: number, mask: number): void;
      stencilFuncSeparate(face: number, func: number, ref: number, mask: number): void;
      stencilMask(mask: number): void;
      stencilMaskSeparate(face: number, mask: number): void;
      stencilOp(fail: number, zfail: number, zpass: number): void;
      stencilOpSeparate(face: number, fail: number, zfail: number, zpass: number): void;
      texImage2D(target: number, level: number, internalformat: number, width: number, height: number, border: number, format: number, type: number, pixels: ArrayBufferView | null): void;
      // texImage2D(target: number, level: number, internalformat: number, format: number, type: number, pixels: ImageBitmap | ImageData | HTMLVideoElement | HTMLImageElement | HTMLCanvasElement): void;
      texParameterf(target: number, pname: number, param: number): void;
      texParameteri(target: number, pname: number, param: number): void;
      texSubImage2D(target: number, level: number, xoffset: number, yoffset: number, width: number, height: number, format: number, type: number, pixels: ArrayBufferView | null): void;
      // texSubImage2D(target: number, level: number, xoffset: number, yoffset: number, format: number, type: number, pixels: ImageBitmap | ImageData | HTMLVideoElement | HTMLImageElement | HTMLCanvasElement): void;
      uniform1f(location: WebGLUniformLocation | null, x: number): void;
      uniform1fv(location: WebGLUniformLocation, v: Float32Array | number[]): void;
      uniform1i(location: WebGLUniformLocation | null, x: number): void;
      uniform1iv(location: WebGLUniformLocation, v: Int32Array | number[]): void;
      uniform2f(location: WebGLUniformLocation | null, x: number, y: number): void;
      uniform2fv(location: WebGLUniformLocation, v: Float32Array | number[]): void;
      uniform2i(location: WebGLUniformLocation | null, x: number, y: number): void;
      uniform2iv(location: WebGLUniformLocation, v: Int32Array | number[]): void;
      uniform3f(location: WebGLUniformLocation | null, x: number, y: number, z: number): void;
      uniform3fv(location: WebGLUniformLocation, v: Float32Array | number[]): void;
      uniform3i(location: WebGLUniformLocation | null, x: number, y: number, z: number): void;
      uniform3iv(location: WebGLUniformLocation, v: Int32Array | number[]): void;
      uniform4f(location: WebGLUniformLocation | null, x: number, y: number, z: number, w: number): void;
      uniform4fv(location: WebGLUniformLocation, v: Float32Array | number[]): void;
      uniform4i(location: WebGLUniformLocation | null, x: number, y: number, z: number, w: number): void;
      uniform4iv(location: WebGLUniformLocation, v: Int32Array | number[]): void;
      uniformMatrix2fv(location: WebGLUniformLocation, transpose: boolean, value: Float32Array | number[]): void;
      uniformMatrix3fv(location: WebGLUniformLocation, transpose: boolean, value: Float32Array | number[]): void;
      uniformMatrix4fv(location: WebGLUniformLocation, transpose: boolean, value: Float32Array | number[]): void;
      useProgram(program: WebGLProgram | null): void;
      validateProgram(program: WebGLProgram | null): void;
      vertexAttrib1f(indx: number, x: number): void;
      vertexAttrib1fv(indx: number, values: Float32Array | number[]): void;
      vertexAttrib2f(indx: number, x: number, y: number): void;
      vertexAttrib2fv(indx: number, values: Float32Array | number[]): void;
      vertexAttrib3f(indx: number, x: number, y: number, z: number): void;
      vertexAttrib3fv(indx: number, values: Float32Array | number[]): void;
      vertexAttrib4f(indx: number, x: number, y: number, z: number, w: number): void;
      vertexAttrib4fv(indx: number, values: Float32Array | number[]): void;
      vertexAttribPointer(indx: number, size: number, type: number, normalized: boolean, stride: number, offset: number): void;
      viewport(x: number, y: number, width: number, height: number): void;
      readonly ACTIVE_ATTRIBUTES: number;
      readonly ACTIVE_TEXTURE: number;
      readonly ACTIVE_UNIFORMS: number;
      readonly ALIASED_LINE_WIDTH_RANGE: number;
      readonly ALIASED_POINT_SIZE_RANGE: number;
      readonly ALPHA: number;
      readonly ALPHA_BITS: number;
      readonly ALWAYS: number;
      readonly ARRAY_BUFFER: number;
      readonly ARRAY_BUFFER_BINDING: number;
      readonly ATTACHED_SHADERS: number;
      readonly BACK: number;
      readonly BLEND: number;
      readonly BLEND_COLOR: number;
      readonly BLEND_DST_ALPHA: number;
      readonly BLEND_DST_RGB: number;
      readonly BLEND_EQUATION: number;
      readonly BLEND_EQUATION_ALPHA: number;
      readonly BLEND_EQUATION_RGB: number;
      readonly BLEND_SRC_ALPHA: number;
      readonly BLEND_SRC_RGB: number;
      readonly BLUE_BITS: number;
      readonly BOOL: number;
      readonly BOOL_VEC2: number;
      readonly BOOL_VEC3: number;
      readonly BOOL_VEC4: number;
      readonly BROWSER_DEFAULT_WEBGL: number;
      readonly BUFFER_SIZE: number;
      readonly BUFFER_USAGE: number;
      readonly BYTE: number;
      readonly CCW: number;
      readonly CLAMP_TO_EDGE: number;
      readonly COLOR_ATTACHMENT0: number;
      readonly COLOR_BUFFER_BIT: number;
      readonly COLOR_CLEAR_VALUE: number;
      readonly COLOR_WRITEMASK: number;
      readonly COMPILE_STATUS: number;
      readonly COMPRESSED_TEXTURE_FORMATS: number;
      readonly CONSTANT_ALPHA: number;
      readonly CONSTANT_COLOR: number;
      readonly CONTEXT_LOST_WEBGL: number;
      readonly CULL_FACE: number;
      readonly CULL_FACE_MODE: number;
      readonly CURRENT_PROGRAM: number;
      readonly CURRENT_VERTEX_ATTRIB: number;
      readonly CW: number;
      readonly DECR: number;
      readonly DECR_WRAP: number;
      readonly DELETE_STATUS: number;
      readonly DEPTH_ATTACHMENT: number;
      readonly DEPTH_BITS: number;
      readonly DEPTH_BUFFER_BIT: number;
      readonly DEPTH_CLEAR_VALUE: number;
      readonly DEPTH_COMPONENT: number;
      readonly DEPTH_COMPONENT16: number;
      readonly DEPTH_FUNC: number;
      readonly DEPTH_RANGE: number;
      readonly DEPTH_STENCIL: number;
      readonly DEPTH_STENCIL_ATTACHMENT: number;
      readonly DEPTH_TEST: number;
      readonly DEPTH_WRITEMASK: number;
      readonly DITHER: number;
      readonly DONT_CARE: number;
      readonly DST_ALPHA: number;
      readonly DST_COLOR: number;
      readonly DYNAMIC_DRAW: number;
      readonly ELEMENT_ARRAY_BUFFER: number;
      readonly ELEMENT_ARRAY_BUFFER_BINDING: number;
      readonly EQUAL: number;
      readonly FASTEST: number;
      readonly FLOAT: number;
      readonly FLOAT_MAT2: number;
      readonly FLOAT_MAT3: number;
      readonly FLOAT_MAT4: number;
      readonly FLOAT_VEC2: number;
      readonly FLOAT_VEC3: number;
      readonly FLOAT_VEC4: number;
      readonly FRAGMENT_SHADER: number;
      readonly FRAMEBUFFER: number;
      readonly FRAMEBUFFER_ATTACHMENT_OBJECT_NAME: number;
      readonly FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE: number;
      readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE: number;
      readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL: number;
      readonly FRAMEBUFFER_BINDING: number;
      readonly FRAMEBUFFER_COMPLETE: number;
      readonly FRAMEBUFFER_INCOMPLETE_ATTACHMENT: number;
      readonly FRAMEBUFFER_INCOMPLETE_DIMENSIONS: number;
      readonly FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: number;
      readonly FRAMEBUFFER_UNSUPPORTED: number;
      readonly FRONT: number;
      readonly FRONT_AND_BACK: number;
      readonly FRONT_FACE: number;
      readonly FUNC_ADD: number;
      readonly FUNC_REVERSE_SUBTRACT: number;
      readonly FUNC_SUBTRACT: number;
      readonly GENERATE_MIPMAP_HINT: number;
      readonly GEQUAL: number;
      readonly GREATER: number;
      readonly GREEN_BITS: number;
      readonly HIGH_FLOAT: number;
      readonly HIGH_INT: number;
      readonly IMPLEMENTATION_COLOR_READ_FORMAT: number;
      readonly IMPLEMENTATION_COLOR_READ_TYPE: number;
      readonly INCR: number;
      readonly INCR_WRAP: number;
      readonly INT: number;
      readonly INT_VEC2: number;
      readonly INT_VEC3: number;
      readonly INT_VEC4: number;
      readonly INVALID_ENUM: number;
      readonly INVALID_FRAMEBUFFER_OPERATION: number;
      readonly INVALID_OPERATION: number;
      readonly INVALID_VALUE: number;
      readonly INVERT: number;
      readonly KEEP: number;
      readonly LEQUAL: number;
      readonly LESS: number;
      readonly LINE_LOOP: number;
      readonly LINE_STRIP: number;
      readonly LINE_WIDTH: number;
      readonly LINEAR: number;
      readonly LINEAR_MIPMAP_LINEAR: number;
      readonly LINEAR_MIPMAP_NEAREST: number;
      readonly LINES: number;
      readonly LINK_STATUS: number;
      readonly LOW_FLOAT: number;
      readonly LOW_INT: number;
      readonly LUMINANCE: number;
      readonly LUMINANCE_ALPHA: number;
      readonly MAX_COMBINED_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_CUBE_MAP_TEXTURE_SIZE: number;
      readonly MAX_FRAGMENT_UNIFORM_VECTORS: number;
      readonly MAX_RENDERBUFFER_SIZE: number;
      readonly MAX_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_TEXTURE_SIZE: number;
      readonly MAX_VARYING_VECTORS: number;
      readonly MAX_VERTEX_ATTRIBS: number;
      readonly MAX_VERTEX_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_VERTEX_UNIFORM_VECTORS: number;
      readonly MAX_VIEWPORT_DIMS: number;
      readonly MEDIUM_FLOAT: number;
      readonly MEDIUM_INT: number;
      readonly MIRRORED_REPEAT: number;
      readonly NEAREST: number;
      readonly NEAREST_MIPMAP_LINEAR: number;
      readonly NEAREST_MIPMAP_NEAREST: number;
      readonly NEVER: number;
      readonly NICEST: number;
      readonly NO_ERROR: number;
      readonly NONE: number;
      readonly NOTEQUAL: number;
      readonly ONE: number;
      readonly ONE_MINUS_CONSTANT_ALPHA: number;
      readonly ONE_MINUS_CONSTANT_COLOR: number;
      readonly ONE_MINUS_DST_ALPHA: number;
      readonly ONE_MINUS_DST_COLOR: number;
      readonly ONE_MINUS_SRC_ALPHA: number;
      readonly ONE_MINUS_SRC_COLOR: number;
      readonly OUT_OF_MEMORY: number;
      readonly PACK_ALIGNMENT: number;
      readonly POINTS: number;
      readonly POLYGON_OFFSET_FACTOR: number;
      readonly POLYGON_OFFSET_FILL: number;
      readonly POLYGON_OFFSET_UNITS: number;
      readonly RED_BITS: number;
      readonly RENDERBUFFER: number;
      readonly RENDERBUFFER_ALPHA_SIZE: number;
      readonly RENDERBUFFER_BINDING: number;
      readonly RENDERBUFFER_BLUE_SIZE: number;
      readonly RENDERBUFFER_DEPTH_SIZE: number;
      readonly RENDERBUFFER_GREEN_SIZE: number;
      readonly RENDERBUFFER_HEIGHT: number;
      readonly RENDERBUFFER_INTERNAL_FORMAT: number;
      readonly RENDERBUFFER_RED_SIZE: number;
      readonly RENDERBUFFER_STENCIL_SIZE: number;
      readonly RENDERBUFFER_WIDTH: number;
      readonly RENDERER: number;
      readonly REPEAT: number;
      readonly REPLACE: number;
      readonly RGB: number;
      readonly RGB5_A1: number;
      readonly RGB565: number;
      readonly RGBA: number;
      readonly RGBA4: number;
      readonly SAMPLE_ALPHA_TO_COVERAGE: number;
      readonly SAMPLE_BUFFERS: number;
      readonly SAMPLE_COVERAGE: number;
      readonly SAMPLE_COVERAGE_INVERT: number;
      readonly SAMPLE_COVERAGE_VALUE: number;
      readonly SAMPLER_2D: number;
      readonly SAMPLER_CUBE: number;
      readonly SAMPLES: number;
      readonly SCISSOR_BOX: number;
      readonly SCISSOR_TEST: number;
      readonly SHADER_TYPE: number;
      readonly SHADING_LANGUAGE_VERSION: number;
      readonly SHORT: number;
      readonly SRC_ALPHA: number;
      readonly SRC_ALPHA_SATURATE: number;
      readonly SRC_COLOR: number;
      readonly STATIC_DRAW: number;
      readonly STENCIL_ATTACHMENT: number;
      readonly STENCIL_BACK_FAIL: number;
      readonly STENCIL_BACK_FUNC: number;
      readonly STENCIL_BACK_PASS_DEPTH_FAIL: number;
      readonly STENCIL_BACK_PASS_DEPTH_PASS: number;
      readonly STENCIL_BACK_REF: number;
      readonly STENCIL_BACK_VALUE_MASK: number;
      readonly STENCIL_BACK_WRITEMASK: number;
      readonly STENCIL_BITS: number;
      readonly STENCIL_BUFFER_BIT: number;
      readonly STENCIL_CLEAR_VALUE: number;
      readonly STENCIL_FAIL: number;
      readonly STENCIL_FUNC: number;
      readonly STENCIL_INDEX: number;
      readonly STENCIL_INDEX8: number;
      readonly STENCIL_PASS_DEPTH_FAIL: number;
      readonly STENCIL_PASS_DEPTH_PASS: number;
      readonly STENCIL_REF: number;
      readonly STENCIL_TEST: number;
      readonly STENCIL_VALUE_MASK: number;
      readonly STENCIL_WRITEMASK: number;
      readonly STREAM_DRAW: number;
      readonly SUBPIXEL_BITS: number;
      readonly TEXTURE: number;
      readonly TEXTURE_2D: number;
      readonly TEXTURE_BINDING_2D: number;
      readonly TEXTURE_BINDING_CUBE_MAP: number;
      readonly TEXTURE_CUBE_MAP: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_X: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_Y: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_Z: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_X: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_Y: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_Z: number;
      readonly TEXTURE_MAG_FILTER: number;
      readonly TEXTURE_MIN_FILTER: number;
      readonly TEXTURE_WRAP_S: number;
      readonly TEXTURE_WRAP_T: number;
      readonly TEXTURE0: number;
      readonly TEXTURE1: number;
      readonly TEXTURE10: number;
      readonly TEXTURE11: number;
      readonly TEXTURE12: number;
      readonly TEXTURE13: number;
      readonly TEXTURE14: number;
      readonly TEXTURE15: number;
      readonly TEXTURE16: number;
      readonly TEXTURE17: number;
      readonly TEXTURE18: number;
      readonly TEXTURE19: number;
      readonly TEXTURE2: number;
      readonly TEXTURE20: number;
      readonly TEXTURE21: number;
      readonly TEXTURE22: number;
      readonly TEXTURE23: number;
      readonly TEXTURE24: number;
      readonly TEXTURE25: number;
      readonly TEXTURE26: number;
      readonly TEXTURE27: number;
      readonly TEXTURE28: number;
      readonly TEXTURE29: number;
      readonly TEXTURE3: number;
      readonly TEXTURE30: number;
      readonly TEXTURE31: number;
      readonly TEXTURE4: number;
      readonly TEXTURE5: number;
      readonly TEXTURE6: number;
      readonly TEXTURE7: number;
      readonly TEXTURE8: number;
      readonly TEXTURE9: number;
      readonly TRIANGLE_FAN: number;
      readonly TRIANGLE_STRIP: number;
      readonly TRIANGLES: number;
      readonly UNPACK_ALIGNMENT: number;
      readonly UNPACK_COLORSPACE_CONVERSION_WEBGL: number;
      readonly UNPACK_FLIP_Y_WEBGL: number;
      readonly UNPACK_PREMULTIPLY_ALPHA_WEBGL: number;
      readonly UNSIGNED_BYTE: number;
      readonly UNSIGNED_INT: number;
      readonly UNSIGNED_SHORT: number;
      readonly UNSIGNED_SHORT_4_4_4_4: number;
      readonly UNSIGNED_SHORT_5_5_5_1: number;
      readonly UNSIGNED_SHORT_5_6_5: number;
      readonly VALIDATE_STATUS: number;
      readonly VENDOR: number;
      readonly VERSION: number;
      readonly VERTEX_ATTRIB_ARRAY_BUFFER_BINDING: number;
      readonly VERTEX_ATTRIB_ARRAY_ENABLED: number;
      readonly VERTEX_ATTRIB_ARRAY_NORMALIZED: number;
      readonly VERTEX_ATTRIB_ARRAY_POINTER: number;
      readonly VERTEX_ATTRIB_ARRAY_SIZE: number;
      readonly VERTEX_ATTRIB_ARRAY_STRIDE: number;
      readonly VERTEX_ATTRIB_ARRAY_TYPE: number;
      readonly VERTEX_SHADER: number;
      readonly VIEWPORT: number;
      readonly ZERO: number;
  }

  declare var WebGLRenderingContext: {
      prototype: WebGLRenderingContext;
      new(): WebGLRenderingContext;
      readonly ACTIVE_ATTRIBUTES: number;
      readonly ACTIVE_TEXTURE: number;
      readonly ACTIVE_UNIFORMS: number;
      readonly ALIASED_LINE_WIDTH_RANGE: number;
      readonly ALIASED_POINT_SIZE_RANGE: number;
      readonly ALPHA: number;
      readonly ALPHA_BITS: number;
      readonly ALWAYS: number;
      readonly ARRAY_BUFFER: number;
      readonly ARRAY_BUFFER_BINDING: number;
      readonly ATTACHED_SHADERS: number;
      readonly BACK: number;
      readonly BLEND: number;
      readonly BLEND_COLOR: number;
      readonly BLEND_DST_ALPHA: number;
      readonly BLEND_DST_RGB: number;
      readonly BLEND_EQUATION: number;
      readonly BLEND_EQUATION_ALPHA: number;
      readonly BLEND_EQUATION_RGB: number;
      readonly BLEND_SRC_ALPHA: number;
      readonly BLEND_SRC_RGB: number;
      readonly BLUE_BITS: number;
      readonly BOOL: number;
      readonly BOOL_VEC2: number;
      readonly BOOL_VEC3: number;
      readonly BOOL_VEC4: number;
      readonly BROWSER_DEFAULT_WEBGL: number;
      readonly BUFFER_SIZE: number;
      readonly BUFFER_USAGE: number;
      readonly BYTE: number;
      readonly CCW: number;
      readonly CLAMP_TO_EDGE: number;
      readonly COLOR_ATTACHMENT0: number;
      readonly COLOR_BUFFER_BIT: number;
      readonly COLOR_CLEAR_VALUE: number;
      readonly COLOR_WRITEMASK: number;
      readonly COMPILE_STATUS: number;
      readonly COMPRESSED_TEXTURE_FORMATS: number;
      readonly CONSTANT_ALPHA: number;
      readonly CONSTANT_COLOR: number;
      readonly CONTEXT_LOST_WEBGL: number;
      readonly CULL_FACE: number;
      readonly CULL_FACE_MODE: number;
      readonly CURRENT_PROGRAM: number;
      readonly CURRENT_VERTEX_ATTRIB: number;
      readonly CW: number;
      readonly DECR: number;
      readonly DECR_WRAP: number;
      readonly DELETE_STATUS: number;
      readonly DEPTH_ATTACHMENT: number;
      readonly DEPTH_BITS: number;
      readonly DEPTH_BUFFER_BIT: number;
      readonly DEPTH_CLEAR_VALUE: number;
      readonly DEPTH_COMPONENT: number;
      readonly DEPTH_COMPONENT16: number;
      readonly DEPTH_FUNC: number;
      readonly DEPTH_RANGE: number;
      readonly DEPTH_STENCIL: number;
      readonly DEPTH_STENCIL_ATTACHMENT: number;
      readonly DEPTH_TEST: number;
      readonly DEPTH_WRITEMASK: number;
      readonly DITHER: number;
      readonly DONT_CARE: number;
      readonly DST_ALPHA: number;
      readonly DST_COLOR: number;
      readonly DYNAMIC_DRAW: number;
      readonly ELEMENT_ARRAY_BUFFER: number;
      readonly ELEMENT_ARRAY_BUFFER_BINDING: number;
      readonly EQUAL: number;
      readonly FASTEST: number;
      readonly FLOAT: number;
      readonly FLOAT_MAT2: number;
      readonly FLOAT_MAT3: number;
      readonly FLOAT_MAT4: number;
      readonly FLOAT_VEC2: number;
      readonly FLOAT_VEC3: number;
      readonly FLOAT_VEC4: number;
      readonly FRAGMENT_SHADER: number;
      readonly FRAMEBUFFER: number;
      readonly FRAMEBUFFER_ATTACHMENT_OBJECT_NAME: number;
      readonly FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE: number;
      readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_CUBE_MAP_FACE: number;
      readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_LEVEL: number;
      readonly FRAMEBUFFER_BINDING: number;
      readonly FRAMEBUFFER_COMPLETE: number;
      readonly FRAMEBUFFER_INCOMPLETE_ATTACHMENT: number;
      readonly FRAMEBUFFER_INCOMPLETE_DIMENSIONS: number;
      readonly FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: number;
      readonly FRAMEBUFFER_UNSUPPORTED: number;
      readonly FRONT: number;
      readonly FRONT_AND_BACK: number;
      readonly FRONT_FACE: number;
      readonly FUNC_ADD: number;
      readonly FUNC_REVERSE_SUBTRACT: number;
      readonly FUNC_SUBTRACT: number;
      readonly GENERATE_MIPMAP_HINT: number;
      readonly GEQUAL: number;
      readonly GREATER: number;
      readonly GREEN_BITS: number;
      readonly HIGH_FLOAT: number;
      readonly HIGH_INT: number;
      readonly IMPLEMENTATION_COLOR_READ_FORMAT: number;
      readonly IMPLEMENTATION_COLOR_READ_TYPE: number;
      readonly INCR: number;
      readonly INCR_WRAP: number;
      readonly INT: number;
      readonly INT_VEC2: number;
      readonly INT_VEC3: number;
      readonly INT_VEC4: number;
      readonly INVALID_ENUM: number;
      readonly INVALID_FRAMEBUFFER_OPERATION: number;
      readonly INVALID_OPERATION: number;
      readonly INVALID_VALUE: number;
      readonly INVERT: number;
      readonly KEEP: number;
      readonly LEQUAL: number;
      readonly LESS: number;
      readonly LINE_LOOP: number;
      readonly LINE_STRIP: number;
      readonly LINE_WIDTH: number;
      readonly LINEAR: number;
      readonly LINEAR_MIPMAP_LINEAR: number;
      readonly LINEAR_MIPMAP_NEAREST: number;
      readonly LINES: number;
      readonly LINK_STATUS: number;
      readonly LOW_FLOAT: number;
      readonly LOW_INT: number;
      readonly LUMINANCE: number;
      readonly LUMINANCE_ALPHA: number;
      readonly MAX_COMBINED_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_CUBE_MAP_TEXTURE_SIZE: number;
      readonly MAX_FRAGMENT_UNIFORM_VECTORS: number;
      readonly MAX_RENDERBUFFER_SIZE: number;
      readonly MAX_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_TEXTURE_SIZE: number;
      readonly MAX_VARYING_VECTORS: number;
      readonly MAX_VERTEX_ATTRIBS: number;
      readonly MAX_VERTEX_TEXTURE_IMAGE_UNITS: number;
      readonly MAX_VERTEX_UNIFORM_VECTORS: number;
      readonly MAX_VIEWPORT_DIMS: number;
      readonly MEDIUM_FLOAT: number;
      readonly MEDIUM_INT: number;
      readonly MIRRORED_REPEAT: number;
      readonly NEAREST: number;
      readonly NEAREST_MIPMAP_LINEAR: number;
      readonly NEAREST_MIPMAP_NEAREST: number;
      readonly NEVER: number;
      readonly NICEST: number;
      readonly NO_ERROR: number;
      readonly NONE: number;
      readonly NOTEQUAL: number;
      readonly ONE: number;
      readonly ONE_MINUS_CONSTANT_ALPHA: number;
      readonly ONE_MINUS_CONSTANT_COLOR: number;
      readonly ONE_MINUS_DST_ALPHA: number;
      readonly ONE_MINUS_DST_COLOR: number;
      readonly ONE_MINUS_SRC_ALPHA: number;
      readonly ONE_MINUS_SRC_COLOR: number;
      readonly OUT_OF_MEMORY: number;
      readonly PACK_ALIGNMENT: number;
      readonly POINTS: number;
      readonly POLYGON_OFFSET_FACTOR: number;
      readonly POLYGON_OFFSET_FILL: number;
      readonly POLYGON_OFFSET_UNITS: number;
      readonly RED_BITS: number;
      readonly RENDERBUFFER: number;
      readonly RENDERBUFFER_ALPHA_SIZE: number;
      readonly RENDERBUFFER_BINDING: number;
      readonly RENDERBUFFER_BLUE_SIZE: number;
      readonly RENDERBUFFER_DEPTH_SIZE: number;
      readonly RENDERBUFFER_GREEN_SIZE: number;
      readonly RENDERBUFFER_HEIGHT: number;
      readonly RENDERBUFFER_INTERNAL_FORMAT: number;
      readonly RENDERBUFFER_RED_SIZE: number;
      readonly RENDERBUFFER_STENCIL_SIZE: number;
      readonly RENDERBUFFER_WIDTH: number;
      readonly RENDERER: number;
      readonly REPEAT: number;
      readonly REPLACE: number;
      readonly RGB: number;
      readonly RGB5_A1: number;
      readonly RGB565: number;
      readonly RGBA: number;
      readonly RGBA4: number;
      readonly SAMPLE_ALPHA_TO_COVERAGE: number;
      readonly SAMPLE_BUFFERS: number;
      readonly SAMPLE_COVERAGE: number;
      readonly SAMPLE_COVERAGE_INVERT: number;
      readonly SAMPLE_COVERAGE_VALUE: number;
      readonly SAMPLER_2D: number;
      readonly SAMPLER_CUBE: number;
      readonly SAMPLES: number;
      readonly SCISSOR_BOX: number;
      readonly SCISSOR_TEST: number;
      readonly SHADER_TYPE: number;
      readonly SHADING_LANGUAGE_VERSION: number;
      readonly SHORT: number;
      readonly SRC_ALPHA: number;
      readonly SRC_ALPHA_SATURATE: number;
      readonly SRC_COLOR: number;
      readonly STATIC_DRAW: number;
      readonly STENCIL_ATTACHMENT: number;
      readonly STENCIL_BACK_FAIL: number;
      readonly STENCIL_BACK_FUNC: number;
      readonly STENCIL_BACK_PASS_DEPTH_FAIL: number;
      readonly STENCIL_BACK_PASS_DEPTH_PASS: number;
      readonly STENCIL_BACK_REF: number;
      readonly STENCIL_BACK_VALUE_MASK: number;
      readonly STENCIL_BACK_WRITEMASK: number;
      readonly STENCIL_BITS: number;
      readonly STENCIL_BUFFER_BIT: number;
      readonly STENCIL_CLEAR_VALUE: number;
      readonly STENCIL_FAIL: number;
      readonly STENCIL_FUNC: number;
      readonly STENCIL_INDEX: number;
      readonly STENCIL_INDEX8: number;
      readonly STENCIL_PASS_DEPTH_FAIL: number;
      readonly STENCIL_PASS_DEPTH_PASS: number;
      readonly STENCIL_REF: number;
      readonly STENCIL_TEST: number;
      readonly STENCIL_VALUE_MASK: number;
      readonly STENCIL_WRITEMASK: number;
      readonly STREAM_DRAW: number;
      readonly SUBPIXEL_BITS: number;
      readonly TEXTURE: number;
      readonly TEXTURE_2D: number;
      readonly TEXTURE_BINDING_2D: number;
      readonly TEXTURE_BINDING_CUBE_MAP: number;
      readonly TEXTURE_CUBE_MAP: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_X: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_Y: number;
      readonly TEXTURE_CUBE_MAP_NEGATIVE_Z: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_X: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_Y: number;
      readonly TEXTURE_CUBE_MAP_POSITIVE_Z: number;
      readonly TEXTURE_MAG_FILTER: number;
      readonly TEXTURE_MIN_FILTER: number;
      readonly TEXTURE_WRAP_S: number;
      readonly TEXTURE_WRAP_T: number;
      readonly TEXTURE0: number;
      readonly TEXTURE1: number;
      readonly TEXTURE10: number;
      readonly TEXTURE11: number;
      readonly TEXTURE12: number;
      readonly TEXTURE13: number;
      readonly TEXTURE14: number;
      readonly TEXTURE15: number;
      readonly TEXTURE16: number;
      readonly TEXTURE17: number;
      readonly TEXTURE18: number;
      readonly TEXTURE19: number;
      readonly TEXTURE2: number;
      readonly TEXTURE20: number;
      readonly TEXTURE21: number;
      readonly TEXTURE22: number;
      readonly TEXTURE23: number;
      readonly TEXTURE24: number;
      readonly TEXTURE25: number;
      readonly TEXTURE26: number;
      readonly TEXTURE27: number;
      readonly TEXTURE28: number;
      readonly TEXTURE29: number;
      readonly TEXTURE3: number;
      readonly TEXTURE30: number;
      readonly TEXTURE31: number;
      readonly TEXTURE4: number;
      readonly TEXTURE5: number;
      readonly TEXTURE6: number;
      readonly TEXTURE7: number;
      readonly TEXTURE8: number;
      readonly TEXTURE9: number;
      readonly TRIANGLE_FAN: number;
      readonly TRIANGLE_STRIP: number;
      readonly TRIANGLES: number;
      readonly UNPACK_ALIGNMENT: number;
      readonly UNPACK_COLORSPACE_CONVERSION_WEBGL: number;
      readonly UNPACK_FLIP_Y_WEBGL: number;
      readonly UNPACK_PREMULTIPLY_ALPHA_WEBGL: number;
      readonly UNSIGNED_BYTE: number;
      readonly UNSIGNED_INT: number;
      readonly UNSIGNED_SHORT: number;
      readonly UNSIGNED_SHORT_4_4_4_4: number;
      readonly UNSIGNED_SHORT_5_5_5_1: number;
      readonly UNSIGNED_SHORT_5_6_5: number;
      readonly VALIDATE_STATUS: number;
      readonly VENDOR: number;
      readonly VERSION: number;
      readonly VERTEX_ATTRIB_ARRAY_BUFFER_BINDING: number;
      readonly VERTEX_ATTRIB_ARRAY_ENABLED: number;
      readonly VERTEX_ATTRIB_ARRAY_NORMALIZED: number;
      readonly VERTEX_ATTRIB_ARRAY_POINTER: number;
      readonly VERTEX_ATTRIB_ARRAY_SIZE: number;
      readonly VERTEX_ATTRIB_ARRAY_STRIDE: number;
      readonly VERTEX_ATTRIB_ARRAY_TYPE: number;
      readonly VERTEX_SHADER: number;
      readonly VIEWPORT: number;
      readonly ZERO: number;
  };

  interface WebGLShader extends WebGLObject {
  }

  declare var WebGLShader: {
      prototype: WebGLShader;
      new(): WebGLShader;
  };

  interface WebGLShaderPrecisionFormat {
      readonly precision: number;
      readonly rangeMax: number;
      readonly rangeMin: number;
  }

  declare var WebGLShaderPrecisionFormat: {
      prototype: WebGLShaderPrecisionFormat;
      new(): WebGLShaderPrecisionFormat;
  };

  interface WebGLTexture extends WebGLObject {
  }

  declare var WebGLTexture: {
      prototype: WebGLTexture;
      new(): WebGLTexture;
  };

  interface WebGLUniformLocation {
  }

  declare var WebGLUniformLocation: {
      prototype: WebGLUniformLocation;
      new(): WebGLUniformLocation;
  };