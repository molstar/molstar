/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import {
  bindForAttachment,
  isAttached,
  setAttached,
  setAttachmentStorage,
  setFramebufferHeight,
} from './framebuffer.js';
import type { AttachmentRef, Framebuffer } from './framebuffer.js';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from '../renderbuffer.js';
import { endRenderPass } from './render-pass.js';
import type { WebGLState } from './state.js';

const getNextRenderbufferId = idFactory();

export function getFormat(gl: GLRenderingContext, format: RenderbufferFormat) {
  switch (format) {
    case 'depth16':
      return gl.DEPTH_COMPONENT16;
    case 'stencil8':
      return gl.STENCIL_INDEX8;
    case 'rgba4':
      return gl.RGBA4;
    case 'depth-stencil':
      return gl.DEPTH_STENCIL;
    case 'depth24':
      if (isWebGL2(gl)) return gl.DEPTH_COMPONENT24;
      else throw new Error('WebGL2 needed for `depth24` renderbuffer format');
    case 'depth32f':
      if (isWebGL2(gl)) return gl.DEPTH_COMPONENT32F;
      else throw new Error('WebGL2 needed for `depth32f` renderbuffer format');
    case 'depth24-stencil8':
      if (isWebGL2(gl)) return gl.DEPTH24_STENCIL8;
      else throw new Error('WebGL2 needed for `depth24-stencil8` renderbuffer format');
    case 'depth32f-stencil8':
      if (isWebGL2(gl)) return gl.DEPTH32F_STENCIL8;
      else throw new Error('WebGL2 needed for `depth32f-stencil8` renderbuffer format');
  }
}

function getFormatSize(format: RenderbufferFormat): number {
  switch (format) {
    case 'depth16':
      return 4;
    case 'stencil8':
      return 2;
    case 'rgba4':
      return 4;
    case 'depth-stencil':
      return 4;
    case 'depth24':
      return 3;
    case 'depth32f':
      return 4;
    case 'depth24-stencil8':
      return 4;
    case 'depth32f-stencil8':
      return 5;
  }
}

export function getAttachment(gl: GLRenderingContext, attachment: RenderbufferAttachment) {
  switch (attachment) {
    case 'depth':
      return gl.DEPTH_ATTACHMENT;
    case 'stencil':
      return gl.STENCIL_ATTACHMENT;
    case 'depth-stencil':
      return gl.DEPTH_STENCIL_ATTACHMENT;
    case 'color0':
      return gl.COLOR_ATTACHMENT0;
  }
}

function getRenderbuffer(gl: GLRenderingContext) {
  const renderbuffer = gl.createRenderbuffer();
  if (renderbuffer === null) {
    throw new Error('Could not create WebGL renderbuffer');
  }
  return renderbuffer;
}

/** depth formats are precision preferences, WebGL1 falls back to `depth16` and `depth-stencil` */
function getSupportedFormat(gl: GLRenderingContext, format: RenderbufferFormat): RenderbufferFormat {
  if (isWebGL2(gl)) return format;
  switch (format) {
    case 'depth24':
    case 'depth32f':
      return 'depth16';
    case 'depth24-stencil8':
    case 'depth32f-stencil8':
      return 'depth-stencil';
  }
  return format;
}

/** A renderbuffer of the WebGL backend, attached to its framebuffers internally */
export interface GLRenderbuffer extends Renderbuffer {
  attachFramebuffer: (framebuffer: Framebuffer) => void;
  detachFramebuffer: (framebuffer: Framebuffer) => void;
}

export function asGLRenderbuffer(renderbuffer: Renderbuffer): GLRenderbuffer {
  if (typeof (renderbuffer as Partial<GLRenderbuffer>).attachFramebuffer !== 'function') {
    throw new Error('expected a WebGL renderbuffer');
  }
  return renderbuffer as GLRenderbuffer;
}

export function createRenderbuffer(
  gl: GLRenderingContext,
  format: RenderbufferFormat,
  attachment: RenderbufferAttachment,
  _width: number,
  _height: number,
  state: WebGLState,
): GLRenderbuffer {
  const _renderbuffer = getRenderbuffer(gl);
  format = getSupportedFormat(gl, format);

  const bind = () => gl.bindRenderbuffer(gl.RENDERBUFFER, _renderbuffer);
  const _format = getFormat(gl, format);
  const _attachment = getAttachment(gl, attachment);

  /** the size of the storage, the defined size once allocated, see `ensureStorage` */
  let storedWidth = 0;
  let storedHeight = 0;
  function specify(w: number, h: number) {
    bind();
    gl.renderbufferStorage(gl.RENDERBUFFER, _format, w, h);
    (storedWidth = w), (storedHeight = h);
  }
  // 1×1 until a framebuffer it is attached to is used (lazy: unused attachments take no memory)
  specify(1, 1);

  let destroyed = false;

  function ensureStorage() {
    if (destroyed) return;
    if (storedWidth !== _width || storedHeight !== _height) specify(_width, _height);
  }

  const id = getNextRenderbufferId();
  /** the height of the attachments of framebuffers it is attached to (`setFramebufferHeight`), made once */
  const getHeight = () => _height;
  /** what it is as an attachment (`Framebuffer.getAttachments`), made once */
  const ref: AttachmentRef = { kind: 'renderbuffer', id };
  return {
    id,

    getByteCount: () => getFormatSize(format) * storedWidth * storedHeight,

    bind,
    attachFramebuffer: (framebuffer: Framebuffer) => {
      // attached already (passes attach their targets every frame): nothing changes
      if (isAttached(framebuffer, _attachment, _renderbuffer)) {
        setFramebufferHeight(framebuffer, getHeight);
        return;
      }
      endRenderPass(state);
      bindForAttachment(gl, framebuffer);
      setAttachmentStorage(gl, framebuffer, _attachment, ref, ensureStorage);
      setAttached(gl, framebuffer, _attachment, _renderbuffer);
      bind();
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, _attachment, gl.RENDERBUFFER, _renderbuffer);
      setFramebufferHeight(framebuffer, getHeight);
    },
    detachFramebuffer: (framebuffer: Framebuffer) => {
      endRenderPass(state);
      bindForAttachment(gl, framebuffer);
      setAttachmentStorage(gl, framebuffer, _attachment, undefined);
      setAttached(gl, framebuffer, _attachment, undefined);
      bind();
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, _attachment, gl.RENDERBUFFER, null);
    },
    releaseStorage: () => {
      if (!destroyed && (storedWidth !== 1 || storedHeight !== 1)) specify(1, 1);
    },
    setSize: (width: number, height: number) => {
      if (width === _width && height === _height) return;
      _width = width;
      _height = height;
      // storage of another size released, allocated again when used
      if (storedWidth !== 1 || storedHeight !== 1) specify(1, 1);
    },
    destroy: () => {
      if (destroyed) return;
      gl.deleteRenderbuffer(_renderbuffer);
      destroyed = true;
    },
  };
}
