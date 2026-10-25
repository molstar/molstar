/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { AttachmentRef, Framebuffer } from '../framebuffer.js';
import type { WebGLExtensions } from './extensions.js';
import { endRenderPass } from './render-pass.js';
import type { WebGLState } from './state.js';

const getNextFramebufferId = idFactory();

function getFramebufferStatusDescription(gl: GLRenderingContext, status: number) {
  switch (status) {
    case gl.FRAMEBUFFER_COMPLETE:
      return 'complete';
    case gl.FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
      return 'incomplete attachment';
    case gl.FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
      return 'incomplete missing attachment';
    case gl.FRAMEBUFFER_INCOMPLETE_DIMENSIONS:
      return 'incomplete dimensions';
    case gl.FRAMEBUFFER_UNSUPPORTED:
      return 'unsupported';
  }
  if (isWebGL2(gl)) {
    switch (status) {
      case gl.FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:
        return 'incomplete multisample';
      case gl.RENDERBUFFER_SAMPLES:
        return 'renderbuffer samples';
    }
  }
  return 'unknown error';
}

export function checkFramebufferStatus(gl: GLRenderingContext, message?: string) {
  const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
  if (status !== gl.FRAMEBUFFER_COMPLETE) {
    const description = getFramebufferStatusDescription(gl, status);
    throw new Error(`Framebuffer status: ${description}${message ? ` (${message})` : ''}`);
  }
}

function getFramebuffer(gl: GLRenderingContext) {
  const framebuffer = gl.createFramebuffer();
  if (framebuffer === null) {
    throw new Error('Could not create WebGL framebuffer');
  }
  return framebuffer;
}

/**
 * Allocates the storage of each attachment of framebuffers (by attachment point): textures and renderbuffers get their
 * storage when a framebuffer they are attached to is bound for use (`Framebuffer.bind`: drawing, reading, copying), not
 * when attached, so that targets of passes that never draw take no memory
 */
const AttachmentStorage = new WeakMap<WebGLFramebuffer, Map<number, { resource: AttachmentRef; ensure: () => void }>>();

/** A framebuffer of `createFramebuffer`: its GL object, also of copies (the resources wrap what they create) */
type GLFramebufferObject = Framebuffer & { getGLFramebuffer?: () => WebGLFramebuffer };

/**
 * Sets `resource` as `attachment` of `framebuffer` with what allocates its storage (`ensure`), or with `resource`
 * undefined removes the attachment
 */
export function setAttachmentStorage(
  gl: GLRenderingContext,
  framebuffer: Framebuffer,
  attachment: number,
  resource: AttachmentRef | undefined,
  ensure?: () => void,
) {
  const glFramebuffer = (framebuffer as GLFramebufferObject).getGLFramebuffer?.();
  // not an object with stand-ins of `gl` in tests
  if (typeof glFramebuffer !== 'object' || glFramebuffer === null) return;
  let storage = AttachmentStorage.get(glFramebuffer);
  if (!storage) AttachmentStorage.set(glFramebuffer, (storage = new Map()));
  // the combined point replaces the depth and stencil ones and the other way round
  if (attachment === gl.DEPTH_STENCIL_ATTACHMENT) {
    storage.delete(gl.DEPTH_ATTACHMENT);
    storage.delete(gl.STENCIL_ATTACHMENT);
  } else if (attachment === gl.DEPTH_ATTACHMENT || attachment === gl.STENCIL_ATTACHMENT) {
    storage.delete(gl.DEPTH_STENCIL_ATTACHMENT);
  }
  if (resource && ensure) storage.set(attachment, { resource, ensure });
  else storage.delete(attachment);
}

/** Binds `framebuffer` to change its attachments, without allocating their storage (see `setAttachmentStorage`) */
export function bindForAttachment(gl: GLRenderingContext, framebuffer: Framebuffer) {
  const glFramebuffer = (framebuffer as GLFramebufferObject).getGLFramebuffer?.();
  if (glFramebuffer) gl.bindFramebuffer(gl.FRAMEBUFFER, glFramebuffer);
  else framebuffer.bind();
}

/** The height of the attachments of framebuffers, to flip rectangles to GL's bottom-left origin */
const AttachmentHeights = new WeakMap<Framebuffer, () => number>();

/** Records the height of an attachment of `framebuffer` (all attachments have the same size) */
export function setFramebufferHeight(framebuffer: Framebuffer, getHeight: () => number) {
  AttachmentHeights.set(framebuffer, getHeight);
}

/** The height of the attachments of `framebuffer`, 0 without attachments */
export function getFramebufferHeight(framebuffer: Framebuffer) {
  return AttachmentHeights.get(framebuffer)?.() ?? 0;
}

/** GL framebuffers with an integer texture as first color attachment, cleared with `clearBufferiv` */
const IntegerColorFramebuffers = new WeakSet<WebGLFramebuffer>();

/** Records whether the first color attachment of the bound framebuffer is an integer texture */
export function setIntegerColorAttachment(gl: GLRenderingContext, isInteger: boolean) {
  const framebuffer = gl.getParameter(gl.FRAMEBUFFER_BINDING) as WebGLFramebuffer | null;
  if (!framebuffer) return;
  if (isInteger) IntegerColorFramebuffers.add(framebuffer);
  else IntegerColorFramebuffers.delete(framebuffer);
}

/** Whether the first color attachment of the bound framebuffer is an integer texture */
export function hasIntegerColorAttachment(gl: GLRenderingContext) {
  const framebuffer = gl.getParameter(gl.FRAMEBUFFER_BINDING) as WebGLFramebuffer | null;
  return !!framebuffer && IntegerColorFramebuffers.has(framebuffer);
}

function setDrawBuffers(extensions: WebGLExtensions, count: number) {
  const { drawBuffers } = extensions;
  if (!drawBuffers) {
    if (count > 1) throw new Error('multiple draw buffers not supported');
    return;
  }
  const buffers: number[] = [];
  for (let i = 0; i < count; ++i) buffers.push(drawBuffers.COLOR_ATTACHMENT0 + i);
  drawBuffers.drawBuffers(buffers);
}

export function createFramebuffer(gl: GLRenderingContext, extensions: WebGLExtensions, state: WebGLState): Framebuffer {
  const _framebuffer = getFramebuffer(gl);

  let destroyed = false;

  const framebuffer: GLFramebufferObject = {
    id: getNextFramebufferId(),
    getGLFramebuffer: () => _framebuffer,

    bind: () => {
      // the attachments' storage before the framebuffer is used (it binds textures and renderbuffers)
      const storage = typeof _framebuffer === 'object' ? AttachmentStorage.get(_framebuffer) : undefined;
      if (storage) for (const { ensure } of storage.values()) ensure();
      gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer);
    },
    setDrawBuffers: (count: number) => {
      endRenderPass(state);
      gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer);
      setDrawBuffers(extensions, count);
    },
    getAttachments: () => {
      const storage = typeof _framebuffer === 'object' ? AttachmentStorage.get(_framebuffer) : undefined;
      return storage ? [...storage.values()].map((s) => s.resource) : [];
    },
    destroy: () => {
      if (destroyed) return;
      gl.deleteFramebuffer(_framebuffer);
      destroyed = true;
    },
  };
  return framebuffer;
}

//

export function createNullFramebuffer(): Framebuffer {
  return {
    id: getNextFramebufferId(),

    bind: () => {},
    setDrawBuffers: () => {},
    getAttachments: () => [],
    destroy: () => {},
  };
}
