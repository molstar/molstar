/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { Framebuffer } from '../framebuffer.js';
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
  let _framebuffer = getFramebuffer(gl);
  let drawBufferCount = 1;

  let destroyed = false;

  return {
    id: getNextFramebufferId(),

    bind: () => gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer),
    setDrawBuffers: (count: number) => {
      endRenderPass(state);
      gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer);
      setDrawBuffers(extensions, count);
      drawBufferCount = count;
    },
    reset: () => {
      _framebuffer = getFramebuffer(gl);
      if (drawBufferCount !== 1) {
        gl.bindFramebuffer(gl.FRAMEBUFFER, _framebuffer);
        setDrawBuffers(extensions, drawBufferCount);
      }
    },
    destroy: () => {
      if (destroyed) return;
      gl.deleteFramebuffer(_framebuffer);
      destroyed = true;
    },
  };
}

//

export function createNullFramebuffer(): Framebuffer {
  return {
    id: getNextFramebufferId(),

    bind: () => { },
    setDrawBuffers: () => { },
    reset: () => { },
    destroy: () => { },
  };
}
