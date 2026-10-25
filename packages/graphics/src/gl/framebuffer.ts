/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** A texture or renderbuffer by its id (the resources may wrap the objects they create, ids stay) */
export interface AttachmentRef {
  readonly kind: 'texture' | 'renderbuffer';
  readonly id: number;
}

export interface Framebuffer {
  readonly id: number;

  bind: () => void;
  /** Number of color attachments written by draws; default 1 */
  setDrawBuffers: (count: number) => void;
  /** The textures and renderbuffers attached (layers of 3d textures as their texture) */
  getAttachments: () => readonly AttachmentRef[];
  destroy: () => void;
}
