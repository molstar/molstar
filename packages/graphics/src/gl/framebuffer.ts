/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Framebuffer {
  readonly id: number;

  bind: () => void;
  /** Number of color attachments written by draws, kept on reset; default 1 */
  setDrawBuffers: (count: number) => void;
  reset: () => void;
  destroy: () => void;
}
