/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Texture } from './texture.js';

export type TextureReadStatus = 'pending' | 'resolved' | 'failed';

/** `rgba8`: four bytes per texel (`Uint8Array`); `r32f`: the red channel of float textures (`Float32Array`) */
export type TextureReadFormat = 'rgba8' | 'r32f';
export type TextureReadArray<F extends TextureReadFormat> = F extends 'r32f' ? Float32Array : Uint8Array;

/**
 * Repeated reads of the same rectangle of 2d textures (e.g. the pick targets) into CPU arrays, one per texture (`width *
 * height` texels in the reader's format), rows top to bottom (the rectangle and the rows have a top-left origin).
 * Backends read synchronously, asynchronously, or both; see `GraphicsContext.createTextureReader`.
 */
export interface TextureReader<A extends Uint8Array | Float32Array = Uint8Array> {
  /** whether `read` reads synchronously */
  readonly synchronous: boolean;
  /** whether `start` reads asynchronously */
  readonly asynchronous: boolean;
  /** reads now, into `arrays` (only where `synchronous`) */
  read(textures: readonly Texture[], x: number, y: number, width: number, height: number, arrays: A[]): void;
  /** starts reading (only where `asynchronous`), see `check` */
  start(textures: readonly Texture[], x: number, y: number, width: number, height: number): void;
  /** the state of the read started last, its data copied into `arrays` once resolved; a read is resolved only once */
  check(arrays: A[]): TextureReadStatus;
  /** drops a started read (e.g. one pending for too long) */
  reset(): void;
  /** bytes of GPU memory held for reads */
  getByteCount(): number;
  dispose(): void;
}

export interface TextureReaderOptions<F extends TextureReadFormat = TextureReadFormat> {
  /** `rgba8` by default */
  readonly format?: F;
}
