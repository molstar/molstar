/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { beforeEach, describe, expect, it, vi } from 'vitest';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TextureReader, TextureReadStatus } from '@molstar/graphics/gl/texture-reader';
import { AsyncPickStatus, createPickBuffers, type PickPass } from '../pick.js';

let time = 0;
vi.mock('@molstar/core/util/now', () => ({ now: () => time }));

function setup(maxAsyncReadLag: number) {
  let status: TextureReadStatus = 'pending';
  const reader: TextureReader = {
    synchronous: false,
    asynchronous: true,
    read: () => {},
    start: () => {
      status = 'pending';
    },
    check: () => status,
    reset: () => {},
    getByteCount: () => 0,
    dispose: () => {},
  };
  const graphics = { createTextureReader: () => reader } as unknown as GraphicsContext;
  const pickPass = { getTextures: () => [] } as unknown as PickPass;
  const buffers = createPickBuffers(graphics, pickPass, maxAsyncReadLag);
  buffers.setViewport(0, 0, 1, 1);
  return { buffers, resolve: () => (status = 'resolved') };
}

describe('PickBuffers', () => {
  beforeEach(() => {
    time = 0;
  });

  it('counts the lag of a read in frames, not checks', () => {
    const { buffers, resolve } = setup(2);
    buffers.asyncRead();
    // a click and the hover polling the same read in the same frames
    for (let frame = 1; frame <= 2; ++frame) {
      time = frame * 17;
      for (let i = 0; i < 4; ++i) expect(buffers.check()).toBe(AsyncPickStatus.Pending);
    }
    resolve();
    expect(buffers.check()).toBe(AsyncPickStatus.Resolved);
  });

  it('gives a read up after the lag', () => {
    const { buffers } = setup(2);
    buffers.asyncRead();
    time = 17;
    expect(buffers.check()).toBe(AsyncPickStatus.Pending);
    time = 34;
    expect(buffers.check()).toBe(AsyncPickStatus.Pending);
    time = 51;
    expect(buffers.check()).toBe(AsyncPickStatus.Failed);
  });
});
