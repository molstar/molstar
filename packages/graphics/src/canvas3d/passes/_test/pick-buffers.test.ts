/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { beforeEach, describe, expect, it, vi } from 'vitest';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { TextureReader, TextureReadStatus } from '@molstar/graphics/gl/texture-reader';
import { PickingId } from '@molstar/graphics/geo/geometry/picking';
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

  it('decodes object, instance and group ids in 32 bits, the white clear as the null', () => {
    // as `packPickId`: id + 1, its 24 low bits in rgb, the top byte in alpha
    const pack = (id: number) => {
      const v = id + 1;
      return [(v >>> 16) & 255, (v >>> 8) & 255, v & 255, v >>> 24];
    };
    let ids: number[][] = [];
    const reader = {
      synchronous: true,
      asynchronous: false,
      read: (_t: unknown, _x: number, _y: number, _w: number, _h: number, arrays: Uint8Array[]) => {
        for (let k = 0; k < 3; ++k) arrays[k].set(ids[k]);
      },
      start: () => {},
      check: () => 'resolved',
      reset: () => {},
      getByteCount: () => 0,
      dispose: () => {},
    } as unknown as TextureReader;
    const graphics = { createTextureReader: () => reader } as unknown as GraphicsContext;
    const buffers = createPickBuffers(graphics, { getTextures: () => [] } as unknown as PickPass, 2);
    buffers.setViewport(0, 0, 1, 1);

    ids = [pack(7), pack(2 ** 24 + 1), pack(2 ** 24 + 3)];
    buffers.read();
    expect(buffers.getPickingId(0, 0)).toEqual({ objectId: 7, instanceId: 2 ** 24 + 1, groupId: 2 ** 24 + 3 });

    // the group target's clear: the whole instance (e.g. a mesh's cap), as the visuals read PickingId.Null
    const white = [255, 255, 255, 255];
    ids = [pack(7), pack(0), white];
    buffers.read();
    expect(buffers.getPickingId(0, 0)).toEqual({ objectId: 7, instanceId: 0, groupId: PickingId.Null });
    // no object
    ids = [white, pack(0), pack(0)];
    buffers.read();
    expect(buffers.getPickingId(0, 0)).toBeUndefined();
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
