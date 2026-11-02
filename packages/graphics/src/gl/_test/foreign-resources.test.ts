/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { DirectVolume } from '@molstar/graphics/geo/geometry/direct-volume/direct-volume';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import type { StorageBuffer } from '../compute.js';
import type { GraphicsContext } from '../context.js';
import {
  createForeignResourcesCheck,
  createRenderObject,
  hasForeignResources,
  transferTextures,
} from '../render-object.js';
import { createNullTexture, type Texture, transferTexture } from '../texture.js';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

function createTextureMesh() {
  const textureMesh = TextureMesh.createEmpty();
  const props = PD.getDefaultValues(TextureMesh.Params);
  const values = TextureMesh.Utils.createValuesSimple(textureMesh, props, ColorNames.orange, 1);
  const state = TextureMesh.Utils.createRenderableState(props);
  return createRenderObject('texture-mesh', values, state, -1);
}

function fakeBuffer(): StorageBuffer {
  return {
    id: 0,
    byteLength: 16,
    write: () => {},
    clear: () => {},
    copyTo: () => {},
    readAsync: async () => new ArrayBuffer(16),
    getByteCount: () => 16,
    destroy: () => {},
  };
}

/** A context owning every texture, and of the buffers `owned` (no compute context without) */
function fakeGraphics(owned?: StorageBuffer[]) {
  return {
    resources: { ownsTexture: () => true },
    compute: owned && { ownsBuffer: (b: StorageBuffer) => owned.includes(b) },
  } as unknown as GraphicsContext;
}

// GPU surfaces in the vertex arena hold storage buffers, not textures: after a restore or a backend switch they are
// of the previous context, computed again by their visuals
describe('foreign resources of render objects', () => {
  it('leaves a 3d texture to its producer on a context without them', () => {
    const data = { array: new Uint8Array(8), width: 2, height: 2, depth: 2 };
    const source = (kind: string) =>
      ({
        getSource: () => ({ kind, format: 'r', type: 'ubyte', filter: 'linear', data, mipmap: false }),
      }) as unknown as Texture;
    const load = vi.fn();
    const graphics = (texture3d: boolean) =>
      ({ capabilities: { texture3d }, resources: { texture: () => ({ load }) } }) as unknown as GraphicsContext;
    expect(transferTexture(source('volume-uint8'), graphics(false))).toBeUndefined();
    expect(load).not.toHaveBeenCalled();
    expect(transferTexture(source('volume-uint8'), graphics(true))).toBeDefined();
    expect(transferTexture(source('image-uint8'), graphics(false))).toBeDefined();
  });

  it('leaves the 2d grid of a volume to its visual on a context with 3d textures', () => {
    const directVolume = DirectVolume.createEmpty();
    const props = PD.getDefaultValues(DirectVolume.Params);
    const values = DirectVolume.Utils.createValuesSimple(directVolume, props, ColorNames.orange, 1);
    const o = createRenderObject('direct-volume', values, DirectVolume.Utils.createRenderableState(props), -1);
    const grid = {
      ...createNullTexture(),
      id: 1,
      getDepth: () => 0,
      getSource: () => ({
        kind: 'image-uint8',
        format: 'rgba',
        type: 'ubyte',
        filter: 'linear',
        data: {},
        mipmap: false,
      }),
    } as unknown as Texture;
    const transferred = { ...createNullTexture(), id: 2, load: () => {} } as unknown as Texture;
    const graphics = (texture3d: boolean) =>
      ({
        capabilities: { texture3d },
        resources: { ownsTexture: (t: Texture) => t === transferred, texture: () => transferred },
      }) as unknown as GraphicsContext;
    ValueCell.update(values.tGridTex, grid);
    expect(hasForeignResources(o, graphics(true))).toBe(true);
    // stays foreign: the visual creates it again
    expect(transferTextures(o, graphics(true))).toBe(false);
    expect(values.tGridTex.ref.value).toBe(grid);
    expect(hasForeignResources(o, graphics(true))).toBe(true);
    // still 2d on WebGL1
    expect(transferTextures(o, graphics(false))).toBe(true);
    expect(values.tGridTex.ref.value).toBe(transferred);
  });

  it('finds storage buffers of another compute context, also of draw arguments', () => {
    const o = createTextureMesh();
    const values = o.values as unknown as {
      sIsoTable: ValueCell<StorageBuffer | undefined>;
      drawIndirect: ValueCell<{ buffer: StorageBuffer; byteOffset: number } | undefined>;
    };
    const table = fakeBuffer();
    expect(hasForeignResources(o, fakeGraphics([]))).toBe(false);

    ValueCell.update(values.sIsoTable, table);
    expect(hasForeignResources(o, fakeGraphics([table]))).toBe(false);
    expect(hasForeignResources(o, fakeGraphics([]))).toBe(true);
    // WebGL: no compute context, no buffer is its
    expect(hasForeignResources(o, fakeGraphics())).toBe(true);
    // not transferred: computed again
    expect(transferTextures(o, fakeGraphics([]))).toBe(false);
    expect(transferTextures(o, fakeGraphics([table]))).toBe(true);

    const args = fakeBuffer();
    ValueCell.update(values.drawIndirect, { buffer: args, byteOffset: 0 });
    expect(hasForeignResources(o, fakeGraphics([table]))).toBe(true);
    expect(hasForeignResources(o, fakeGraphics([table, args]))).toBe(false);
  });

  it('fills a texture mesh last filled into the vertex arena into textures', () => {
    const mesh = TextureMesh.createEmpty();
    ValueCell.update(mesh.isoChunks, true);
    ValueCell.update(mesh.isoTable, fakeBuffer());
    ValueCell.update(mesh.isoVertices, fakeBuffer());
    ValueCell.update(mesh.drawIndirect, { buffer: fakeBuffer(), byteOffset: 0 });
    const texture = createNullTexture();
    TextureMesh.create(3, 1, texture, texture, texture, mesh.boundingSphere, mesh);
    expect(mesh.isoChunks.ref.value).toBe(false);
    expect(mesh.isoTable.ref.value).toBeUndefined();
    expect(mesh.isoVertices.ref.value).toBeUndefined();
    expect(mesh.drawIndirect.ref.value).toBeUndefined();
  });
});

// a visual checks its render object on every update: a clean one is remembered for the context
describe('the foreign resources check of a visual', () => {
  it('walks the values once while the render object and the context stay', () => {
    const o = createTextureMesh();
    const table = fakeBuffer();
    (o.values as unknown as { sIsoTable: ValueCell<StorageBuffer | undefined> }).sIsoTable.ref.value = table;
    const graphics = fakeGraphics([table]);
    const owns = vi.spyOn(graphics.compute!, 'ownsBuffer');
    const check = createForeignResourcesCheck();
    expect(check(o, graphics)).toBe(false);
    expect(owns).toHaveBeenCalledTimes(1);
    expect(check(o, graphics)).toBe(false);
    expect(check(o, graphics)).toBe(false);
    expect(owns).toHaveBeenCalledTimes(1);
  });

  it('checks again for another context (a restore, a backend switch) or another render object', () => {
    const o = createTextureMesh();
    const table = fakeBuffer();
    const check = createForeignResourcesCheck();
    expect(check(o, fakeGraphics([table]))).toBe(false);
    (o.values as unknown as { sIsoTable: ValueCell<StorageBuffer | undefined> }).sIsoTable.ref.value = table;
    // the buffer is of the previous context
    expect(check(o, fakeGraphics([]))).toBe(true);
    // a foreign result is not remembered: checked again until it is clean
    const next = fakeGraphics([table]);
    expect(check(o, next)).toBe(false);
    expect(check(createTextureMesh(), fakeGraphics([]))).toBe(false);
  });
});
