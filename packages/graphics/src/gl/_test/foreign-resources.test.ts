/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { TextureMesh } from '@molstar/graphics/geo/geometry/texture-mesh/texture-mesh';
import type { StorageBuffer } from '../compute.js';
import type { GraphicsContext } from '../context.js';
import {
  createForeignResourcesCheck,
  createRenderObject,
  hasForeignResources,
  transferTextures,
} from '../render-object.js';
import { createNullTexture } from '../texture.js';
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
