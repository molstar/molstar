/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { MeshSchema } from '../renderable/mesh.js';
import {
  DataPageName,
  getDataLayerAccessors,
  getDataLayerSchema,
  getDataOffsetName,
  usesDataLayers,
} from '../webgpu/data-layers.js';

describe('data layers', () => {
  it('extend a schema with layers by their offsets and the page, once', () => {
    const schema = getDataLayerSchema(MeshSchema);
    expect(getDataLayerSchema(MeshSchema)).toBe(schema);
    expect(schema[getDataOffsetName('tColor')]).toMatchObject({ type: 'uniform', kind: 'i', block: 'object' });
    expect(schema[getDataOffsetName('tMarker')]).toBeDefined();
    expect(schema[DataPageName]).toMatchObject({ type: 'storage-buffer', element: 'u32' });
    // the layers stay in the schema (WebGL, compute and the values' textures)
    expect(schema.tColor).toBe(MeshSchema.tColor);
  });

  it('leave schemas without layers as they are', () => {
    const schema = { uFoo: MeshSchema.uAlpha };
    expect(getDataLayerSchema(schema)).toBe(schema);
  });

  it('generate the accessors of referenced layers of the schema only', () => {
    const schema = getDataLayerSchema(MeshSchema);
    const referenced = new Set(['read_tColor', 'read_tMarker']);
    const isReferenced = (n: string) => referenced.has(n);
    const code = getDataLayerAccessors(schema, isReferenced);
    expect(code).toContain('fn read_tColor(i: i32) -> vec4<f32>');
    expect(code).toContain(`unpack4x8unorm(${DataPageName}[u32(uObject.uColorOffset) + u32(i)])`);
    // a byte per element, four per word
    expect(code).toContain(`${DataPageName}[u32(uObject.uMarkerOffset) + (u32(i) >> 2u)]`);
    expect(code).not.toContain('read_tOverpaint');
    expect(usesDataLayers(schema, isReferenced)).toBe(true);
    // a program of another schema with a function of that name of its own (a full-screen pass reading `tColor`)
    expect(getDataLayerAccessors(MeshSchema, isReferenced)).toBe('');
    expect(usesDataLayers(MeshSchema, isReferenced)).toBe(false);
  });
});
