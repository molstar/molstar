/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getAttributeLayout } from '../attribute-layout.js';
import { CylindersSchema } from '../renderable/cylinders.js';
import { DirectVolumeSchema } from '../renderable/direct-volume.js';
import { ImageSchema } from '../renderable/image.js';
import { LinesSchema } from '../renderable/lines.js';
import { MeshSchema } from '../renderable/mesh.js';
import { PointsSchema } from '../renderable/points.js';
import {
  AttributeSpec,
  BaseSchema,
  getDefineLookup,
  getTextureDimension,
  GlobalTextureSchema,
  isUsed,
  type RenderableSchema,
  SegmentSchema,
  TextureSpec,
} from '../renderable/schema.js';
import { SpheresSchema } from '../renderable/spheres.js';
import { TextSchema } from '../renderable/text.js';
import { TextureMeshSchema } from '../renderable/texture-mesh.js';
import { getBindingLayout } from '../uniform-layout.js';

describe('texture spec', () => {
  it('defaults', () => {
    const t = (spec: TextureSpec<any>) => [spec.dimension, spec.sampleType];
    expect(t(TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'))).toEqual(['2d', 'float']);
    expect(t(TextureSpec('volume-float32', 'r', 'float', 'linear'))).toEqual(['3d', 'float']);
    expect(t(TextureSpec('image-float32', 'rgba', 'float', 'nearest'))).toEqual(['2d', 'unfilterable-float']);
    expect(t(TextureSpec('image-int32', 'r', 'int', 'nearest'))).toEqual(['2d', 'sint']);
    expect(t(TextureSpec('texture', 'depth', 'ushort', 'nearest'))).toEqual(['2d', 'depth']);
    expect(t(TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, { dimension: 'cube' }))).toEqual([
      'cube',
      'float',
    ]);
  });

  it('define-dependent dimension', () => {
    const spec = DirectVolumeSchema.tGridTex;
    expect(getTextureDimension(spec, getDefineLookup(['dGridTexType_3d']))).toBe('3d');
    expect(getTextureDimension(spec, getDefineLookup(['dGridTexType_2d']))).toBe('2d');
    expect(() => getTextureDimension(spec)).toThrow();
    expect(() => getTextureDimension(spec, getDefineLookup([]))).toThrow();
  });
});

describe('used if', () => {
  it('conditions', () => {
    const spec = AttributeSpec('float32', 1, 0, ['dA', '!dB']);
    expect(isUsed(AttributeSpec('float32', 1, 0), getDefineLookup([]))).toBe(true);
    expect(isUsed(spec, getDefineLookup(['dA', 'dB']))).toBe(true);
    expect(isUsed(spec, getDefineLookup([]))).toBe(true);
    expect(isUsed(spec, getDefineLookup(['dB']))).toBe(false);
    // unknown names may hold
    expect(isUsed(spec, (name) => (name === 'dA' ? false : undefined))).toBe(true);
    expect(isUsed(spec, (name) => (name === 'dA' ? false : name === 'dB'))).toBe(false);
  });

  it('vertex id emulation', () => {
    expect(isUsed(BaseSchema.aVertex, getDefineLookup(['emulateVertexId']))).toBe(true);
    expect(isUsed(BaseSchema.aVertex, getDefineLookup([]))).toBe(false);
  });
});

describe('binding layout', () => {
  it('resolves textures with defines', () => {
    const schema = { ...GlobalTextureSchema, ...DirectVolumeSchema };
    const all = getBindingLayout({ ...GlobalTextureSchema, ...BaseSchema });
    expect(all.global.textures.map((t) => t.name)).toEqual(Object.keys(GlobalTextureSchema));
    expect(all.object.textures.map((t) => t.binding)).toEqual(all.object.textures.map((_, i) => i + 1));

    const layout = getBindingLayout(schema, getDefineLookup(['dGridTexType_3d', 'dUsePalette']));
    expect(layout.global.textures.map((t) => t.name)).toEqual(['tDepth']);
    const object = layout.object.textures.map((t) => t.name);
    expect(object).toContain('tPalette');
    expect(object).not.toContain('tColorGrid');
    expect(layout.object.textures.find((t) => t.name === 'tGridTex')).toEqual({
      name: 'tGridTex',
      binding: object.indexOf('tGridTex') + 1,
      dimension: '3d',
      sampleType: 'float',
      filter: 'linear',
    });
  });
});

describe('attribute layout', () => {
  it('locations', () => {
    const layout = getAttributeLayout({
      aInstance: AttributeSpec('float32', 1, 1),
      aTransform: AttributeSpec('float32', 16, 1),
      aVertex: AttributeSpec('float32', 1, 0, ['emulateVertexId']),
      aPosition: AttributeSpec('float32', 3, 0),
      aGroup: AttributeSpec('float32', 1, 0),
    });
    expect(layout.entries.map((e) => [e.name, e.location])).toEqual([
      ['aPosition', 0],
      ['aInstance', 1],
      ['aTransform', 2],
      ['aVertex', 6],
      ['aGroup', 7],
    ]);
    expect(layout.byName.aTransform.locationCount).toBe(4);
    expect(layout.size).toBe(8);
  });

  it('location 0 falls back to instanced attributes', () => {
    const layout = getAttributeLayout({ ...BaseSchema, ...TextureMeshSchema });
    expect(layout.entries[0].name).toBe('aInstance');
  });

  // WebGL guarantees at least 16 vertex attributes, as does WebGPU's default `maxVertexAttributes`
  const schemas: [string, RenderableSchema][] = [
    ['mesh', MeshSchema],
    ['texture-mesh', TextureMeshSchema],
    ['spheres', SpheresSchema],
    ['cylinders', CylindersSchema],
    ['points', PointsSchema],
    ['lines', LinesSchema],
    ['text', TextSchema],
    ['image', ImageSchema],
    ['direct-volume', DirectVolumeSchema],
  ];
  for (const [name, schema] of schemas) {
    it(`fits 16 locations: ${name}`, () => {
      const layout = getAttributeLayout({ ...schema, ...SegmentSchema });
      expect(layout.size).toBeLessThanOrEqual(16);
      const first = schema[layout.entries[0].name] as AttributeSpec<'float32'>;
      expect(first.usedIf).toBeUndefined();
    });
  }
});
