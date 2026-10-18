/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { GraphicsContext } from '../context.js';
import { createSignatureProgramItem, getSchemaLinearTextures, type ObjectSignature } from '../precompile.js';
import { getRenderableProgramSource } from '../render-object.js';
import type { DefineValues } from '../shader-code.js';

/** A context that only records the define values program items are created with */
function createRecordingContext() {
  const created: DefineValues[] = [];
  const graphics = {
    createProgramItem: (_drawMode: unknown, _code: unknown, _schema: unknown, defineValues: DefineValues) => {
      created.push(defineValues);
      return { id: 0, getProgram: () => undefined, destroy: () => {} };
    },
  } as unknown as GraphicsContext;
  return { graphics, created };
}

describe('precompile signatures', () => {
  it('create program items with the define values in schema order, the context defines from the globals', () => {
    const { schema } = getRenderableProgramSource({ type: 'spheres', merged: false })!;
    const defineKeys = Object.keys(schema).filter((k) => schema[k].type === 'define' && k !== 'dRenderVariant');
    // the signature's defines in another order (e.g. read from JSON), without the context defines
    const defines: { [k: string]: string | number | boolean } = {};
    for (const k of [...defineKeys].reverse()) {
      if (k === 'dLightCount' || k === 'dColorMarker') continue;
      defines[k] = k === 'dGeometryType' ? 'spheres' : false;
    }
    const signature: ObjectSignature = {
      type: 'spheres',
      merged: false,
      defines,
      linearTextures: [],
      doubleSided: false,
      hasReflection: false,
    };
    const { graphics, created } = createRecordingContext();
    createSignatureProgramItem(graphics, signature, 'wboit', { dLightCount: 2, dColorMarker: true });
    expect(created.length).toBe(1);
    // WebGL keys programs by the define values in order: as render items have them (`splitValues`)
    expect(Object.keys(created[0])).toEqual(defineKeys);
    expect(created[0].dLightCount.ref.value).toBe(2);
    expect(created[0].dColorMarker.ref.value).toBe(true);
  });

  it('take the linear textures from the schema, as far as their conditions hold', () => {
    const text = getRenderableProgramSource({ type: 'text', merged: false })!.schema;
    expect(getSchemaLinearTextures(text, { dColorType: 'texture' })).toEqual(['tFont']);
    const mesh = getRenderableProgramSource({ type: 'mesh', merged: false })!.schema;
    expect(getSchemaLinearTextures(mesh, { dColorType: 'texture' })).toEqual([]);
    expect(getSchemaLinearTextures(mesh, { dColorType: 'grid' })).toContain('tColorGrid');
  });
});
