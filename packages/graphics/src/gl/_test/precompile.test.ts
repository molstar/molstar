/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { GraphicsContext } from '../context.js';
import {
  createSignatureProgramItem,
  getSchemaLinearTextures,
  type ObjectSignature,
  PrecompileRecorder,
} from '../precompile.js';
import { createRenderObject, getRenderableProgramSource } from '../render-object.js';
import { Scene } from '../scene.js';
import { ValueCell } from '@molstar/core/util/value-cell';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import type { DefineValues } from '../shader-code.js';
import { createContext } from '../webgl/context.js';
import { getDataTextureSchema } from '../webgl/data-textures.js';
import { createGl } from './gl.shim.js';

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

describe('precompile signatures on WebGL', () => {
  it('make the programs render items make: data textures, and loop bounds as defines on WebGL1', () => {
    const ctx = createContext(createGl(32, 32));
    const { schema } = getRenderableProgramSource({ type: 'mesh', merged: false })!;
    const defines: { [k: string]: string | number | boolean } = {};
    for (const k of Object.keys(schema)) {
      if (schema[k].type !== 'define' || k === 'dRenderVariant' || k === 'dLightCount' || k === 'dColorMarker')
        continue;
      const kind = (schema[k] as { kind: string }).kind;
      defines[k] = k === 'dGeometryType' ? 'mesh' : k === 'dColorType' ? 'texture' : kind === 'number' ? 0 : false;
    }
    const signature: ObjectSignature = {
      type: 'mesh',
      merged: false,
      defines,
      linearTextures: [],
      doubleSided: false,
      hasReflection: false,
    };
    const schemas: unknown[] = [];
    const program = ctx.resources.program;
    ctx.resources.program = (defineValues, shaderCode, s) => {
      schemas.push(s);
      return program(defineValues, shaderCode, s);
    };
    const item = createSignatureProgramItem(ctx, signature, 'blended', { dLightCount: 1, dColorMarker: true })!;
    item.getProgram('color');
    expect(schemas.length).toBe(1);
    const programSchema = schemas[0] as typeof schema;
    // the data arrays as data textures, with their accessors
    expect(programSchema.tColor).toEqual(getDataTextureSchema(schema).tColor);
    // WebGL1: loop bounds as number defines, valued from the defines mirroring them
    expect(ctx.isWebGL2).toBe(false);
    expect(programSchema.uClipObjectCount.type).toBe('define');
    item.destroy();
    ctx.destroy();
  });

  it('record a renderable again only where what its signature reads changed, and after a reset', () => {
    const ctx = createContext(createGl(32, 32));
    const props = PD.getDefaultValues(Lines.Params);
    const values = Lines.Utils.createValuesSimple(Lines.createEmpty(), props, ColorNames.orange, 1);
    const object = createRenderObject('lines', values, Lines.Utils.createRenderableState(props), -1);
    const scene = Scene.create(ctx);
    scene.add(object);
    scene.commit();
    const renderable = scene.renderables[0];
    const recorder = new PrecompileRecorder();
    expect(recorder.recordAgain([renderable])).toBe(true);
    // nothing it reads changed (a state set elsewhere): the same, not signed again
    ValueCell.update(values.uColor, values.uColor.ref.value);
    expect(recorder.recordAgain([renderable])).toBe(false);
    expect(recorder.size).toBe(1);
    // a define changed: another signature
    ValueCell.update(values.dOverpaint, true);
    expect(recorder.recordAgain([renderable])).toBe(true);
    expect(recorder.size).toBe(2);
    // forgotten: recorded again, its signature the same
    recorder.reset();
    expect(recorder.recordAgain([renderable])).toBe(false);
    expect(recorder.size).toBe(1);
    scene.clear();
    ctx.destroy();
  });
});
