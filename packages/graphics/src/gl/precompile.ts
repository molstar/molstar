/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import type { GraphicsContext } from './context.js';
import {
  getRenderableKind,
  getRenderableProgramSource,
  type GraphicsRenderObject,
  isMergedRenderObject,
  type RenderableKind,
} from './render-object.js';
import { getDefineLookup, isUsed, type RenderableSchema } from './renderable/schema.js';
import type { ProgramItem, Transparency } from './render-item.js';
import type { GraphicsRenderable } from './renderable.js';
import type { GlobalDefines } from './renderable/schema.js';
import type { DefineType, DefineValues } from './shader-code.js';

/**
 * What the programs of a render object depend on, without its data: its kind (schema and shader code), its define
 * values and what selects its render states (`renderer.collect`); no render variant, transparency mode, render state or
 * target, which the context and the props that draw it give. Plain data, the entries of a precompile list.
 */
export interface ObjectSignature extends RenderableKind {
  /** the values of the schema's defines, without the render variant and the context's (`ContextDefines`) */
  readonly defines: { readonly [k: string]: string | number | boolean };
  /** the textures read with linear filtering (`RenderItem.getLinearTextures`) */
  readonly linearTextures: readonly string[];
  readonly doubleSided: boolean;
  readonly hasReflection: boolean;
  /**
   * How the object may draw, which selects the passes it is drawn in (`renderer.collectSignatures`): opaque (alpha 1,
   * not fully transparent), transparent (alpha below 1 or a transparency layer), emissive; not part of the key,
   * signatures of the same key merge them (`mergeSignatures`)
   */
  readonly opaque: boolean;
  readonly transparent: boolean;
  readonly emissive: boolean;
  /**
   * Culled on the GPU (a merged object's `gpuCulling`): its programs are the same, the culling's kernels prepared with
   * them; not part of the key. Unset in lists derived from params (whether an object is culled there depends on its
   * instances)
   */
  readonly gpuCulled?: boolean;
}

/** How an object with `values` (its final alpha, `uAlpha`) draws now, see `ObjectSignature.opaque` */
function getDrawFlags(values: { [k: string]: ValueCell<unknown> | undefined }) {
  const alpha = (values.uAlpha?.ref.value as number) ?? 1;
  const transparency = (values.transparencyAverage?.ref.value as number) ?? 0;
  const emissive =
    ((values.emissiveAverage?.ref.value as number) ?? 0) + ((values.uEmissive?.ref.value as number) ?? 0);
  return {
    opaque: alpha === 1 && transparency !== 1,
    transparent: alpha < 1 || transparency > 0,
    emissive: emissive > 0,
  };
}

/** `a` with the draw flags of `b` added (signatures of the same key) */
export function mergeSignatures(a: ObjectSignature, b: ObjectSignature): ObjectSignature {
  if (
    (a.opaque || !b.opaque) &&
    (a.transparent || !b.transparent) &&
    (a.emissive || !b.emissive) &&
    (a.gpuCulled || !b.gpuCulled)
  )
    return a;
  return {
    ...a,
    opaque: a.opaque || b.opaque,
    transparent: a.transparent || b.transparent,
    emissive: a.emissive || b.emissive,
    gpuCulled: a.gpuCulled || b.gpuCulled,
  };
}

/**
 * The defines a context sets for all objects (`GlobalDefines`): not in signatures, taken from the context that prepares
 * one
 */
export const ContextDefines = new Set<string>(['dLightCount', 'dColorMarker']);

/** The signature of `renderable`, `undefined` if it wasn't created from a render object (see `createRenderable`) */
export function getObjectSignature(renderable: GraphicsRenderable): ObjectSignature | undefined {
  const kind = getRenderableKind(renderable);
  const source = kind && getRenderableProgramSource(kind);
  if (!kind || !source) return;
  const defines = getDefines(source.schema, renderable.values as { [k: string]: ValueCell<unknown> | undefined });
  return {
    type: kind.type,
    merged: kind.merged,
    defines,
    linearTextures: [...(renderable.renderItem.getLinearTextures?.() ?? [])].sort(),
    doubleSided: !!renderable.values.uDoubleSided?.ref.value,
    hasReflection: !!renderable.values.hasReflection?.ref.value,
    ...getDrawFlags(renderable.values as { [k: string]: ValueCell<unknown> | undefined }),
    gpuCulled: !!renderable.values.gpuCulling?.ref.value,
  };
}

type DefineValue = string | number | boolean;

/** The textures read with linear filtering: linear in `schema` and used with `defines` (see `isUsed`) */
export function getSchemaLinearTextures(schema: RenderableSchema, defines: { readonly [k: string]: DefineValue }) {
  const names: string[] = [];
  for (const [k, v] of Object.entries(defines)) {
    if (typeof v === 'string') names.push(`${k}_${v}`);
    else if (v !== false) names.push(k);
  }
  const lookup = getDefineLookup(names);
  const linear: string[] = [];
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type === 'texture' && spec.filter === 'linear' && isUsed(spec, lookup)) linear.push(k);
  }
  return linear.sort();
}

function getDefines(schema: RenderableSchema, values: { [k: string]: ValueCell<unknown> | undefined }) {
  const defines: { [k: string]: DefineValue } = {};
  for (const k of Object.keys(schema)) {
    if (schema[k].type !== 'define' || k === 'dRenderVariant' || ContextDefines.has(k)) continue;
    const value = values[k]?.ref.value;
    if (typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean') defines[k] = value;
  }
  return defines;
}

/**
 * The signature of a render object from its values, without a renderable (and so without a graphics context): its
 * linear textures by the schema (`getSchemaLinearTextures`), not by the textures a render item binds
 */
export function getRenderObjectSignature(object: GraphicsRenderObject): ObjectSignature | undefined {
  const kind = { type: object.type, merged: isMergedRenderObject(object) };
  const source = getRenderableProgramSource(kind);
  if (!source) return;
  const values = object.values as unknown as { [k: string]: ValueCell<unknown> | undefined };
  const defines = getDefines(source.schema, values);
  return {
    ...kind,
    defines,
    linearTextures: getSchemaLinearTextures(source.schema, defines),
    doubleSided: !!values.uDoubleSided?.ref.value,
    hasReflection: !!values.hasReflection?.ref.value,
    ...getDrawFlags(values),
    gpuCulled: !!values.gpuCulling?.ref.value,
  };
}

/** A key equal for equal signatures, to deduplicate them */
export function getSignatureKey(signature: ObjectSignature) {
  const defines = Object.keys(signature.defines)
    .sort()
    .map((k) => `${k}=${signature.defines[k]}`)
    .join(',');
  return [
    signature.type,
    signature.merged ? 'merged' : '',
    defines,
    signature.linearTextures.join(','),
    signature.doubleSided ? 'double' : '',
    signature.hasReflection ? 'reflection' : '',
  ].join('|');
}

/**
 * The programs of `signature` in `graphics` (see `GraphicsContext.createProgramItem`), for objects drawn with
 * `transparency`; the context's own defines (light count, color marker) from `globals`. `undefined` for a kind the
 * renderables don't know
 */
export function createSignatureProgramItem(
  graphics: GraphicsContext,
  signature: ObjectSignature,
  transparency: Transparency,
  globals: GlobalDefines,
): ProgramItem | undefined {
  const source = getRenderableProgramSource(signature);
  if (!source) return;
  // in schema order, as render items have them: WebGL keys its programs by the defines in order
  const defineValues: DefineValues = {};
  for (const k of Object.keys(source.schema)) {
    if (source.schema[k].type !== 'define' || k === 'dRenderVariant') continue;
    const value = ContextDefines.has(k) ? globals[k as keyof GlobalDefines] : signature.defines[k];
    // number defines (e.g. counts) are define values as well, see `getWgslDefines`
    if (value !== undefined) defineValues[k] = ValueCell.create(value as DefineType);
  }
  return graphics.createProgramItem(
    'triangles',
    source.shaderCode,
    source.schema,
    defineValues,
    signature.linearTextures,
    transparency,
  );
}

/**
 * The objects a scene drew (or will draw), as signatures: what `Canvas3D.precompile` prepares the programs and pipelines
 * of, for given props, before the objects exist (e.g. the snapshots of a state, before they are applied)
 */
export interface PrecompileList {
  readonly objects: readonly ObjectSignature[];
}

/** Records the signatures of renderables, deduplicated, since the last reset */
export class PrecompileRecorder {
  private readonly signatures = new Map<string, ObjectSignature>();
  private recorded = new WeakSet<GraphicsRenderable>();

  /** records the renderables not recorded before, or since `invalidate` */
  record(renderables: Iterable<GraphicsRenderable>) {
    for (const r of renderables) {
      if (this.recorded.has(r)) continue;
      this.recorded.add(r);
      const signature = getObjectSignature(r);
      if (!signature) continue;
      const key = getSignatureKey(signature);
      const recorded = this.signatures.get(key);
      this.signatures.set(key, recorded ? mergeSignatures(recorded, signature) : signature);
    }
  }

  /** the next `record` takes the renderables again, e.g. after their values were updated (defines may have changed) */
  invalidate() {
    this.recorded = new WeakSet();
  }

  /** forgets what was recorded */
  reset() {
    this.signatures.clear();
    this.invalidate();
  }

  get size() {
    return this.signatures.size;
  }

  getList(): PrecompileList {
    return { objects: [...this.signatures.values()] };
  }
}
