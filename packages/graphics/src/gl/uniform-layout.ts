/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import {
  type DefineLookup,
  getTextureDimension,
  GlobalTextureSchema,
  GlobalUniformSchema,
  isUsed,
  type RenderableSchema,
  type TextureDimension,
  type TextureSampleType,
  type UniformBlock,
  type UniformScope,
} from './renderable/schema.js';
import type { TextureFilter } from './texture.js';
import type { UniformKind, UniformsList } from './uniform.js';

/**
 * Backend-neutral description of the uniform blocks and texture bindings of a schema.
 * The block layout follows std140, which equals the WGSL uniform address space layout.
 * The block tiers map to fixed binding indices (WebGL2) or bind group indices (WebGPU).
 * Shader declarations are generated in the backend folders (`webgl/`, `webgpu/`).
 */

/** Block (struct) names and binding indices, shared across all programs of a context. */
export const UniformBlockInfo: { [k in UniformBlock]: { name: string; bindingIndex: number } } = {
  global: { name: 'Global', bindingIndex: 0 },
  material: { name: 'Material', bindingIndex: 1 },
  object: { name: 'Object', bindingIndex: 2 },
};

export const UniformBlocks: UniformBlock[] = ['global', 'material', 'object'];

export type UniformBlockMember = {
  name: string;
  kind: UniformKind;
  offset: number;
  /** array length, 0 for non-array members */
  length: number;
  /** std140 array stride in bytes, 0 for non-array members */
  stride: number;
};

export interface UniformBlockLayout {
  readonly block: UniformBlock;
  readonly members: UniformBlockMember[];
  /** std140 data size in bytes */
  readonly size: number;
  readonly hash: string;
}

export function alignTo(offset: number, alignment: number) {
  return Math.ceil(offset / alignment) * alignment;
}

export function getElementKind(kind: UniformKind): UniformKind {
  return (kind.endsWith('[]') ? kind.slice(0, -2) : kind) as UniformKind;
}

/** [baseAlignment, byteSize] according to std140 rules */
function getStd140AlignmentSize(kind: UniformKind): [number, number] {
  switch (kind) {
    case 'f':
    case 'i':
    case 'b':
      return [4, 4];
    case 'v2':
    case 'iv2':
      return [8, 8];
    case 'v3':
    case 'iv3':
      return [16, 12];
    case 'v4':
    case 'iv4':
    case 'q':
      return [16, 16];
    case 'm3':
      return [16, 48];
    case 'm4':
      return [16, 64];
    default:
      throw new Error(`unsupported uniform block member kind '${kind}'`);
  }
}

/** number of values per array element in the flat source array */
function getElementItemSize(kind: UniformKind) {
  switch (kind) {
    case 'f':
    case 'i':
    case 'b':
      return 1;
    case 'v2':
    case 'iv2':
      return 2;
    case 'v3':
    case 'iv3':
      return 3;
    case 'v4':
    case 'iv4':
    case 'q':
      return 4;
    case 'm3':
      return 9;
    case 'm4':
      return 16;
    default:
      throw new Error(`unsupported uniform block member kind '${kind}'`);
  }
}

/**
 * Compute the std140 layout of the uniforms of `schema` assigned to `block`.
 * The layout is fully determined by the (ordered) member list.
 */
export function getUniformBlockLayout(block: UniformBlock, schema: RenderableSchema): UniformBlockLayout | undefined {
  // assumes uniform specs are not added to a schema after first use
  let layouts = layoutCache.get(schema);
  if (layouts === undefined) {
    layouts = new Map();
    layoutCache.set(schema, layouts);
  }
  if (layouts.has(block)) return layouts.get(block);
  const layout = computeUniformBlockLayout(block, schema);
  layouts.set(block, layout);
  return layout;
}

const layoutCache = new WeakMap<RenderableSchema, Map<UniformBlock, UniformBlockLayout | undefined>>();

function computeUniformBlockLayout(block: UniformBlock, schema: RenderableSchema): UniformBlockLayout | undefined {
  // the global block must have an identical layout across all programs, so it
  // is always derived from `GlobalUniformSchema` (if the schema includes it)
  if (block === 'global' && schema !== GlobalUniformSchema) {
    const hasGlobals = Object.keys(schema).some((k) => {
      const spec = schema[k];
      return spec.type === 'uniform' && spec.block === 'global';
    });
    if (!hasGlobals) return;
    schema = GlobalUniformSchema;
  }

  const std140 = getStd140Members(block, schema);
  if (std140 === undefined) return;
  const { members, size } = std140;

  return {
    block,
    members,
    size,
    hash: `${block}|${members.map((m) => `${m.name}:${m.kind}${m.length ? m.length : ''}`).join(',')}`,
  };
}

/** std140 members and data size of the uniforms of `schema` in `scope`, in schema order */
export function getStd140Members(
  scope: UniformScope,
  schema: RenderableSchema,
): { members: UniformBlockMember[]; size: number } | undefined {
  const members: UniformBlockMember[] = [];
  let offset = 0;
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'uniform' || spec.block !== scope) continue;
    const [alignment, size] = getStd140AlignmentSize(getElementKind(spec.kind));
    if (spec.length !== undefined) {
      // std140: array elements are aligned and padded to vec4
      const stride = alignTo(size, 16);
      offset = alignTo(offset, 16);
      members.push({ name: k, kind: spec.kind, offset, length: spec.length, stride });
      offset += stride * spec.length;
    } else {
      offset = alignTo(offset, alignment);
      members.push({ name: k, kind: spec.kind, offset, length: 0, stride: 0 });
      offset += size;
    }
  }
  if (members.length === 0) return;
  return { members, size: alignTo(offset, 16) };
}

/** Block tier of `spec` when uniform buffers are used, `undefined` for classic uniforms */
export function getUniformBlock(spec: RenderableSchema[string]): UniformBlock | undefined {
  if (spec.type !== 'uniform' || spec.block === undefined || spec.block === 'draw') return;
  return spec.block;
}

//

export type TextureBinding = {
  readonly name: string;
  /** within the tier, 0 is the uniform block */
  readonly binding: number;
  readonly dimension: TextureDimension;
  readonly sampleType: TextureSampleType;
  readonly filter: TextureFilter;
};
export type BindingTier = { readonly uniforms?: UniformBlockLayout; readonly textures: readonly TextureBinding[] };
export type BindingLayout = { readonly [k in UniformBlock]: BindingTier };

/**
 * Uniform block and textures per tier, i.e. the future bind group layouts:
 * 'global' textures are those of `GlobalTextureSchema`, 'material' textures
 * have the 'material' variant, all other textures are per 'object'.
 * With `defines`, textures known to be unused are skipped and define-dependent
 * dimensions are resolved; without, all textures are included.
 */
export function getBindingLayout(schema: RenderableSchema, defines?: DefineLookup): BindingLayout {
  const textures: { [k in UniformBlock]: TextureBinding[] } = { global: [], material: [], object: [] };
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'texture' || (defines && !isUsed(spec, defines))) continue;
    // global by spec, not by name: compute schemas may have their own textures named like global ones (e.g. `tDepth`)
    const isGlobal = (GlobalTextureSchema as RenderableSchema)[k] === spec;
    const tier = isGlobal ? textures.global : spec.variant === 'material' ? textures.material : textures.object;
    tier.push({
      name: k,
      binding: tier.length + 1,
      dimension: getTextureDimension(spec, defines),
      sampleType: spec.sampleType,
      filter: spec.filter,
    });
  }
  return {
    global: { uniforms: getUniformBlockLayout('global', schema), textures: textures.global },
    material: { uniforms: getUniformBlockLayout('material', schema), textures: textures.material },
    object: { uniforms: getUniformBlockLayout('object', schema), textures: textures.object },
  };
}

//

/** views on the std140 data of one block instance */
export interface UniformBlockData {
  readonly f32: Float32Array;
  readonly i32: Int32Array;
}

export function createUniformBlockData(buffer: ArrayBuffer, byteOffset: number, byteSize: number): UniformBlockData {
  return {
    f32: new Float32Array(buffer, byteOffset, byteSize / 4),
    i32: new Int32Array(buffer, byteOffset, byteSize / 4),
  };
}

/** write array element `v[s...]` of `kind` at float offset `o` */
function packElement(f32: Float32Array, i32: Int32Array, kind: UniformKind, o: number, v: ArrayLike<any>, s: number) {
  switch (kind) {
    case 'f':
      f32[o] = v[s];
      break;
    case 'i':
      i32[o] = v[s];
      break;
    case 'b':
      i32[o] = v[s] ? 1 : 0;
      break;
    case 'v2':
      f32[o] = v[s];
      f32[o + 1] = v[s + 1];
      break;
    case 'v3':
      f32[o] = v[s];
      f32[o + 1] = v[s + 1];
      f32[o + 2] = v[s + 2];
      break;
    case 'v4':
    case 'q':
      f32[o] = v[s];
      f32[o + 1] = v[s + 1];
      f32[o + 2] = v[s + 2];
      f32[o + 3] = v[s + 3];
      break;
    case 'iv2':
      i32[o] = v[s];
      i32[o + 1] = v[s + 1];
      break;
    case 'iv3':
      i32[o] = v[s];
      i32[o + 1] = v[s + 1];
      i32[o + 2] = v[s + 2];
      break;
    case 'iv4':
      i32[o] = v[s];
      i32[o + 1] = v[s + 1];
      i32[o + 2] = v[s + 2];
      i32[o + 3] = v[s + 3];
      break;
    case 'm3':
      for (let c = 0; c < 3; ++c) {
        f32[o + c * 4] = v[s + c * 3];
        f32[o + c * 4 + 1] = v[s + c * 3 + 1];
        f32[o + c * 4 + 2] = v[s + c * 3 + 2];
      }
      break;
    case 'm4':
      for (let i = 0; i < 16; ++i) f32[o + i] = v[s + i];
      break;
  }
}

function getMemberWriter(m: UniformBlockMember, { f32, i32 }: UniformBlockData): (value: any) => void {
  const o = m.offset >> 2;
  if (m.length) {
    const kind = getElementKind(m.kind);
    const itemSize = getElementItemSize(kind);
    const stride = m.stride >> 2;
    const { length } = m;
    return (v) => {
      const n = Math.min(length, Math.floor(v.length / itemSize));
      for (let i = 0; i < n; ++i) packElement(f32, i32, kind, o + i * stride, v, i * itemSize);
    };
  }
  switch (m.kind) {
    case 'f':
      return (v) => {
        f32[o] = v;
      };
    case 'i':
      return (v) => {
        i32[o] = v;
      };
    case 'b':
      return (v) => {
        i32[o] = v ? 1 : 0;
      };
    default:
      return (v) => packElement(f32, i32, m.kind, o, v, 0);
  }
}

/**
 * Create a writer that copies the given uniform values into `data`.
 * Each call only copies values whose version changed and calls `onWrite` if anything was written.
 * Cell instances are resolved once, keep using the same writer for them.
 */
/** Writes the changed values into the block data; `dispose` stops observing the cells */
export type UniformBlockWriter = (() => void) & { dispose: () => void };

export function createUniformBlockWriter(
  layout: UniformBlockLayout,
  data: UniformBlockData,
  uniformValues: UniformsList,
  onWrite: () => void,
): UniformBlockWriter {
  const memberByName = new Map<string, UniformBlockMember>();
  for (const m of layout.members) memberByName.set(m.name, m);

  // resolve names, cells and member writers once so the returned
  // writer only does version checks in the per-draw hot path
  const cells: ValueCell<any>[] = [];
  const writers: ((value: any) => void)[] = [];
  for (let i = 0, il = uniformValues.length; i < il; ++i) {
    const [k, v] = uniformValues[i];
    if (!v) continue;
    const m = memberByName.get(k);
    if (m === undefined) continue;
    cells.push(v);
    writers.push(getMemberWriter(m, data));
  }
  const count = cells.length;
  // -1 forces the initial write, cell versions are non-negative
  const versions = new Int32Array(count).fill(-1);
  // the versions are checked only after a cell changed (most draws have none), or every time with volatile cells
  const observer: ValueCell.Observer = { changed: true };
  for (const cell of cells) ValueCell.observe(cell, observer);
  const volatile = cells.some((c) => c.volatile);
  const write = () => {
    if (!observer.changed && !volatile) return;
    observer.changed = false;
    let written = false;
    for (let i = 0; i < count; ++i) {
      const ref = cells[i].ref;
      if (ref.version !== versions[i]) {
        writers[i](ref.value);
        versions[i] = ref.version;
        written = true;
      }
    }
    if (written) onWrite();
  };
  return Object.assign(write, {
    dispose: () => {
      for (const cell of cells) ValueCell.unobserve(cell, observer);
    },
  });
}
