/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { type ComputeContext, type StorageBuffer, StorageBufferSpec } from './compute.js';
import type { GraphicsContext } from './context.js';
import { CulledData, type CulledSegments } from './culled-data.js';
import { type AttributeSpec, DefineSpec, type RenderableSchema } from './renderable/schema.js';
import { getTransformParts, hasInstancePoses, type InstancePosesValue } from './transform-parts.js';

/**
 * Objects of poses drawn without GPU culling on WebGPU (`dInstanceStorage`): their instanced attributes read from
 * storage (`sInstanceData`, `CulledData`) at `instance_index` (the draw's first instance included: the grid's order,
 * as the attributes), not packed as rows into vertex buffers (`affineRowsCell`); the reader of the culled programs
 * (`webgpu/culled.ts`) without their indices
 */
export const InstanceStorageDefine = 'dInstanceStorage';

export const InstanceStorageSchema = {
  sInstanceData: StorageBufferSpec('vec4<f32>', 'read', [InstanceStorageDefine]),
  dInstanceStorage: DefineSpec('boolean'),
} as const;

/** The instanced attributes read from `sInstanceData` where it is used */
const StorageAttributes = ['aTransform', 'aInstance'] as const;

/**
 * `schema` with `InstanceStorageSchema`, its instanced attributes (`StorageAttributes`, `extra` besides) unused with
 * `dInstanceStorage` (and with `and`, e.g. a culled variant's own define: unused with either)
 */
export function withInstanceStorage<S extends RenderableSchema>(
  schema: S,
  extra: readonly string[] = [],
  and?: string,
): S & typeof InstanceStorageSchema {
  const out: { [k: string]: RenderableSchema[string] } = { ...schema, ...InstanceStorageSchema };
  const condition = and ? `!${and}&!${InstanceStorageDefine}` : `!${InstanceStorageDefine}`;
  for (const name of [...StorageAttributes, ...extra]) {
    const spec = out[name] as AttributeSpec<never> | undefined;
    if (spec) out[name] = { ...spec, usedIf: [condition] };
  }
  return out as S & typeof InstanceStorageSchema;
}

/** Whether `ctx` reads the instances of objects with poses from storage (WebGPU) */
export function canUseInstanceStorage(ctx: GraphicsContext | undefined) {
  return ctx?.capabilities.backend === 'webgpu' && !!ctx.compute;
}

/** Whether an object of `values` in `ctx` reads its instances from storage: it has poses (`instancePoses`) */
export function usesInstanceStorage(
  ctx: GraphicsContext | undefined,
  values: { readonly instancePoses?: ValueCell<InstancePosesValue> },
) {
  return canUseInstanceStorage(ctx) && hasInstancePoses(values);
}

export interface InstanceStorageValues {
  readonly instanceCount: ValueCell<number>;
  readonly aTransform: ValueCell<unknown>;
  readonly aInstance: ValueCell<unknown>;
  readonly instancePoses?: ValueCell<InstancePosesValue>;
}

export interface InstanceStorage {
  /** the values the programs read (`InstanceStorageSchema`) */
  readonly values: { sInstanceData: ValueCell<StorageBuffer>; dInstanceStorage: ValueCell<boolean> };
  /** writes the instances where they changed */
  sync(): void;
  destroy(): void;
}

/**
 * The instances' data of `values` (`CulledData`: matrices as rows, poses by index) in a buffer of its own, rewritten
 * where the transforms or `aInstance` changed; `getSegments` a merged object's (their groups), else a segment a part
 */
export function createInstanceStorage(
  compute: ComputeContext,
  values: InstanceStorageValues,
  getSegments?: () => { segments: CulledSegments; groups: Float32Array },
): InstanceStorage {
  const data = ValueCell.create(compute.createBuffer(16, 'instance data'));
  let versions = '';
  function sync() {
    const key = `${values.aTransform.ref.version},${values.instancePoses?.ref.version},${values.aInstance.ref.version},${values.instanceCount.ref.value}`;
    if (key === versions) return;
    versions = key;
    const count = values.instanceCount.ref.value;
    const parts = getTransformParts(values, count);
    const merged = getSegments?.();
    const segments = merged?.segments ?? {
      count: parts.length,
      instanceBases: parts.map((p) => p.first),
      instanceCounts: parts.map((p) => p.transforms.count),
    };
    const layout = CulledData.layout(
      parts,
      values.aInstance.ref.value as Uint32Array,
      segments,
      merged?.groups,
      count,
      1 << 18,
    );
    if (data.ref.value.byteLength < layout.vec4s * 16) {
      data.ref.value.destroy();
      ValueCell.update(data, compute.createBuffer(layout.vec4s * 16, 'instance data'));
    }
    const buffer = data.ref.value;
    layout.write((d, o) => buffer.write(d, o));
  }
  sync();
  return {
    values: { sInstanceData: data, dInstanceStorage: ValueCell.create(true) },
    sync,
    destroy: () => data.ref.value.destroy(),
  };
}
