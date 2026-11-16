/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { InstanceTransforms, type InstancePoses } from '@molstar/core/math/geometry/instance-transforms';
import { ValueCell } from '@molstar/core/util';
import { fillSerial } from '@molstar/core/util/array';
import { calcPendingInstanceGrid, type InstanceGrid } from '@molstar/core/math/geometry/instance-grid';
import { SegmentedArray } from './segmented-array.js';

/** A merged object's instances by their members' encodings: each member's at its first instance */
export interface SegmentedTransforms {
  readonly kind: 'segmented';
  readonly count: number;
  readonly parts: readonly TransformPart[];
}

/** Instances `first` to `first + transforms.count` of an object */
export interface TransformPart {
  readonly transforms: InstanceTransforms;
  readonly first: number;
}

/**
 * The value of `instancePoses`: an object's instances as poses (as drawn: the grid's order, the representation's
 * matrix applied), a merged object's members' parts where one of them has poses, none where `aTransform` holds them
 */
export type InstancePosesValue = InstancePoses | SegmentedTransforms | undefined;

/**
 * Whether an object's instances are culled and their levels of detail decimated by cells:
 * a cell size above 1 (`instanceCellSize`); 0 or 1 culls them as a whole and draws their levels at full detail. A
 * property of the object, not of its grid (which may be built later, or not at all where the GPU culls)
 */
export function hasInstanceCells(values: { readonly instanceCellSize: ValueCell<number> }) {
  return values.instanceCellSize.ref.value > 1;
}

/** The values of an object with an instance grid, maybe pending (`InstanceGrid.pending`) */
export interface InstanceGridValues extends TransformValues {
  readonly instanceGrid: ValueCell<InstanceGrid>;
  readonly aInstance: ValueCell<Uint32Array>;
}

const NoMatrices = new Float32Array(0);

/**
 * The grid a pending grid of `values` stands for (`InstanceGrid.pending`), built from their instances as they are,
 * their data not reordered: its order (`cellInstance`) for a reader drawing them by index
 */
export function buildPendingGrid(values: InstanceGridValues): InstanceGrid {
  const t = getObjectTransforms(values);
  const n = t.count;
  return calcPendingInstanceGrid(values.instanceGrid.ref.value, {
    instanceCount: n,
    instance: fillSerial(new Uint32Array(n)),
    transform: t.kind === 'matrices' ? t.matrices : NoMatrices,
    poses: t.kind === 'poses' ? t : undefined,
  });
}

/**
 * A pending grid of `values` built and applied as a built one is (`updatePosesData`): the instances in its order,
 * `aInstance` and the poses' `order` its `cellInstance`; for what draws the instances as stored (no culler reading
 * them by index: an object not merged, WebGL). Poses only (a pending grid's), else nothing
 */
export function applyPendingGrid(values: InstanceGridValues): boolean {
  const poses = values.instancePoses?.ref.value;
  if (!values.instanceGrid.ref.value.pending || poses?.kind !== 'poses') return false;
  const grid = buildPendingGrid(values);
  ValueCell.update(values.instanceGrid, grid);
  ValueCell.update(values.aInstance, grid.cellInstance);
  ValueCell.update(values.instancePoses!, { ...poses, order: grid.cellInstance });
  // a new value: what packs rows reads `aTransform`'s version
  ValueCell.update(values.aTransform, new Float32Array(0));
  return true;
}

/** The values an object's transforms are read from */
export interface TransformValues {
  readonly aTransform: ValueCell<unknown>;
  readonly instanceCount: ValueCell<number>;
  readonly instancePoses?: ValueCell<InstancePosesValue>;
  /**
   * the instances' bounding sphere the box sphere at any count (`calculateTransformDataBoundingSphere`), not the
   * boundary helper's tighter one of a few: an object of a large whole (a particle target of a list of millions,
   * thousands of targets) whose tight spheres cost more than they cull; seeded without a pass where known
   */
  readonly instanceSphereFast?: ValueCell<boolean>;
}

/**
 * The transforms of the first `count` instances of `values` as parts by their encoding: the
 * poses of `instancePoses`, else the matrices of `aTransform` (one array, or a merged object's members' as parts of a
 * `SegmentedArray`); the parts not composed into one array (millions of instances)
 */
export function getTransformParts(values: TransformValues, count: number): TransformPart[] {
  const poses = values.instancePoses?.ref.value;
  const parts: TransformPart[] = [];
  if (poses?.kind === 'segmented') {
    for (const p of poses.parts) pushPart(parts, p.transforms, p.first, count);
  } else if (poses) {
    pushPart(parts, poses, 0, count);
  } else {
    const t = values.aTransform.ref.value as Float32Array | SegmentedArray<Float32Array>;
    if (SegmentedArray.is(t)) {
      for (const p of t.parts) {
        pushPart(parts, InstanceTransforms.ofMatrices(p.array, p.count / 16), p.offset / 16, count);
      }
    } else {
      pushPart(parts, InstanceTransforms.ofMatrices(t, t.length / 16), 0, count);
    }
  }
  return parts;
}

function pushPart(parts: TransformPart[], transforms: InstanceTransforms, first: number, count: number) {
  const n = Math.min(transforms.count, count - first);
  if (n <= 0) return;
  parts.push({ transforms: n === transforms.count ? transforms : { ...transforms, count: n }, first });
}

/**
 * The first three rows (`affine-rows.ts`) of the first `count` instances of `values` (`getTransformParts`) into `out`
 * from `outOffset`, an instance each `outStride` floats
 */
export function writeAffineRowsOf(
  values: TransformValues,
  count: number,
  out: Float32Array,
  outOffset = 0,
  outStride = 12,
) {
  return writeAffineRowsRange(values, 0, count, out, outOffset, outStride);
}

/**
 * The first three rows (`affine-rows.ts`) of the instances `first` to `first + count` of `values` into `out` from
 * `outOffset`, an instance each `outStride` floats
 */
export function writeAffineRowsRange(
  values: TransformValues,
  first: number,
  count: number,
  out: Float32Array,
  outOffset = 0,
  outStride = 12,
) {
  const end = first + count;
  for (const p of getTransformParts(values, end)) {
    const a = Math.max(first, p.first);
    const b = Math.min(end, p.first + p.transforms.count);
    if (a >= b) continue;
    InstanceTransforms.writeAffineRows(
      p.transforms,
      a - p.first,
      b - a,
      out,
      outOffset + (a - first) * outStride,
      outStride,
    );
  }
  return out;
}

/** Whether `values` hold poses (`instancePoses`), not matrices alone */
export function hasInstancePoses(values: { readonly instancePoses?: ValueCell<InstancePosesValue> }) {
  return !!values.instancePoses?.ref.value;
}

/**
 * An unmerged object's instances by their encoding: its poses as drawn (`instancePoses`), else the matrices of
 * `aTransform`
 */
export function getObjectTransforms(values: TransformValues): InstanceTransforms {
  const poses = values.instancePoses?.ref.value;
  if (poses?.kind === 'poses') return poses;
  return InstanceTransforms.ofMatrices(values.aTransform.ref.value as Float32Array, values.instanceCount.ref.value);
}
