/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { type State, type StateObjectCell, StateSelection } from '@molstar/core/state';
import { Task } from '@molstar/core/task';
import type { Representation } from '@molstar/graphics/repr/representation';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import type { PluginContext } from '@molstar/plugin/context';
import {
  getTransformFromParams,
  type TransformParam,
  transformParamsNeedCentroid,
} from '@molstar/plugin/state/transforms/helpers';
import { getShapeProviderCenter, TransformShapeProvider } from '@molstar/plugin/state/transforms/shape/transform';
import type { ShapeEntry } from './shapes.js';

// A shape moved by hand (`TransformShapeProvider`, a decorator on the shape's provider, as a structure's conformation
// transform): its decorator and the matrix it applies, from where the shape was loaded. A shape's center is known only
// once it is made: the decorator makes it for a turn about it, the handles take its representations' bounds.

/** A shape's transform decorator (`TransformShapeProvider` below it), if any */
export function shapeTransformCell(state: State, shapeRef: string): StateObjectCell | undefined {
  return state.select(StateSelection.Generators.byRef(shapeRef).subtree().withTransformer(TransformShapeProvider))[0];
}

/** A shape's center as loaded, as the decorator turns about it (the shape made once) */
export async function shapeCenter(plugin: PluginContext, provider: ShapeProvider<any, any, any>): Promise<Vec3> {
  return plugin.runTask(Task.create('Shape center', (ctx) => getShapeProviderCenter(ctx, provider)));
}

/**
 * The matrix of a transform param, turning about the shape's center where it says so: `center` (where it was loaded)
 * if given, else as the decorator (`shapeCenter`)
 */
export async function shapeMatrixOf(
  plugin: PluginContext,
  t: TransformParam,
  provider: ShapeProvider<any, any, any>,
  center?: Vec3,
): Promise<Mat4> {
  if (!transformParamsNeedCentroid(t)) return getTransformFromParams(t, Vec3.unit);
  return getTransformFromParams(t, center ?? (await shapeCenter(plugin, provider)));
}

/**
 * The center a shape turns about, where it was loaded: as drawn where it is (`shapeDrawnCenter`, the handles' pivot),
 * else as its provider makes it (`shapeCenter`)
 */
export async function shapePivot(
  plugin: PluginContext,
  entry: ShapeEntry,
  provider: ShapeProvider<any, any, any>,
  matrix: Mat4,
): Promise<Vec3> {
  return shapeDrawnCenter(entry, matrix) ?? (await shapeCenter(plugin, provider));
}

/** A shape's transform now (identity without one); `provider` is the shape as loaded (before the decorator) */
export async function currentShapeTransform(
  plugin: PluginContext,
  shapeRef: string,
  provider: ShapeProvider<any, any, any>,
) {
  const t = shapeTransformCell(plugin.state.data, shapeRef)?.transform.params?.transform as TransformParam | undefined;
  return t ? shapeMatrixOf(plugin, t, provider) : Mat4.identity();
}

/**
 * A shape's transform now, where the handles need it at once: `shape.transform` writes a matrix, its reset a move by
 * nothing (no turn: the center not needed)
 */
export function currentShapeMatrix(state: State, shapeRef: string): Mat4 {
  const t = shapeTransformCell(state, shapeRef)?.transform.params?.transform as TransformParam | undefined;
  return t ? getTransformFromParams(t, Vec3.unit) : Mat4.identity();
}

/** A shape's center as drawn, from its representations' bounds (the union's center), back to where it was loaded */
export function shapeDrawnCenter(entry: ShapeEntry, matrix: Mat4): Vec3 | undefined {
  let center: Vec3 | undefined;
  let weight = 0;
  for (const cell of entry.representations) {
    const repr = (cell.obj?.data as { repr?: Representation.Any } | undefined)?.repr;
    for (const ro of repr?.renderObjects ?? []) {
      const { center: c, radius } = ro.values.boundingSphere.ref.value;
      const w = Math.max(radius, 1e-3);
      center = Vec3.scaleAndAdd(center ?? Vec3(), center ?? Vec3(), c, w);
      weight += w;
    }
  }
  if (!center || !weight) return undefined;
  Vec3.scale(center, center, 1 / weight);
  return Vec3.transformMat4(center, center, Mat4.invert(Mat4(), matrix));
}
