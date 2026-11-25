/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Geometry, GeometryUtils } from '../geometry/geometry.js';
import type { RuntimeContext } from '@molstar/core/task';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Shape } from '@molstar/model/model/shape';
import type { GraphicsContext } from '../../gl/context.js';

export interface ShapeProvider<D, G extends Geometry, P extends Geometry.Params<G>> {
  label: string;
  data: D;
  params: P;
  getShape: ShapeGetter<D, G, P>;
  geometryUtils: GeometryUtils<G>;
  /**
   * A group's transparency (0 opaque, 1 invisible), for a shape whose parts differ (a translucent wall, its solid
   * edges): applied by the shape representation as a transparency layer on the groups that have one
   */
  transparency?: (group: number) => number;
  /** Releases what a shape's geometry owns on the GPU (e.g. SOG's textures), when its representation is destroyed */
  dispose?: (geometry: G) => void;
}

/**
 * The shape of `data`; `graphics` the context it is drawn with, for geometry made on it (GPU textures: the getter is
 * called again with a new context where `GeometryUtils.isCompatible` says the shape's geometry is not of it)
 */
export type ShapeGetter<D, G extends Geometry, P extends Geometry.Params<G>> = (
  ctx: RuntimeContext,
  data: D,
  props: PD.Values<P>,
  shape?: Shape<G>,
  graphics?: GraphicsContext,
) => Shape<G> | Promise<Shape<G>>;
