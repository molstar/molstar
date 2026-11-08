/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { type Mat4, Vec4 } from '@molstar/core/math/linear-algebra';
import { transformPositionArray, type GroupMapping } from '../../util.js';
import type { GeometryUtils } from '../geometry.js';
import { createColors } from '../color-data.js';
import { createMarkers } from '../marker-data.js';
import { createSizes } from '../size-data.js';
import type { TransformData } from '../transform-data.js';
import type { LocationIterator } from '../../util/location-iterator.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type DataArray, calculateTransformBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { PointsValues } from '@molstar/graphics/gl/renderable/points';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import type { Color } from '@molstar/core/util/color';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import { createEmptyClipping } from '../clipping-data.js';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';
import {
  createCenterGroupData,
  createPositionGroupIterator,
  interleavePositionGroup,
  PositionGroupStride,
  setPositionGroup,
} from '../center-group.js';
import * as Lod from '../lod-levels.js';

/** Point cloud */
export interface Points {
  readonly kind: 'points';

  /** Number of vertices in the point cloud */
  pointCount: number;

  /** Center xyz and group per point (`PositionGroupStride` floats), the only copy (the shaders read it as `tPositionGroup`) */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the points */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to point indices times `PositionGroupStride` (offsets into `positionGroup`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;

  shaderData: Points.ShaderData;
}

export namespace Points {
  /** Data read by the shaders, points are drawn as quads of 6 vertices that pull their center and group */
  export interface ShaderData {
    /** center xyz and group per point, in level of detail order where levels are set */
    readonly positionGroup: ValueCell<DataArray<Float32Array>>;
    readonly lodLevels: ValueCell<Lod.LodLevelsValue>;
    /** orders the points into the levels of detail `lodLevels` (the current ones by default) */
    update(props?: { lodLevels: Lod.LodLevels }): void;
  }

  /** From separate `centers` (xyz) and `groups`, interleaved (`fromPositionGroup` takes them as they are stored) */
  export function create(centers: Float32Array, groups: Float32Array, pointCount: number, points?: Points): Points {
    return fromPositionGroup(interleavePositionGroup(centers, groups, pointCount), pointCount, points);
  }

  /** From center xyz and group per point (`PositionGroupStride` floats), kept as the geometry's array */
  export function fromPositionGroup(positionGroup: Float32Array, pointCount: number, points?: Points): Points {
    return points ? update(positionGroup, pointCount, points) : fromArray(positionGroup, pointCount);
  }

  export function createEmpty(points?: Points): Points {
    return fromPositionGroup(points ? points.positionGroup.ref.value.array : new Float32Array(0), 0, points);
  }

  function fromArray(positionGroup: Float32Array, pointCount: number): Points {
    const data = createCenterGroupData(positionGroup, pointCount, () => points.pointCount);
    const lodLevels = ValueCell.create([] as Lod.LodLevelsValue);

    const points = {
      kind: 'points' as const,
      pointCount,
      positionGroup: data.positionGroup,
      get boundingSphere() {
        return data.getBoundingSphere();
      },
      get groupMapping() {
        return data.getGroupMapping();
      },
      setBoundingSphere: data.setBoundingSphere,
      hasBoundingSphere: data.hasBoundingSphere,
      shaderData: {
        positionGroup: data.positionGroup,
        lodLevels,
        update(props?: { lodLevels: Lod.LodLevels }) {
          const lodLevelsProp = props?.lodLevels ?? Lod.getLodLevels(lodLevels.ref.value);
          // strides as given: a point's size is in pixels, not adjusted for `sizeFactor` as spheres' strides are
          ValueCell.update(lodLevels, Lod.orderLodLevels(data.positionGroup, points.pointCount, lodLevelsProp, 1));
        },
      },
    };
    return points;
  }

  function update(positionGroup: Float32Array, pointCount: number, points: Points) {
    points.pointCount = pointCount;
    setPositionGroup(points.positionGroup, positionGroup, pointCount);
    points.shaderData.update();
    return points;
  }

  export function transform(points: Points, t: Mat4) {
    const pg = points.positionGroup.ref.value;
    transformPositionArray(t, pg.array, 0, points.pointCount, PositionGroupStride);
    ValueCell.update(points.positionGroup, pg);
  }

  //

  export const StyleTypes = {
    square: 'Square',
    circle: 'Circle',
    fuzzy: 'Fuzzy',
    ring: 'Ring',
  };
  export type StyleTypes = keyof typeof StyleTypes;
  export const StyleTypeNames = Object.keys(StyleTypes) as StyleTypes[];

  export const Params = {
    ...BaseGeometry.Params,
    sizeFactor: PD.Numeric(3, { min: 0, max: 10, step: 0.1 }),
    pointSizeAttenuation: PD.Boolean(false),
    pointStyle: PD.Select('square', PD.objectToOptions(StyleTypes)),
    /** A ring point's ring: its width as a fraction of the radius (1 the full disk) */
    pointRingWidth: PD.Numeric(
      0.25,
      { min: 0, max: 1, step: 0.01 },
      {
        description: "The ring's width as a fraction of the point's radius (1: the full disk).",
        hideIf: (p) => p.pointStyle !== 'ring',
      },
    ),
    /**
     * Each point twice its radius toward the camera, along the view ray (as large and where it was on screen): in front
     * of what meets it, drawn and picked (a gizmo's knob over the edges ending in it)
     */
    pointFront: PD.Boolean(false, { isHidden: true }),
    /**
     * CSS pixels toward the camera (negative: away), in depth only (as large and where it was on screen): a length of
     * that many pixels at the point's depth, so it scales with a part at a set size on screen (a gizmo's hit disk in
     * front of, or behind, its other parts)
     */
    pointDepthOffset: PD.Numeric(0, { min: -1000, max: 1000, step: 1 }, { isHidden: true }),
    animation: getAnimationParam(),
    lodLevels: Lod.getLodLevelsParam(),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<Points, Params> = {
    Params,
    createEmpty,
    createValues,
    createValuesSimple,
    updateValues,
    updateBoundingSphere,
    createRenderableState,
    updateRenderableState,
    createPositionIterator,
  };

  function createPositionIterator(points: Points, transform: TransformData): LocationIterator {
    return createPositionGroupIterator(points.positionGroup, points.pointCount, transform);
  }

  function createValues(
    points: Points,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): PointsValues {
    const { instanceCount, groupCount } = locationIt;
    // first: what is computed per vertex reads the points in this order
    points.shaderData.update({ lodLevels: props.lodLevels });
    const positionIt = createPositionIterator(points, transform);

    const color = createColors(locationIt, positionIt, theme.color);
    const size = createSizes(locationIt, positionIt, theme.size);
    const marker = resolveInstanceGranularity(props.instanceGranularity, groupCount, instanceCount)
      ? createMarkers(instanceCount, 'instance')
      : createMarkers(instanceCount * groupCount, 'groupInstance');
    const overpaint = createEmptyOverpaint();
    const transparency = createEmptyTransparency();
    const emissive = createEmptyEmissive();
    const material = createEmptySubstance();
    const clipping = createEmptyClipping();
    const wiggle = createEmptyWiggle();

    // each point is a quad of two triangles, the vertex count is the number of points (data per point)
    const counts = { drawCount: points.pointCount * 2 * 3, vertexCount: points.pointCount, groupCount, instanceCount };

    const invariantBoundingSphere = Sphere3D.clone(points.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('points'),

      tPositionGroup: points.shaderData.positionGroup,
      lodLevels: points.shaderData.lodLevels,
      boundingSphere: ValueCell.create(boundingSphere),
      invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
      uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
      ...color,
      ...size,
      ...marker,
      ...overpaint,
      ...transparency,
      ...emissive,
      ...material,
      ...clipping,
      ...wiggle,
      ...transform,

      ...BaseGeometry.createValues(props, counts),
      uSizeFactor: ValueCell.create(props.sizeFactor),
      dPointSizeAttenuation: ValueCell.create(props.pointSizeAttenuation),
      dPointStyle: ValueCell.create(props.pointStyle),
      uPointRingWidth: ValueCell.create(props.pointRingWidth),
      uPointDepthOffset: ValueCell.create(props.pointDepthOffset),
      dPointFront: ValueCell.create(props.pointFront),
      ...createAnimationValues(props.animation),
    };
  }

  function createValuesSimple(
    points: Points,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(points, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: PointsValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
    ValueCell.updateIfChanged(values.dPointSizeAttenuation, props.pointSizeAttenuation);
    ValueCell.updateIfChanged(values.dPointStyle, props.pointStyle);
    ValueCell.updateIfChanged(values.dPointFront, props.pointFront);
    ValueCell.updateIfChanged(values.uPointRingWidth, props.pointRingWidth);
    ValueCell.updateIfChanged(values.uPointDepthOffset, props.pointDepthOffset);
    updateAnimationValues(values, props.animation);

    const lodLevels = values.lodLevels.ref.value as Lod.LodLevelsValue;
    if (!Lod.areLodLevelsEqual(props.lodLevels, Lod.getLodLevels(lodLevels))) {
      const count = values.uVertexCount.ref.value;
      ValueCell.update(values.lodLevels, Lod.orderLodLevels(values.tPositionGroup, count, props.lodLevels, 1));
    }
  }

  function updateBoundingSphere(values: PointsValues, points: Points) {
    const invariantBoundingSphere = Sphere3D.clone(points.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      values.aTransform.ref.value,
      values.instanceCount.ref.value,
      0,
    );

    if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
      ValueCell.update(values.boundingSphere, boundingSphere);
    }
    if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
      ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
      ValueCell.update(
        values.uInvariantBoundingSphere,
        Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere),
      );
    }
  }

  function createRenderableState(props: PD.Values<Params>): RenderableState {
    const state = BaseGeometry.createRenderableState(props);
    updateRenderableState(state, props);
    return state;
  }

  function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
    BaseGeometry.updateRenderableState(state, props);
    state.opaque = state.opaque && props.pointStyle !== 'fuzzy';
    state.writeDepth = state.opaque;
  }
}
