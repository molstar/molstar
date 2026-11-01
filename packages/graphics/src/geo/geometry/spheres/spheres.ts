/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from '@molstar/core/util';
import type { GeometryUtils } from '../geometry.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { TransformData } from '../transform-data.js';
import type { LocationIterator } from '../../util/location-iterator.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import { createSphereResolveValues, type SpheresValues } from '@molstar/graphics/gl/renderable/spheres';
import { createColors } from '../color-data.js';
import {
  createCenterGroupData,
  createPositionGroupIterator,
  interleavePositionGroup,
  setPositionGroup,
} from '../center-group.js';
import { createMarkers } from '../marker-data.js';
import { type DataArray, calculateTransformBoundingSphere } from '@molstar/graphics/gl/renderable/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { createSizes, getMaxSize, type SizeData } from '../size-data.js';
import type { Color } from '@molstar/core/util/color';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import type { GroupMapping } from '../../util.js';
import { createEmptyClipping } from '../clipping-data.js';
import { Vec4 } from '@molstar/core/math/linear-algebra';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { createInteriorValues, getInteriorParam, updateInteriorValues } from '../interior.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';
import * as Lod from '../lod-levels.js';

export interface Spheres {
  readonly kind: 'spheres';

  /** Number of spheres */
  sphereCount: number;

  /**
   * Center xyz and group per sphere (`PositionGroupStride` floats), the only copy (the shaders read it as
   * `tPositionGroup`): in level of detail order where levels are set, the order has no meaning beyond the geometry
   */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the spheres */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to sphere indices times `PositionGroupStride` (offsets into `positionGroup`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;

  shaderData: Spheres.ShaderData;
}

export namespace Spheres {
  export interface ShaderData {
    readonly positionGroup: ValueCell<DataArray<Float32Array>>;
    readonly lodLevels: ValueCell<LodLevelsValue>;
    readonly sizeFactor: ValueCell<number>;
    update(props?: LodProps): void;
  }

  /** From separate `centers` (xyz) and `groups`, interleaved (`fromPositionGroup` takes them as they are stored) */
  export function create(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
    return fromPositionGroup(interleavePositionGroup(centers, groups, sphereCount), sphereCount, spheres);
  }

  /** From center xyz and group per sphere (`PositionGroupStride` floats), kept as the geometry's array */
  export function fromPositionGroup(positionGroup: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
    return spheres ? update(positionGroup, sphereCount, spheres) : fromArray(positionGroup, sphereCount);
  }

  export function createEmpty(spheres?: Spheres): Spheres {
    return fromPositionGroup(spheres ? spheres.positionGroup.ref.value.array : new Float32Array(0), 0, spheres);
  }

  function fromArray(positionGroup: Float32Array, sphereCount: number): Spheres {
    const data = createCenterGroupData(positionGroup, sphereCount, () => spheres.sphereCount);
    const lodLevels = ValueCell.create([] as LodLevelsValue);
    const sizeFactor = ValueCell.create(0);

    const spheres = {
      kind: 'spheres' as const,
      sphereCount,
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
        sizeFactor,
        update(props?: LodProps) {
          const lodLevelsProp = props?.lodLevels ?? getLodLevels(lodLevels.ref.value);
          const sizeFactorProp = props?.sizeFactor ?? sizeFactor.ref.value;
          const pixelSizeProp = props?.lodPixelSize ?? lodLevels.ref.value.pixelSize ?? 0;
          const radiusProp = props?.lodRadius ?? lodLevels.ref.value.radius ?? LodNominalRadius;

          const newLodLevels = Lod.orderLodLevels(
            data.positionGroup,
            spheres.sphereCount,
            lodLevelsProp,
            sizeFactorProp,
            pixelSizeProp,
            radiusProp,
          );

          ValueCell.update(lodLevels, newLodLevels);
          ValueCell.update(sizeFactor, sizeFactorProp);
        },
      },
    };
    // no fill here: `createValues` fills the position-group texture with the real LOD and size (perf-audit-2 08/P11)
    return spheres;
  }

  function update(positionGroup: Float32Array, sphereCount: number, spheres: Spheres) {
    spheres.sphereCount = sphereCount;
    setPositionGroup(spheres.positionGroup, positionGroup, sphereCount);
    spheres.shaderData.update();
    return spheres;
  }

  // the levels of detail (`lod-levels.ts`), shared with points
  export type LodLevels = Lod.LodLevels;
  export type LodLevelsValue = Lod.LodLevelsValue;
  export const LodLevelsPresets = Lod.LodLevelsPresets;
  export const getLodPixelBoundary = Lod.getLodPixelBoundary;
  const { LodNominalRadius, areLodLevelsEqual, getLodLevels } = Lod;

  /** `lodRadius` the radius the pixel boundaries assume (`getLodRadius`) */
  type LodProps = { lodLevels: LodLevels; sizeFactor: number; lodPixelSize?: number; lodRadius?: number };

  /**
   * The radius the pixel boundaries assume at size factor 1: `lodRadius`, or with 0 the largest of the size data (dots,
   * sized by their theme)
   */
  function getLodRadius(lodRadius: number, size: SizeData) {
    return lodRadius > 0 ? lodRadius : getMaxSize(size) || LodNominalRadius;
  }

  /** The `lodRadius` param: hidden, the nominal atom radius for spheres and 0 (from the size data) for dots */
  export function getLodRadiusParam(defaultValue = LodNominalRadius) {
    return PD.Numeric(
      defaultValue,
      { min: 0, max: 100, step: 0.1 },
      {
        ...BaseGeometry.CullingLodCategory,
        isHidden: true,
        description:
          'The sphere radius levels of detail by size on screen are chosen for, at size factor 1; 0 for the largest size.',
      },
    );
  }

  /**
   * Strides growing by 10 (scale biases 1, 3, 2.7, 2.3): levels by size on screen that reach a whole mesoscale scene
   * with four levels and look as `quality` at 3 px (`screen-space-lod.md`)
   */
  const LodPixelChain = [
    [1, 1],
    [10, 3],
    [100, 2.7],
    [1000, 2.3],
  ] as const;

  /** As distances at 1280x720 and a field of view of 45°, the levels by size on screen at `pixelSize` */
  const LodPixelReferenceDepthScale = (2 * Math.tan(Math.PI / 8)) / 720;

  /** The chain's levels for `pixelSize`, their distances those at the reference view (for `lodPixelSize` 0) */
  function getLodPixelLevels(pixelSize: number): LodLevels {
    let minDistance = 1;
    return LodPixelChain.map(([stride, scaleBias], j) => {
      const level = { minDistance, maxDistance: 10000000, overlap: 0, stride, scaleBias };
      if (j < LodPixelChain.length - 1) {
        level.maxDistance = Math.round(getLodPixelBoundary(level, pixelSize) / LodPixelReferenceDepthScale / 10) * 10;
      }
      minDistance = level.maxDistance;
      return level;
    });
  }

  /** Levels by size on screen (`lodPixelSize`) for the graphics modes of `LodLevelsPresets` */
  export const LodPixelPresets: {
    [key in keyof typeof LodLevelsPresets]: { lodLevels: LodLevels; lodPixelSize: number };
  } = {
    performance: { lodLevels: getLodPixelLevels(6), lodPixelSize: 6 },
    balanced: { lodLevels: getLodPixelLevels(4), lodPixelSize: 4 },
    quality: { lodLevels: getLodPixelLevels(3), lodPixelSize: 3 },
    ultra: { lodLevels: getLodPixelLevels(1.5), lodPixelSize: 1.5 },
  };

  /** Levels of detail by size on screen (`screen-space-lod.md`): the smallest sphere drawn, in CSS pixels; 0 for distances */
  export function getLodPixelSizeParam(defaultValue = 0) {
    return PD.Numeric(
      defaultValue,
      { min: 0, max: 10, step: 0.5 },
      {
        ...BaseGeometry.CullingLodCategory,
        description:
          'If above zero, the levels of detail are chosen by the size of their spheres on screen: the smallest sphere drawn, in CSS pixels. The distances of the levels are ignored.',
      },
    );
  }

  export const Params = {
    ...BaseGeometry.Params,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    doubleSided: PD.Boolean(false, BaseGeometry.CustomQualityParamInfo),
    ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
    celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
    xrayShaded: PD.Select<boolean | 'inverted'>(
      false,
      [
        [false, 'Off'],
        [true, 'On'],
        ['inverted', 'Inverted'],
      ],
      BaseGeometry.ShadingCategory,
    ),
    transparentBackfaces: PD.Select(
      'off',
      PD.arrayToOptions(['off', 'on', 'opaque'] as const),
      BaseGeometry.ShadingCategory,
    ),
    solidInterior: PD.Boolean(true, BaseGeometry.ShadingCategory),
    clipPrimitive: PD.Boolean(false, {
      ...BaseGeometry.ShadingCategory,
      description: 'Clip whole sphere instead of cutting it.',
    }),
    approximate: PD.Boolean(false, {
      ...BaseGeometry.ShadingCategory,
      description: 'Faster rendering, but has artifacts.',
    }),
    alphaThickness: PD.Numeric(
      0,
      { min: 0, max: 20, step: 1 },
      { ...BaseGeometry.ShadingCategory, description: 'If not zero, adjusts alpha for radius.' },
    ),
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    interior: getInteriorParam(),
    animation: getAnimationParam(),
    lodLevels: Lod.getLodLevelsParam(),
    lodPixelSize: getLodPixelSizeParam(),
    lodRadius: getLodRadiusParam(),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<Spheres, Params> = {
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

  function createPositionIterator(spheres: Spheres, transform: TransformData): LocationIterator {
    return createPositionGroupIterator(spheres.positionGroup, spheres.sphereCount, transform);
  }

  function createValues(
    spheres: Spheres,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): SpheresValues {
    const { instanceCount, groupCount } = locationIt;
    // the sizes' largest for the levels of detail, which order the spheres: what is computed per vertex reads that order
    const size = createSizes(locationIt, createPositionIterator(spheres, transform), theme.size);
    spheres.shaderData.update({
      lodLevels: props.lodLevels,
      sizeFactor: props.sizeFactor,
      lodPixelSize: props.lodPixelSize,
      lodRadius: getLodRadius(props.lodRadius, size),
    });
    const positionIt = createPositionIterator(spheres, transform);
    if (theme.size.granularity.startsWith('vertex')) createSizes(locationIt, positionIt, theme.size, size);
    const color = createColors(locationIt, positionIt, theme.color);
    const marker = resolveInstanceGranularity(props.instanceGranularity, groupCount, instanceCount)
      ? createMarkers(instanceCount, 'instance')
      : createMarkers(instanceCount * groupCount, 'groupInstance');
    const overpaint = createEmptyOverpaint();
    const transparency = createEmptyTransparency();
    const emissive = createEmptyEmissive();
    const material = createEmptySubstance();
    const clipping = createEmptyClipping();
    const wiggle = createEmptyWiggle();

    const counts = {
      drawCount: spheres.sphereCount * 2 * 3,
      vertexCount: spheres.sphereCount,
      groupCount,
      instanceCount,
    };

    // not by the radius: that of a single sphere is 0, or a rounding error
    const padding = spheres.sphereCount > 0 ? getMaxSize(size) * props.sizeFactor : 0;
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('spheres'),

      tPositionGroup: spheres.shaderData.positionGroup,

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

      padding: ValueCell.create(padding),

      ...BaseGeometry.createValues(props, counts),
      uSizeFactor: spheres.shaderData.sizeFactor,
      uDoubleSided: ValueCell.create(props.doubleSided),
      dIgnoreLight: ValueCell.create(props.ignoreLight),
      dCelShaded: ValueCell.create(props.celShaded),
      dXrayShaded: ValueCell.create(
        props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
      ),
      dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
      dSolidInterior: ValueCell.create(props.solidInterior),
      dClipPrimitive: ValueCell.create(props.clipPrimitive),
      dApproximate: ValueCell.create(props.approximate),
      uAlphaThickness: ValueCell.create(props.alphaThickness),
      dBump: ValueCell.create(BaseGeometry.hasBump(props)),
      uBumpFrequency: ValueCell.create(props.bumpFrequency),
      uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

      lodLevels: spheres.shaderData.lodLevels,

      ...createInteriorValues(props.interior),
      ...createAnimationValues(props.animation),
      ...createSphereResolveValues(),
    };
  }

  function createValuesSimple(
    spheres: Spheres,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(spheres, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: SpheresValues, props: PD.Values<Params>) {
    const prevSizeFactor = values.uSizeFactor.ref.value;

    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
    ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
    ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
    ValueCell.updateIfChanged(
      values.dXrayShaded,
      props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
    );
    ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
    ValueCell.updateIfChanged(values.dSolidInterior, props.solidInterior);
    ValueCell.updateIfChanged(values.dClipPrimitive, props.clipPrimitive);
    ValueCell.updateIfChanged(values.dApproximate, props.approximate);
    ValueCell.updateIfChanged(values.uAlphaThickness, props.alphaThickness);
    ValueCell.updateIfChanged(values.dBump, BaseGeometry.hasBump(props));
    ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
    ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);
    updateInteriorValues(values, props.interior);
    updateAnimationValues(values, props.animation);

    const lodValue = values.lodLevels.ref.value as LodLevelsValue;
    const lodLevels = getLodLevels(lodValue);
    const lodRadius = getLodRadius(props.lodRadius, values);
    if (
      !areLodLevelsEqual(props.lodLevels, lodLevels) ||
      props.sizeFactor !== prevSizeFactor ||
      props.lodPixelSize !== (lodValue.pixelSize ?? 0) ||
      (props.lodPixelSize > 0 && lodRadius !== lodValue.radius)
    ) {
      const count = values.uVertexCount.ref.value;
      const lodLevels = Lod.orderLodLevels(
        values.tPositionGroup,
        count,
        props.lodLevels,
        props.sizeFactor,
        props.lodPixelSize,
        lodRadius,
      );
      ValueCell.update(values.lodLevels, lodLevels);
    }
  }

  function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
    const padding = spheres.sphereCount > 0 ? getMaxSize(values) * values.uSizeFactor.ref.value : 0;
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
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
    ValueCell.update(values.padding, padding);
  }

  function createRenderableState(props: PD.Values<Params>): RenderableState {
    const state = BaseGeometry.createRenderableState(props);
    updateRenderableState(state, props);
    return state;
  }

  function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
    BaseGeometry.updateRenderableState(state, props);
    state.opaque = state.opaque && !props.xrayShaded;
    state.writeDepth = state.opaque;
  }
}
