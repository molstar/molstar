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
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { SpheresValues } from '@molstar/graphics/gl/renderable/spheres';
import { createColors } from '../color-data.js';
import { createCenterGroupData, setPositionGroup } from '../center-group.js';
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
import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { createInteriorValues, getInteriorParam, updateInteriorValues } from '../interior.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';

export interface Spheres {
  readonly kind: 'spheres';

  /** Number of spheres */
  sphereCount: number;

  /** Center buffer as array of xyz values wrapped in a value cell */
  readonly centerBuffer: ValueCell<Float32Array>;
  /** Group buffer as array of group ids for each vertex wrapped in a value cell */
  readonly groupBuffer: ValueCell<Float32Array>;

  /** Bounding sphere of the spheres */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to sphere indices */
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

  export function create(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
    return spheres ? update(centers, groups, sphereCount, spheres) : fromArrays(centers, groups, sphereCount);
  }

  export function createEmpty(spheres?: Spheres): Spheres {
    const cb = spheres ? spheres.centerBuffer.ref.value : new Float32Array(0);
    const gb = spheres ? spheres.groupBuffer.ref.value : new Float32Array(0);
    return create(cb, gb, 0, spheres);
  }

  function fromArrays(centers: Float32Array, groups: Float32Array, sphereCount: number): Spheres {
    const data = createCenterGroupData(centers, groups, () => spheres.sphereCount);
    const lodLevels = ValueCell.create([] as LodLevelsValue);
    const sizeFactor = ValueCell.create(0);

    const spheres = {
      kind: 'spheres' as const,
      sphereCount,
      centerBuffer: data.centerBuffer,
      groupBuffer: data.groupBuffer,
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

          const strides = getStrides(lodLevelsProp, sizeFactorProp);
          const offsets = data.updatePositionGroup(strides);
          const newLodLevels = offsets
            ? getLodLevelsValue(lodLevelsProp, sizeFactorProp, offsets, spheres.sphereCount, pixelSizeProp, radiusProp)
            : ([] as LodLevelsValue);

          ValueCell.update(lodLevels, newLodLevels);
          ValueCell.update(sizeFactor, sizeFactorProp);
        },
      },
    };
    // no fill here: `createValues` fills the position-group texture with the real LOD and size (perf-audit-2 08/P11)
    return spheres;
  }

  function update(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres: Spheres) {
    spheres.sphereCount = sphereCount;
    ValueCell.update(spheres.centerBuffer, centers);
    ValueCell.update(spheres.groupBuffer, groups);
    spheres.shaderData.update();
    return spheres;
  }

  type LodLevels = {
    minDistance: number;
    maxDistance: number;
    overlap: number;
    stride: number;
    scaleBias: number;
  }[];

  function areLodLevelsEqual(a: LodLevels, b: LodLevels) {
    if (a.length !== b.length) return false;
    for (let i = 0, il = a.length; i < il; ++i) {
      if (a[i].maxDistance !== b[i].maxDistance) return false;
      if (a[i].minDistance !== b[i].minDistance) return false;
      if (a[i].overlap !== b[i].overlap) return false;
      if (a[i].stride !== b[i].stride) return false;
      if (a[i].scaleBias !== b[i].scaleBias) return false;
    }
    return true;
  }

  /** `lodRadius` the radius the pixel boundaries assume (`getLodRadius`) */
  type LodProps = { lodLevels: LodLevels; sizeFactor: number; lodPixelSize?: number; lodRadius?: number };

  /**
   * Per level its range (distances in Å, or with `pixelSize` above 0 the Å a CSS pixel covers), fade, vertex count and
   * sphere scale (`renderable.ts` `LodLevelsValue`), then the level's props
   */
  type LodLevelsValue = [
    min: number,
    max: number,
    overlap: number,
    count: number,
    scale: number,
    stride: number,
    scaleBias: number,
    minDistance: number,
    maxDistance: number,
  ][] & { readonly pixelSize?: number; readonly radius?: number };

  /** A carbon's van der Waals radius: the nominal sphere radius levels by size on screen are chosen for (atoms) */
  const LodNominalRadius = 1.7;

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
   * Where level `j`'s spheres get smaller than `pixelSize` CSS pixels across, in Å a pixel covers: its nominal spheres
   * (`LodNominalRadius` times the stride's scale) `pixelSize` pixels across (`screen-space-lod.md`). The preset's stride:
   * the stride adjusted for `sizeFactor` draws spheres the size of those of the preset's stride at size factor 1
   */
  export function getLodPixelBoundary(level: LodLevels[0], pixelSize: number, radius = LodNominalRadius) {
    return (2 * radius * Math.pow(Math.max(1, level.stride), 1 / level.scaleBias)) / pixelSize;
  }

  function getLodLevelsValue(
    prop: LodLevels,
    sizeFactor: number,
    offsets: number[],
    count: number,
    pixelSize: number,
    radius: number,
  ): LodLevelsValue {
    const pixel = pixelSize > 0;
    const value = prop.map((l, i): LodLevelsValue[0] => {
      const stride = getAdjustedStride(l, sizeFactor);
      return [
        pixel ? (i === 0 ? 0 : getLodPixelBoundary(prop[i - 1], pixelSize, radius)) : l.minDistance,
        pixel ? (i === prop.length - 1 ? Infinity : getLodPixelBoundary(l, pixelSize, radius)) : l.maxDistance,
        l.overlap,
        offsets[offsets.length - 1 - i],
        Math.pow(Math.min(count, stride), 1 / l.scaleBias),
        l.stride,
        l.scaleBias,
        l.minDistance,
        l.maxDistance,
      ];
    });
    return pixel ? Object.assign(value, { pixelSize, radius }) : value;
  }

  function getLodLevels(lodLevelsValue: LodLevelsValue): LodLevels {
    return lodLevelsValue.map((l) => ({
      minDistance: l[7],
      maxDistance: l[8],
      overlap: l[2],
      stride: l[5],
      scaleBias: l[6],
    }));
  }

  function getAdjustedStride(lodLevel: LodLevels[0], sizeFactor: number) {
    return Math.max(1, Math.round(lodLevel.stride / Math.pow(sizeFactor, lodLevel.scaleBias)));
  }

  function getStrides(lodLevels: LodLevels, sizeFactor: number) {
    return lodLevels.map((l) => getAdjustedStride(l, sizeFactor)).reverse();
  }

  export const LodLevelsPresets: { [key in 'performance' | 'balanced' | 'quality' | 'ultra']: LodLevels } = {
    performance: [
      { minDistance: 1, maxDistance: 300, overlap: 0, stride: 1, scaleBias: 1 },
      { minDistance: 300, maxDistance: 2000, overlap: 0, stride: 40, scaleBias: 3 },
      { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 150, scaleBias: 3 },
      { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 300, scaleBias: 2.5 },
    ],
    balanced: [
      { minDistance: 1, maxDistance: 500, overlap: 0, stride: 1, scaleBias: 1 },
      { minDistance: 500, maxDistance: 2000, overlap: 0, stride: 15, scaleBias: 3 },
      { minDistance: 2000, maxDistance: 6000, overlap: 0, stride: 70, scaleBias: 2.7 },
      { minDistance: 6000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.5 },
    ],
    quality: [
      { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1, scaleBias: 1 },
      { minDistance: 1000, maxDistance: 4000, overlap: 0, stride: 10, scaleBias: 3 },
      { minDistance: 4000, maxDistance: 10000, overlap: 0, stride: 50, scaleBias: 2.7 },
      { minDistance: 10000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2.3 },
    ],
    ultra: [
      { minDistance: 1, maxDistance: 5000, overlap: 0, stride: 1, scaleBias: 1 },
      { minDistance: 5000, maxDistance: 10000, overlap: 0, stride: 10, scaleBias: 3 },
      { minDistance: 10000, maxDistance: 30000, overlap: 0, stride: 50, scaleBias: 2.5 },
      { minDistance: 30000, maxDistance: 10000000, overlap: 0, stride: 200, scaleBias: 2 },
    ],
  };

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
    lodLevels: PD.ObjectList(
      {
        minDistance: PD.Numeric(0),
        maxDistance: PD.Numeric(0),
        overlap: PD.Numeric(0),
        stride: PD.Numeric(0),
        scaleBias: PD.Numeric(3, { min: 0.1, max: 10, step: 0.1 }),
      },
      (o) => `${o.stride}`,
      {
        ...BaseGeometry.CullingLodCategory,
        defaultValue: [] as LodLevels,
        presets: Object.entries(LodLevelsPresets).map(([k, v]) => [v, k]),
      },
    ),
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
    const groupCount = spheres.sphereCount;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    const v = spheres.centerBuffer.ref.value;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      if (instanceIndex < 0) {
        Vec3.fromArray(p, v, groupIndex * 3);
      } else {
        Vec3.transformMat4Offset(p, v, m, 0, groupIndex * 3, instanceIndex * 16);
      }
      return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
  }

  function createValues(
    spheres: Spheres,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): SpheresValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(spheres, transform);

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

    const counts = {
      drawCount: spheres.sphereCount * 2 * 3,
      vertexCount: spheres.sphereCount * 6,
      groupCount,
      instanceCount,
    };

    const padding = spheres.boundingSphere.radius ? getMaxSize(size) * props.sizeFactor : 0;
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    spheres.shaderData.update({
      lodLevels: props.lodLevels,
      sizeFactor: props.sizeFactor,
      lodPixelSize: props.lodPixelSize,
      lodRadius: getLodRadius(props.lodRadius, size),
    });

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
      uBumpFrequency: ValueCell.create(props.bumpFrequency),
      uBumpAmplitude: ValueCell.create(props.bumpAmplitude),

      lodLevels: spheres.shaderData.lodLevels,
      centerBuffer: spheres.centerBuffer,
      groupBuffer: spheres.groupBuffer,

      ...createInteriorValues(props.interior),
      ...createAnimationValues(props.animation),
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
      const count = values.uVertexCount.ref.value / 6;
      const strides = getStrides(props.lodLevels, props.sizeFactor);
      const offsets = setPositionGroup(
        values.tPositionGroup.ref.value,
        values.centerBuffer.ref.value,
        values.groupBuffer.ref.value,
        count,
        strides,
      );
      const lodLevels = offsets
        ? getLodLevelsValue(props.lodLevels, props.sizeFactor, offsets, count, props.lodPixelSize, lodRadius)
        : [];
      ValueCell.update(values.tPositionGroup, values.tPositionGroup.ref.value);
      ValueCell.update(values.lodLevels, lodLevels);
    }
  }

  function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
    const padding = spheres.boundingSphere.radius ? getMaxSize(values) * values.uSizeFactor.ref.value : 0;
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
