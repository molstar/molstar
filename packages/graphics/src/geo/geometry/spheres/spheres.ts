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
import {
  createSiteColorValues,
  createSphereResolveValues,
  type SpheresValues,
} from '@molstar/graphics/gl/renderable/spheres';
import { createColors } from '../color-data.js';
import {
  createCenterGroupData,
  createPositionGroupIterator,
  interleavePositionGroup,
  setPositionGroup,
} from '../center-group.js';
import { createMarkers } from '../marker-data.js';
import {
  type DataArray,
  calculateTransformDataBoundingSphere,
  createDataArray,
} from '@molstar/graphics/gl/renderable/util';
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
import * as Lod from '../lod-levels.js';
import {
  createSitesData,
  createSitesIterator,
  decodeSitesPositionGroup,
  getSitesDataArray,
  type SphereSites,
} from './sphere-sites.js';
import { PositionGroupStride } from '../center-group.js';
import type { GraphicsContext } from '@molstar/graphics/gl/context';

export interface Spheres {
  readonly kind: 'spheres';

  /** Number of spheres */
  sphereCount: number;

  /**
   * Center xyz and group per sphere (`PositionGroupStride` floats), the only copy (the shaders read it as
   * `tPositionGroup`), in item order: a level of detail draws every stride-th. A placeholder of one where the centers are
   * quantized (`sites`)
   */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;
  /** The quantized centers (`fromSites`, sphere `i` of group `i`), undefined where they are `positionGroup`'s */
  readonly sites: SphereSites | undefined;

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
    /** the quantized centers' data (`tSite`, a placeholder without), origin and scale; `dSites` whether drawn from them */
    readonly sites: ValueCell<DataArray<Int16Array>>;
    readonly uSiteOrigin: ValueCell<Vec3>;
    readonly uSiteScale: ValueCell<number>;
    readonly dSites: ValueCell<boolean>;
    readonly lodLevels: ValueCell<LodLevelsValue>;
    readonly sizeFactor: ValueCell<number>;
    update(props?: LodProps): void;
  }

  /** From separate `centers` (xyz) and `groups`, interleaved (`fromPositionGroup` takes them as they are stored) */
  export function create(centers: Float32Array, groups: Float32Array, sphereCount: number, spheres?: Spheres): Spheres {
    return fromPositionGroup(interleavePositionGroup(centers, groups, sphereCount), sphereCount, spheres);
  }

  /**
   * From center xyz and group per sphere (`PositionGroupStride` floats), kept as the geometry's array; `base`: the
   * spheres every `base`-th of a whole's (coarse sites, `ParticleSites.stride`), their levels of detail of the whole, as
   * quantized sites' (`SphereSites.stride`)
   */
  export function fromPositionGroup(
    positionGroup: Float32Array,
    sphereCount: number,
    spheres?: Spheres,
    base = 1,
  ): Spheres {
    return spheres && !spheres.sites
      ? update(positionGroup, sphereCount, spheres, base)
      : fromData(positionGroup, undefined, sphereCount, base);
  }

  /**
   * From quantized centers (`SphereSites`: sphere `i` of group `i`), drawn from their array as it is (no copy; WebGL2 and
   * WebGPU, `isSitesSupported`)
   */
  export function fromSites(sites: SphereSites, sphereCount: number, spheres?: Spheres): Spheres {
    return spheres?.sites ? updateSites(sites, sphereCount, spheres) : fromData(undefined, sites, sphereCount);
  }

  /** Whether `graphics` draws quantized centers (`fromSites`): integer data, not WebGL1 */
  export function isSitesSupported(graphics: GraphicsContext | undefined) {
    return !!graphics && graphics.capabilities.backend !== 'webgl';
  }

  export function createEmpty(spheres?: Spheres): Spheres {
    const array = spheres && !spheres.sites ? spheres.positionGroup.ref.value.array : new Float32Array(0);
    return fromPositionGroup(array, 0, spheres);
  }

  function fromData(
    positionGroup: Float32Array | undefined,
    sites: SphereSites | undefined,
    sphereCount: number,
    base = 1,
  ): Spheres {
    const sitesCell = sites ? ValueCell.create(sites) : undefined;
    // the levels' base of centers not quantized (`fromPositionGroup`'s), quantized ones' their stride
    const baseCell = ValueCell.create(base);
    const getCount = () => spheres.sphereCount;
    const centerGroup = positionGroup ? createCenterGroupData(positionGroup, sphereCount, getCount) : undefined;
    const data = centerGroup ?? createSitesData(sitesCell!, getCount);
    const positionGroupCell: ValueCell<DataArray<Float32Array>> =
      centerGroup?.positionGroup ?? ValueCell.create(createDataArray(1, PositionGroupStride, Float32Array));
    const lodLevels = ValueCell.create([] as LodLevelsValue);
    const sizeFactor = ValueCell.create(0);

    const spheres: Spheres & { sitesCell?: ValueCell<SphereSites>; baseCell: ValueCell<number> } = {
      kind: 'spheres' as const,
      sphereCount,
      positionGroup: positionGroupCell,
      get sites() {
        return sitesCell?.ref.value;
      },
      get boundingSphere() {
        return data.getBoundingSphere();
      },
      get groupMapping() {
        return data.getGroupMapping();
      },
      setBoundingSphere: data.setBoundingSphere,
      hasBoundingSphere: data.hasBoundingSphere,
      sitesCell,
      baseCell,
      shaderData: {
        positionGroup: positionGroupCell,
        sites: ValueCell.create(sites ? getSitesDataArray(sites, sphereCount) : createDataArray(1, 4, Int16Array)),
        uSiteOrigin: ValueCell.create(sites ? Vec3.clone(sites.origin) : Vec3()),
        uSiteScale: ValueCell.create(sites?.scale ?? 1),
        dSites: ValueCell.create(!!sites),
        lodLevels,
        sizeFactor,
        update(props?: LodProps) {
          const sizeFactorProp = props?.sizeFactor ?? sizeFactor.ref.value;
          const pixelSizeProp = props?.lodPixelSize ?? lodLevels.ref.value.pixelSize ?? 0;
          const radiusProp = props?.lodRadius ?? lodLevels.ref.value.radius ?? LodNominalRadius;
          const base = sitesCell ? (sitesCell.ref.value.stride ?? 1) : baseCell.ref.value;
          const value = Lod.getLodPixelLevelsValue(
            spheres.sphereCount,
            pixelSizeProp,
            sizeFactorProp,
            radiusProp,
            base,
          );
          ValueCell.update(lodLevels, value);
          ValueCell.update(sizeFactor, sizeFactorProp);
        },
      },
    };
    // no fill here: `createValues` fills the position-group texture with the real LOD and size (perf-audit-2 08/P11)
    return spheres;
  }

  function update(positionGroup: Float32Array, sphereCount: number, spheres: Spheres, base: number) {
    spheres.sphereCount = sphereCount;
    ValueCell.updateIfChanged((spheres as unknown as { baseCell: ValueCell<number> }).baseCell, base);
    setPositionGroup(spheres.positionGroup, positionGroup, sphereCount);
    spheres.shaderData.update();
    return spheres;
  }

  function updateSites(sites: SphereSites, sphereCount: number, spheres: Spheres) {
    spheres.sphereCount = sphereCount;
    const { shaderData } = spheres;
    ValueCell.update((spheres as unknown as { sitesCell: ValueCell<SphereSites> }).sitesCell, sites);
    ValueCell.update(shaderData.sites, getSitesDataArray(sites, sphereCount));
    ValueCell.update(shaderData.uSiteOrigin, Vec3.copy(shaderData.uSiteOrigin.ref.value, sites.origin));
    ValueCell.updateIfChanged(shaderData.uSiteScale, sites.scale);
    shaderData.update();
    return spheres;
  }

  // the levels of detail (`lod-levels.ts`), shared with points
  export type LodLevels = Lod.LodLevels;
  export type LodLevelsValue = Lod.LodLevelsValue;
  export const getLodPixelBoundary = Lod.getLodPixelBoundary;
  export const LodPixelPresets = Lod.LodPixelPresets;
  export const getLodPixelSizeParam = Lod.getLodPixelSizeParam;
  export const getLodPixelLevelCount = Lod.getLodPixelLevelCount;
  export type LodPixelPreset = Lod.LodPixelPreset;
  const { LodNominalRadius } = Lod;

  /** `lodRadius` the radius the pixel boundaries assume (`getLodRadius`) */
  type LodProps = { sizeFactor: number; lodPixelSize?: number; lodRadius?: number };

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

  /**
   * The centers and groups of `values`' spheres as `tPositionGroup` holds them (`PositionGroupStride` floats): its
   * array, or quantized centers (`dSites`) decoded (a transient copy)
   */
  export function getPositionGroupArray(values: SpheresValues): Float32Array {
    if (!values.dSites.ref.value) return values.tPositionGroup.ref.value.array;
    const count = values.uVertexCount.ref.value;
    return decodeSitesPositionGroup(
      values.tSite.ref.value.array,
      count,
      values.uSiteOrigin.ref.value,
      values.uSiteScale.ref.value,
    );
  }

  function createPositionIterator(spheres: Spheres, transform: TransformData): LocationIterator {
    const { sites } = spheres;
    if (sites) return createSitesIterator(sites, spheres.sphereCount, transform);
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
    const positionIt = createPositionIterator(spheres, transform);
    // the sizes' largest for the levels of detail
    const size = createSizes(locationIt, positionIt, theme.size);
    spheres.shaderData.update({
      sizeFactor: props.sizeFactor,
      lodPixelSize: props.lodPixelSize,
      lodRadius: getLodRadius(props.lodRadius, size),
    });
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
    const boundingSphere = calculateTransformDataBoundingSphere(invariantBoundingSphere, transform);

    return {
      dGeometryType: ValueCell.create('spheres'),

      tPositionGroup: spheres.shaderData.positionGroup,
      tSite: spheres.shaderData.sites,
      uSiteOrigin: spheres.shaderData.uSiteOrigin,
      uSiteScale: spheres.shaderData.uSiteScale,
      dSites: spheres.shaderData.dSites,

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
      uLodStride: ValueCell.create(1),

      ...createInteriorValues(props.interior),
      ...createAnimationValues(props.animation),
      ...createSphereResolveValues(),
      // set by a particle target of sites (`setSitesTargetSiteColor`)
      ...createSiteColorValues(),
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
    const lodRadius = getLodRadius(props.lodRadius, values);
    if (
      props.sizeFactor !== prevSizeFactor ||
      props.lodPixelSize !== (lodValue.pixelSize ?? 0) ||
      (props.lodPixelSize > 0 && lodRadius !== lodValue.radius)
    ) {
      const count = values.uVertexCount.ref.value;
      const value = Lod.getLodPixelLevelsValue(
        count,
        props.lodPixelSize,
        props.sizeFactor,
        lodRadius,
        lodValue.base ?? 1,
      );
      ValueCell.update(values.lodLevels, value);
    }
  }

  function updateBoundingSphere(values: SpheresValues, spheres: Spheres) {
    const padding = spheres.sphereCount > 0 ? getMaxSize(values) * values.uSizeFactor.ref.value : 0;
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), spheres.boundingSphere, padding);
    // the instances as before where only the padding changed (`calculateTransformDataBoundingSphere`)
    const boundingSphere = calculateTransformDataBoundingSphere(invariantBoundingSphere, values);

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
