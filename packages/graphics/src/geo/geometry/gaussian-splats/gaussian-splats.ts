/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { hashFnv32a } from '@molstar/core/data/util';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Color } from '@molstar/core/util/color';
import {
  type DataArray,
  calculateTransformBoundingSphere,
  createDataArray,
  type TextureImage,
} from '@molstar/graphics/gl/renderable/util';
import type { GaussianSplatsValues } from '@molstar/graphics/gl/renderable/gaussian-splats';
import { createSplatOrder } from '@molstar/graphics/gl/renderable/gaussian-splats-order';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { StorageBuffer } from '@molstar/graphics/gl/compute';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import type { GroupMapping } from '../../util.js';
import type { LocationIterator } from '../../util/location-iterator.js';
import type { GeometryUtils } from '../geometry.js';
import type { TransformData } from '../transform-data.js';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createColors } from '../color-data.js';
import { createMarkers } from '../marker-data.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import { createEmptyClipping } from '../clipping-data.js';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import {
  int32ToFloatBits,
  createCenterGroupData,
  createPositionGroupIterator,
  PositionGroupStride,
  setPositionGroup,
} from '../center-group.js';
import { getMaxSplatScale, NoSplatSh, type SplatSh, SplatStride, SplatTexels } from './splat-data.js';
import { ShBasisCounts } from '@molstar/model/formats/shape/gaussian-splat/data';

/**
 * Gaussian splats: anisotropic 3D Gaussians with an opacity, colored by group (3D Gaussian splatting scenes), drawn as
 * one quad of 6 vertices per splat
 */
export interface GaussianSplats {
  readonly kind: 'gaussian-splats';

  /** Number of splats */
  splatCount: number;

  /** Center xyz and group per splat (`PositionGroupStride` floats), the shaders' `tPositionGroup` */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;
  /**
   * Rotation, scales and opacity per splat (`SplatStride` int16, `SplatTexels` elements, see `setSplat`), the shaders'
   * `tSplat`
   */
  readonly splatData: ValueCell<DataArray<Int16Array>>;
  /** The view-dependent color (`SplatSh`, degree 0 without), the shaders' `tSplatSh` */
  sh: SplatSh;
  readonly shData: ValueCell<DataArray<Uint8Array>>;
  /** its degree and scale, the shaders' `uSplatShDegree` and `uSplatShScale` (shared with the values: updated in place) */
  readonly shDegree: ValueCell<number>;
  readonly shScale: ValueCell<number>;

  /** `packed` (the arrays above on the GPU) or `sog` (SOG's images read in place, `sog`), the shaders' `dSplatEncoding` */
  readonly encoding: ValueCell<GaussianSplats.Encoding>;
  /** SOG's images, where the encoding is `sog` (null textures otherwise) */
  readonly sog: GaussianSplats.Sog;
  /** the GPU's `tPositionGroup`: `positionGroup` (`packed`), a placeholder (`sog`: the shaders decode the means) */
  readonly gpuPositionGroup: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the splats, their centers' padded by `BoundsSigmas` of the largest scale */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to splat indices times `PositionGroupStride` (offsets into `positionGroup`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
}

export namespace GaussianSplats {
  /**
   * Standard deviations the bounds cover past a center (along the largest axis): `exp(-½ 3²)`, 1.1% of the peak, at
   * the edge
   */
  export const BoundsSigmas = 3;

  export type Encoding = 'packed' | 'sog';

  /** SOG's images read in place (the `sog` encoding, see `createSog`), shared with the values */
  export interface Sog {
    /** the context its textures are of (undefined: none) */
    graphics: GraphicsContext | undefined;
    readonly meansLow: ValueCell<Texture>;
    readonly meansHigh: ValueCell<Texture>;
    readonly scales: ValueCell<Texture>;
    readonly quats: ValueCell<Texture>;
    readonly sh0: ValueCell<Texture>;
    readonly shLabels: ValueCell<Texture>;
    readonly shCentroids: ValueCell<Texture>;
    /** the codebooks: scales (r), base color (g), view-dependent color (b), 256 entries */
    readonly codebook: ValueCell<TextureImage<Float32Array>>;
    readonly imageDim: ValueCell<Vec2>;
    readonly meansMin: ValueCell<Vec3>;
    readonly meansMax: ValueCell<Vec3>;
    readonly version3: ValueCell<boolean>;
    readonly hasSh: ValueCell<boolean>;
    readonly shCoefficients: ValueCell<number>;
    readonly shCount: ValueCell<number>;
    readonly shDim: ValueCell<Vec2>;
    /** the largest scale of the codebook: the bounds' padding (no scales a splat on the CPU) */
    maxScale: number;
    /**
     * the CPU's centers still being read (`setSogCenters`; the GPU decodes its own): `positionGroup` the groups alone,
     * the bounds `centersBounds`'
     */
    centersPending: boolean;
    /** the bounds of the centers from the means' ranges (`meta.json`) */
    readonly centersBounds: Sphere3D;
  }

  /** The input of `createSog`: SOG's images uploaded to `graphics`, the centers decoded from the means */
  export interface SogInput {
    readonly graphics: GraphicsContext;
    readonly count: number;
    /**
     * the centers (`PositionGroupStride` floats a splat, the group the splat's index); none: read later
     * (`setSogCenters`), the bounds `centersBounds` until then
     */
    readonly positionGroup?: Float32Array;
    /** the bounds of the centers from the means' ranges */
    readonly centersBounds: Sphere3D;
    readonly textures: {
      readonly meansLow: Texture;
      readonly meansHigh: Texture;
      readonly scales: Texture;
      readonly quats: Texture;
      readonly sh0: Texture;
      readonly shLabels?: Texture;
      readonly shCentroids?: Texture;
    };
    readonly meansMin: Vec3;
    readonly meansMax: Vec3;
    /** 256 entries of 4 floats, see `Sog.codebook` */
    readonly codebook: Float32Array;
    readonly version: 2 | 3;
    /** the palette: the degree drawn (at most its bands), its coefficients an entry and its entries */
    readonly sh?: { readonly degree: number; readonly coefficients: number; readonly count: number };
    readonly maxScale: number;
  }

  /**
   * From center xyz and group per splat (`PositionGroupStride` floats, see `interleavePositionGroup`) and their shapes
   * (`SplatStride` int16, see `setSplat`), kept as the geometry's arrays (the `packed` encoding)
   */
  export function create(
    positionGroup: Float32Array,
    splatData: Int16Array,
    splatCount: number,
    splats?: GaussianSplats,
    sh: SplatSh = NoSplatSh,
  ): GaussianSplats {
    const s = splats ?? createObject();
    setArrays(s, positionGroup, splatData, splatCount, sh);
    releaseSog(s.sog);
    ValueCell.updateIfChanged(s.encoding, 'packed');
    ValueCell.update(s.gpuPositionGroup, s.positionGroup.ref.value);
    return s;
  }

  /**
   * From SOG's images uploaded to `input.graphics` (the `sog` encoding): the textures owned by the geometry from then
   * on (destroyed when replaced, `dispose`), the centers on the CPU for the sorts and the bounds
   */
  export function createSog(input: SogInput, splats?: GaussianSplats): GaussianSplats {
    const s = splats ?? createObject();
    const { sog } = s;
    const sh = input.sh && input.sh.degree > 0 ? input.sh : undefined;
    const t = input.textures;
    const hasSh = !!sh && !!t.shLabels && !!t.shCentroids;
    Sphere3D.copy(sog.centersBounds, input.centersBounds);
    const positionGroup = input.positionGroup ?? createGroupsOnly(input.count);
    // the packed layers empty
    setArrays(s, positionGroup, new Int16Array(0), input.count, NoSplatSh, 0);
    // the degree drawn from the palette (the packed layer empty)
    ValueCell.updateIfChanged(s.shDegree, hasSh ? Math.min(3, sh!.degree) : 0);
    const keep = new Set<Texture>([t.meansLow, t.meansHigh, t.scales, t.quats, t.sh0]);
    if (t.shLabels) keep.add(t.shLabels);
    if (t.shCentroids) keep.add(t.shCentroids);
    releaseSog(sog, keep);
    sog.graphics = input.graphics;
    sog.centersPending = !input.positionGroup;
    ValueCell.update(sog.meansLow, t.meansLow);
    ValueCell.update(sog.meansHigh, t.meansHigh);
    ValueCell.update(sog.scales, t.scales);
    ValueCell.update(sog.quats, t.quats);
    ValueCell.update(sog.sh0, t.sh0);
    ValueCell.update(sog.shLabels, hasSh ? t.shLabels! : NullTexture);
    ValueCell.update(sog.shCentroids, hasSh ? t.shCentroids! : NullTexture);
    // unused palette images: the geometry's to destroy all the same
    if (!hasSh) {
      t.shLabels?.destroy();
      t.shCentroids?.destroy();
    }
    ValueCell.update(sog.codebook, { array: input.codebook, width: 256, height: 1 });
    // the images' size the textures' own (a bitmap's is 0 once closed)
    ValueCell.update(sog.imageDim, Vec2.create(t.meansLow.getWidth(), t.meansLow.getHeight()));
    ValueCell.update(sog.meansMin, Vec3.clone(input.meansMin));
    ValueCell.update(sog.meansMax, Vec3.clone(input.meansMax));
    ValueCell.updateIfChanged(sog.version3, input.version === 3);
    ValueCell.updateIfChanged(sog.hasSh, hasSh);
    ValueCell.updateIfChanged(sog.shCoefficients, hasSh ? sh!.coefficients : 0);
    ValueCell.updateIfChanged(sog.shCount, hasSh ? sh!.count : 0);
    ValueCell.update(
      sog.shDim,
      hasSh ? Vec2.create(t.shCentroids!.getWidth(), t.shCentroids!.getHeight()) : Vec2.create(1, 1),
    );
    sog.maxScale = input.maxScale;
    ValueCell.updateIfChanged(s.encoding, 'sog');
    ValueCell.update(s.gpuPositionGroup, NoPositionGroup);
    return s;
  }

  /** The centers of the `sog` encoding read later (`SogInput.positionGroup` none), the CPU's for the sorts and bounds */
  export function setSogCenters(splats: GaussianSplats, positionGroup: Float32Array) {
    if (splats.encoding.ref.value !== 'sog' || !splats.sog.centersPending) return;
    splats.sog.centersPending = false;
    setPositionGroup(splats.positionGroup, positionGroup, splats.splatCount);
  }

  /** The groups of `count` splats (each its index) without centers */
  function createGroupsOnly(count: number) {
    const array = new Float32Array(count * PositionGroupStride);
    for (let i = 0; i < count; ++i) array[i * PositionGroupStride + 3] = int32ToFloatBits(i);
    return array;
  }

  /** Releases what the geometry owns on the GPU (SOG's textures) */
  export function dispose(splats: GaussianSplats) {
    releaseSog(splats.sog);
  }

  /** The largest scale of the splats (the bounds' padding) */
  export function getMaxScale(splats: GaussianSplats) {
    return splats.encoding.ref.value === 'sog'
      ? splats.sog.maxScale
      : getMaxSplatScale(splats.splatData.ref.value.array, splats.splatCount);
  }

  /** a texture of none: a SOG image of a geometry not of that encoding (never bound: `usedIf`) */
  const NullTexture = createNullTexture();
  /** the GPU's centers of the `sog` encoding: one element (a data array is never empty) */
  const NoPositionGroup = createDataArray(1, PositionGroupStride, Float32Array);

  /** SOG's textures destroyed (those in `keep` apart: the new ones) and unset */
  function releaseSog(sog: Sog, keep?: Set<Texture>) {
    for (const cell of [sog.meansLow, sog.meansHigh, sog.scales, sog.quats, sog.sh0, sog.shLabels, sog.shCentroids]) {
      const texture = cell.ref.value;
      if (texture === NullTexture) continue;
      if (!keep?.has(texture)) texture.destroy();
      ValueCell.update(cell, NullTexture);
    }
    sog.graphics = undefined;
    sog.centersPending = false;
  }

  function createShDataArray(sh: SplatSh, count: number, array?: Uint8Array) {
    const n = ShBasisCounts[sh.degree];
    if (sh.data.length < count * n * 4)
      throw new Error(`view-dependent color of ${sh.data.length} bytes for ${count} splats`);
    // at least one element: an empty layer
    return n === 0 ? createDataArray(1, 4, Uint8Array, array) : createDataArray(count * n, 4, Uint8Array, sh.data);
  }

  export function createEmpty(splats?: GaussianSplats): GaussianSplats {
    const pg = splats ? splats.positionGroup.ref.value.array : new Float32Array(0);
    const sd = splats ? splats.splatData.ref.value.array : new Int16Array(0);
    return create(pg, sd, 0, splats);
  }

  function createSogCells(): Sog {
    const texture = () => ValueCell.create(NullTexture);
    return {
      graphics: undefined,
      meansLow: texture(),
      meansHigh: texture(),
      scales: texture(),
      quats: texture(),
      sh0: texture(),
      shLabels: texture(),
      shCentroids: texture(),
      codebook: ValueCell.create({ array: new Float32Array(256 * 4), width: 256, height: 1 }),
      imageDim: ValueCell.create(Vec2.create(1, 1)),
      meansMin: ValueCell.create(Vec3()),
      meansMax: ValueCell.create(Vec3()),
      version3: ValueCell.create(false),
      hasSh: ValueCell.create(false),
      shCoefficients: ValueCell.create(0),
      shCount: ValueCell.create(0),
      shDim: ValueCell.create(Vec2.create(1, 1)),
      maxScale: 0,
      centersPending: false,
      centersBounds: Sphere3D(),
    };
  }

  /** An empty geometry of the `packed` encoding, its arrays set by `create` or `createSog` */
  function createObject(): GaussianSplats {
    const data = createCenterGroupData(new Float32Array(0), 0, () => splats.splatCount);
    const boundingSphere = Sphere3D();
    let currentHash = -1;
    const hashCode = () =>
      hashFnv32a([
        splats.splatCount,
        splats.positionGroup.ref.version,
        splats.splatData.ref.version,
        splats.sog.codebook.ref.version,
        splats.sog.centersPending ? 1 : 0,
      ]);

    const splats: GaussianSplats = {
      kind: 'gaussian-splats' as const,
      splatCount: 0,
      positionGroup: data.positionGroup,
      splatData: ValueCell.create(createSplatDataArray(new Int16Array(0), 0)),
      sh: NoSplatSh,
      shData: ValueCell.create(createShDataArray(NoSplatSh, 0)),
      shDegree: ValueCell.create<number>(0),
      shScale: ValueCell.create(NoSplatSh.scale),
      encoding: ValueCell.create<Encoding>('packed'),
      sog: createSogCells(),
      gpuPositionGroup: ValueCell.create(data.positionGroup.ref.value),
      get boundingSphere() {
        const newHash = hashCode();
        if (newHash !== currentHash) {
          // SOG's centers still being read: their bounds from the means' ranges
          const pending = splats.encoding.ref.value === 'sog' && splats.sog.centersPending;
          const centers = pending ? splats.sog.centersBounds : data.getBoundingSphere();
          Sphere3D.expand(boundingSphere, centers, BoundsSigmas * getMaxScale(splats));
          currentHash = newHash;
        }
        return boundingSphere;
      },
      get groupMapping() {
        return data.getGroupMapping();
      },
      setBoundingSphere(sphere: Sphere3D) {
        Sphere3D.copy(boundingSphere, sphere);
        currentHash = hashCode();
      },
      hasBoundingSphere() {
        return currentHash === hashCode();
      },
    };
    return splats;
  }

  function setArrays(
    splats: GaussianSplats,
    positionGroup: Float32Array,
    splatData: Int16Array,
    splatCount: number,
    sh: SplatSh,
    /** the splats in `splatData` and `sh` (none: the `sog` encoding) */
    packedCount = splatCount,
  ) {
    splats.splatCount = splatCount;
    setPositionGroup(splats.positionGroup, positionGroup, splatCount);
    ValueCell.update(splats.splatData, createSplatDataArray(splatData, packedCount));
    splats.sh = sh;
    ValueCell.update(splats.shData, createShDataArray(sh, packedCount, splats.shData.ref.value.array));
    ValueCell.updateIfChanged(splats.shDegree, sh.degree);
    ValueCell.updateIfChanged(splats.shScale, sh.scale);
  }

  function createSplatDataArray(array: Int16Array, count: number) {
    if (array.length < count * SplatStride) {
      throw new Error(`splat data of ${array.length} int16 for ${count} splats (${SplatStride} each)`);
    }
    return createDataArray(count * SplatTexels, 4, Int16Array, array);
  }

  /** Arrays for `count` splats: center xyz and group, and their shapes (`create`'s inputs) */
  export function createArrays(count: number) {
    return {
      positionGroup: new Float32Array(count * PositionGroupStride),
      splatData: new Int16Array(count * SplatStride),
    };
  }

  export const Params = {
    ...BaseGeometry.Params,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.01 }, { description: 'Scales the size of every splat.' }),
    kernelSize: PD.Numeric(
      3,
      { min: 1, max: 4, step: 0.1 },
      {
        description:
          'The extent of a splat in standard deviations (at 3 its edge is 1.1% of its center, at 2 13.5%: visibly cut).',
      },
    ),
    alphaCutoff: PD.Numeric(
      1 / 255,
      { min: 0, max: 0.5, step: 0.001 },
      { description: 'Fragments of a lower coverage (the falloff times the opacity) are not drawn.' },
    ),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<GaussianSplats, Params> = {
    Params,
    createEmpty,
    createValues,
    createValuesSimple,
    updateValues,
    updateBoundingSphere,
    createRenderableState,
    updateRenderableState,
    createPositionIterator,
    // SOG's textures are of the context they were uploaded to: the shape's getter uploads them again to another
    isCompatible: (splats, graphics) => splats.encoding.ref.value !== 'sog' || splats.sog.graphics === graphics,
  };

  function createPositionIterator(splats: GaussianSplats, transform: TransformData): LocationIterator {
    return createPositionGroupIterator(splats.positionGroup, splats.splatCount, transform);
  }

  /** The bounds at `sizeFactor` and `kernelSize`: the geometry's are at 1 and `BoundsSigmas` */
  function getInvariantBoundingSphere(splats: GaussianSplats, sizeFactor: number, kernelSize: number) {
    const sphere = Sphere3D.clone(splats.boundingSphere);
    const scale = (sizeFactor * kernelSize) / BoundsSigmas;
    if (scale > 1) {
      const padding = BoundsSigmas * getMaxScale(splats);
      Sphere3D.expand(sphere, sphere, padding * (scale - 1));
    }
    return sphere;
  }

  function createValues(
    splats: GaussianSplats,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): GaussianSplatsValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(splats, transform);

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

    // each splat is a quad of two triangles, the vertex count is the number of splats (data per splat)
    const counts = { drawCount: splats.splatCount * 2 * 3, vertexCount: splats.splatCount, groupCount, instanceCount };

    const invariantBoundingSphere = getInvariantBoundingSphere(splats, props.sizeFactor, props.kernelSize);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('gaussianSplats'),

      dSplatEncoding: splats.encoding,
      splatCenters: splats.positionGroup,
      tPositionGroup: splats.gpuPositionGroup,
      tSplat: splats.splatData,
      tSplatSh: splats.shData,
      uSplatShDegree: splats.shDegree,
      uSplatShScale: splats.shScale,
      tSogMeansLow: splats.sog.meansLow,
      tSogMeansHigh: splats.sog.meansHigh,
      tSogScales: splats.sog.scales,
      tSogQuats: splats.sog.quats,
      tSogSh0: splats.sog.sh0,
      uSogImageDim: splats.sog.imageDim,
      uSogMeansMin: splats.sog.meansMin,
      uSogMeansMax: splats.sog.meansMax,
      tSogCodebook: splats.sog.codebook,
      dSogVersion3: splats.sog.version3,
      dSogSh: splats.sog.hasSh,
      tSogShLabels: splats.sog.shLabels,
      tSogShCentroids: splats.sog.shCentroids,
      uSogShCoefficients: splats.sog.shCoefficients,
      uSogShCount: splats.sog.shCount,
      uSogShDim: splats.sog.shDim,
      // the splats' own order until sorted for a view (by the renderer)
      tSplatOrder: ValueCell.create(createSplatOrder(splats.splatCount, instanceCount)),
      // sorted on the GPU (WebGPU) once the renderer has, the order its buffer
      dSplatGpuOrder: ValueCell.create(false),
      sSplatOrder: ValueCell.create(undefined as unknown as StorageBuffer),
      boundingSphere: ValueCell.create(boundingSphere),
      invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
      uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
      ...color,
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
      uKernelSize: ValueCell.create(props.kernelSize),
      uSplatAlphaCutoff: ValueCell.create(props.alphaCutoff),
    };
  }

  function createValuesSimple(
    splats: GaussianSplats,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(splats, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: GaussianSplatsValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
    ValueCell.updateIfChanged(values.uKernelSize, props.kernelSize);
    ValueCell.updateIfChanged(values.uSplatAlphaCutoff, props.alphaCutoff);
  }

  function updateBoundingSphere(values: GaussianSplatsValues, splats: GaussianSplats) {
    const invariantBoundingSphere = getInvariantBoundingSphere(
      splats,
      values.uSizeFactor.ref.value,
      values.uKernelSize.ref.value,
    );
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
    // blended back to front, never writing depth
    state.opaque = false;
    state.writeDepth = false;
  }
}
