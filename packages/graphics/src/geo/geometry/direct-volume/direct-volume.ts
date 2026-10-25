/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { hashFnv32a } from '@molstar/core/data/util';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import type { DirectVolumeValues } from '@molstar/graphics/gl/renderable/direct-volume';
import {
  calculateTransformBoundingSphere,
  createTextureImage,
  type TextureImage,
} from '@molstar/graphics/gl/renderable/util';
import { clamp } from '@molstar/core/math/interpolate';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import type { Theme } from '@molstar/graphics/theme/theme';
import { ValueCell } from '@molstar/core/util';
import type { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Box } from '../../primitive/box.js';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createColors } from '../color-data.js';
import type { GeometryUtils } from '../geometry.js';
import { createMarkers } from '../marker-data.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import type { TransformData } from '../transform-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import { createTransferFunctionTexture, getControlPointsFromVec2Array } from './transfer-function.js';
import { createEmptyClipping } from '../clipping-data.js';
import { Grid } from '@molstar/model/model/volume';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';

const VolumeBox = Box();

/** the four diagonals of a unit cube, by their signs */
const BoxDiagonals = [
  [1, 1, 1],
  [-1, 1, 1],
  [1, -1, 1],
  [1, 1, -1],
] as const;

/**
 * The value range of the grid's bricks of `BrickSize`³ cells, for skipping bricks the transfer function makes
 * transparent: per brick the lowest and highest value (normalized, times 255, rounded outwards), of its cells and the
 * next cell along each axis (trilinear samples in a brick read those); `dim` the bricks per axis, 0 when not known (no
 * skipping)
 */
export interface BrickRange {
  readonly array: Uint8Array;
  readonly dim: Vec3;
}

/** The edge of a brick in cells, see `BrickRange` */
export const BrickSize = 8;

const EmptyBrickRange: BrickRange = { array: new Uint8Array(0), dim: Vec3.create(0, 0, 0) };

export interface DirectVolume {
  readonly kind: 'direct-volume';

  readonly gridTexture: ValueCell<Texture>;
  readonly gridTextureDim: ValueCell<Vec3>;
  readonly gridDimension: ValueCell<Vec3>;
  readonly gridStats: ValueCell<Vec4>; // [min, max, mean, sigma]
  readonly bboxSize: ValueCell<Vec3>;
  readonly bboxMin: ValueCell<Vec3>;
  readonly bboxMax: ValueCell<Vec3>;
  readonly transform: ValueCell<Mat4>;

  readonly cellDim: ValueCell<Vec3>;
  readonly unitToCartn: ValueCell<Mat4>;
  readonly cartnToUnit: ValueCell<Mat4>;
  readonly packedGroup: ValueCell<boolean>;
  readonly axisOrder: ValueCell<Vec3>;
  readonly dataType: ValueCell<'byte' | 'float' | 'halfFloat'>;
  readonly brickRange: ValueCell<BrickRange>;

  /** Bounding sphere of the volume */
  readonly boundingSphere: Sphere3D;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;

  readonly meta: { [k: string]: unknown };
}

export namespace DirectVolume {
  export function create(
    bbox: Box3D,
    gridDimension: Vec3,
    transform: Mat4,
    unitToCartn: Mat4,
    cellDim: Vec3,
    texture: Texture,
    stats: Grid['stats'],
    packedGroup: boolean,
    axisOrder: Vec3,
    dataType: 'byte' | 'float' | 'halfFloat',
    directVolume?: DirectVolume,
    brickRange: BrickRange = EmptyBrickRange,
  ): DirectVolume {
    return directVolume
      ? update(
          bbox,
          gridDimension,
          transform,
          unitToCartn,
          cellDim,
          texture,
          stats,
          packedGroup,
          axisOrder,
          dataType,
          brickRange,
          directVolume,
        )
      : fromData(
          bbox,
          gridDimension,
          transform,
          unitToCartn,
          cellDim,
          texture,
          stats,
          packedGroup,
          axisOrder,
          dataType,
          brickRange,
        );
  }

  function hashCode(directVolume: DirectVolume) {
    return hashFnv32a([
      directVolume.bboxSize.ref.version,
      directVolume.gridDimension.ref.version,
      directVolume.gridTexture.ref.version,
      directVolume.transform.ref.version,
      directVolume.gridStats.ref.version,
    ]);
  }

  function fromData(
    bbox: Box3D,
    gridDimension: Vec3,
    transform: Mat4,
    unitToCartn: Mat4,
    cellDim: Vec3,
    texture: Texture,
    stats: Grid['stats'],
    packedGroup: boolean,
    axisOrder: Vec3,
    dataType: 'byte' | 'float' | 'halfFloat',
    brickRange: BrickRange,
  ): DirectVolume {
    const boundingSphere = Sphere3D();
    let currentHash = -1;

    const width = texture.getWidth();
    const height = texture.getHeight();
    const depth = texture.getDepth();

    const directVolume = {
      kind: 'direct-volume' as const,
      gridDimension: ValueCell.create(gridDimension),
      gridTexture: ValueCell.create(texture),
      gridTextureDim: ValueCell.create(Vec3.create(width, height, depth)),
      gridStats: ValueCell.create(Vec4.create(stats.min, stats.max, stats.mean, stats.sigma)),
      bboxMin: ValueCell.create(bbox.min),
      bboxMax: ValueCell.create(bbox.max),
      bboxSize: ValueCell.create(Vec3.sub(Vec3(), bbox.max, bbox.min)),
      transform: ValueCell.create(transform),
      cellDim: ValueCell.create(cellDim),
      unitToCartn: ValueCell.create(unitToCartn),
      cartnToUnit: ValueCell.create(Mat4.invert(Mat4(), unitToCartn)),
      get boundingSphere() {
        const newHash = hashCode(directVolume);
        if (newHash !== currentHash) {
          const b = getBoundingSphere(directVolume.gridDimension.ref.value, directVolume.transform.ref.value);
          Sphere3D.copy(boundingSphere, b);
          currentHash = newHash;
        }
        return boundingSphere;
      },
      packedGroup: ValueCell.create(packedGroup),
      axisOrder: ValueCell.create(axisOrder),
      dataType: ValueCell.create(dataType),
      brickRange: ValueCell.create(brickRange),
      setBoundingSphere(sphere: Sphere3D) {
        Sphere3D.copy(boundingSphere, sphere);
        currentHash = hashCode(directVolume);
      },
      hasBoundingSphere() {
        return currentHash === hashCode(directVolume);
      },
      meta: {},
    };
    return directVolume;
  }

  function update(
    bbox: Box3D,
    gridDimension: Vec3,
    transform: Mat4,
    unitToCartn: Mat4,
    cellDim: Vec3,
    texture: Texture,
    stats: Grid['stats'],
    packedGroup: boolean,
    axisOrder: Vec3,
    dataType: 'byte' | 'float' | 'halfFloat',
    brickRange: BrickRange,
    directVolume: DirectVolume,
  ): DirectVolume {
    const width = texture.getWidth();
    const height = texture.getHeight();
    const depth = texture.getDepth();

    ValueCell.update(directVolume.gridDimension, gridDimension);
    ValueCell.update(directVolume.gridTexture, texture);
    ValueCell.update(
      directVolume.gridTextureDim,
      Vec3.set(directVolume.gridTextureDim.ref.value, width, height, depth),
    );
    ValueCell.update(
      directVolume.gridStats,
      Vec4.set(directVolume.gridStats.ref.value, stats.min, stats.max, stats.mean, stats.sigma),
    );
    ValueCell.update(directVolume.bboxMin, bbox.min);
    ValueCell.update(directVolume.bboxMax, bbox.max);
    ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min));
    ValueCell.update(directVolume.transform, transform);
    ValueCell.update(directVolume.cellDim, cellDim);
    ValueCell.update(directVolume.unitToCartn, unitToCartn);
    ValueCell.update(directVolume.cartnToUnit, Mat4.invert(Mat4(), unitToCartn));
    ValueCell.updateIfChanged(directVolume.packedGroup, packedGroup);
    ValueCell.updateIfChanged(directVolume.axisOrder, Vec3.fromArray(directVolume.axisOrder.ref.value, axisOrder, 0));
    ValueCell.updateIfChanged(directVolume.dataType, dataType);
    ValueCell.update(directVolume.brickRange, brickRange);
    return directVolume;
  }

  export function createEmpty(directVolume?: DirectVolume): DirectVolume {
    const bbox = Box3D();
    const gridDimension = Vec3();
    const transform = Mat4.identity();
    const unitToCartn = Mat4.identity();
    const cellDim = Vec3();
    const texture = createNullTexture();
    const stats = Grid.One.stats;
    const packedGroup = false;
    const axisOrder = Vec3.create(0, 1, 2);
    const dataType = 'byte';
    return create(
      bbox,
      gridDimension,
      transform,
      unitToCartn,
      cellDim,
      texture,
      stats,
      packedGroup,
      axisOrder,
      dataType,
      directVolume,
    );
  }

  export const Params = {
    ...BaseGeometry.Params,
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
    controlPoints: PD.LineGraph(
      [
        Vec2.create(0.19, 0.0),
        Vec2.create(0.2, 0.05),
        Vec2.create(0.25, 0.05),
        Vec2.create(0.26, 0.0),
        Vec2.create(0.79, 0.0),
        Vec2.create(0.8, 0.05),
        Vec2.create(0.85, 0.05),
        Vec2.create(0.86, 0.0),
      ],
      { isEssential: true, yAxis: { scale: 'log' } },
    ),
    stepsPerCell: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
    surfaceOpacity: PD.Numeric(
      0.3,
      { min: 0.05, max: 0.95, step: 0.05 },
      {
        description:
          'Opacity along a ray at which the volume counts as solid: where highlights and selections draw their edges ' +
          'and picking finds it. Lower for faint transfer functions.',
      },
    ),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<DirectVolume, Params> = {
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

  function createPositionIterator(directVolume: DirectVolume, transform: TransformData): LocationIterator {
    const t = directVolume.transform.ref.value;
    const [x, y, z] = directVolume.gridDimension.ref.value;
    const groupCount = x * y * z;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      const k = Math.floor(groupIndex / z);
      p[0] = Math.floor(k / y);
      p[1] = k % y;
      p[2] = groupIndex % z;
      Vec3.transformMat4(p, p, t);
      if (instanceIndex >= 0) {
        Vec3.transformMat4Offset(p, p, m, 0, 0, instanceIndex * 16);
      }
      return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
  }

  /** the steps along the longest diagonal of the grid's box (`unitToCartn`'s unit cube), the ray marcher's loop bound */
  function getMaxSteps(unitToCartn: Mat4, cellDim: Vec3, stepsPerCell: number) {
    const stepScale = getStepScale(cellDim, stepsPerCell);
    if (!(stepScale > 0)) return 0;
    const m = unitToCartn;
    let diagonal = 0;
    for (const [sx, sy, sz] of BoxDiagonals) {
      const x = sx * m[0] + sy * m[4] + sz * m[8];
      const y = sx * m[1] + sy * m[5] + sz * m[9];
      const z = sx * m[2] + sy * m[6] + sz * m[10];
      diagonal = Math.max(diagonal, Math.sqrt(x * x + y * y + z * z));
    }
    return Math.ceil(diagonal / stepScale) + 1;
  }

  function getStepScale(cellDim: Vec3, stepsPerCell: number) {
    return Math.min(...cellDim) * (1 / stepsPerCell);
  }

  function getTransferScale(stepsPerCell: number) {
    return 1 / stepsPerCell;
  }

  /**
   * The `jumpLength` of the Gaussian volume's params (no CPU grid for bricks: it probes a step that far ahead in empty
   * space instead), 0 for volumes (they skip transparent bricks)
   */
  function getJumpLength(props: PD.Values<Params>) {
    return (props as { jumpLength?: number }).jumpLength ?? 0;
  }

  /**
   * Per brick 255 where the transfer function is above 0 somewhere in the brick's value range (`BrickRange`), else 0:
   * the ray marcher skips the bricks of 0
   */
  function getBrickVisibility(brickRange: BrickRange, transfer: Uint8Array, image?: TextureImage<Uint8Array>) {
    const [bx, by, bz] = brickRange.dim;
    const count = bx * by * bz;
    const out = createTextureImage(count, 1, Uint8Array, image?.array);
    const n = transfer.length;
    // the transfer entries above 0 before each entry
    const above = new Uint32Array(n + 1);
    for (let i = 0; i < n; ++i) above[i + 1] = above[i] + (transfer[i] > 0 ? 1 : 0);
    const { array } = brickRange;
    for (let i = 0; i < count; ++i) {
      // the transfer texels that values in [lo, hi] / 255 filter from (linear filtering of `n` texels, clamped)
      const t0 = clamp(Math.floor((array[2 * i] / 255) * n - 0.5), 0, n - 1);
      const t1 = clamp(Math.ceil((array[2 * i + 1] / 255) * n - 0.5), 0, n - 1);
      out.array[i] = above[t1 + 1] - above[t0] > 0 ? 255 : 0;
    }
    return out;
  }

  function createValues(
    directVolume: DirectVolume,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): DirectVolumeValues {
    const { gridTexture, gridTextureDim, gridStats } = directVolume;
    const { bboxSize, bboxMin, bboxMax, gridDimension, transform: gridTransform } = directVolume;

    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(directVolume, transform);

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

    const [x, y, z] = gridDimension.ref.value;
    const counts = { drawCount: VolumeBox.indices.length, vertexCount: x * y * z, groupCount, instanceCount };

    const invariantBoundingSphere = Sphere3D.clone(directVolume.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    const controlPoints = getControlPointsFromVec2Array(props.controlPoints);
    const transferTex = createTransferFunctionTexture(controlPoints);
    const brickRange = directVolume.brickRange.ref.value;
    const brickImage = getBrickVisibility(brickRange, transferTex.ref.value.array);

    return {
      dGeometryType: ValueCell.create('directVolume'),

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

      aPosition: ValueCell.create(VolumeBox.vertices as Float32Array),
      elements: ValueCell.create(VolumeBox.indices as Uint32Array),
      boundingSphere: ValueCell.create(boundingSphere),
      invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
      uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),

      uBboxMin: bboxMin,
      uBboxMax: bboxMax,
      uBboxSize: bboxSize,
      uMaxSteps: ValueCell.create(
        getMaxSteps(directVolume.unitToCartn.ref.value, directVolume.cellDim.ref.value, props.stepsPerCell),
      ),
      uStepScale: ValueCell.create(getStepScale(directVolume.cellDim.ref.value, props.stepsPerCell)),
      uJumpLength: ValueCell.create(getJumpLength(props)),
      uTransform: gridTransform,
      uGridDim: gridDimension,
      tTransferTex: transferTex,
      uTransferScale: ValueCell.create(getTransferScale(props.stepsPerCell)),
      uSurfaceOpacity: ValueCell.create(props.surfaceOpacity),
      brickRange: directVolume.brickRange,
      tBrickTex: ValueCell.create(brickImage),
      uBrickTexDim: ValueCell.create(Vec2.create(brickImage.width, brickImage.height)),
      uBrickDim: ValueCell.create(Vec3.clone(brickRange.dim)),

      dGridTexType: ValueCell.create(gridTexture.ref.value.getDepth() > 0 ? '3d' : '2d'),
      uGridTexDim: gridTextureDim,
      tGridTex: gridTexture,
      uGridStats: gridStats,

      uCellDim: directVolume.cellDim,
      uCartnToUnit: directVolume.cartnToUnit,
      uUnitToCartn: directVolume.unitToCartn,
      dPackedGroup: directVolume.packedGroup,
      dAxisOrder: ValueCell.create(directVolume.axisOrder.ref.value.join('')),

      dIgnoreLight: ValueCell.create(props.ignoreLight),
      dCelShaded: ValueCell.create(props.celShaded),
      dXrayShaded: ValueCell.create(
        props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
      ),

      meta: ValueCell.create(directVolume.meta),
    };
  }

  function createValuesSimple(
    directVolume: DirectVolume,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(directVolume, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: DirectVolumeValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
    ValueCell.updateIfChanged(
      values.dXrayShaded,
      props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
    );

    const controlPoints = getControlPointsFromVec2Array(props.controlPoints);
    createTransferFunctionTexture(controlPoints, values.tTransferTex);
    // the transfer function or the bricks (a new geometry) may have changed
    const brickRange = values.brickRange.ref.value as BrickRange;
    const brickImage = getBrickVisibility(brickRange, values.tTransferTex.ref.value.array, values.tBrickTex.ref.value);
    ValueCell.update(values.tBrickTex, brickImage);
    ValueCell.update(values.uBrickTexDim, Vec2.set(values.uBrickTexDim.ref.value, brickImage.width, brickImage.height));
    if (!Vec3.exactEquals(values.uBrickDim.ref.value, brickRange.dim)) {
      ValueCell.update(values.uBrickDim, Vec3.copy(values.uBrickDim.ref.value, brickRange.dim));
    }

    ValueCell.updateIfChanged(
      values.uMaxSteps,
      getMaxSteps(values.uUnitToCartn.ref.value, values.uCellDim.ref.value, props.stepsPerCell),
    );
    ValueCell.updateIfChanged(values.uStepScale, getStepScale(values.uCellDim.ref.value, props.stepsPerCell));
    ValueCell.updateIfChanged(values.uTransferScale, getTransferScale(props.stepsPerCell));
    ValueCell.updateIfChanged(values.uSurfaceOpacity, props.surfaceOpacity);
    ValueCell.updateIfChanged(values.uJumpLength, getJumpLength(props));
  }

  function updateBoundingSphere(values: DirectVolumeValues, directVolume: DirectVolume) {
    const invariantBoundingSphere = Sphere3D.clone(directVolume.boundingSphere);
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
    state.opaque = false;
    state.writeDepth = false;
    return state;
  }

  function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
    BaseGeometry.updateRenderableState(state, props);
    state.opaque = false;
    state.writeDepth = false;
  }
}

//

function getBoundingSphere(gridDimension: Vec3, gridTransform: Mat4) {
  return Sphere3D.fromDimensionsAndTransform(Sphere3D(), gridDimension, gridTransform);
}
