/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { type Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { transformPositionArray, type GroupMapping, createGroupMapping } from '../../util.js';
import type { GeometryUtils } from '../geometry.js';
import { createColors } from '../color-data.js';
import { createMarkers } from '../marker-data.js';
import { createSizes } from '../size-data.js';
import type { TransformData } from '../transform-data.js';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import type { LinesValues } from '@molstar/graphics/gl/renderable/lines';
import type { Mesh } from '../mesh/mesh.js';
import { LinesBuilder } from './lines-builder.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  calculateInvariantBoundingSphere,
  calculateTransformBoundingSphere,
  createDataArray,
  type DataArray,
} from '@molstar/graphics/gl/renderable/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { Theme } from '@molstar/graphics/theme/theme';
import type { Color } from '@molstar/core/util/color';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import { hashFnv32a } from '@molstar/core/data/util';
import { createEmptyClipping } from '../clipping-data.js';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';

/**
 * Floats per line of `Lines.lineData` (2 texels of `tLine`): start xyz and group, end xyz and whether the line continues
 * the one before (1: a strip's line after its first, starting where the line before ended)
 */
export const LineStride = 8;

/** Wide line */
export interface Lines {
  readonly kind: 'lines';

  /** Number of lines */
  lineCount: number;

  /**
   * Per line (`LineStride` floats, the only copy; the shaders read it as `tLine`): start xyz and group, end xyz and
   * whether it continues the line before, as the lines of a strip after its first: its start corners then are the end
   * corners of the line before (a connected strip, drawn by `StripLinesBuilder`)
   */
  readonly lineData: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the lines */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to line indices times `LineStride` (offsets into `lineData`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
}

export namespace Lines {
  /** From the data of `lineCount` lines (`LineStride` floats each), kept as the geometry's array */
  export function create(lineData: Float32Array, lineCount: number, lines?: Lines): Lines {
    const data = createDataArray(lineCount * 2, 4, Float32Array, lineData);
    if (lines) {
      lines.lineCount = lineCount;
      ValueCell.update(lines.lineData, data);
      return lines;
    }
    return fromData(data, lineCount);
  }

  export function createEmpty(lines?: Lines): Lines {
    return create(lines ? lines.lineData.ref.value.array : new Float32Array(0), 0, lines);
  }

  export function fromMesh(mesh: Mesh, lines?: Lines) {
    const vb = mesh.vertexBuffer.ref.value;
    const ib = mesh.indexBuffer.ref.value;
    const gb = mesh.groupBuffer.ref.value;

    const builder = LinesBuilder.create(mesh.triangleCount * 3, mesh.triangleCount / 10, lines);

    // TODO avoid duplicate lines
    for (let i = 0, il = mesh.triangleCount * 3; i < il; i += 3) {
      const i0 = ib[i],
        i1 = ib[i + 1],
        i2 = ib[i + 2];
      const x0 = vb[i0 * 3],
        y0 = vb[i0 * 3 + 1],
        z0 = vb[i0 * 3 + 2];
      const x1 = vb[i1 * 3],
        y1 = vb[i1 * 3 + 1],
        z1 = vb[i1 * 3 + 2];
      const x2 = vb[i2 * 3],
        y2 = vb[i2 * 3 + 1],
        z2 = vb[i2 * 3 + 2];
      builder.add(x0, y0, z0, x1, y1, z1, gb[i0]);
      builder.add(x0, y0, z0, x2, y2, z2, gb[i0]);
      builder.add(x1, y1, z1, x2, y2, z2, gb[i1]);
    }

    return builder.getLines();
  }

  function hashCode(lines: Lines) {
    return hashFnv32a([lines.lineCount, lines.lineData.ref.version]);
  }

  function fromData(data: DataArray<Float32Array>, lineCount: number): Lines {
    const boundingSphere = Sphere3D();
    let groupMapping: GroupMapping;

    let currentHash = -1;
    let currentGroup = -1;

    const lines = {
      kind: 'lines' as const,
      lineCount,
      lineData: ValueCell.create(data),
      get boundingSphere() {
        const newHash = hashCode(lines);
        if (newHash !== currentHash) {
          const { array } = lines.lineData.ref.value;
          const count = lines.lineCount;
          const s = calculateInvariantBoundingSphere(array, count, 1, LineStride);
          const e = calculateInvariantBoundingSphere(array.subarray(4), count, 1, LineStride);

          Sphere3D.expandBySphere(boundingSphere, s, e);
          currentHash = newHash;
        }
        return boundingSphere;
      },
      get groupMapping() {
        if (lines.lineData.ref.version !== currentGroup) {
          groupMapping = createGroupMapping(lines.lineData.ref.value.array, lines.lineCount, LineStride, 3);
          currentGroup = lines.lineData.ref.version;
        }
        return groupMapping;
      },
      setBoundingSphere(sphere: Sphere3D) {
        Sphere3D.copy(boundingSphere, sphere);
        currentHash = hashCode(lines);
      },
      hasBoundingSphere() {
        return currentHash === hashCode(lines);
      },
    };
    return lines;
  }

  export function transform(lines: Lines, t: Mat4) {
    const data = lines.lineData.ref.value;
    transformPositionArray(t, data.array, 0, lines.lineCount, LineStride);
    transformPositionArray(t, data.array, 4, lines.lineCount, LineStride);
    ValueCell.update(lines.lineData, data);
  }

  //

  export const Params = {
    ...BaseGeometry.Params,
    sizeFactor: PD.Numeric(2, { min: 0, max: 10, step: 0.1 }),
    lineSizeAttenuation: PD.Boolean(false),
    animation: getAnimationParam(),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<Lines, Params> = {
    Params,
    createEmpty,
    createValues,
    createValuesSimple,
    updateValues,
    updateBoundingSphere,
    createRenderableState: BaseGeometry.createRenderableState,
    updateRenderableState: BaseGeometry.updateRenderableState,
    createPositionIterator,
  };

  /** per data vertex (2 per line, see `Geometry.getVertexCount`): the line's start, then its end */
  function createPositionIterator(lines: Lines, transform: TransformData): LocationIterator {
    const groupCount = lines.lineCount * 2;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    const d = lines.lineData.ref.value.array;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      const o = (groupIndex >> 1) * LineStride + (groupIndex & 1) * 4;
      if (instanceIndex < 0) {
        Vec3.fromArray(p, d, o);
      } else {
        Vec3.transformMat4Offset(p, d, m, 0, o, instanceIndex * 16);
      }
      return location;
    };
    return LocationIterator(groupCount, instanceCount, 1, getLocation);
  }

  function createValues(
    lines: Lines,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): LinesValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(lines, transform);

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

    const counts = { drawCount: lines.lineCount * 2 * 3, vertexCount: lines.lineCount * 2, groupCount, instanceCount };

    const invariantBoundingSphere = Sphere3D.clone(lines.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('lines'),

      tLine: lines.lineData,
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
      dLineSizeAttenuation: ValueCell.create(props.lineSizeAttenuation),
      uDoubleSided: ValueCell.create(true),
      dFlipSided: ValueCell.create(false),
      ...createAnimationValues(props.animation),
    };
  }

  function createValuesSimple(
    lines: Lines,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(lines, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: LinesValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
    ValueCell.updateIfChanged(values.dLineSizeAttenuation, props.lineSizeAttenuation);
    updateAnimationValues(values, props.animation);
  }

  function updateBoundingSphere(values: LinesValues, lines: Lines) {
    const invariantBoundingSphere = Sphere3D.clone(lines.boundingSphere);
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
}
