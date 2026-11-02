/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '@molstar/core/util';
import { type Mat4, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { transformPositionArray, type GroupMapping, createGroupMapping } from '../../util.js';
import type { GeometryUtils } from '../geometry.js';
import { createColors } from '../color-data.js';
import { createMarkers } from '../marker-data.js';
import { createSizes, getMaxSize } from '../size-data.js';
import type { TransformData } from '../transform-data.js';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
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
import type { CylindersValues } from '@molstar/graphics/gl/renderable/cylinders';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';
import { getInteriorParam, updateInteriorValues, createInteriorValues } from '../interior.js';
import { getAnimationParam, createAnimationValues, updateAnimationValues } from '../animation.js';

/**
 * Floats per cylinder of `Cylinders.cylinderData` (3 texels of `tCylinder`): start xyz and group, end xyz and radius
 * scale, cap and color mode (2 unused)
 */
export const CylinderStride = 12;

export interface Cylinders {
  readonly kind: 'cylinders';

  /** Number of cylinders */
  cylinderCount: number;

  /**
   * Per cylinder (`CylinderStride` floats, the only copy; the shaders read it as `tCylinder`): start xyz and group, end
   * xyz and radius scale, cap (1 top, 2 bottom, 3 both) and color mode:
   * - for colorMode between 0 and 1 use colorMode to interpolate
   * - for colorMode == 2 do nothing, i.e., use given theme color
   * - for colorMode == 3 use position on cylinder axis to interpolate
   */
  readonly cylinderData: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the cylinders */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to cylinder indices times `CylinderStride` (offsets into `cylinderData`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
}

export namespace Cylinders {
  /** From the data of `cylinderCount` cylinders (`CylinderStride` floats each), kept as the geometry's array */
  export function create(cylinderData: Float32Array, cylinderCount: number, cylinders?: Cylinders): Cylinders {
    const data = createDataArray(cylinderCount * 3, 4, Float32Array, cylinderData);
    if (cylinders) {
      cylinders.cylinderCount = cylinderCount;
      ValueCell.update(cylinders.cylinderData, data);
      return cylinders;
    }
    return fromData(data, cylinderCount);
  }

  export function createEmpty(cylinders?: Cylinders): Cylinders {
    return create(cylinders ? cylinders.cylinderData.ref.value.array : new Float32Array(0), 0, cylinders);
  }

  function hashCode(cylinders: Cylinders) {
    return hashFnv32a([cylinders.cylinderCount, cylinders.cylinderData.ref.version]);
  }

  function fromData(data: DataArray<Float32Array>, cylinderCount: number): Cylinders {
    const boundingSphere = Sphere3D();
    let groupMapping: GroupMapping;

    let currentHash = -1;
    let currentGroup = -1;

    const cylinders = {
      kind: 'cylinders' as const,
      cylinderCount,
      cylinderData: ValueCell.create(data),
      get boundingSphere() {
        const newHash = hashCode(cylinders);
        if (newHash !== currentHash) {
          const { array } = cylinders.cylinderData.ref.value;
          const count = cylinders.cylinderCount;
          const s = calculateInvariantBoundingSphere(array, count, 1, CylinderStride);
          const e = calculateInvariantBoundingSphere(array.subarray(4), count, 1, CylinderStride);

          Sphere3D.expandBySphere(boundingSphere, s, e);
          currentHash = newHash;
        }
        return boundingSphere;
      },
      get groupMapping() {
        if (cylinders.cylinderData.ref.version !== currentGroup) {
          const { array } = cylinders.cylinderData.ref.value;
          groupMapping = createGroupMapping(array, cylinders.cylinderCount, CylinderStride, 3);
          currentGroup = cylinders.cylinderData.ref.version;
        }
        return groupMapping;
      },
      setBoundingSphere(sphere: Sphere3D) {
        Sphere3D.copy(boundingSphere, sphere);
        currentHash = hashCode(cylinders);
      },
      hasBoundingSphere() {
        return currentHash === hashCode(cylinders);
      },
    };
    return cylinders;
  }

  export function transform(cylinders: Cylinders, t: Mat4) {
    const data = cylinders.cylinderData.ref.value;
    transformPositionArray(t, data.array, 0, cylinders.cylinderCount, CylinderStride);
    transformPositionArray(t, data.array, 4, cylinders.cylinderCount, CylinderStride);
    ValueCell.update(cylinders.cylinderData, data);
  }

  //

  export const Params = {
    ...BaseGeometry.Params,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
    sizeAspectRatio: PD.Numeric(1, { min: 0, max: 3, step: 0.01 }),
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
    bumpFrequency: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    bumpAmplitude: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }, BaseGeometry.ShadingCategory),
    interior: getInteriorParam(),
    animation: getAnimationParam(),
    colorMode: PD.Select(
      'default',
      PD.arrayToOptions(['default', 'interpolate'] as const),
      BaseGeometry.ShadingCategory,
    ),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<Cylinders, Params> = {
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

  /** per data vertex (2 per cylinder, see `Geometry.getVertexCount`): the cylinder's start, then its end */
  function createPositionIterator(cylinders: Cylinders, transform: TransformData): LocationIterator {
    const groupCount = cylinders.cylinderCount * 2;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    const d = cylinders.cylinderData.ref.value.array;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      const o = (groupIndex >> 1) * CylinderStride + (groupIndex & 1) * 4;
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
    cylinders: Cylinders,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): CylindersValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(cylinders, transform);

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
      drawCount: cylinders.cylinderCount * 4 * 3,
      vertexCount: cylinders.cylinderCount * 2,
      groupCount,
      instanceCount,
    };

    const padding = getMaxSize(size) * props.sizeFactor;
    const invariantBoundingSphere = Sphere3D.clone(cylinders.boundingSphere);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );
    return {
      dGeometryType: ValueCell.create('cylinders'),

      tCylinder: cylinders.cylinderData,
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
      uSizeFactor: ValueCell.create(props.sizeFactor * props.sizeAspectRatio),
      uDoubleSided: ValueCell.create(props.doubleSided),
      dIgnoreLight: ValueCell.create(props.ignoreLight),
      dCelShaded: ValueCell.create(props.celShaded),
      dXrayShaded: ValueCell.create(
        props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
      ),
      dTransparentBackfaces: ValueCell.create(props.transparentBackfaces),
      dSolidInterior: ValueCell.create(props.solidInterior),
      dBump: ValueCell.create(BaseGeometry.hasBump(props)),
      uBumpFrequency: ValueCell.create(props.bumpFrequency),
      uBumpAmplitude: ValueCell.create(props.bumpAmplitude),
      dDualColor: ValueCell.create(props.colorMode === 'interpolate'),

      ...createInteriorValues(props.interior),
      ...createAnimationValues(props.animation),
    };
  }

  function createValuesSimple(
    cylinders: Cylinders,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(cylinders, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: CylindersValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor * props.sizeAspectRatio);
    ValueCell.updateIfChanged(values.uDoubleSided, props.doubleSided);
    ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
    ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
    ValueCell.updateIfChanged(
      values.dXrayShaded,
      props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off',
    );
    ValueCell.updateIfChanged(values.dTransparentBackfaces, props.transparentBackfaces);
    ValueCell.updateIfChanged(values.dSolidInterior, props.solidInterior);
    ValueCell.updateIfChanged(values.dBump, BaseGeometry.hasBump(props));
    ValueCell.updateIfChanged(values.uBumpFrequency, props.bumpFrequency);
    ValueCell.updateIfChanged(values.uBumpAmplitude, props.bumpAmplitude);
    ValueCell.updateIfChanged(values.dDualColor, props.colorMode === 'interpolate');
    updateInteriorValues(values, props.interior);
    updateAnimationValues(values, props.animation);
  }

  function updateBoundingSphere(values: CylindersValues, cylinders: Cylinders) {
    const invariantBoundingSphere = Sphere3D.clone(cylinders.boundingSphere);
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
    state.opaque = state.opaque && !props.xrayShaded;
    state.writeDepth = state.opaque;
  }
}
