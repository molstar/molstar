/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util';
import type { GeometryUtils } from '../geometry.js';
import { LocationIterator, PositionLocation } from '../../util/location-iterator.js';
import type { TransformData } from '../transform-data.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import { createColors } from '../color-data.js';
import { createSizes, getMaxSize } from '../size-data.js';
import { createMarkers } from '../marker-data.js';
import { ColorNames } from '@molstar/core/util/color/names';
import { Sphere3D } from '@molstar/core/math/geometry';
import {
  type TextureImage,
  createTextureImage,
  calculateInvariantBoundingSphere,
  calculateTransformBoundingSphere,
  createDataArray,
  type DataArray,
} from '@molstar/graphics/gl/renderable/util';
import type { TextValues } from '@molstar/graphics/gl/renderable/text';
import { Color } from '@molstar/core/util/color';
import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { FontAtlasParams } from './font-atlas.js';
import type { RenderableState } from '@molstar/graphics/gl/renderable';
import { clamp } from '@molstar/core/math/interpolate';
import { createRenderObject as _createRenderObject } from '@molstar/graphics/gl/render-object';
import { BaseGeometry, resolveInstanceGranularity } from '../base.js';
import { createEmptyOverpaint } from '../overpaint-data.js';
import { createEmptyTransparency } from '../transparency-data.js';
import { hashFnv32a } from '@molstar/core/data/util';
import { type GroupMapping, createGroupMapping } from '../../util.js';
import { createEmptyClipping } from '../clipping-data.js';
import { createEmptySubstance } from '../substance-data.js';
import { createEmptyEmissive } from '../emissive-data.js';
import { createEmptyWiggle } from '../wiggle-data.js';

type TextAttachment =
  | 'bottom-left'
  | 'bottom-center'
  | 'bottom-right'
  | 'middle-left'
  | 'middle-center'
  | 'middle-right'
  | 'top-left'
  | 'top-center'
  | 'top-right';

/**
 * Floats per item of `Text.textData` (4 texels of `tText`): the label's position xyz and group; its depth and whether the
 * item is a tether (2 unused); then for a glyph or a background its corners (left, top, right, bottom) and texture
 * coordinates (left, top, right, bottom; 10 for the background), for a tether its four corners (base B, tip, base
 * center, base A: the two triangles of the shared indices)
 */
export const TextStride = 16;

/** Text */
export interface Text {
  readonly kind: 'text';

  /** Number of items of the text: its glyphs, backgrounds and tethers */
  charCount: number;

  /** Font Atlas */
  readonly fontTexture: ValueCell<TextureImage<Uint8Array>>;

  /** Per item (`TextStride` floats, the only copy; the shaders read it as `tText`), see `TextStride` */
  readonly textData: ValueCell<DataArray<Float32Array>>;

  /** Bounding sphere of the text */
  readonly boundingSphere: Sphere3D;
  /** Maps group ids to item indices times `TextStride` (offsets into `textData`'s array) */
  readonly groupMapping: GroupMapping;

  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
}

export namespace Text {
  /** From the data of `charCount` items (`TextStride` floats each), kept as the geometry's array */
  export function create(
    fontTexture: TextureImage<Uint8Array>,
    textData: Float32Array,
    charCount: number,
    text?: Text,
  ): Text {
    const data = createDataArray(charCount * 4, 4, Float32Array, textData);
    if (text) {
      text.charCount = charCount;
      ValueCell.update(text.fontTexture, fontTexture);
      ValueCell.update(text.textData, data);
      return text;
    }
    return fromData(fontTexture, data, charCount);
  }

  export function createEmpty(text?: Text): Text {
    const ft = text ? text.fontTexture.ref.value : createTextureImage(0, 1, Uint8Array);
    return create(ft, text ? text.textData.ref.value.array : new Float32Array(0), 0, text);
  }

  function hashCode(text: Text) {
    return hashFnv32a([text.charCount, text.fontTexture.ref.version, text.textData.ref.version]);
  }

  function fromData(fontTexture: TextureImage<Uint8Array>, data: DataArray<Float32Array>, charCount: number): Text {
    const boundingSphere = Sphere3D();
    let groupMapping: GroupMapping;

    let currentHash = -1;
    let currentGroup = -1;

    const text = {
      kind: 'text' as const,
      charCount,
      fontTexture: ValueCell.create(fontTexture),
      textData: ValueCell.create(data),
      get boundingSphere() {
        const newHash = hashCode(text);
        if (newHash !== currentHash) {
          const b = calculateInvariantBoundingSphere(text.textData.ref.value.array, text.charCount, 1, TextStride);
          Sphere3D.copy(boundingSphere, b);
          currentHash = newHash;
        }
        return boundingSphere;
      },
      get groupMapping() {
        if (text.textData.ref.version !== currentGroup) {
          groupMapping = createGroupMapping(text.textData.ref.value.array, text.charCount, TextStride, 3);
          currentGroup = text.textData.ref.version;
        }
        return groupMapping;
      },
      setBoundingSphere(sphere: Sphere3D) {
        Sphere3D.copy(boundingSphere, sphere);
        currentHash = hashCode(text);
      },
      hasBoundingSphere() {
        return currentHash === hashCode(text);
      },
    };
    return text;
  }

  export const Params = {
    ...BaseGeometry.Params,
    ...FontAtlasParams,
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),

    borderWidth: PD.Numeric(0, { min: 0, max: 0.5, step: 0.01 }),
    borderColor: PD.Color(ColorNames.grey),
    offsetX: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    offsetY: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    offsetZ: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    background: PD.Boolean(false),
    backgroundMargin: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
    backgroundColor: PD.Color(ColorNames.grey),
    backgroundOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    tether: PD.Boolean(false),
    tetherLength: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }),
    tetherBaseWidth: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }),

    attachment: PD.Select('middle-center', [
      ['bottom-left', 'bottom-left'],
      ['bottom-center', 'bottom-center'],
      ['bottom-right', 'bottom-right'],
      ['middle-left', 'middle-left'],
      ['middle-center', 'middle-center'],
      ['middle-right', 'middle-right'],
      ['top-left', 'top-left'],
      ['top-center', 'top-center'],
      ['top-right', 'top-right'],
    ] as [TextAttachment, string][]),
  };
  export type Params = typeof Params;

  export const Utils: GeometryUtils<Text, Params> = {
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

  function createPositionIterator(text: Text, transform: TransformData): LocationIterator {
    const groupCount = text.charCount;
    const instanceCount = transform.instanceCount.ref.value;
    const location = PositionLocation();
    const p = location.position;
    // per data vertex (1 per item, see `Geometry.getVertexCount`): the label's position
    const d = text.textData.ref.value.array;
    const m = transform.aTransform.ref.value;
    const getLocation = (groupIndex: number, instanceIndex: number) => {
      const o = groupIndex * TextStride;
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
    text: Text,
    transform: TransformData,
    locationIt: LocationIterator,
    theme: Theme,
    props: PD.Values<Params>,
  ): TextValues {
    const { instanceCount, groupCount } = locationIt;
    const positionIt = createPositionIterator(text, transform);

    const color = createColors(locationIt, positionIt, theme.color);
    const size = createSizes(locationIt, positionIt, theme.size);
    const marker = resolveInstanceGranularity(props.instanceGranularity, groupCount, instanceCount)
      ? createMarkers(instanceCount, 'instance')
      : createMarkers(instanceCount * groupCount, 'groupInstance');
    const overpaint = createEmptyOverpaint();
    const transparency = createEmptyTransparency();
    const emissive = createEmptyEmissive();
    const substance = createEmptySubstance();
    const clipping = createEmptyClipping();
    const wiggle = createEmptyWiggle();

    const counts = { drawCount: text.charCount * 2 * 3, vertexCount: text.charCount, groupCount, instanceCount };

    const scale = getMaxSize(size) * props.sizeFactor;
    const padding = getPadding(text.textData.ref.value.array, text.charCount, scale);
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), text.boundingSphere, padding);
    const boundingSphere = calculateTransformBoundingSphere(
      invariantBoundingSphere,
      transform.aTransform.ref.value,
      instanceCount,
      0,
    );

    return {
      dGeometryType: ValueCell.create('text'),

      tText: text.textData,
      boundingSphere: ValueCell.create(boundingSphere),
      invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
      uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
      ...color,
      ...size,
      ...marker,
      ...overpaint,
      ...transparency,
      ...emissive,
      ...substance,
      ...clipping,
      ...wiggle,
      ...transform,

      tFont: text.fontTexture,
      padding: ValueCell.create(padding),

      ...BaseGeometry.createValues(props, counts),
      uSizeFactor: ValueCell.create(props.sizeFactor),

      uBorderWidth: ValueCell.create(clamp(props.borderWidth, 0, 0.5)),
      uBorderColor: ValueCell.create(Color.toArrayNormalized(props.borderColor, Vec3.zero(), 0)),
      uOffsetX: ValueCell.create(props.offsetX),
      uOffsetY: ValueCell.create(props.offsetY),
      uOffsetZ: ValueCell.create(props.offsetZ),
      uBackgroundColor: ValueCell.create(Color.toArrayNormalized(props.backgroundColor, Vec3.zero(), 0)),
      uBackgroundOpacity: ValueCell.create(props.backgroundOpacity),
    };
  }

  function createValuesSimple(
    text: Text,
    props: Partial<PD.Values<Params>>,
    colorValue: Color,
    sizeValue: number,
    transform?: TransformData,
  ) {
    const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
    const p = { ...PD.getDefaultValues(Params), ...props };
    return createValues(text, s.transform, s.locationIterator, s.theme, p);
  }

  function updateValues(values: TextValues, props: PD.Values<Params>) {
    BaseGeometry.updateValues(values, props);
    ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);

    ValueCell.updateIfChanged(values.uBorderWidth, props.borderWidth);
    if (Color.fromNormalizedArray(values.uBorderColor.ref.value, 0) !== props.borderColor) {
      Color.toArrayNormalized(props.borderColor, values.uBorderColor.ref.value, 0);
      ValueCell.update(values.uBorderColor, values.uBorderColor.ref.value);
    }
    ValueCell.updateIfChanged(values.uOffsetX, props.offsetX);
    ValueCell.updateIfChanged(values.uOffsetY, props.offsetY);
    ValueCell.updateIfChanged(values.uOffsetZ, props.offsetZ);
    if (Color.fromNormalizedArray(values.uBackgroundColor.ref.value, 0) !== props.backgroundColor) {
      Color.toArrayNormalized(props.backgroundColor, values.uBackgroundColor.ref.value, 0);
      ValueCell.update(values.uBackgroundColor, values.uBackgroundColor.ref.value);
    }
    ValueCell.updateIfChanged(values.uBackgroundOpacity, props.backgroundOpacity);
  }

  function updateBoundingSphere(values: TextValues, text: Text) {
    const scale = getMaxSize(values) * values.uSizeFactor.ref.value;
    const padding = getPadding(values.tText.ref.value.array, text.charCount, scale);
    const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), text.boundingSphere, padding);
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
    state.pickable = false;
    state.opaque = false;
    state.writeDepth = true;
  }
}

/** the largest corner offset (times `scale`) or depth of the `charCount` items of `data` (`TextStride`) */
function getPadding(data: Float32Array, charCount: number, scale: number) {
  let maxOffsetX = 0;
  let maxOffsetY = 0;
  let maxDepth = 0;
  for (let i = 0; i < charCount; ++i) {
    const o = i * TextStride;
    const d = Math.abs(data[o + 4]);
    if (d > maxDepth) maxDepth = d;
    // the corners' x at 8, 10, 12, 14 and y at 9, 11, 13, 15 (a glyph's or background's texture coordinates are not
    // corners: its rectangle is left, top, right, bottom)
    const tether = data[o + 5] > 0.5;
    for (let k = 8; k < (tether ? 16 : 12); k += 2) {
      const ox = Math.abs(data[o + k]);
      if (ox > maxOffsetX) maxOffsetX = ox;
      const oy = Math.abs(data[o + k + 1]);
      if (oy > maxOffsetY) maxOffsetY = oy;
    }
  }
  return Math.max(maxDepth, scale * Math.sqrt(maxOffsetX * maxOffsetX + maxOffsetY * maxOffsetY));
}
