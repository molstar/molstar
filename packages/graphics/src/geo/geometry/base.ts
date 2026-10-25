/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RenderableState } from '@molstar/graphics/gl/renderable';
import { ValueCell } from '@molstar/core/util';
import { clamp } from '@molstar/core/math/interpolate';
import type { BaseValues } from '@molstar/graphics/gl/renderable/schema';
import type { DrawIndirect } from '@molstar/graphics/gl/compute';
import { LocationIterator } from '../util/location-iterator.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type TransformData, createIdentityTransform } from './transform-data.js';
import type { Theme } from '@molstar/graphics/theme/theme';
import { ColorNames } from '@molstar/core/util/color/names';
import { NullLocation } from '@molstar/model/model/location';
import { UniformColorTheme } from '@molstar/graphics/theme/color/uniform';
import { UniformSizeTheme } from '@molstar/graphics/theme/size/uniform';
import { smoothstep } from '@molstar/core/math/interpolate';
import { Material } from '@molstar/core/util/material';
import { Clip } from '@molstar/core/util/clip';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';

export const VisualQualityInfo = {
  custom: {},
  auto: {},
  highest: {},
  higher: {},
  high: {},
  medium: {},
  low: {},
  lower: {},
  lowest: {},
};
export type VisualQuality = keyof typeof VisualQualityInfo;
export const VisualQualityNames = Object.keys(VisualQualityInfo) as VisualQuality[];
export const VisualQualityOptions = PD.arrayToOptions(VisualQualityNames);

//

export const ColorSmoothingParams = {
  smoothColors: PD.MappedStatic('auto', {
    auto: PD.Group({}),
    on: PD.Group({
      resolutionFactor: PD.Numeric(2, { min: 0.5, max: 6, step: 0.1 }),
      sampleStride: PD.Numeric(3, { min: 1, max: 12, step: 1 }),
    }),
    off: PD.Group({}),
  }),
};
export type ColorSmoothingParams = typeof ColorSmoothingParams;

export function hasColorSmoothingProp(props: PD.Values<any>): props is PD.Values<ColorSmoothingParams> {
  return !!props.smoothColors;
}

export function getColorSmoothingProps(
  smoothColors: PD.Values<ColorSmoothingParams>['smoothColors'],
  preferSmoothing?: boolean,
  resolution?: number,
) {
  if (
    (smoothColors.name === 'on' || (smoothColors.name === 'auto' && preferSmoothing)) &&
    resolution &&
    resolution < 3
  ) {
    let stride = 3;
    if (smoothColors.name === 'on') {
      resolution *= smoothColors.params.resolutionFactor;
      stride = smoothColors.params.sampleStride;
    } else {
      // https://graphtoy.com/?f1(x,t)=(2-smoothstep(0,1.1,x))*x&coords=0.7,0.6,1.8
      resolution *= 2 - smoothstep(0, 1.1, resolution);
      resolution = Math.max(0.5, resolution);
      if (resolution > 1.2) stride = 2;
    }
    return { resolution, stride };
  }
}

//

export type InstanceGranularityValue = true | false | 'auto';
export const InstanceGranularityOptions: [InstanceGranularityValue, string][] = [
  [true, 'On'],
  [false, 'Off'],
  ['auto', 'Auto'],
];

/**
 * Threshold (in `groupCount * instanceCount`, e.g. number of marker-texture
 * slots) above which `instanceGranularity: 'auto'` resolves to `true`.
 */
export const AutoInstanceGranularityThreshold = 50_000_000;

/**
 * Resolves the `instanceGranularity` param value to a boolean.
 */
export function resolveInstanceGranularity(
  value: InstanceGranularityValue,
  groupCount: number,
  instanceCount: number,
): boolean {
  if (value === 'auto') return groupCount * instanceCount > AutoInstanceGranularityThreshold;
  return value;
}

//

export type GpuCullingValue = true | false | 'auto';
export const GpuCullingOptions: [GpuCullingValue, string][] = [
  [true, 'On'],
  [false, 'Off'],
  ['auto', 'Auto'],
];

/**
 * Instances per draw (a segment and level) from which `gpuCulling: 'auto'` culls on the GPU: it replays a draw per
 * segment and level every frame (about 0.28 µs each), the CPU culling costs about 10 ns per instance in frames whose
 * view changed (`merged-render-objects.md`, step 6)
 */
export const AutoGpuCullingThreshold = 30;

/**
 * Instances from which `gpuCulling: 'auto'` culls on the GPU: a GPU-culled object costs about 0.1 ms per frame of its
 * own (its dispatches, the Hi-Z levels it asks for), the CPU's 10 ns per instance (`merging/instances-*`)
 */
export const AutoGpuCullingMinInstances = 20000;

/** Resolves the `gpuCulling` param value for `instanceCount` instances drawn by `drawCount` segments times levels */
export function resolveGpuCulling(value: GpuCullingValue, instanceCount: number, drawCount: number): boolean {
  if (value === 'auto') {
    return instanceCount >= Math.max(AutoGpuCullingMinInstances, AutoGpuCullingThreshold * Math.max(1, drawCount));
  }
  return value;
}

export type RenderMergedValue = true | false | 'auto';
export const RenderMergedOptions: [RenderMergedValue, string][] = [
  [true, 'On'],
  [false, 'Off'],
  ['auto', 'Auto'],
];

/**
 * Threshold (in number of mergeable units) above which `renderMerged: 'auto'`
 * resolves to `true`.
 */
export const AutoRenderMergedThreshold = 30;

/**
 * Resolves the `renderMerged` param value to a boolean.
 */
export function resolveRenderMerged(value: RenderMergedValue, unitCount: number): boolean {
  if (value === 'auto') return unitCount > AutoRenderMergedThreshold;
  return value;
}

//

export namespace BaseGeometry {
  export const MaterialCategory: PD.Info = { category: 'Material' };
  export const ShadingCategory: PD.Info = { category: 'Shading' };
  export const CullingLodCategory: PD.Info = { category: 'Culling & LOD' };
  export const CustomQualityParamInfo: PD.Info = {
    category: 'Custom Quality',
    hideIf: (params: PD.Values<Params>) => typeof params.quality !== 'undefined' && params.quality !== 'custom',
  };

  export const Params = {
    alpha: PD.Numeric(
      1,
      { min: 0, max: 1, step: 0.01 },
      { label: 'Opacity', isEssential: true, description: 'How opaque/transparent the representation is rendered.' },
    ),
    quality: PD.Select<VisualQuality>('auto', VisualQualityOptions, {
      isEssential: true,
      description: 'Visual/rendering quality of the representation.',
    }),
    material: Material.getParam(),
    clip: PD.Group(Clip.Params),
    emissive: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    density: PD.Numeric(
      0.2,
      { min: 0, max: 1, step: 0.01 },
      { description: 'Density value to estimate object thickness.' },
    ),
    instanceGranularity: PD.Select<InstanceGranularityValue>('auto', InstanceGranularityOptions, {
      description:
        'Use instance granularity for marker, transparency, clipping, overpaint, substance data to save memory. When set to `auto`, granularity is enabled if `groupCount * instanceCount` exceeds `AutoInstanceGranularityThreshold`.',
    }),
    renderMerged: PD.Select<RenderMergedValue>('auto', RenderMergedOptions, {
      description:
        'Render all units of the representation through a single shared render item per geometry kind to reduce draw calls. When set to `auto`, enabled if the number of mergeable units exceeds `AutoRenderMergedThreshold`. Supported for spheres, mesh, cylinders, lines, and points geometry.',
    }),
    gpuCulling: PD.Select<GpuCullingValue>('auto', GpuCullingOptions, {
      ...CullingLodCategory,
      description:
        'Cull and choose the level of detail per cell of instances on the GPU where the device supports it (WebGPU with compute and indirect draws by first instance), merging the render objects for it. When set to `auto`, enabled with many instances (`AutoGpuCullingMinInstances`), and many per draw (`AutoGpuCullingThreshold`).',
    }),
    lod: PD.Vec3(Vec3(), undefined, {
      ...CullingLodCategory,
      description: 'Level of detail.',
      fieldLabels: { x: 'Min Distance', y: 'Max Distance', z: 'Overlap (Shader)' },
    }),
    cellSize: PD.Numeric(
      200,
      { min: 0, max: 5000, step: 100 },
      { ...CullingLodCategory, description: 'Instance grid cell size.' },
    ),
    batchSize: PD.Numeric(
      2000,
      { min: 0, max: 50000, step: 500 },
      { ...CullingLodCategory, description: 'Instance grid batch size.' },
    ),
  };
  export type Params = typeof Params;

  export type Counts = { drawCount: number; vertexCount: number; groupCount: number; instanceCount: number };

  export function createSimple(colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData) {
    if (!transform) transform = createIdentityTransform();
    const locationIterator = LocationIterator(
      1,
      transform.instanceCount.ref.value,
      1,
      () => NullLocation,
      false,
      () => false,
    );
    const theme: Theme = {
      color: UniformColorTheme({}, { value: colorValue, lightness: 0, saturation: 0 }),
      size: UniformSizeTheme({}, { value: sizeValue }),
    };
    return { transform, locationIterator, theme };
  }

  export function createValues(props: PD.Values<Params>, counts: Counts) {
    const clip = Clip.getClip(props.clip);
    return {
      alpha: ValueCell.create(props.alpha),
      uAlpha: ValueCell.create(clamp(props.alpha, 0, 1)),
      uVertexCount: ValueCell.create(counts.vertexCount),
      // filled by the render item when used, see `BaseSchema.aVertex`
      aVertex: ValueCell.create(new Float32Array(0)),
      uGroupCount: ValueCell.create(counts.groupCount),
      drawCount: ValueCell.create(counts.drawCount),
      drawIndirect: ValueCell.create<DrawIndirect | undefined>(undefined),
      uMetalness: ValueCell.create(props.material.metalness),
      uRoughness: ValueCell.create(props.material.roughness),
      uBumpiness: ValueCell.create(props.material.bumpiness),
      uEmissive: ValueCell.create(props.emissive),
      uDensity: ValueCell.create(props.density),

      dClipObjects: ValueCell.create(clip.objects.count > 0),
      dClipObjectCount: ValueCell.create(clip.objects.count),
      uClipObjectCount: ValueCell.create(clip.objects.count),
      dClipVariant: ValueCell.create(clip.variant),
      clipObjects: ValueCell.create(clip.objects),
      uClipObjectFrame: ValueCell.create(clip.objects.frame),
      uClipObjectShape: ValueCell.create(clip.objects.shape),

      uInstanceGranularity: ValueCell.create(
        resolveInstanceGranularity(props.instanceGranularity, counts.groupCount, counts.instanceCount),
      ),
      uLod: ValueCell.create(Vec4.create(props.lod[0], props.lod[1], props.lod[2], 0)),
    };
  }

  /** The clip props each values were updated with last (`updateValues`) */
  const clipProps = new WeakMap<BaseValues, Clip.Props>();

  export function updateValues(values: BaseValues, props: PD.Values<Params>) {
    ValueCell.updateIfChanged(values.alpha, props.alpha);
    ValueCell.updateIfChanged(values.uAlpha, clamp(props.alpha, 0, 1));
    ValueCell.updateIfChanged(values.uMetalness, props.material.metalness);
    ValueCell.updateIfChanged(values.uRoughness, props.material.roughness);
    ValueCell.updateIfChanged(values.uBumpiness, props.material.bumpiness);
    ValueCell.updateIfChanged(values.uEmissive, props.emissive);
    ValueCell.updateIfChanged(values.uDensity, props.density);

    // the clip of the props applied last: updates keep the props' clip unless it changed (e.g. a recolor of many units)
    if (clipProps.get(values) !== props.clip) {
      clipProps.set(values, props.clip);
      const clip = Clip.getClip(props.clip);
      ValueCell.updateIfChanged(values.dClipObjects, clip.objects.count > 0);
      ValueCell.updateIfChanged(values.dClipObjectCount, clip.objects.count);
      ValueCell.updateIfChanged(values.uClipObjectCount, clip.objects.count);
      ValueCell.updateIfChanged(values.dClipVariant, clip.variant);
      ValueCell.update(values.clipObjects, clip.objects);
      ValueCell.update(values.uClipObjectFrame, clip.objects.frame);
      ValueCell.update(values.uClipObjectShape, clip.objects.shape);
    }

    ValueCell.updateIfChanged(
      values.uInstanceGranularity,
      resolveInstanceGranularity(
        props.instanceGranularity,
        values.uGroupCount.ref.value,
        values.instanceCount.ref.value,
      ),
    );
    ValueCell.update(values.uLod, Vec4.set(values.uLod.ref.value, props.lod[0], props.lod[1], props.lod[2], 0));
  }

  export function createRenderableState(props: Partial<PD.Values<Params>> = {}): RenderableState {
    const opaque = props.alpha === undefined ? true : props.alpha === 1;
    return RenderableState.create({
      disposed: false,
      visible: true,
      pickable: true,
      colorOnly: false,
      opaque,
      writeDepth: opaque,
    });
  }

  export function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
    state.opaque = props.alpha >= 1;
    state.writeDepth = state.opaque;
  }
}
