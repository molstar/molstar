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
import type { InteriorProps } from './interior.js';

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
 * Instances from which `gpuCulling: 'auto'` culls on the GPU: a GPU-culled object costs about 0.1-0.3 ms per frame of
 * its own (its dispatches, the Hi-Z levels it asks for, a draw per segment and level), the CPU's 10 ns per instance in
 * frames whose view changed. 100 chains as spacefill rotating, CPU / GPU: 30k instances 0.52 / 0.78 ms, 50k 0.75 /
 * 0.75, 100k 1.40 / 0.84; as cartoon the CPU ahead up to 100k (1.25 / 1.32) (`merging/instances-*`)
 */
export const AutoGpuCullingMinInstances = 50000;

/**
 * Spheres from which `gpuCulling: 'auto'` culls on the GPU where the device has the sphere raster, which takes merged
 * objects culled there only: faster at 300k spheres and up on both GPUs measured (3J3Q 1.48 -> 0.63 ms, 61.9 -> 11.5 ms on
 * an integrated GPU), about 0.1 ms slower at a few hundred (`sphere-compute-raster.md`). Such an object is merged for
 * it, the renderer culling it on the GPU only where the raster takes it (`Renderable.rasterSpheres`,
 * `sphereRasterObjects`), else on the CPU
 */
export const AutoGpuCullingMinRasterSpheres = 100_000;

/**
 * Resolves the `gpuCulling` param value for `instanceCount` instances drawn by `drawCount` segments times levels;
 * `rasterSpheres`: the spheres of those instances the sphere raster could draw, 0 where the device has no raster
 */
export function resolveGpuCulling(
  value: GpuCullingValue,
  instanceCount: number,
  drawCount: number,
  rasterSpheres = 0,
): boolean {
  if (value === 'auto') {
    if (rasterSpheres >= AutoGpuCullingMinRasterSpheres) return true;
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
  /** Whether objects fade at their lod distances (`dLodFade`): `lod` with a min or max distance */
  export function hasLodFade(lod: Vec3): boolean {
    return lod[0] !== 0 || lod[1] !== 0;
  }

  /**
   * Whether bumps show by the props (`dBump`): frequency, amplitude and bumpiness above 0, the bumpiness of the material
   * or of the interior's substance; a substance layer (`dSubstance`) can add bumpiness per group as well
   */
  export function hasBump(props: {
    bumpFrequency: number;
    bumpAmplitude: number;
    material: Material;
    interior: InteriorProps;
  }): boolean {
    if (!(props.bumpFrequency > 0 && props.bumpAmplitude > 0)) return false;
    const { substance, substanceStrength } = props.interior;
    return props.material.bumpiness > 0 || (substance.bumpiness > 0 && substanceStrength > 0);
  }

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
    occludedAlpha: PD.Numeric(
      1,
      { min: 0, max: 1, step: 0.01 },
      {
        isHidden: true,
        description:
          'In the overlay layer: the opacity of the parts behind the scene (1 drawn on top, 0 hidden by the scene, ' +
          'between x-ray).',
      },
    ),
    clip: PD.Group(Clip.Params),
    clipGhost: PD.Boolean(true, {
      description: 'What the clip objects remove drawn as a ghost, where the canvas draws ghosts (`clipGhost`).',
    }),
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
        "Cull and choose the level of detail per cell of instances on the GPU where the device supports it (WebGPU with compute and indirect draws by first instance), merging the render objects for it. When set to `auto`, enabled with many instances (`AutoGpuCullingMinInstances`), and many per draw (`AutoGpuCullingThreshold`), or with many spheres where the device has the sphere raster, which takes such objects only (`AutoGpuCullingMinRasterSpheres`): those culled on the GPU only where the raster takes them (the renderer's `sphereRaster` and `sphereRasterObjects`), else on the CPU.",
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
    const clip = getSharedClip(props.clip);
    return {
      alpha: ValueCell.create(props.alpha),
      uAlpha: ValueCell.create(clamp(props.alpha, 0, 1)),
      // set by the canvas for its overlay layer
      dOverlay: ValueCell.create(false),
      uOccludedAlpha: ValueCell.create(clamp(props.occludedAlpha, 0, 1)),
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
      clipGhost: ValueCell.create(props.clipGhost),

      uInstanceGranularity: ValueCell.create(
        resolveInstanceGranularity(props.instanceGranularity, counts.groupCount, counts.instanceCount),
      ),
      uLod: ValueCell.create(Vec4.create(props.lod[0], props.lod[1], props.lod[2], 0)),
      dLodFade: ValueCell.create(hasLodFade(props.lod)),
    };
  }

  /** The clip props each values were updated with last (`updateValues`) */
  const clipProps = new WeakMap<BaseValues, Clip.Props>();
  /**
   * The clip of a clip props object, made once: the visuals of a representation share their props' clip (each target of
   * a particle list, each unit), so a moved clip object is computed once, not per visual. Read only.
   */
  const clips = new WeakMap<Clip.Props, Clip>();
  function getSharedClip(props: Clip.Props) {
    let clip = clips.get(props);
    if (!clip) {
      clip = Clip.getClip(props);
      clips.set(props, clip);
    }
    return clip;
  }

  /** Whether `a` and `b` set the same clip defines (`dClipObjects`, `dClipObjectCount`, `dClipVariant`) */
  export function isSameClipDefines(a: Clip.Props, b: Clip.Props) {
    const count = (p: Clip.Props) => Math.min(p.objects.length, Clip.MaxObjectCount);
    return a.variant === b.variant && count(a) === count(b);
  }

  /** The clip values of `values` for `clip` (as `updateValues` sets them), nothing else */
  export function updateClip(values: BaseValues, clip: Clip.Props) {
    if (clipProps.get(values) === clip) return;
    clipProps.set(values, clip);
    const c = getSharedClip(clip);
    ValueCell.updateIfChanged(values.dClipObjects, c.objects.count > 0);
    ValueCell.updateIfChanged(values.dClipObjectCount, c.objects.count);
    ValueCell.updateIfChanged(values.uClipObjectCount, c.objects.count);
    ValueCell.updateIfChanged(values.dClipVariant, c.variant);
    ValueCell.update(values.clipObjects, c.objects);
    ValueCell.update(values.uClipObjectFrame, c.objects.frame);
    ValueCell.update(values.uClipObjectShape, c.objects.shape);
  }

  export function updateValues(values: BaseValues, props: PD.Values<Params>) {
    ValueCell.updateIfChanged(values.alpha, props.alpha);
    ValueCell.updateIfChanged(values.uAlpha, clamp(props.alpha, 0, 1));
    ValueCell.updateIfChanged(values.uOccludedAlpha, clamp(props.occludedAlpha, 0, 1));
    ValueCell.updateIfChanged(values.uMetalness, props.material.metalness);
    ValueCell.updateIfChanged(values.uRoughness, props.material.roughness);
    ValueCell.updateIfChanged(values.uBumpiness, props.material.bumpiness);
    ValueCell.updateIfChanged(values.uEmissive, props.emissive);
    ValueCell.updateIfChanged(values.uDensity, props.density);
    ValueCell.updateIfChanged(values.clipGhost, props.clipGhost);

    // the clip of the props applied last: updates keep the props' clip unless it changed (e.g. a recolor of many units)
    updateClip(values, props.clip);

    ValueCell.updateIfChanged(
      values.uInstanceGranularity,
      resolveInstanceGranularity(
        props.instanceGranularity,
        values.uGroupCount.ref.value,
        values.instanceCount.ref.value,
      ),
    );
    ValueCell.update(values.uLod, Vec4.set(values.uLod.ref.value, props.lod[0], props.lod[1], props.lod[2], 0));
    ValueCell.updateIfChanged(values.dLodFade, hasLodFade(props.lod));
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
