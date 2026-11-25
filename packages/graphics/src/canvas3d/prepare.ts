/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// A canvas's preparation: the pipelines of its frames and first interaction compiled ahead (`prepare`, the warm-up),
// those of objects not yet built (`precompile`, `prepareObjects`), and the pass targets no frame draws released

import { getDenoiseKernels, isComputeDenoiseSupported } from '@molstar/graphics/gl/compute/denoise';
import { prepareLayerSampling } from '@molstar/graphics/gl/compute/layer-sampling';
import { SphereRasterKernels } from '@molstar/graphics/gl/compute/sphere-raster';
import { isSplatSortSupported, SplatSortKernels } from '@molstar/graphics/gl/compute/splat-sort';
import { Object3D } from '@molstar/graphics/gl/object3d';
import {
  createSignatureProgramItem,
  getSignatureKey,
  type ObjectSignature,
  type PrecompileList,
} from '@molstar/graphics/gl/precompile';
import type { ProgramItem } from '@molstar/graphics/gl/render-item';
import { AutoGpuCullingMinRasterSpheres } from '@molstar/graphics/gl/renderable/merged';
import { MergedCullKernels } from '@molstar/graphics/gl/renderable/merged-culler';
import { type GlobalDefines, MaxLightCount } from '@molstar/graphics/gl/renderable/schema';
import type { Scene } from '@molstar/graphics/gl/scene';
import { isDebugMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import type { Canvas3DPrecompileProps, Canvas3DPrepareOptions, Canvas3DProps } from './canvas3d.js';
import { ShaderManager } from './helper/shader-manager.js';
import { isMaterialColorMarker, MarkingPass } from './passes/marking.js';
import { MultiSamplePass } from './passes/multi-sample.js';
import { createPrepareBuilder, type PassBuilder, type PrepareCollect } from './passes/pass-builder.js';
import { createTargetCollector } from './passes/pass-targets.js';
import type { Canvas3DParts } from './parts.js';

export interface Canvas3DPreparer {
  /** Compiles the pipelines frames with the props draw (`Canvas3D.prepare`) */
  prepare: (options?: Canvas3DPrepareOptions) => Promise<void>;
  /** Compiles the pipelines of a list's objects before they are built (`Canvas3D.precompile`) */
  precompile: (list: PrecompileList, propsList?: readonly Canvas3DPrecompileProps[]) => Promise<void>;
  /** Plans the pipelines of objects being built by `owner` (`Canvas3D.prepareObjects`) */
  prepareObjects: (owner: object, list: PrecompileList) => void;
  /** `owner`'s objects were added (or removed): its prepared items may go once committed */
  markPreparedObjectsAdded: (owner: object) => void;
  /** releases the prepared objects' items whose owners were all added and are committed, once the work is ready */
  releasePreparedObjects: () => void;
  /** Releases the pass targets frames with the current props and scene don't draw, a while after they last changed */
  scheduleTargetRelease: () => void;
  /** prepares the frame and the first interaction in idle time, after the scene changed */
  scheduleWarmUp: () => void;
  /** the scheduled warm-up cancelled (it would plan a frame on a destroyed context) */
  cancelWarmUp: () => void;
  /** the program items of the last `precompile` and of prepared objects destroyed (a context replaced, disposed) */
  destroyItems: () => void;
}

/** what `planSignatures` plans: with the interaction (marking and picking), the stand-in scenes' averages */
interface PlanSignaturesOptions {
  readonly interaction: boolean;
  /** the opacity and emissive averages the stand-in scenes take: each branch drawing objects by default */
  readonly opacity?: readonly number[];
  readonly emissive?: readonly number[];
}

/** how long prepared objects' items wait for their owners (a build aborted or failed adds nothing) */
const PreparedObjectsMaxMs = 60_000;

export function createCanvas3DPreparer(parts: Canvas3DParts): Canvas3DPreparer {
  /**
   * the camera frames with `props` draw with: the stereo camera's eyes updated (as `render` does each frame), before
   * the first stereo frame zero matrices (non-invertible)
   */
  function getCamera(props: Canvas3DProps) {
    if (props.camera.stereo.name !== 'on') return parts.camera;
    parts.stereoCamera.update();
    return parts.stereoCamera;
  }

  /** Plans the passes of a frame with `props` (as `render` chooses them, without XR); with `interaction`, hovering */
  function planFrame(b: PassBuilder, props: Canvas3DProps, interaction: boolean) {
    const { renderer, scene, helper, frame, passes, camera, hiZ, pickHelper, multiSampleHelper } = parts;
    const cam = getCamera(props);
    const ctx = { renderer, camera: cam, scene, helper, frame };
    const illuminationEnabled = passes.illumination.supported && props.illumination.enabled;
    if (illuminationEnabled) {
      passes.illumination.plan(b, { ...ctx, camera }, props, true);
    } else {
      if (MultiSamplePass.isEnabled(props.multiSample)) {
        multiSampleHelper.plan(b, ctx, props, true);
      } else {
        passes.draw.plan(b, ctx, props, true);
      }
      hiZ.plan(b, camera, props.hiZ);
    }
    if (interaction) {
      // picking, and marking, which draws only marked objects
      pickHelper.plan(b, cam, frame);
      if (MarkingPass.isEnabled(props.marking)) passes.draw.plan(b, ctx, props, true, false, true);
    }
  }

  function prepare(options?: Canvas3DPrepareOptions): Promise<void> {
    if (parts.contextLost) return Promise.resolve();
    const { graphics, renderer, scene, helper } = parts;
    const props: Canvas3DProps = options?.props ? { ...parts.props, ...options.props } : parts.props;
    // the variants the props require are linked for every render item when set (`ShaderManager.updateRequired`)
    graphics.resources.linkPrograms(ShaderManager.getRequired(graphics, scene, helper.overlay, props));
    if (options?.interaction) {
      // the first marking links the marking variant of every render item (see `ShaderManager.updateRequired`)
      graphics.resources.linkPrograms(MarkingPass.isEnabled(props.marking) ? ['pick', 'marking'] : ['pick']);
    }
    const b = createPrepareBuilder(graphics, renderer);
    planFrame(b, props, !!options?.interaction);
    return b.whenReady();
  }

  /**
   * Releases the storage of the pass targets that frames with the current props and scene don't draw into,
   * hovering included (e.g. of a transparency mode or an effect turned off), a while after the props or the scene
   * last changed (see `PassTargets`)
   */
  function scheduleTargetRelease() {
    const { targets } = parts.passes;
    targets.scheduleRelease(() => {
      if (parts.disposed || parts.contextLost || targets !== parts.passes.targets) return;
      const b = createTargetCollector();
      planFrame(b, parts.props, true);
      targets.releaseUnused(b.used);
    });
  }

  /** the program items of the last `precompile`, by transparency, globals and signature */
  let precompileItems = new Map<string, ProgramItem>();
  function destroyPrecompileItems() {
    for (const item of precompileItems.values()) item.destroy();
    precompileItems.clear();
  }

  /** the scene's globals with `props` (as `setProps` sets them) */
  function getGlobals(props: Canvas3DProps): GlobalDefines {
    const { renderer } = parts;
    return {
      dColorMarker: isMaterialColorMarker(
        { colorMarker: props.renderer?.colorMarker ?? renderer.props.colorMarker },
        props.marking,
      ),
      dLightCount: Math.min(props.renderer?.light?.length ?? renderer.props.light.length, MaxLightCount),
    };
  }

  /**
   * Plans the passes frames with each of `propsList` draw `list`'s objects with, against stand-in scenes, requesting
   * their pipelines: the program items in `items` by transparency, globals and signature (taken from `reuse` where
   * there); whether objects are culled on the GPU
   */
  function planSignatures(
    list: PrecompileList,
    propsList: readonly Canvas3DPrecompileProps[],
    items: Map<string, ProgramItem>,
    reuse: ReadonlyMap<string, ProgramItem> | undefined,
    options: PlanSignaturesOptions,
  ) {
    const { graphics, renderer, scene, helper, frame, passes, camera, hiZ, pickHelper, multiSampleHelper } = parts;
    // objects with a culler (its kernels compiled as it is made, also where it culls on the CPU this frame)
    const gpuCulled = list.objects.some((s) => s.gpuCulled);
    // the passes plan against scenes of the list's objects: groups with their counts (the requests come from the
    // signatures, see `collect`) and averages that take the branches drawing objects (opaque and transparent or
    // only opaque, emissive or not); marking as for the first interaction
    const group = (n: number): Scene.Group => ({ ...Object3D.create(), renderables: new Array(n) });
    for (const { transparency: mode, ...overrides } of propsList) {
      const props: Canvas3DProps = { ...parts.props, ...overrides };
      // the Hi-Z's GPU test where objects cull there with these props (for the raster alone only where it takes them)
      const culled = list.objects.some((s) => isCulledOnGpu(s, props.renderer.sphereRaster));
      const transparency = passes.draw.resolveTransparency(mode ?? passes.draw.transparency);
      const globals = getGlobals(props);
      const primitives: { signature: ObjectSignature; item: ProgramItem }[] = [];
      const volumes: { signature: ObjectSignature; item: ProgramItem }[] = [];
      for (const signature of list.objects) {
        // splats blended whatever the scene's transparency (their own pass, `GaussianSplatsRenderable`)
        const mode = signature.type === 'gaussian-splats' ? 'blended' : transparency;
        const key = `${mode}|${globals.dLightCount}|${globals.dColorMarker}|${getSignatureKey(signature)}`;
        let item = items.get(key) ?? reuse?.get(key);
        if (!item) {
          try {
            item = createSignatureProgramItem(graphics, signature, mode, globals);
          } catch (e) {
            if (isDebugMode) console.warn('Skipped a precompile list entry', signature, e);
          }
          if (!item) continue;
        }
        items.set(key, item);
        (signature.type === 'direct-volume' ? volumes : primitives).push({ signature, item });
      }

      const listPrimitives = group(primitives.length);
      const hasGaussianSplats = primitives.some((p) => p.signature.type === 'gaussian-splats');
      const listVolumes = group(volumes.length);
      // the averages decide between exclusive branches (e.g. bloom composites the background over emissive objects
      // only, postprocessing blends transparency only with transparent objects): planned with none and some
      const listScenes: Scene[] = [];
      for (const opacityAverage of options.opacity ?? [1, 0.5]) {
        for (const emissiveAverage of options.emissive ?? [0, 0.5]) {
          listScenes.push(
            Object.create(scene, {
              primitives: { value: listPrimitives },
              volumes: { value: listVolumes },
              hasOpaque: { value: true },
              hasGaussianSplats: { value: hasGaussianSplats },
              opacityAverage: { value: opacityAverage },
              transparencyMin: { value: 0 },
              emissiveAverage: { value: emissiveAverage },
              markerAverage: { value: 0 },
            }),
          );
        }
      }
      // marking as for the first interaction, with something marked: marking dims the rest (`MarkingPass`, dim)
      const markedScenes = listScenes.map((s) => Object.create(s, { markerAverage: { value: 0.5 } }) as Scene);
      const allItems = [...primitives, ...volumes];
      const collect: PrepareCollect = (pass, g, target, all) => {
        // the scene's own groups too: the pick helper plans with them
        const items =
          g === listPrimitives || g === scene.primitives
            ? primitives
            : g === listVolumes || g === scene.volumes
              ? volumes
              : g === scene || listScenes.includes(g as Scene) || markedScenes.includes(g as Scene)
                ? allItems
                : undefined;
        return items && renderer.collectSignatures(pass, items, target, { all });
      };

      const b = createPrepareBuilder(graphics, renderer, collect);
      const cam = getCamera(props);
      passes.draw.withTransparency(transparency, () => {
        for (let i = 0; i < listScenes.length; ++i) {
          const ctx = { renderer, camera: cam, scene: listScenes[i], helper, frame };
          // as `prepare` with the interaction
          if (passes.illumination.supported && props.illumination.enabled) {
            passes.illumination.plan(b, { ...ctx, camera }, props, true);
          } else {
            if (MultiSamplePass.isEnabled(props.multiSample)) {
              multiSampleHelper.plan(b, ctx, props, true);
            } else {
              passes.draw.plan(b, ctx, props, true);
            }
            // as `planFrame`: with the props' Hi-Z (e.g. the Mesoscale Explorer turns it on), built for the GPU's test
            // where the list has objects culled there
            hiZ.plan(b, camera, props.hiZ, culled);
          }
          if (options.interaction && MarkingPass.isEnabled(props.marking)) {
            passes.draw.plan(b, { ...ctx, scene: markedScenes[i] }, props, true, false, true);
          }
        }
        if (options.interaction) pickHelper.plan(b, cam, frame);
      });
    }
    return gpuCulled;
  }

  /**
   * whether `s`'s objects are culled on the GPU with the raster's mode `mode`: those culled for the raster alone
   * (`rasterSpheres`) only where it takes them (`'on'`, `'only'`; `'auto'` from `AutoGpuCullingMinRasterSpheres`), on
   * the CPU otherwise
   */
  function isCulledOnGpu(s: ObjectSignature, mode: Canvas3DProps['renderer']['sphereRaster']) {
    if (!s.gpuCulled) return false;
    if (s.rasterSpheres === undefined) return true;
    const { graphics } = parts;
    if (mode === 'off' || !graphics.compute || !graphics.capabilities.atomicMinMax64) return false;
    return mode !== 'auto' || s.rasterSpheres >= AutoGpuCullingMinRasterSpheres;
  }

  /**
   * whether the sphere raster may draw spheres of `list` (`RendererProps.sphereRaster`, a device with it; `auto`: of
   * objects culled on the GPU)
   */
  function hasSphereRaster(list: PrecompileList) {
    const { graphics } = parts;
    const mode = parts.renderer.props.sphereRaster;
    if (mode === 'off' || !graphics.compute || !graphics.capabilities.atomicMinMax64) return false;
    return list.objects.some(
      (s) =>
        s.type === 'spheres' &&
        s.opaque &&
        (mode !== 'auto' ||
          (s.gpuCulled === true &&
            (s.rasterSpheres === undefined || s.rasterSpheres >= AutoGpuCullingMinRasterSpheres))),
    );
  }

  /** the GPU sort's kernels where `list` has splats and compute shaders run, else compiled at their first frame */
  function prepareSplatSort(list: PrecompileList) {
    const { compute } = parts.graphics;
    if (!compute || !list.objects.some((s) => s.type === 'gaussian-splats') || !isSplatSortSupported(compute)) return;
    return compute.prepareKernels(SplatSortKernels);
  }

  /**
   * the per-vertex sampling's kernels where `list` has merged meshes, which may smooth a layer per vertex and instance
   * on the GPU where compute shaders run (their visual compiles them once it smooths; whether it does is not in the
   * signature), else compiled once such a visual is made
   */
  function prepareLayerSamplingOf(list: PrecompileList) {
    const merged = list.objects.some((s) => s.type === 'mesh' && s.merged);
    return merged ? prepareLayerSampling(parts.graphics) : undefined;
  }

  function precompile(list: PrecompileList, propsList: readonly Canvas3DPrecompileProps[] = [{}]): Promise<void> {
    if (parts.contextLost) return Promise.resolve();
    const { graphics, passes } = parts;
    const previous = precompileItems;
    precompileItems = new Map();
    const gpuCulled = planSignatures(list, propsList, precompileItems, previous, { interaction: true });
    for (const [key, item] of previous) {
      if (!precompileItems.has(key)) item.destroy();
    }
    // the culling's kernels of objects culled on the GPU, else compiled once such an object is created
    const kernels = graphics.compute && gpuCulled ? graphics.compute.prepareKernels(MergedCullKernels) : undefined;
    // the sphere raster's kernels where it may draw the list's spheres, else compiled at its first frame
    const raster = hasSphereRaster(list) ? graphics.compute!.prepareKernels(SphereRasterKernels) : undefined;
    // the illumination's denoiser, a kernel where compute shaders run
    const denoise =
      passes.illumination.supported &&
      isComputeDenoiseSupported(graphics) &&
      propsList.some((o) => {
        const { enabled, denoise } = { ...parts.props.illumination, ...o.illumination };
        return enabled && denoise;
      })
        ? graphics.compute.prepareKernels(getDenoiseKernels())
        : undefined;
    const splatSort = prepareSplatSort(list);
    const layerSampling = prepareLayerSamplingOf(list);
    return Promise.all([
      graphics.resources.pending.whenReady(),
      kernels,
      raster,
      denoise,
      splatSort,
      layerSampling,
    ]).then(() => {});
  }

  /**
   * Objects being built (`prepareObjects`): the requests until the caller awaits, planned together then, and
   * the program items of those planned, kept until each owner's objects were added and committed and are ready
   */
  let preparingObjects: Map<object, PrecompileList> | undefined;
  const preparedObjects: {
    owners: Set<object>;
    added: Set<object>;
    items: Map<string, ProgramItem>;
    time: number;
  }[] = [];
  let preparedObjectsTimer: ReturnType<typeof setTimeout> | undefined;

  function prepareObjects(owner: object, list: PrecompileList) {
    if (!parts.graphics.capabilities.parallelCompile || parts.contextLost || parts.disposed) return;
    if (list.objects.length === 0) return;
    if (!preparingObjects) {
      preparingObjects = new Map();
      // at the builder's first await: a timer could wait until a build that never yields to it is committed and drawn
      queueMicrotask(planPreparedObjects);
    }
    preparingObjects.set(owner, list);
  }

  function planPreparedObjects() {
    const requests = preparingObjects;
    preparingObjects = undefined;
    if (!requests || parts.disposed || parts.contextLost) return;
    const { graphics, scene } = parts;
    const objects: ObjectSignature[] = [];
    requests.forEach((l) => objects.push(...l.objects));
    // the branches the objects and the scene they join take (not every one, as `precompile`): the frame's passes
    // with the current props; the interaction is the warm-up's, once they are drawn
    const opacity = new Set<number>();
    const emissive = new Set<number>([0]);
    for (const o of objects) {
      if (o.opaque) opacity.add(1);
      if (o.transparent) opacity.add(0.5);
      if (o.emissive) emissive.add(0.5);
    }
    if (scene.opacityAverage < 1) opacity.add(0.5);
    if (scene.emissiveAverage > 0) emissive.add(0.5);
    if (opacity.size === 0) opacity.add(1);
    const reuse = new Map<string, ProgramItem>();
    for (const entry of preparedObjects) entry.items.forEach((item, key) => reuse.set(key, item));
    const items = new Map<string, ProgramItem>();
    const gpuCulled = planSignatures({ objects }, [{}], items, reuse, {
      interaction: false,
      opacity: [...opacity],
      emissive: [...emissive],
    });
    // the culling's kernels of objects culled on the GPU, as `precompile`
    if (graphics.compute && gpuCulled) graphics.compute.prepareKernels(MergedCullKernels);
    if (hasSphereRaster({ objects })) graphics.compute!.prepareKernels(SphereRasterKernels);
    prepareSplatSort({ objects });
    prepareLayerSamplingOf({ objects });
    // the items this request uses too move to it: with an earlier one's they could go before its owners are added
    for (const entry of preparedObjects) for (const key of items.keys()) entry.items.delete(key);
    preparedObjects.push({ owners: new Set(requests.keys()), added: new Set(), items, time: now() });
    // a build aborted or failed commits nothing
    preparedObjectsTimer ??= setTimeout(onPreparedObjectsTimer, PreparedObjectsMaxMs);
  }

  function onPreparedObjectsTimer() {
    preparedObjectsTimer = undefined;
    releasePreparedObjects();
    if (preparedObjects.length > 0) preparedObjectsTimer = setTimeout(onPreparedObjectsTimer, PreparedObjectsMaxMs);
  }

  function releasePreparedObjects() {
    const done: Map<string, ProgramItem>[] = [];
    for (let i = preparedObjects.length - 1; i >= 0; --i) {
      const entry = preparedObjects[i];
      if (entry.added.size < entry.owners.size && now() - entry.time < PreparedObjectsMaxMs) continue;
      done.push(entry.items);
      preparedObjects.splice(i, 1);
    }
    if (done.length === 0) return;
    // the committed objects hold the programs by then
    parts.graphics.resources.pending.whenReady().then(() => {
      for (const items of done) for (const item of items.values()) item.destroy();
    });
  }

  function markPreparedObjectsAdded(owner: object) {
    for (const entry of preparedObjects) {
      if (entry.owners.has(owner)) entry.added.add(owner);
    }
  }

  function destroyPreparedObjects() {
    preparingObjects = undefined;
    clearTimeout(preparedObjectsTimer);
    preparedObjectsTimer = undefined;
    for (const entry of preparedObjects) for (const item of entry.items.values()) item.destroy();
    preparedObjects.length = 0;
  }

  /** the scheduled warm-up, cancelled by `dispose` (it would plan a frame on the destroyed context) */
  let warmUpCancel: (() => void) | undefined;
  function scheduleWarmUp() {
    if (warmUpCancel || !parts.attribs.warmUp || parts.disposed) return;
    const warmUp = () => {
      warmUpCancel = undefined;
      if (parts.disposed) return;
      // the frame's pipelines first, the interaction's once they are ready: compiling together, they share the
      // compiler's threads, and the frame waits for its own
      prepare().then(() => {
        if (!parts.disposed) prepare({ interaction: true });
      });
      scheduleTargetRelease();
    };
    if (typeof requestIdleCallback !== 'undefined') {
      const handle = requestIdleCallback(warmUp, { timeout: 1000 });
      warmUpCancel = () => cancelIdleCallback(handle);
    } else {
      const handle = setTimeout(warmUp, 0);
      warmUpCancel = () => clearTimeout(handle);
    }
  }

  return {
    prepare,
    precompile,
    prepareObjects,
    markPreparedObjectsAdded,
    releasePreparedObjects,
    scheduleTargetRelease,
    scheduleWarmUp,
    cancelWarmUp: () => {
      warmUpCancel?.();
      warmUpCancel = undefined;
    },
    destroyItems: () => {
      destroyPrecompileItems();
      destroyPreparedObjects();
    },
  };
}
