/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Herman Bergwerf <post@hbergwerf.nl>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { BehaviorSubject, type Observable, Subject, type Subscription, debounceTime, merge } from 'rxjs';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { type GraphicsRenderObject, hasForeignResources } from '@molstar/graphics/gl/render-object';
import { type Frame, type GraphicsRenderable, createFrame } from '@molstar/graphics/gl/renderable';
import { Renderer, RendererParams, type RendererStats } from '@molstar/graphics/gl/renderer';
import { MaxLightCount } from '@molstar/graphics/gl/renderable/schema';
import { getVertexArenaStats } from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { Scene } from '@molstar/graphics/gl/scene';
import {
  createSignatureProgramItem,
  getSignatureKey,
  type ObjectSignature,
  type PrecompileList,
  PrecompileRecorder,
} from '@molstar/graphics/gl/precompile';
import type { ProgramItem } from '@molstar/graphics/gl/render-item';
import type { GlobalDefines } from '@molstar/graphics/gl/renderable/schema';
import { MergedCullKernels } from '@molstar/graphics/gl/renderable/merged-culler';
import { getDenoiseKernels, isComputeDenoiseSupported } from '@molstar/graphics/gl/compute/denoise';
import { Object3D } from '@molstar/graphics/gl/object3d';
import {
  createGraphicsContext,
  type GraphicsBackendDescriptor,
  type GraphicsContextAttribs,
} from '@molstar/graphics/gl/backend';
import type { GraphicsContext, GraphicsDestroyOptions } from '@molstar/graphics/gl/context';
import { type EasingFunction, EasingParamDefinition } from '@molstar/core/math/easing';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { Ray3D } from '@molstar/core/math/geometry/primitives/ray3d';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { degToRad, radToDeg } from '@molstar/core/math/misc';
import { EmptyLoci, Loci, isEmptyLoci } from '@molstar/model/model/loci';
import { Representation } from '@molstar/graphics/repr/representation';
import { shallowEqual } from '@molstar/core/util';
import type { AssetManager } from '@molstar/core/util/assets';
import { isMobileBrowser } from '@molstar/core/util/browser';
import {
  addConsoleStatsProvider,
  isDebugMode,
  isTimingMode,
  removeConsoleStatsProvider,
} from '@molstar/core/util/debug';
import { type ButtonsType, InputObserver, type ModifiersKeys } from '@molstar/core/util/input/input-observer';
import type { MarkerAction } from '@molstar/core/util/marker-action';
import { now } from '@molstar/core/util/now';
import { toFixed } from '@molstar/core/util/number';
import { deepClone } from '@molstar/core/util/object';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { produce } from '@molstar/core/util/produce';
import { SetUtils } from '@molstar/core/util/set';
import { Camera } from './camera.js';
import { StereoCamera, StereoCameraParams } from './camera/stereo.js';
import type { CameraTransitionManager } from './camera/transition.js';
import { TransitionTrajectoryParamDefinition, type TransitionTrajectory } from './camera/transition-functions.js';
import { Viewport } from './camera/util.js';
import { MarkIndex } from './mark-index.js';
import { DefaultTrackballControlsAttribs, TrackballControls, TrackballControlsParams } from './controls/trackball.js';
import { CameraHelperParams } from './helper/camera-helper.js';
import type { DebugRegistry } from './helper/debug-registry.js';
import { HandleHelperParams } from './helper/handle-helper.js';
import { Helper } from './helper/helper.js';
import { Canvas3dInteractionHelper, Canvas3dInteractionHelperParams } from './helper/interaction-events.js';
import { PickHelper } from './helper/pick-helper.js';
import { PointerHelperParams } from './helper/pointer-helper.js';
import { RayHelper } from './helper/ray-helper.js';
import { ShaderManager } from './helper/shader-manager.js';
import { DefaultXRManagerAttribs, XRManager, XRManagerParams } from './helper/xr-manager.js';
import { HiZParams, HiZPass } from './passes/hi-z.js';
import { IlluminationParams } from './passes/illumination.js';
import { ImagePass, type ImageProps } from './passes/image.js';
import { isMaterialColorMarker, MarkingParams, MarkingPass, SingleSample } from './passes/marking.js';
import { createPrepareBuilder, type PassBuilder, type PrepareCollect } from './passes/pass-builder.js';
import { createTargetCollector } from './passes/pass-targets.js';
import type { TransparencyMode } from './passes/draw.js';
import { createSupportedPropsLimiter } from './passes/supported-props.js';
import { MultiSampleHelper, MultiSampleParams, MultiSamplePass } from './passes/multi-sample.js';
import { getJitterOffsets, getTemporalSamplesPerFrame } from './passes/jitter.js';
import { Passes } from './passes/passes.js';
import { type AsyncPickData, DefaultPickOptions, type PickResult } from './passes/pick.js';
import { PostprocessingParams } from './passes/postprocessing.js';

export const CameraFogParams = {
  intensity: PD.Numeric(15, { min: 1, max: 100, step: 1 }),
};
export const Canvas3DParams = {
  camera: PD.Group(
    {
      mode: PD.Select('perspective', PD.arrayToOptions(['perspective', 'orthographic'] as const), { label: 'Camera' }),
      helper: PD.Group(CameraHelperParams, { isFlat: true }),
      stereo: PD.MappedStatic(
        'off',
        {
          on: PD.Group(StereoCameraParams),
          off: PD.Group({}),
        },
        { cycle: true, hideIf: (p) => p?.mode !== 'perspective' },
      ),
      fov: PD.Numeric(45, { min: 10, max: 130, step: 1 }, { label: 'Field of View' }),
      manualReset: PD.Boolean(false, { isHidden: true }),
    },
    { pivot: 'mode' },
  ),
  cameraFog: PD.MappedStatic(
    'on',
    {
      on: PD.Group(CameraFogParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Show fog in the distance' },
  ),
  cameraClipping: PD.Group(
    {
      radius: PD.Numeric(
        100,
        { min: 0, max: 99, step: 1 },
        { label: 'Clipping', description: 'How much of the scene to show.' },
      ),
      far: PD.Boolean(true, { description: 'Hide scene in the distance' }),
      minNear: PD.Numeric(
        1,
        { min: 0.1, max: 100, step: 0.1 },
        {
          description:
            'Minimal allowed distance of near clipping plane from the camera. Note, may cause performance issues rendering impostors when set too small and cause issues with outline rendering when too close to 0.',
        },
      ),
    },
    { pivot: 'radius' },
  ),
  viewport: PD.MappedStatic('canvas', {
    canvas: PD.Group({}),
    'static-frame': PD.Group({
      x: PD.Numeric(0),
      y: PD.Numeric(0),
      width: PD.Numeric(128),
      height: PD.Numeric(128),
    }),
    'relative-frame': PD.Group({
      x: PD.Numeric(0.33, { min: 0, max: 1, step: 0.01 }),
      y: PD.Numeric(0.33, { min: 0, max: 1, step: 0.01 }),
      width: PD.Numeric(0.5, { min: 0.01, max: 1, step: 0.01 }),
      height: PD.Numeric(0.5, { min: 0.01, max: 1, step: 0.01 }),
    }),
  }),

  cameraResetDurationMs: PD.Numeric(
    250,
    { min: 0, max: 1000, step: 1 },
    { description: 'The time it takes to reset the camera.' },
  ),
  cameraResetEasing: EasingParamDefinition('linear'),
  cameraResetTrajectory: TransitionTrajectoryParamDefinition('linear'),
  sceneRadiusFactor: PD.Numeric(1, { min: 1, max: 10, step: 0.1 }),
  transparentBackground: PD.Boolean(false),
  checkeredTransparentBackground: PD.Boolean(false),
  dpoitIterations: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
  volumeResolutionScale: PD.Numeric(
    0.5,
    { min: 0.25, max: 1, step: 0.05 },
    {
      description:
        'Resolution of direct volumes relative to the drawing buffer, upsampled with the depth (volumes are smooth, ' +
        'their cost is per pixel)',
    },
  ),
  pickPadding: PD.Numeric(
    3,
    { min: 0, max: 10, step: 1 },
    { description: 'Extra pixels to around target to check in case target is empty.' },
  ),
  userInteractionReleaseMs: PD.Numeric(
    250,
    { min: 0, max: 1000, step: 1 },
    { description: 'The time before the user is not considered interacting anymore.' },
  ),

  multiSample: PD.Group(MultiSampleParams),
  postprocessing: PD.Group(PostprocessingParams),
  marking: PD.Group(MarkingParams),
  illumination: PD.Group(IlluminationParams),
  hiZ: PD.Group(HiZParams),
  renderer: PD.Group(RendererParams),
  trackball: PD.Group(TrackballControlsParams),
  interaction: PD.Group(Canvas3dInteractionHelperParams),
  handle: PD.Group(HandleHelperParams),
  pointer: PD.Group(PointerHelperParams),
  xr: PD.Group(XRManagerParams, { label: 'XR' }),
};
export const DefaultCanvas3DParams = PD.getDefaultValues(Canvas3DParams);
export type Canvas3DProps = PD.Values<typeof Canvas3DParams>;
export type PartialCanvas3DProps = {
  [K in keyof Canvas3DProps]?: Canvas3DProps[K] extends { name: string; params: any }
    ? Canvas3DProps[K]
    : Partial<Canvas3DProps[K]>;
};

export type Canvas3DPrepareOptions = {
  /** props other than the current ones, e.g. before setting them */
  props?: Partial<Canvas3DProps>;
  /** also what the first interaction draws: picking, and marking as if objects were marked */
  interaction?: boolean;
};

/** Props of frames `Canvas3D.precompile` prepares for, with the context's transparency mode (the current by default) */
export type Canvas3DPrecompileProps = Partial<Canvas3DProps> & { transparency?: TransparencyMode };

export const DefaultCanvas3DAttribs = {
  trackball: DefaultTrackballControlsAttribs,
  xr: DefaultXRManagerAttribs,
  /**
   * prepare the frame and the first interaction (picking, marking) in idle time after the scene changed, so that the
   * first hover doesn't wait for compiling (`Canvas3D.prepare`)
   */
  warmUp: true,
};
export type Canvas3DAttribs = typeof DefaultCanvas3DAttribs;
export type PartialCanvas3DAttribs = {
  [K in keyof Canvas3DAttribs]?: Canvas3DAttribs[K] extends { name: string; params: any }
    ? Canvas3DAttribs[K]
    : Partial<Canvas3DAttribs[K]>;
};

export { Canvas3DContext };

/** Can be used to create multiple Canvas3D objects */
interface Canvas3DContext {
  /** replaced by `setBackend` */
  readonly canvas?: HTMLCanvasElement;
  /** replaced by `setBackend` */
  readonly graphics: GraphicsContext;
  readonly input: InputObserver;
  /** replaced by `setBackend` */
  readonly passes: Passes;
  readonly attribs: Readonly<Canvas3DContext.Attribs>;
  readonly props: Readonly<Canvas3DContext.Props>;
  readonly contextLost?: Subject<now.Timestamp>;
  readonly contextRestored?: Subject<now.Timestamp>;
  /**
   * After `setBackend` replaced `graphics`, `canvas` and `passes`, while the previous context is still alive (it is
   * destroyed once the observers returned): consumers release what they created with it and create it again
   */
  readonly contextChanged?: Subject<Canvas3DContext.ContextChange>;
  readonly assetManager: AssetManager;
  readonly changed?: BehaviorSubject<undefined>;
  readonly pixelScale: number;

  syncPixelScale(): void;
  setProps: (props?: Partial<Canvas3DContext.Props>) => void;
  /**
   * Replaces the graphics context with that of the first of `backends` that can be created, on a new canvas element
   * taking the place of the current one (a canvas keeps the context type it gave out first), see `contextChanged`.
   * Resolves to false, keeping the current context, if none can be created, while presenting XR, without a canvas in
   * a document or while another change is in progress. While released (`release`), the backends are kept for `restore`,
   * which creates the context from them (resolves to true then, the canvas stays released).
   */
  setBackend?: (backends: readonly GraphicsBackendDescriptor[]) => Promise<boolean>;
  /**
   * Releases the graphics context on purpose, freeing its GPU resources (`GraphicsContext.release`, e.g. a canvas out of
   * view): reported as a loss (`contextLost`), but not restored until `restore`. False where it can't (see
   * `GraphicsContext.release`), while presenting XR or changing the backend.
   */
  release?: () => boolean;
  /**
   * Restores a released context (`release`) as after a loss, once the loss was reported (`contextLost`): resolves to
   * whether it was restored (`contextRestored`) within `ReleasedRestoreTimeoutMs`; not restored, it stays released
   */
  restore?: () => Promise<boolean>;
  /** Released (`release`) and not restored yet */
  readonly isReleased?: boolean;
  dispose: (options?: GraphicsDestroyOptions) => void;
}

namespace Canvas3DContext {
  export const DefaultAttribs = {
    powerPreference: 'high-performance' as WebGLContextAttributes['powerPreference'],
    failIfMajorPerformanceCaveat: false,
    /** true by default to avoid issues with Safari (Jan 2021) */
    antialias: true,
    /** true to support multiple Canvas3D objects with a single context */
    preserveDrawingBuffer: true,
    /**
     * Falls back to the next of the backends the context was created from when the context is lost for good (not
     * restored within `ContextRestoreTimeoutMs` while the page is visible, a restore that failed, or more than
     * `MaxContextLosses` losses within `ContextLossWindowMs`); each backend at most once, never back up the list
     */
    contextFallback: true,

    handleResize: () => {},
  };
  export type Attribs = typeof DefaultAttribs;

  export const Params = {
    resolutionMode: PD.Select('auto', PD.arrayToOptions(['auto', 'scaled', 'native'] as const)),
    pixelScale: PD.Numeric(1, { min: 0.1, max: 2, step: 0.05 }),
    pickScale: PD.Numeric(0.25, { min: 0.1, max: 1, step: 0.05 }),
    transparency: PD.Select('wboit', [
      ['blended', 'Blended'],
      ['wboit', 'Weighted, Blended'],
      ['dpoit', 'Depth Peeling'],
    ] as const),
  };
  export const DefaultProps = PD.getDefaultValues(Params);
  export type Props = PD.Values<typeof Params>;

  export interface ContextChange {
    readonly previous: GraphicsContext;
    /**
     * `switch`: `setBackend`; `fallback`: the previous context was lost for good (see `Attribs.contextFallback`);
     * `restore`: the previous context was lost and is replaced by a new one on the same canvas (`GraphicsContext.recreate`)
     */
    readonly reason: 'switch' | 'fallback' | 'restore';
  }

  /** time a lost context has to be restored while the page is visible, before falling back */
  export const ContextRestoreTimeoutMs = 3000;
  /**
   * time a released context has to be restored by `restore` while the page is visible: no fallback waits for it, and a
   * new device may take a while (a busy GPU process took over 6 s)
   */
  export const ReleasedRestoreTimeoutMs = 30_000;
  /** losses within `ContextLossWindowMs` beyond which the context counts as lost for good, even if restored */
  export const MaxContextLosses = 3;
  export const ContextLossWindowMs = 60_000;

  /** Calls `onTimeout` after `ms` of visible page (paused while hidden: browsers delay restores of background tabs) */
  function visibleTimeout(ms: number, onTimeout: () => void) {
    let timer: ReturnType<typeof setTimeout> | undefined;
    let remaining = ms;
    let started = 0;
    const hasDocument = typeof document !== 'undefined';
    const run = () => {
      if (timer !== undefined || (hasDocument && document.hidden)) return;
      started = Date.now();
      timer = setTimeout(() => {
        cancel();
        onTimeout();
      }, remaining);
    };
    const pause = () => {
      if (timer === undefined) return;
      clearTimeout(timer);
      timer = undefined;
      remaining = Math.max(0, remaining - (Date.now() - started));
    };
    const onVisibility = () => (document.hidden ? pause() : run());
    const cancel = () => {
      pause();
      if (hasDocument) document.removeEventListener('visibilitychange', onVisibility);
    };
    if (hasDocument) document.addEventListener('visibilitychange', onVisibility);
    run();
    return { cancel };
  }

  /** The pixel scale of the drawing buffer for `p`, see `resolutionMode` */
  export function getPixelScale(p: Props) {
    const scaled = p.pixelScale / (typeof window !== 'undefined' ? window?.devicePixelRatio || 1 : 1);
    if (p.resolutionMode === 'auto') {
      return isMobileBrowser() ? scaled : p.pixelScale;
    }
    return p.resolutionMode === 'native' ? p.pixelScale : scaled;
  }

  function getContextAttribs(a: Attribs, p: Props): GraphicsContextAttribs {
    return {
      powerPreference: a.powerPreference ?? 'default',
      failIfMajorPerformanceCaveat: a.failIfMajorPerformanceCaveat,
      antialias: a.antialias,
      preserveDrawingBuffer: a.preserveDrawingBuffer,
      pixelScale: getPixelScale(p),
    };
  }

  /** A canvas element to take the place of `canvas`, with its attributes (id, class, style) and size */
  function cloneCanvas(canvas: HTMLCanvasElement) {
    const next = canvas.ownerDocument.createElement('canvas');
    for (const { name, value } of Array.from(canvas.attributes)) next.setAttribute(name, value);
    next.width = canvas.width;
    next.height = canvas.height;
    return next;
  }

  /**
   * Creates the graphics context with the first available of `backends`, e.g.
   * `[WebGPUBackend, WebGL2Backend, WebGLBackend]`, see `createGraphicsContext`
   */
  export async function create(
    canvas: HTMLCanvasElement,
    assetManager: AssetManager,
    backends: readonly GraphicsBackendDescriptor[],
    attribs: Partial<Attribs> = {},
    props: Partial<Props> = {},
  ): Promise<Canvas3DContext> {
    const a = { ...DefaultAttribs, ...attribs };
    const p = { ...DefaultProps, ...props };
    const created = await createGraphicsContext(canvas, backends, getContextAttribs(a, p));
    if (!created) {
      throw new Error(`Could not create a graphics context with any of: ${backends.map((b) => b.id).join(', ')}`);
    }
    return fromGraphics(canvas, created.context, assetManager, a, p, backends);
  }

  /**
   * The context for `graphics`, a context created on `canvas` (e.g. by a backend descriptor, see `create`, or the WebGL
   * convenience `createWebGLCanvas3DContext` of `@molstar/graphics/canvas3d/webgl`); `backends`, the list it was
   * created from, for falling back when it is lost for good (see `Attribs.contextFallback`)
   */
  export function fromGraphics(
    _canvas: HTMLCanvasElement,
    _graphics: GraphicsContext,
    assetManager: AssetManager,
    attribs: Partial<Attribs> = {},
    props: Partial<Props> = {},
    backends: readonly GraphicsBackendDescriptor[] = [],
  ): Canvas3DContext {
    const a = { ...DefaultAttribs, ...attribs };
    const p = { ...DefaultProps, ...props };
    let canvas = _canvas;
    let graphics = _graphics;
    const syncPixelScale = () => {
      const pixelScale = getPixelScale(p);
      input.setPixelScale(pixelScale);
      graphics.setPixelScale(pixelScale);
    };

    const pixelScale = getPixelScale(p);
    const input = InputObserver.fromElement(canvas, { pixelScale, preventGestures: true });
    const createPasses = () =>
      new Passes(graphics, assetManager, { pickScale: p.pickScale, transparency: p.transparency });
    let passes = createPasses();

    // debugging: hold down shift+ctrl+alt and press any mouse button to lose the context, which is then restored
    const simulateContextLoss = (e: MouseEvent) => {
      if (!e.shiftKey || !e.ctrlKey || !e.altKey) return;
      if (graphics.simulateContextLoss()) console.log('lose context');
    };
    if (isDebugMode) canvas.addEventListener('mousedown', simulateContextLoss, false);

    // the losses the context reports (on WebGL from the events of the canvas); a lost context is replaced by a new one
    // on the same canvas once that can be created, and everything is built again on it as after a backend switch
    // (`contextChanged` with `restore`)
    const contextLost = new Subject<now.Timestamp>();
    const contextRestored = new Subject<now.Timestamp>();
    // fallback: the backend in use, those lost for good or that couldn't be created (never tried again), the times
    // of recent losses
    let current = backends.find((b) => b.id === _graphics.capabilities.backend);
    const failed = new Set<GraphicsBackendDescriptor>();
    let losses: number[] = [];
    let pendingRestore: { cancel: () => void } | undefined;
    /** lost on purpose (`release`): restored by `restore`, not on the loss */
    let released = false;
    /** backends set while released (`setBackend`), for `restore` */
    let backendsOnRestore: readonly GraphicsBackendDescriptor[] | undefined;

    const fallback = (why: string) => {
      if (!a.contextFallback || changing) return;
      if (current) failed.add(current);
      // only further down the list, each backend once: no loop back to a lost one
      const start = current ? backends.indexOf(current) + 1 : 0;
      const candidates = backends.slice(start).filter((b) => !failed.has(b));
      if (candidates.length === 0) {
        console.error(`Graphics context ${why}, no backend left to fall back to`);
        return;
      }
      console.warn(`Graphics context ${why}, falling back to ${candidates.map((b) => b.id).join(', ')}`);
      changeBackend(candidates, 'fallback').then((ok) => {
        if (!ok) console.error(`Graphics context ${why}, falling back failed`);
      });
    };

    const onContextLost = (t: now.Timestamp) => {
      if (isDebugMode) console.log(released ? 'context released' : 'context lost');
      contextLost.next(t);
      pendingRestore?.cancel();
      if (released) return;
      const time = Date.now();
      losses = losses.filter((l) => time - l < ContextLossWindowMs);
      losses.push(time);
      if (losses.length > MaxContextLosses) {
        fallback(`lost ${losses.length} times within ${ContextLossWindowMs / 1000} s`);
        return;
      }
      recreateLost();
    };

    /**
     * Replaces the lost context with a new one on the same canvas, falling back when that fails or takes too long; a
     * released context (`restore`) stays released instead, to be restored later
     */
    const recreateLost = (wasReleased = false): Promise<boolean> => {
      let settled = false;
      let done!: (restored: boolean) => void;
      const result = new Promise<boolean>((resolve) => {
        done = resolve;
      });
      const failed = (why: string) => {
        if (wasReleased) released = true;
        else fallback(why);
        done(false);
      };
      const timeoutMs = wasReleased ? ReleasedRestoreTimeoutMs : ContextRestoreTimeoutMs;
      const timeout = visibleTimeout(timeoutMs, () => {
        if (settled) return;
        settled = true;
        failed(`not restored within ${timeoutMs / 1000} s`);
      });
      pendingRestore = {
        cancel: () => {
          // also by `replaceContext` once restored (settled then)
          if (!settled) done(false);
          settled = true;
          timeout.cancel();
        },
      };
      const lost = graphics;
      graphics.recreate().then((next) => {
        if (settled || lost !== graphics || !next) {
          next?.destroy({ doNotForceWebGLContextLoss: true });
          if (!settled && lost === graphics) {
            settled = true;
            timeout.cancel();
            failed('not restored');
          }
          done(false);
          return;
        }
        settled = true;
        timeout.cancel();
        const restored = replaceContext({ canvas, graphics: next, backend: current }, 'restore');
        if (isDebugMode) console.log(restored ? 'context restored' : 'context not restored');
        if (restored) contextRestored.next(now());
        else failed('not restored');
        done(restored);
      });
      return result;
    };

    const subscribeContext = () => [graphics.contextLost.subscribe(onContextLost)];
    let contextSubs = subscribeContext();

    const changed = new BehaviorSubject<undefined>(undefined);
    const contextChanged = new Subject<ContextChange>();

    let changing = false;
    type NextContext = { canvas: HTMLCanvasElement; graphics: GraphicsContext; backend?: GraphicsBackendDescriptor };

    /**
     * Puts `next` in place of the current context and notifies `contextChanged`, then destroys the previous one (a
     * restore keeps the canvas: the lost context is dropped without GPU calls)
     */
    const replaceContext = (next: NextContext, reason: ContextChange['reason']): boolean => {
      // a new canvas takes the place of the current one in the document; a restore keeps the canvas (no parent needed)
      if ((next.canvas !== canvas && !canvas.parentNode) || graphics.xr.isPresenting()) {
        next.graphics.destroy(reason === 'restore' ? { doNotForceWebGLContextLoss: true } : undefined);
        return false;
      }
      const previous = graphics;
      for (const s of contextSubs) s.unsubscribe();
      canvas.removeEventListener('mousedown', simulateContextLoss, false);
      if (next.canvas !== canvas) canvas.replaceWith(next.canvas);

      canvas = next.canvas;
      graphics = next.graphics;
      input.setElement(canvas);
      if (isDebugMode) canvas.addEventListener('mousedown', simulateContextLoss, false);
      contextSubs = subscribeContext();
      passes = createPasses();
      syncPixelScale();
      pendingRestore?.cancel();
      pendingRestore = undefined;
      released = false;
      // a new backend starts without losses, a restore keeps counting them
      if (reason !== 'restore') losses = [];
      if (next.backend) {
        const { backend } = next;
        current = backends.includes(backend) ? backend : backends.find((b) => b.id === backend.id);
      }

      contextChanged.next({ previous, reason });
      // also a released context replaced by `restore` with other backends
      previous.destroy(reason === 'restore' || previous.isContextLost ? { lost: true } : undefined);
      return true;
    };

    /** set by `dispose`: a backend change finishing afterwards frees what it made and changes nothing */
    let disposed = false;

    const changeBackend = async (
      candidates: readonly GraphicsBackendDescriptor[],
      reason: ContextChange['reason'],
    ): Promise<boolean> => {
      if (disposed || changing || !canvas.parentNode || graphics.xr.isPresenting()) return false;
      changing = true;
      try {
        // each backend on a fresh canvas: a failed backend may leave its context type on the canvas
        let next: NextContext | undefined;
        for (const backend of candidates) {
          const nextCanvas = cloneCanvas(canvas);
          const created = await createGraphicsContext(nextCanvas, [backend], getContextAttribs(a, p));
          if (disposed) {
            created?.context.destroy();
            return false;
          }
          if (created) {
            next = { canvas: nextCanvas, graphics: created.context, backend };
            break;
          }
          // not tried again by a fallback
          if (reason === 'fallback') failed.add(backend);
        }
        return next ? replaceContext(next, reason) : false;
      } finally {
        changing = false;
      }
    };
    const setBackend = async (candidates: readonly GraphicsBackendDescriptor[]) => {
      // a released canvas (e.g. out of view) gets no live context before `restore`
      if (released && !changing) {
        backendsOnRestore = candidates;
        return true;
      }
      return changeBackend(candidates, 'switch');
    };

    return {
      get canvas() {
        return canvas;
      },
      get graphics() {
        return graphics;
      },
      input,
      get passes() {
        return passes;
      },
      attribs: a,
      get props() {
        return { ...p };
      },
      contextLost,
      contextRestored,
      contextChanged,
      assetManager,
      changed,
      get pixelScale() {
        return getPixelScale(p);
      },

      syncPixelScale,
      setProps: (props?: Partial<Props>) => {
        if (!props) return;

        let hasChanged = false;
        let pixelScaleNeedsUpdate = false;

        if (props.resolutionMode !== undefined && props.resolutionMode !== p.resolutionMode) {
          p.resolutionMode = props.resolutionMode;
          pixelScaleNeedsUpdate = true;
        }

        if (props.pixelScale !== undefined && props.pixelScale !== p.pixelScale) {
          p.pixelScale = props.pixelScale;
          pixelScaleNeedsUpdate = true;
        }

        if (pixelScaleNeedsUpdate) {
          syncPixelScale();
          a.handleResize();
          hasChanged = true;
        }

        if (props.pickScale !== undefined && props.pickScale !== p.pickScale) {
          p.pickScale = props.pickScale;
          passes.setPickScale(props.pickScale);
          hasChanged = true;
        }

        if (props.transparency !== undefined && props.transparency !== p.transparency) {
          p.transparency = props.transparency;
          passes.setTransparency(props.transparency);
          hasChanged = true;
        }

        if (hasChanged) changed.next(undefined);
      },
      setBackend,
      release: () => {
        if (released || changing || graphics.isContextLost || graphics.xr.isPresenting()) return false;
        released = true;
        if (graphics.release()) return true;
        released = false;
        return false;
      },
      restore: async () => {
        if (!released || changing || !graphics.isContextLost) return false;
        const backends = backendsOnRestore;
        backendsOnRestore = undefined;
        if (backends) {
          // the backends set while released; the current one restored if none of them can be created
          if (await changeBackend(backends, 'switch')) {
            contextRestored.next(now());
            return true;
          }
          if (!released || changing || !graphics.isContextLost) return false;
        }
        released = false;
        return recreateLost(true);
      },
      get isReleased() {
        return released;
      },
      dispose: (options?: GraphicsDestroyOptions) => {
        disposed = true;
        input.dispose();

        canvas.removeEventListener('mousedown', simulateContextLoss, false);
        for (const s of contextSubs) s.unsubscribe();
        pendingRestore?.cancel();
        graphics.destroy(options);

        contextLost.complete();
        contextRestored.complete();
        contextChanged.complete();
        changed.complete();
      },
    };
  }
}

export { Canvas3D };

export interface Canvas3DCameraResetOptions {
  durationMs?: number;
  snapshot?: Camera.SnapshotProvider;
  keyframes?: CameraTransitionManager.TransitionKeyframes;
  easing?: EasingFunction;
  trajectory?: TransitionTrajectory;
}

interface Canvas3D {
  /** replaced when the context's backend changes, see `contextChanged` */
  readonly graphics: GraphicsContext;
  /**
   * After the graphics context changed (`Canvas3DContext.setBackend`) and the scene was rebuilt on it: render objects
   * are drawn again as they are (textures loaded from data transferred), representations create and update visuals
   * with the new context (`Representation.setGraphics`) and are updated, which creates again what was rendered or
   * computed with the previous context (e.g. GPU surfaces) and geometry relying on what the new context lacks; image
   * passes (`getImagePass`) follow. Render objects of no representation that hold such textures are not drawn anymore
   * (`dropped`).
   */
  readonly contextChanged: Subject<Canvas3D.ContextChange>;
  /** What the graphics context reports as not getting the GPU memory it needed (`GraphicsResources.outOfMemory`) */
  readonly outOfMemory: Subject<string>;

  add(repr: Representation.Any): void;
  remove(repr: Representation.Any): void;
  /**
   * This function must be called if animate() is not set up so that add/remove actions take place.
   */
  commit(isSynchronous?: boolean): void;
  /**
   * Function for external "animation" control
   * Calls commit.
   */
  tick(
    t: now.Timestamp,
    options?: { isSynchronous?: boolean; manualDraw?: boolean; updateControls?: boolean; xrFrame?: XRFrame },
  ): void;
  update(repr?: Representation.Any, keepBoundingSphere?: boolean): void;
  clear(): void;
  syncVisibility(): void;

  requestDraw(): void;
  /**
   * Makes the pipelines and resources frames with the props draw ready without drawing anything (see `PassBuilder`);
   * resolves when they are
   */
  prepare(options?: Canvas3DPrepareOptions): Promise<void>;
  /**
   * Makes the programs and pipelines ready that frames with each of `props` (the current ones by default) draw the
   * objects of `list` with, before the objects exist (e.g. the snapshots of a state, before they are applied), including
   * the first interaction; as this device draws them (e.g. blended transparency where the props' mode is not supported).
   * The programs are kept until the next call; entries of unknown kinds are skipped. Resolves when they are ready
   */
  precompile(list: PrecompileList, props?: readonly Canvas3DPrecompileProps[]): Promise<void>;
  /** The objects drawn since the last reset (`resetPrecompileList`), as signatures, see `precompile` */
  getPrecompileList(): PrecompileList;
  resetPrecompileList(): void;
  /**
   * Whether the last tick found nothing to do: nothing drawn or committed, and nothing under way that needs the next
   * tick (a camera transition, the controls' animation or keys, an interaction, a hover pick, a hi-Z readback). A loop
   * may stop asking for frames while it stays idle (`PluginAnimationLoop`, `PluginConfig.General.PauseLoopWhenIdle`)
   * and ask again on `wake`. Code that changes the scene or the camera past this API must call `requestDraw`.
   */
  readonly isIdle: boolean;
  /** Fires once after an idle tick when something gives the next tick work: a draw requested, input, a scene, camera or props change. */
  readonly wake: Observable<void>;

  /** Reset the timers, used by "animate" */
  resetTime(t: number): void;
  animate(): void;
  /**
   * Pause animation loop and optionally any rendering
   * @param noDraw pause any rendering (drawPaused = true)
   */
  pause(noDraw?: boolean): void;
  /** Sets drawPaused = false without starting the built in animation loop */
  resume(): void;
  /**
   * Frames wait until `until` settles (the last frame stays shown), for at most `maxMs` (1000), as they wait for new
   * props to be prepared: a change made in steps (a state update rebuilding representations, then their own transforms
   * put back) drawn once it is whole. The animation loop and any pause are left as they are.
   */
  holdFrames(until: Promise<unknown>, maxMs?: number): void;

  requestAnimationFrame(callback: FrameRequestCallback | XRFrameRequestCallback): number;
  cancelAnimationFrame(handle: number): void;

  /**
   * The pick data at `target`, `undefined` if nothing is there. Read in the same call where the backend reads
   * synchronously (WebGL); elsewhere (WebGPU) `'pending'` until a read started for the current scene and camera has
   * returned: `asyncIdentify(target)?.wait()` resolves then
   */
  identify(target: Vec2 | Ray3D): PickResult;
  /**
   * Starts reading the pick data at `target`, on every backend: poll the result with `tryGet` (e.g. per frame, as hover
   * and clicks do) or await it with `wait`
   */
  asyncIdentify(target: Vec2 | Ray3D): AsyncPickData | undefined;
  mark(loci: Representation.Loci, action: MarkerAction): void;
  getLoci(pickingId: PickingId | undefined): Representation.Loci;

  notifyDidDraw: boolean;
  readonly didDraw: BehaviorSubject<now.Timestamp>;
  readonly commited: BehaviorSubject<now.Timestamp>;
  readonly commitQueueSize: BehaviorSubject<number>;
  readonly reprCount: BehaviorSubject<number>;
  readonly resized: BehaviorSubject<any>;

  handleResize(): void;
  /** performs handleResize on the next animation frame */
  requestResize(): void;
  /** Focuses camera on scene's bounding sphere, centered and zoomed. */
  requestCameraReset(options?: Canvas3DCameraResetOptions): void;
  readonly camera: Camera;
  readonly boundingSphere: Readonly<Sphere3D>;
  readonly boundingSphereVisible: Readonly<Sphere3D>;
  setProps(
    props: PartialCanvas3DProps | ((old: Canvas3DProps) => Partial<Canvas3DProps> | void),
    doNotRequestDraw?: boolean /* = false */,
  ): void;
  setAttribs(attribs: PartialCanvas3DAttribs): void;
  getImagePass(props: Partial<ImageProps>): ImagePass;
  getRenderObjects(): GraphicsRenderObject[];

  /** Returns a copy of the current Canvas3D instance props */
  readonly props: Readonly<Canvas3DProps>;
  readonly attribs: Readonly<Canvas3DAttribs>;
  readonly input: InputObserver;
  readonly stats: RendererStats;
  readonly interaction: Canvas3dInteractionHelper['events'];

  readonly debugRegistry: DebugRegistry;

  readonly xr: {
    request(): Promise<void>;
    end(): Promise<void>;
    readonly isSupported: BehaviorSubject<boolean>;
    readonly isPresenting: BehaviorSubject<boolean>;
    readonly requestFailed: Subject<string>;
  };

  dispose(): void;
}

const requestAnimationFrame =
  typeof window !== 'undefined'
    ? window.requestAnimationFrame
    : (f: (time: number) => void) => setImmediate(() => f(Date.now())) as unknown as number;
const cancelAnimationFrame =
  typeof window !== 'undefined'
    ? window.cancelAnimationFrame
    : (handle: number) => clearImmediate(handle as unknown as NodeJS.Immediate);

/** Longest time frames wait for the work a props change started, see `holdFramesUntilPrepared` */
const MaxFrameHoldMs = 1000;

function syncCanvasBackground(canvas: HTMLCanvasElement | undefined, canvasProps: Canvas3DProps) {
  // A HeadlessPluginContext has no HTML canvas to style.
  if (!canvas) return;
  if (canvasProps.transparentBackground && canvasProps.checkeredTransparentBackground) {
    Object.assign(canvas.style, {
      'background-image':
        'linear-gradient(45deg, lightgrey 25%, transparent 25%, transparent 75%, lightgrey 75%, lightgrey), linear-gradient(45deg, lightgrey 25%, transparent 25%, transparent 75%, lightgrey 75%, lightgrey)',
      'background-size': '60px 60px',
      'background-position': '0 0, 30px 30px',
    });
  } else {
    Object.assign(canvas.style, {
      'background-image': '',
      'background-size': '',
      'background-position': '',
    });
  }
}

namespace Canvas3D {
  export interface ContextChange {
    readonly dropped: readonly GraphicsRenderObject[];
    /** `switch` (`Canvas3DContext.setBackend`) or `fallback` (the previous context was lost for good) */
    readonly reason: Canvas3DContext.ContextChange['reason'];
  }

  export interface HoverEvent {
    current: Representation.Loci;
    buttons: ButtonsType;
    button: ButtonsType.Flag;
    modifiers: ModifiersKeys;
    page?: Vec2;
    position?: Vec3;
  }
  export interface DragEvent {
    current: Representation.Loci;
    buttons: ButtonsType;
    button: ButtonsType.Flag;
    modifiers: ModifiersKeys;
    pageStart: Vec2;
    pageEnd: Vec2;
  }
  export interface ClickEvent {
    current: Representation.Loci;
    buttons: ButtonsType;
    button: ButtonsType.Flag;
    modifiers: ModifiersKeys;
    page?: Vec2;
    position?: Vec3;
  }

  export function create(
    ctx: Canvas3DContext,
    props: Partial<Canvas3DProps> = {},
    attribs: Partial<Canvas3DAttribs> = {},
  ): Canvas3D {
    const { input, assetManager, contextLost } = ctx;
    let { graphics, passes, canvas } = ctx;
    const p: Canvas3DProps = { ...deepClone(DefaultCanvas3DParams), ...deepClone(props) };
    let limitToSupported = createSupportedPropsLimiter(graphics);
    limitToSupported(p);
    const a = { ...deepClone(DefaultCanvas3DAttribs), ...deepClone(attribs) };

    const reprRenderObjects = new Map<Representation.Any, Set<GraphicsRenderObject>>();
    const reprUpdatedSubscriptions = new Map<Representation.Any, Subscription>();
    /** the representations a loci marked on all of them may change (`applyMark`) */
    const markIndex = new MarkIndex();
    const reprCount = new BehaviorSubject(0);
    const interactionEvent = new Subject<void>();

    let startTime = now();
    const didDraw = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
    const commited = new BehaviorSubject<now.Timestamp>(0 as now.Timestamp);
    const commitQueueSize = new BehaviorSubject<number>(0);

    let x = 0;
    let y = 0;
    let width = 128;
    let height = 128;

    let forceNextRender = false;
    let currentTime = 0;

    // The idle state (`isIdle`): set by a tick that found nothing to do, cleared, with one `wake`, by whatever gives the
    // next tick work.
    let idle = false;
    let drewInTick = false;
    const wake = new Subject<void>();
    function busy() {
      if (!idle) return;
      idle = false;
      wake.next();
    }
    // bumped once per render() call; spans stereo eyes and multi-sample jitter
    // sub-renders of that call, so their cull results can be safely reused
    let frame: Frame = createFrame();

    syncCanvasBackground(canvas, p);
    updateViewport();
    const scene = Scene.create(graphics, passes.draw.transparency, {
      dColorMarker: isMaterialColorMarker(p.renderer, p.marking),
      // `props.renderer` may be partial (shallow merge above); the renderer fills in its defaults
      dLightCount: Math.min(p.renderer.light?.length ?? RendererParams.light.defaultValue.length, MaxLightCount),
    });

    function getSceneRadius() {
      return scene.boundingSphere.radius * p.sceneRadiusFactor;
    }

    const camera = new Camera(
      {
        position: Vec3.create(0, 0, 100),
        mode: p.camera.mode,
        fog: p.cameraFog.name === 'on' ? p.cameraFog.params.intensity : 0,
        clipFar: p.cameraClipping.far,
        minNear: p.cameraClipping.minNear,
        fov: degToRad(p.camera.fov),
      },
      { x, y, width, height },
    );
    const stereoCamera = new StereoCamera(camera, p.camera.stereo.params);

    const controls = TrackballControls.create(input, camera, scene, p.trackball, a.trackball);
    const helper = new Helper(graphics, scene, p);
    let hiZ = new HiZPass(graphics, passes.draw, canvas, p.hiZ);
    // the pass's props are the props (`setProps` updates them), not a copy that goes stale
    p.hiZ = hiZ.props;

    let renderer = Renderer.create(graphics, p.renderer);
    /** the depth pyramid of the previous frame, for GPU culling (`Renderer.setOcclusionPyramid`) */
    const getOcclusionPyramid = () => hiZ.gpuPyramid;
    renderer.setOcclusionTest(hiZ.occlusionTest);
    renderer.setOcclusionPyramid(getOcclusionPyramid);

    let shaderManager = new ShaderManager(graphics, scene);
    shaderManager.updateRequired(p);

    const pickOptions = {
      pickPadding: p.pickPadding,
      maxAsyncReadLag: DefaultPickOptions.maxAsyncReadLag,
    };
    let pickHelper = new PickHelper(
      graphics,
      renderer,
      scene,
      helper,
      passes.pick,
      { x, y, width, height },
      pickOptions,
    );
    let rayHelper = new RayHelper(graphics, renderer, scene, helper, pickOptions);
    const interactionHelper = new Canvas3dInteractionHelper(
      identify,
      asyncIdentify,
      getLoci,
      input,
      camera,
      controls,
      p.interaction,
      pickHelper.synchronous,
    );
    let multiSampleHelper = new MultiSampleHelper(passes.multiSample);
    // the viewport of the renderer, picking and Hi-Z from the start, not only after a resize (a headless canvas may
    // never have one)
    syncViewport();

    passes.draw.postprocessing.background.update(camera, p.postprocessing.background, (changed) => {
      if (changed) requestDraw();
    });

    let cameraResetRequested = false;
    const nextCameraResetOptions: Canvas3DCameraResetOptions = {
      durationMs: undefined,
      snapshot: undefined,
      keyframes: undefined,
      easing: undefined,
      trajectory: undefined,
    };
    let resizeRequested = false;

    //

    function getNonXRProps() {
      return {
        transparency: ctx.props.transparency,
        transparentBackground: p.transparentBackground,
        hiZ: hiZ.props.enabled,
        postprocessing: p.postprocessing.enabled,
        axes: deepClone(helper.camera.props.axes),
      };
    }
    const nonXRProps = getNonXRProps();

    function saveNonXRProps() {
      Object.assign(nonXRProps, getNonXRProps());
    }

    function loadNonXRProps() {
      p.postprocessing.enabled = nonXRProps.postprocessing;
      p.transparentBackground = nonXRProps.transparentBackground;
      ctx.setProps({ transparency: nonXRProps.transparency });
      hiZ.setProps({ enabled: nonXRProps.hiZ });
      helper.camera.setProps({ axes: nonXRProps.axes });
    }

    function setXRProps() {
      p.postprocessing.enabled = !xrManager.props.disablePostprocessing;
      ctx.setProps({ transparency: 'blended' });
      hiZ.setProps({ enabled: false });
      helper.camera.setProps({ axes: { name: 'off', params: {} } });

      if (xrManager.session?.environmentBlendMode === 'alpha-blend') {
        p.transparentBackground = xrPassthrough;
      }
    }

    const xrManager = new XRManager(
      graphics,
      input,
      scene,
      camera,
      stereoCamera,
      helper.pointer,
      interactionHelper,
      p.xr,
      a.xr,
    );

    const xr = {
      request: async () => {
        try {
          await xrManager.request();
        } catch (e) {
          console.error(e);
          xr.requestFailed.next(e);
        }
      },
      end: () => xrManager.end(),
      isSupported: new BehaviorSubject(false),
      isPresenting: new BehaviorSubject(false),
      requestFailed: new Subject<string>(),
    };

    let xrPassthrough = false;

    const xrSubs = [
      xrManager.isSupported.subscribe((e) => xr.isSupported.next(e)),
      xrManager.togglePassthrough.subscribe(() => {
        if (xrManager.session?.environmentBlendMode === 'alpha-blend') {
          xrPassthrough = !p.transparentBackground;
        }
      }),
      xrManager.sessionChanged.subscribe(() => {
        resizeRequested = true;
        if (xrManager.session) {
          saveNonXRProps();
          xrPassthrough = xrManager.session?.environmentBlendMode === 'alpha-blend';
          setXRProps();
        } else {
          loadNonXRProps();
        }
        resume();
        busy();
        xr.isPresenting.next(!!xrManager.session);
      }),
    ];

    //

    let notifyDidDraw = true;

    function getLoci(pickingId: PickingId | undefined) {
      let loci: Loci = EmptyLoci;
      let repr: Representation.Any = Representation.Empty;
      if (pickingId) {
        const cameraHelperLoci = helper.camera.getLoci(pickingId);
        if (cameraHelperLoci !== EmptyLoci) return { loci: cameraHelperLoci, repr };

        loci = helper.handle.getLoci(pickingId);

        reprRenderObjects.forEach((_, _repr) => {
          const _loci = _repr.getLoci(pickingId);
          if (!isEmptyLoci(_loci)) {
            if (!isEmptyLoci(loci)) {
              console.warn('found another loci, this should not happen');
            }
            loci = _loci;
            repr = _repr;
          }
        });
      }
      return { loci, repr };
    }

    let markBuffer: [reprLoci: Representation.Loci, action: MarkerAction][] = [];

    function mark(reprLoci: Representation.Loci, action: MarkerAction) {
      // NOTE: might try to optimize a case with opposite actions for the
      //       same loci. Tho this might end up being more expensive (and error prone)
      //       then just applying everything "naively".
      markBuffer.push([reprLoci, action]);
      busy();
    }

    /** helpers are rendered into the scene color, so their marking can't be redrawn on its own */
    let helperMarkingUpdated = false;

    function resolveMarking() {
      let changed = false;
      helperMarkingUpdated = false;
      for (const [r, l] of markBuffer) {
        changed = applyMark(r, l) || changed;
      }
      markBuffer = [];
      if (changed) {
        // the marked objects' values changed: those updated, not every object of the scene
        scene.updateChanged(true);
        helper.handle.scene.updateChanged(true);
        helper.camera.scene.updateChanged(true);

        shaderManager.updateRequired(p);
        shaderManager.finalizeRequired(true);

        interactionEvent.next();
      }
      return changed;
    }

    function applyMark(reprLoci: Representation.Loci, action: MarkerAction) {
      const { repr, loci } = reprLoci;
      let changed = false;
      if (repr) {
        if (repr.mark(loci, action)) {
          markIndex.note(repr, action);
          changed = true;
        }
      } else {
        markIndex.forEach(loci, action, (_repr) => {
          if (!_repr.mark(loci, action)) return;
          markIndex.note(_repr, action);
          changed = true;
        });
      }
      const handleChanged = helper.handle.mark(loci, action);
      const cameraChanged = helper.camera.mark(loci, action);
      const helperChanged = handleChanged || cameraChanged;
      helperMarkingUpdated = helperMarkingUpdated || helperChanged;
      return changed || helperChanged;
    }

    function render(force: boolean, xrFrame?: XRFrame) {
      if (graphics.isContextLost) return false;
      if (graphics.xr.session && !xrFrame) return false;

      let resized = false;
      if (resizeRequested) {
        handleResize(false);
        resizeRequested = false;
        resized = true;
      }

      const drs = graphics.getDrawingBufferSize();
      if (x > drs.width || x + width < 0 || y > drs.height || y + height < 0) return false;

      if (xrFrame) {
        setXRProps();
        p.transparentBackground = xrPassthrough;
      }

      const markingUpdated = resolveMarking() && (renderer.props.colorMarker || p.marking.enabled);

      let didRender = false;
      controls.update(currentTime);
      const cameraChanged = camera.update();
      const xrChanged = xrManager.update(xrFrame);
      if (!xrChanged && xrFrame) return false;

      const activeAnimation = renderer.props.enableAnimation && scene.hasAnimation;
      const shouldRender = force || cameraChanged || resized || forceNextRender || xrChanged || activeAnimation;
      forceNextRender = false;

      // only a new frameId when something that can affect culling actually changed, so idle
      // ticks (e.g. temporal multi-sample accumulation on a resting camera) reuse the last cull
      if (shouldRender) {
        frame = createFrame();
        if (isDebugMode) console.log('New frame');
      }

      const illuminationEnabled = passes.illumination.supported && p.illumination.enabled && !xrFrame;
      const multiSampleEnabled = MultiSamplePass.isEnabled(p.multiSample) && !xrFrame;

      const markingOffsets = (illuminationEnabled ? p.multiSample.mode === 'on' : multiSampleEnabled)
        ? getJitterOffsets(p.multiSample.sampleLevel)
        : SingleSample;
      const redrawMarking = (restart: boolean, samples: number) => {
        // other passes (e.g. pick) reset this, it has to match the last full render
        renderer.setOcclusionTest(illuminationEnabled ? null : hiZ.occlusionTest);
        renderer.setOcclusionPyramid(illuminationEnabled ? null : getOcclusionPyramid);
        if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
        const redrawn = passes.draw.marking.redraw(
          { renderer, camera, scene, helper, frame },
          p,
          markingOffsets,
          restart,
          samples,
        );
        if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');
        return redrawn;
      };

      // marking-only changes are redrawn over the image of the last render, avoiding a full
      // re-render; the per-object color marker is part of that image, so it needs a re-render
      const canRedrawMarking =
        !helperMarkingUpdated &&
        !shouldRender &&
        !isMaterialColorMarker(renderer.props, p.marking) &&
        !xrFrame &&
        p.camera.stereo.name !== 'on' &&
        (!illuminationEnabled || passes.illumination.iteration > 0);

      if (canRedrawMarking && markingUpdated) {
        // illumination multi-samples are expensive (full depth pass each), so add them over frames
        const progressive = illuminationEnabled || p.multiSample.mode === 'temporal';
        if (redrawMarking(true, progressive ? 1 : markingOffsets.length)) return true;
      }

      if (illuminationEnabled) {
        if (shouldRender || markingUpdated) {
          renderer.setOcclusionTest(null);
          renderer.setOcclusionPyramid(null);
          passes.illumination.restart();
        }

        if (
          passes.illumination.shouldRender(p.illumination) &&
          ((!isActivelyInteracting && scene.count > 0) ||
            passes.illumination.iteration === 0 ||
            p.userInteractionReleaseMs === 0)
        ) {
          if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
          const ctx = { renderer, camera, scene, helper, frame };
          passes.illumination.render(ctx, p, true);
          if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');

          // if only marking has updated, do not set the flag to dirty
          pickHelper.dirty = pickHelper.dirty || shouldRender;
          didRender = true;
        }
      } else {
        const multiSampleChanged = multiSampleHelper.update(markingUpdated || shouldRender, p.multiSample) && !xrFrame;

        if (shouldRender || multiSampleChanged || markingUpdated) {
          renderer.setOcclusionTest(hiZ.occlusionTest);
          renderer.setOcclusionPyramid(getOcclusionPyramid);

          let cam: Camera | StereoCamera = camera;
          if (p.camera.stereo.name === 'on' || xrChanged) {
            if (!xrChanged) stereoCamera.update();
            cam = stereoCamera;
          }

          if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
          const ctx = { renderer, camera: cam, scene, helper, frame };
          if (multiSampleEnabled) {
            const forceOn =
              p.multiSample.reduceFlicker &&
              isMaterialColorMarker(renderer.props, p.marking) &&
              !cameraChanged &&
              markingUpdated &&
              !controls.isAnimating;
            multiSampleHelper.render(ctx, p, true, forceOn);
          } else {
            passes.draw.render(ctx, p, true);
          }
          // the GPU's occlusion test of the objects culled there needs the levels (also with Hi-Z off)
          hiZ.setGpuCulling(renderer.gpuCulledCount > 0);
          hiZ.render(camera);
          if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');

          // if only marking has updated, do not set the flag to dirty
          pickHelper.dirty = pickHelper.dirty || shouldRender;
          didRender = true;
        }
      }

      // keep accumulating marking samples after the scene has settled
      if (!didRender && !markingUpdated && canRedrawMarking && passes.draw.marking.needsMoreSamples) {
        didRender = redrawMarking(false, getTemporalSamplesPerFrame(p.multiSample.sampleLevel));
      }

      return didRender;
    }

    let forceDrawAfterAllCommited = false;
    let drawPaused = false;
    let isContextLost = false;
    /** set by `dispose`: idle callbacks queued before (warm-up, target release) do nothing */
    let canvasDisposed = false;

    /** prepares of new props, and holds (`holdFrames`), that frames wait for, see `holdFramesUntilPrepared` */
    let propsPreparing = 0;
    /**
     * Frames with the current props wait until what they draw with is ready, when that is not the case yet (e.g.
     * WebGPU pipelines of illumination compiling, whose draws would be skipped): the last frame stays shown meanwhile,
     * instead of frames missing those draws (e.g. only the background); for at most `MaxFrameHoldMs`, so that work
     * that does not settle doesn't freeze the canvas
     */
    function holdFrames(until: Promise<unknown>, maxMs = MaxFrameHoldMs) {
      propsPreparing += 1;
      let held = true;
      const release = () => {
        if (!held) return;
        held = false;
        clearTimeout(limit);
        propsPreparing -= 1;
        requestDraw();
      };
      const limit = setTimeout(release, maxMs);
      until.then(release, release);
    }

    function holdFramesUntilPrepared() {
      if (isContextLost || graphics.isContextLost) return;
      const { pending } = graphics.resources;
      const trackedCount = pending.trackedCount;
      const ready = prepare();
      if (pending.trackedCount === trackedCount) return;
      propsPreparing += 1;
      let held = true;
      const release = () => {
        if (!held) return;
        held = false;
        clearTimeout(limit);
        propsPreparing -= 1;
        requestDraw();
      };
      // cleared once ready, so it only outlives the hold when the work never settles
      const limit = setTimeout(release, MaxFrameHoldMs);
      ready.then(release);
    }

    function draw(options?: { force?: boolean; isSynchronous?: boolean; xrFrame?: XRFrame }) {
      if (drawPaused || isContextLost) return;
      if (propsPreparing > 0 || !shaderManager.finalizeRequired(options?.isSynchronous)) {
        forceNextRender = true;
        return;
      }
      if (render(!!options?.force, options?.xrFrame)) {
        drewInTick = true;
        if (notifyDidDraw) didDraw.next((now() - startTime) as now.Timestamp);
      }
    }

    function requestDraw() {
      forceNextRender = true;
      busy();
    }

    let animationFrameHandle = 0;

    function tick(
      t: now.Timestamp,
      options?: { isSynchronous?: boolean; manualDraw?: boolean; updateControls?: boolean; xrFrame?: XRFrame },
    ) {
      if (isContextLost) return;
      if (graphics.xr.session && !options?.xrFrame) return;

      drewInTick = false;
      idle = false;
      currentTime = t;
      renderer.setTime((currentTime - startTime) / 1000);
      commit(options?.isSynchronous);

      // update the controler before the camera transition
      if (options?.updateControls) {
        controls.update(currentTime);
      }

      camera.transition.tick(currentTime);
      hiZ.tick();

      if (options?.manualDraw) {
        return;
      }

      draw({ isSynchronous: options?.isSynchronous, xrFrame: options?.xrFrame });
      if (!camera.transition.inTransition && !graphics.isContextLost) {
        interactionHelper.tick(currentTime);
      }
      idle =
        !drewInTick &&
        !forceNextRender &&
        !forceDrawAfterAllCommited &&
        !cameraResetRequested &&
        !resizeRequested &&
        markBuffer.length === 0 &&
        !scene.needsCommit &&
        !camera.transition.inTransition &&
        !controls.isAnimating &&
        !controls.isMoving &&
        !isActivelyInteracting &&
        !hiZ.isPending &&
        !interactionHelper.isPending &&
        !graphics.xr.session;
    }

    let animationFrameCB: FrameRequestCallback | XRFrameRequestCallback | undefined = undefined;

    function _requestAnimationFrame(callback: FrameRequestCallback | XRFrameRequestCallback): number {
      animationFrameCB = callback;
      return graphics.xr.session
        ? graphics.xr.session.requestAnimationFrame(callback)
        : requestAnimationFrame(callback as FrameRequestCallback);
    }

    function _cancelAnimationFrame(handle: number): void {
      animationFrameCB = undefined;
      graphics.xr.session ? graphics.xr.session.cancelAnimationFrame(handle) : cancelAnimationFrame(handle);
    }

    function _animate(_timestamp: number, xrFrame?: XRFrame) {
      tick(now(), { xrFrame });
      animationFrameHandle = _requestAnimationFrame(_animate);
    }

    function resetTime(t: now.Timestamp) {
      startTime = t;
      controls.start(t);
      interactionHelper.resetTime(t);
    }

    function animate() {
      drawPaused = false;
      const t = now();
      controls.start(t);
      interactionHelper.resetTime(t);
      if (animationFrameHandle === 0) _animate(0);
    }

    function pause(noDraw = false) {
      drawPaused = noDraw;
      if (animationFrameHandle !== 0) {
        _cancelAnimationFrame(animationFrameHandle);
        animationFrameHandle = 0;
      }
    }

    function resume() {
      drawPaused = false;
      if (animationFrameCB) _requestAnimationFrame(animationFrameCB);
    }

    function identify(target: Vec2 | Ray3D): PickResult {
      if (graphics.isContextLost) return undefined;
      shaderManager.finalize(['pick'], true);

      if ('origin' in target) {
        return rayHelper.identify(target, camera);
      } else {
        const cam = p.camera.stereo.name === 'on' ? stereoCamera : camera;
        return pickHelper.identify(target[0], target[1], cam, frame);
      }
    }

    function asyncIdentify(target: Vec2 | Ray3D): AsyncPickData | undefined {
      if (graphics.isContextLost) return undefined;
      shaderManager.finalize(['pick'], true);

      if ('origin' in target) {
        return rayHelper.asyncIdentify(target, camera);
      } else {
        const cam = p.camera.stereo.name === 'on' ? stereoCamera : camera;
        return pickHelper.asyncIdentify(target[0], target[1], cam, frame);
      }
    }

    /** Plans the passes of a frame with `props` (as `render` chooses them, without XR); with `interaction`, hovering */
    function planFrame(b: PassBuilder, props: Canvas3DProps, interaction: boolean) {
      const cam = props.camera.stereo.name === 'on' ? stereoCamera : camera;
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
      if (isContextLost || graphics.isContextLost) return Promise.resolve();
      const props: Canvas3DProps = options?.props ? { ...p, ...options.props } : p;
      // the variants the props require are linked for every render item when set (`ShaderManager.updateRequired`)
      graphics.resources.linkPrograms(ShaderManager.getRequired(graphics, scene, props));
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
      const { targets } = passes;
      targets.scheduleRelease(() => {
        if (canvasDisposed || isContextLost || graphics.isContextLost || targets !== passes.targets) return;
        const b = createTargetCollector();
        planFrame(b, p, true);
        targets.releaseUnused(b.used);
      });
    }

    const precompileRecorder = new PrecompileRecorder();
    /** the program items of the last `precompile`, by transparency, globals and signature */
    let precompileItems = new Map<string, ProgramItem>();
    function destroyPrecompileItems() {
      for (const item of precompileItems.values()) item.destroy();
      precompileItems.clear();
    }

    /** the scene's globals with `props` (as `setProps` sets them) */
    function getGlobals(props: Canvas3DProps): GlobalDefines {
      return {
        dColorMarker: isMaterialColorMarker(
          { colorMarker: props.renderer?.colorMarker ?? renderer.props.colorMarker },
          props.marking,
        ),
        dLightCount: Math.min(props.renderer?.light?.length ?? renderer.props.light.length, MaxLightCount),
      };
    }

    function precompile(list: PrecompileList, propsList: readonly Canvas3DPrecompileProps[] = [{}]): Promise<void> {
      if (isContextLost || graphics.isContextLost) return Promise.resolve();
      const previous = precompileItems;
      precompileItems = new Map();
      const gpuCulled = list.objects.some((s) => s.gpuCulled);
      // the passes plan against scenes of the list's objects: groups with their counts (the requests come from the
      // signatures, see `collect`) and averages that take every branch drawing objects (opaque and transparent or
      // only opaque, emissive or not); marking as for the first interaction
      const group = (n: number): Scene.Group => ({ ...Object3D.create(), renderables: new Array(n) });
      for (const { transparency: mode, ...options } of propsList) {
        const props: Canvas3DProps = { ...p, ...options };
        const transparency = passes.draw.resolveTransparency(mode ?? passes.draw.transparency);
        const globals = getGlobals(props);
        const primitives: { signature: ObjectSignature; item: ProgramItem }[] = [];
        const volumes: { signature: ObjectSignature; item: ProgramItem }[] = [];
        for (const signature of list.objects) {
          const key = `${transparency}|${globals.dLightCount}|${globals.dColorMarker}|${getSignatureKey(signature)}`;
          let item = precompileItems.get(key) ?? previous.get(key);
          if (!item) {
            try {
              item = createSignatureProgramItem(graphics, signature, transparency, globals);
            } catch (e) {
              if (isDebugMode) console.warn('Skipped a precompile list entry', signature, e);
            }
            if (!item) continue;
          }
          precompileItems.set(key, item);
          (signature.type === 'direct-volume' ? volumes : primitives).push({ signature, item });
        }

        const listPrimitives = group(primitives.length);
        const listVolumes = group(volumes.length);
        // the averages decide between exclusive branches (e.g. bloom composites the background over emissive objects
        // only, postprocessing blends transparency only with transparent objects): planned with none and some
        const listScenes: Scene[] = [];
        for (const opacityAverage of [1, 0.5]) {
          for (const emissiveAverage of [0, 0.5]) {
            listScenes.push(
              Object.create(scene, {
                primitives: { value: listPrimitives },
                volumes: { value: listVolumes },
                hasOpaque: { value: true },
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
        const cam = props.camera.stereo.name === 'on' ? stereoCamera : camera;
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
              hiZ.plan(b, camera, props.hiZ, gpuCulled);
            }
            if (MarkingPass.isEnabled(props.marking)) {
              passes.draw.plan(b, { ...ctx, scene: markedScenes[i] }, props, true, false, true);
            }
          }
          pickHelper.plan(b, cam, frame);
        });
      }
      for (const [key, item] of previous) {
        if (!precompileItems.has(key)) item.destroy();
      }
      // the culling's kernels of objects culled on the GPU, else compiled once such an object is created
      const kernels = graphics.compute && gpuCulled ? graphics.compute.prepareKernels(MergedCullKernels) : undefined;
      // the illumination's denoiser, a kernel where compute shaders run
      const denoise =
        passes.illumination.supported &&
        isComputeDenoiseSupported(graphics) &&
        propsList.some((o) => {
          const { enabled, denoise } = { ...p.illumination, ...o.illumination };
          return enabled && denoise;
        })
          ? graphics.compute.prepareKernels(getDenoiseKernels())
          : undefined;
      return Promise.all([graphics.resources.pending.whenReady(), kernels, denoise]).then(() => {});
    }

    /** the scheduled warm-up, cancelled by `dispose` (it would plan a frame on the destroyed context) */
    let warmUpCancel: (() => void) | undefined;
    /** prepares the frame and the first interaction in idle time, after the scene changed */
    function scheduleWarmUp() {
      if (warmUpCancel || !a.warmUp || canvasDisposed) return;
      const warmUp = () => {
        warmUpCancel = undefined;
        if (canvasDisposed) return;
        prepare({ interaction: true });
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

    function commit(isSynchronous: boolean = false) {
      const allCommited = commitScene(isSynchronous);
      shaderManager.updateRequired(p);
      // Only reset the camera after the full scene has been commited.
      if (allCommited) {
        resolveCameraReset();
        if (forceDrawAfterAllCommited) {
          if (helper.debug.isEnabled) helper.debug.update();
          draw({ force: true });
          forceDrawAfterAllCommited = false;
        }
        commited.next(now());
      }
    }

    function resolveCameraReset() {
      if (!cameraResetRequested) return;

      if (!xr.isPresenting.value) {
        xrManager.resetScale();
      }

      const boundingSphere = scene.boundingSphereVisible;
      const { center, radius } = boundingSphere;

      if (radius > 0) {
        const duration = nextCameraResetOptions.durationMs ?? p.cameraResetDurationMs;
        const easing = nextCameraResetOptions.easing ?? p.cameraResetEasing;
        const trajectory = nextCameraResetOptions.trajectory ?? p.cameraResetTrajectory;
        const autoAdjustControls = controls.props.autoAdjustMinMaxDistance;
        if (autoAdjustControls.name === 'on') {
          const minDistance =
            autoAdjustControls.params.minDistanceFactor * radius + autoAdjustControls.params.minDistancePadding;
          const maxDistance = Math.max(
            autoAdjustControls.params.maxDistanceFactor * radius,
            autoAdjustControls.params.maxDistanceMin,
          );
          controls.setProps({ minDistance, maxDistance });
        }
        const focus = camera.getFocus(center, radius);
        const next =
          typeof nextCameraResetOptions.snapshot === 'function'
            ? nextCameraResetOptions.snapshot(scene, camera)
            : nextCameraResetOptions.snapshot;
        const snapshot = next ? { ...focus, ...next } : focus;
        camera.setState({ ...snapshot, radiusMax: getSceneRadius() }, duration, {
          keyframes: nextCameraResetOptions.keyframes,
          easing,
          trajectory,
        });
      }

      nextCameraResetOptions.durationMs = void 0;
      nextCameraResetOptions.snapshot = void 0;
      nextCameraResetOptions.keyframes = void 0;
      nextCameraResetOptions.easing = void 0;
      nextCameraResetOptions.trajectory = void 0;

      cameraResetRequested = false;
    }

    const oldBoundingSphereVisible = Sphere3D();
    const cameraSphere = Sphere3D();

    function shouldResetCamera() {
      if (camera.state.radiusMax === 0) return true;

      if (camera.transition.inTransition || nextCameraResetOptions.snapshot) return false;

      let cameraSphereOverlapsNone = true,
        isEmpty = true;
      Sphere3D.set(cameraSphere, camera.state.target, camera.state.radius);

      // check if any renderable has moved outside of the old bounding sphere
      // and if no renderable is overlapping with the camera sphere
      for (const r of scene.renderables) {
        if (!r.state.visible) continue;

        const b = r.values.boundingSphere.ref.value;
        if (!b.radius) continue;

        isEmpty = false;
        const cameraDist = Vec3.distance(cameraSphere.center, b.center);
        if (
          (cameraDist > cameraSphere.radius || cameraDist > b.radius || b.radius > camera.state.radiusMax) &&
          !Sphere3D.includes(oldBoundingSphereVisible, b)
        )
          return true;
        if (Sphere3D.overlaps(cameraSphere, b)) cameraSphereOverlapsNone = false;
      }

      return cameraSphereOverlapsNone || (!isEmpty && cameraSphere.radius <= 0.1);
    }

    const sceneCommitTimeoutMs = 250;
    function commitScene(isSynchronous: boolean) {
      if (!scene.needsCommit) return true;

      // snapshot the current bounding sphere of visible objects
      Sphere3D.copy(oldBoundingSphereVisible, scene.boundingSphereVisible);

      // clear hi-Z buffer when scene changes
      hiZ.clear();

      if (!scene.commit(isSynchronous ? void 0 : sceneCommitTimeoutMs)) {
        commitQueueSize.next(scene.commitQueueSize);
        return false;
      }
      commitQueueSize.next(0);
      precompileRecorder.record(scene.renderables);
      // frames held for the props (`holdFramesUntilPrepared`) wait for the new scene too: the hold's prepare had the
      // previous one, so the first frame after it would compile (and skip) the draws of what was added meanwhile
      if (propsPreparing > 0) holdFramesUntilPrepared();
      scheduleWarmUp();

      if (helper.debug.isEnabled) helper.debug.update();
      if (!p.camera.manualReset && (reprCount.value === 0 || shouldResetCamera())) {
        cameraResetRequested = true;
      }
      if (oldBoundingSphereVisible.radius === 0) nextCameraResetOptions.durationMs = 0;

      if (!p.camera.manualReset && scene.boundingSphere.radius > 0) {
        camera.setState({ radiusMax: getSceneRadius() }, 0);
      }
      reprCount.next(reprRenderObjects.size);
      if (isDebugMode) consoleStats();

      return true;
    }

    function consoleStats() {
      const items = scene.renderables.map((r) => ({
        drawCount: r.values.drawCount.ref.value,
        instanceCount: r.values.instanceCount.ref.value,
        materialId: r.materialId,
        renderItemId: r.id,
        geometryType: r.values.dGeometryType.ref.value,
        'byteCount [MiB]': toFixed(r.getByteCount() / 1024 / 1024, 3),
        radius: toFixed(r.values.boundingSphere.ref.value.radius, 3),
      }));

      console.groupCollapsed(`${items.length} RenderItems`);

      if (items.length <= 64) {
        console.table(items);
      } else {
        console.log(items);
      }
      console.log(JSON.stringify(graphics.stats, undefined, 4));

      const { texture, cubeTexture, attribute, elements, pixelPack, renderbuffer, storageBuffer } =
        graphics.resources.getByteCounts();
      const arena = graphics.compute ? getVertexArenaStats(graphics.compute) : undefined;
      console.log(
        JSON.stringify(
          {
            texture: `${(texture / 1024 / 1024).toFixed(3)} MiB`,
            cubeTexture: `${(cubeTexture / 1024 / 1024).toFixed(3)} MiB`,
            attribute: `${(attribute / 1024 / 1024).toFixed(3)} MiB`,
            elements: `${(elements / 1024 / 1024).toFixed(3)} MiB`,
            pixelPack: `${(pixelPack / 1024 / 1024).toFixed(3)} MiB`,
            renderbuffer: `${(renderbuffer / 1024 / 1024).toFixed(3)} MiB`,
            storageBuffer: `${(storageBuffer / 1024 / 1024).toFixed(3)} MiB`,
            ...(arena
              ? {
                  vertexArena: `${(arena.usedByteCount / 1024 / 1024).toFixed(3)} of ${(arena.byteCount / 1024 / 1024).toFixed(3)} MiB used (${arena.usedPages} of ${arena.pages} pages, ${arena.arenas} arenas${arena.exhausted ? ', out of memory' : ''})`,
                }
              : {}),
          },
          undefined,
          4,
        ),
      );

      console.log(
        JSON.stringify(
          {
            renderables: `${(scene.renderables.reduce((sum, r) => sum + r.getByteCount(), 0) / 1024 / 1024).toFixed(3)} MiB`,
            passes: {
              draw: `${(passes.draw.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
              illumination: `${(passes.illumination.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
              pick: `${(passes.pick.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
              hiZ: `${(hiZ.getByteCount() / 1024 / 1024).toFixed(3)} MiB`,
            },
          },
          undefined,
          4,
        ),
      );

      console.log(
        JSON.stringify(
          {
            scene: {
              radius: toFixed(scene.boundingSphere.radius, 3),
              radiusVisible: toFixed(scene.boundingSphereVisible.radius, 3),
              transparency: scene.transparency,
              markerAverage: toFixed(scene.markerAverage, 2),
              emissiveAverage: toFixed(scene.emissiveAverage, 2),
              wiggleAverage: toFixed(scene.wiggleAverage, 2),
              opacityAverage: toFixed(scene.opacityAverage, 2),
              transparencyMin: toFixed(scene.transparencyMin, 2),
              hasOpaque: scene.hasOpaque,
              hasAnimation: scene.hasAnimation,
              counts: {
                all: scene.renderables.length,
                primitives: scene.primitives.renderables.length,
                volumes: scene.volumes.renderables.length,
              },
            },
          },
          undefined,
          4,
        ),
      );

      if (isTimingMode) {
        console.log(JSON.stringify(graphics.timer.formatedStats(), undefined, 4));
      }

      console.groupEnd();
    }

    /** render objects left out of the scene by a context change until their representation updated them */
    const droppedRenderObjects = new Set<GraphicsRenderObject>();

    function add(repr: Representation.Any) {
      registerAutoUpdate(repr);
      busy();

      const oldRO = reprRenderObjects.get(repr);
      const newRO = new Set<GraphicsRenderObject>();
      repr.renderObjects.forEach((o) => newRO.add(o));

      if (oldRO) {
        if (!SetUtils.areEqual(newRO, oldRO)) {
          newRO.forEach((o) => {
            if (!oldRO.has(o)) scene.add(o);
          });
          oldRO.forEach((o) => {
            if (!newRO.has(o)) scene.remove(o);
          });
        }
        if (droppedRenderObjects.size > 0) {
          oldRO.forEach((o) => {
            if (!droppedRenderObjects.has(o)) return;
            if (!newRO.has(o)) {
              droppedRenderObjects.delete(o);
            } else if (!hasForeignResources(o, graphics)) {
              droppedRenderObjects.delete(o);
              scene.add(o);
            }
          });
        }
      } else {
        repr.renderObjects.forEach((o) => scene.add(o));
      }
      reprRenderObjects.set(repr, newRO);
      markIndex.set(repr);

      scene.update(repr.renderObjects, false);
      forceDrawAfterAllCommited = true;
      if (isDebugMode) consoleStats();
    }

    function remove(repr: Representation.Any) {
      unregisterAutoUpdate(repr);
      busy();
      markIndex.delete(repr);

      const renderObjects = reprRenderObjects.get(repr);
      if (renderObjects) {
        renderObjects.forEach((o) => {
          scene.remove(o);
          droppedRenderObjects.delete(o);
        });
        reprRenderObjects.delete(repr);
        forceDrawAfterAllCommited = true;
        if (isDebugMode) consoleStats();
      }
    }

    function registerAutoUpdate(repr: Representation.Any) {
      if (reprUpdatedSubscriptions.has(repr)) return;

      reprUpdatedSubscriptions.set(
        repr,
        repr.updated.subscribe((_) => {
          if (!repr.state.syncManually) add(repr);
        }),
      );
    }

    function unregisterAutoUpdate(repr: Representation.Any) {
      const updatedSubscription = reprUpdatedSubscriptions.get(repr);
      if (updatedSubscription) {
        updatedSubscription.unsubscribe();
        reprUpdatedSubscriptions.delete(repr);
      }
    }

    function getProps(): Canvas3DProps {
      const radius =
        scene.boundingSphere.radius > 0
          ? 100 - Math.round((camera.transition.target.radius / getSceneRadius()) * 100)
          : 0;

      return {
        camera: {
          mode: camera.state.mode,
          helper: { ...helper.camera.props },
          stereo: { ...p.camera.stereo },
          fov: Math.round(radToDeg(camera.state.fov)),
          manualReset: !!p.camera.manualReset,
        },
        cameraFog:
          camera.state.fog > 0
            ? { name: 'on' as const, params: { intensity: camera.state.fog } }
            : { name: 'off' as const, params: {} },
        cameraClipping: {
          far: camera.state.clipFar,
          radius,
          minNear: camera.state.minNear,
        },
        cameraResetDurationMs: p.cameraResetDurationMs,
        cameraResetEasing: p.cameraResetEasing,
        cameraResetTrajectory: p.cameraResetTrajectory,
        sceneRadiusFactor: p.sceneRadiusFactor,
        transparentBackground: p.transparentBackground,
        checkeredTransparentBackground: p.checkeredTransparentBackground,
        dpoitIterations: p.dpoitIterations,
        volumeResolutionScale: p.volumeResolutionScale,
        pickPadding: p.pickPadding,
        userInteractionReleaseMs: p.userInteractionReleaseMs,
        viewport: p.viewport,

        postprocessing: { ...p.postprocessing },
        marking: { ...p.marking },
        multiSample: { ...p.multiSample },
        illumination: { ...p.illumination },
        hiZ: { ...hiZ.props },
        renderer: { ...renderer.props },
        trackball: { ...controls.props },
        interaction: { ...interactionHelper.props },
        handle: { ...helper.handle.props },
        pointer: { ...helper.pointer.props },
        xr: { ...xrManager.props },
      };
    }

    const contextLostSub = contextLost?.subscribe(() => {
      isContextLost = true;
      pickHelper.dirty = true;
    });

    const resized = new BehaviorSubject<any>(0);

    function handleResize(draw = true) {
      busy();
      passes.updateSize();
      updateViewport();
      syncViewport();
      if (draw) requestDraw();
      resized.next(+new Date());
    }

    addConsoleStatsProvider(consoleStats);

    const ctxChangedSub = ctx.changed?.subscribe(() => {
      scene.setTransparency(passes.draw.transparency);
      // e.g. the pipelines of another transparency mode: frames wait for them as for a props change
      holdFramesUntilPrepared();
      // e.g. the targets of the transparency mode no longer drawn
      scheduleTargetRelease();
      requestDraw();
    });

    // e.g. images loaded into textures: drawn with placeholders until then
    const onPendingSettled = () => requestDraw();
    const outOfMemory = new Subject<string>();
    const subscribeContext = () => [
      graphics.resources.pending.settled.subscribe(onPendingSettled),
      graphics.resources.outOfMemory.subscribe((message) => outOfMemory.next(message)),
    ];
    let contextSubs = subscribeContext();

    const contextChanged = new Subject<Canvas3D.ContextChange>();
    const imagePasses = new Set<WeakRef<ImagePass>>();

    /**
     * Rebuilds what is bound to the graphics context on the context of `ctx` (after `ctx.setBackend`): the renderer,
     * passes and helpers are created again with their props; the scene and helpers keep their render objects except
     * those holding resources of the previous context, which their owners create again (see `contextChanged`)
     */
    function setContext(change: Canvas3DContext.ContextChange) {
      const { renderer: rendererProps, hiZ: hiZProps } = getProps();
      destroyPrecompileItems();
      renderer.dispose();
      hiZ.dispose();
      pickHelper.dispose();
      rayHelper.dispose();
      for (const s of contextSubs) s.unsubscribe();

      graphics = ctx.graphics;
      passes = ctx.passes;
      canvas = ctx.canvas;

      const dropped = scene.setContext(graphics);
      for (const o of dropped) droppedRenderObjects.add(o);
      helper.debug.setContext(graphics);
      helper.camera.setContext(graphics);
      helper.handle.setContext(graphics);
      helper.pointer.setContext(graphics);
      xrManager.setContext(graphics);
      limitToSupported = createSupportedPropsLimiter(graphics);
      limitToSupported(p);

      hiZ = new HiZPass(graphics, passes.draw, canvas, hiZProps);
      p.hiZ = hiZ.props;
      renderer = Renderer.create(graphics, rendererProps);
      renderer.setOcclusionTest(hiZ.occlusionTest);
      renderer.setOcclusionPyramid(getOcclusionPyramid);
      shaderManager = new ShaderManager(graphics, scene);
      pickHelper = new PickHelper(
        graphics,
        renderer,
        scene,
        helper,
        passes.pick,
        { x, y, width, height },
        { ...pickOptions, pickPadding: p.pickPadding },
      );
      rayHelper = new RayHelper(graphics, renderer, scene, helper, { ...pickOptions, pickPadding: p.pickPadding });
      multiSampleHelper = new MultiSampleHelper(passes.multiSample);
      syncViewport();
      interactionHelper.resetPicks(pickHelper.synchronous);
      passes.draw.postprocessing.background.update(camera, p.postprocessing.background, (changed) => {
        if (changed) requestDraw();
      });
      scene.setTransparency(passes.draw.transparency);
      contextSubs = subscribeContext();

      isContextLost = false;
      syncCanvasBackground(canvas, p);
      handleResize(false);
      shaderManager.updateRequired(p);
      forceDrawAfterAllCommited = true;
      requestDraw();

      for (const ref of [...imagePasses]) {
        const pass = ref.deref();
        if (pass) pass.setContext(graphics, renderer, helper);
        else imagePasses.delete(ref);
      }

      // representations continue with the new context and update: visuals create again what was rendered or computed
      // with the previous context and geometry relying on what the new one lacks (see `GeometryUtils.isCompatible`),
      // the others change nothing
      for (const repr of reprRenderObjects.keys()) {
        repr.setGraphics?.(graphics);
        repr
          .createOrUpdate()
          .run()
          .then(() => {
            if (reprRenderObjects.has(repr)) add(repr);
          })
          .catch((e) => console.error('Updating a representation for the new graphics context failed', e));
      }

      const unowned = dropped.filter((o) => droppedRenderObjects.has(o) && !isOwnedByRepr(o));
      contextChanged.next({ dropped: unowned, reason: change.reason });
    }

    function isOwnedByRepr(o: GraphicsRenderObject) {
      for (const renderObjects of reprRenderObjects.values()) {
        if (renderObjects.has(o)) return true;
      }
      return false;
    }
    const ctxContextChangedSub = ctx.contextChanged?.subscribe(setContext);

    // Monitor user interactions
    let isDragging = false;
    let isActivelyInteracting = false;
    const interactionSubs = [
      input.drag.subscribe(() => {
        isDragging = true;
      }),
      input.interactionEnd.subscribe(() => {
        isDragging = false;
      }),
      merge(input.drag, input.pinch, input.wheel, input.interactionEnd).subscribe(() => {
        interactionEvent.next();
      }),
      interactionEvent.subscribe(() => {
        isActivelyInteracting = true;
      }),
      interactionEvent.pipe(debounceTime(p.userInteractionReleaseMs)).subscribe(() => {
        // no requestDraw here: it would force a full illumination restart (flash) even
        // though the next animation frame already resumes accumulation on its own
        isActivelyInteracting = isDragging;
        busy();
      }),
      // any input may give the next tick work (controls, hover, picking): an idle loop asks for frames again
      merge(
        input.move,
        input.drag,
        input.interactionEnd,
        input.wheel,
        input.pinch,
        input.gesture,
        input.click,
        input.leave,
        input.enter,
        input.resize,
        input.modifiers,
        input.key,
        input.keyDown,
        input.keyUp,
        input.lock,
        input.trackedPointers,
      ).subscribe(busy),
      camera.stateChanged.subscribe(busy),
    ];

    //

    if (isDebugMode && canvas) {
      let occlusionLoci: Loci | undefined = undefined;

      const printOcclusion = (loci: Loci | undefined) => {
        const s = loci && Loci.getBoundingSphere(Loci.normalize(loci, 'residue'));
        hiZ.debugOcclusion(s);
      };

      input.click.subscribe(async (e) => {
        if (!e.modifiers.control || e.button !== 2) return;

        const p = await asyncIdentify(Vec2.create(e.x, e.y))?.wait();
        if (!p) {
          occlusionLoci = undefined;
          printOcclusion(occlusionLoci);
          return;
        }

        const l = getLoci(p.id);
        occlusionLoci = l.loci;
        printOcclusion(occlusionLoci);
      });

      didDraw.subscribe(() => {
        setTimeout(() => {
          printOcclusion(occlusionLoci);
        }, 100);
      });
    }

    //

    return {
      get graphics() {
        return graphics;
      },
      contextChanged,
      outOfMemory,

      add,
      remove,
      commit,
      update: (repr, keepSphere) => {
        // updated values may have changed defines: the updated renderables recorded again
        if (repr) {
          if (!reprRenderObjects.has(repr)) return;
          markIndex.set(repr);
          const { renderObjects } = repr;
          scene.update(renderObjects, !!keepSphere);
          // only the representation's, not the whole scene's: callers update one representation at a time (a drag)
          const renderables: GraphicsRenderable[] = [];
          for (let i = 0, il = renderObjects.length; i < il; ++i) {
            const renderable = scene.get(renderObjects[i]);
            if (renderable) renderables.push(renderable);
          }
          precompileRecorder.recordAgain(renderables);
        } else {
          reprRenderObjects.forEach((_, r) => markIndex.set(r));
          scene.update(void 0, !!keepSphere);
          precompileRecorder.invalidate();
          precompileRecorder.record(scene.renderables);
        }
        forceDrawAfterAllCommited = true;
        busy();
      },
      clear: () => {
        reprUpdatedSubscriptions.forEach((v) => v.unsubscribe());
        reprUpdatedSubscriptions.clear();
        reprRenderObjects.clear();
        markIndex.clear();
        scene.clear();
        helper.debug.clear();
        requestDraw();
        reprCount.next(reprRenderObjects.size);
      },
      syncVisibility: () => {
        if (camera.state.radiusMax === 0) {
          cameraResetRequested = true;
          nextCameraResetOptions.durationMs = 0;
        }

        if (scene.syncVisibility()) {
          if (helper.debug.isEnabled) helper.debug.update();
        }
        requestDraw();
      },

      requestDraw,
      prepare,
      precompile,
      getPrecompileList: () => precompileRecorder.getList(),
      resetPrecompileList: () => precompileRecorder.reset(),
      tick,
      animate,
      resetTime,
      pause,
      resume,
      holdFrames,

      requestAnimationFrame: _requestAnimationFrame,
      cancelAnimationFrame: _cancelAnimationFrame,

      identify,
      asyncIdentify,
      mark,
      getLoci,

      handleResize,
      requestResize: () => {
        resizeRequested = true;
        busy();
      },
      requestCameraReset: (options) => {
        Object.assign(nextCameraResetOptions, options);
        cameraResetRequested = true;
        busy();
      },
      camera,
      boundingSphere: scene.boundingSphere,
      boundingSphereVisible: scene.boundingSphereVisible,
      get notifyDidDraw() {
        return notifyDidDraw;
      },
      set notifyDidDraw(v: boolean) {
        notifyDidDraw = v;
      },
      get isIdle() {
        return idle;
      },
      wake,
      didDraw,
      commited,
      commitQueueSize,
      reprCount,
      resized,
      setProps: (properties, doNotRequestDraw = false) => {
        const props: PartialCanvas3DProps =
          typeof properties === 'function' ? produce(getProps(), properties as any) : properties;
        busy();

        if (props.sceneRadiusFactor !== undefined) {
          p.sceneRadiusFactor = props.sceneRadiusFactor;
          camera.setState({ radiusMax: getSceneRadius() }, 0);
        }

        const cameraState: Partial<Camera.Snapshot> = Object.create(null);
        if (props.camera && props.camera.mode !== undefined && props.camera.mode !== camera.state.mode) {
          cameraState.mode = props.camera.mode;
        }
        const oldFov = Math.round(radToDeg(camera.state.fov));
        if (props.camera && props.camera.fov !== undefined && props.camera.fov !== oldFov) {
          cameraState.fov = degToRad(props.camera.fov);
        }
        if (props.cameraFog !== undefined && props.cameraFog.params) {
          const newFog = props.cameraFog.name === 'on' ? props.cameraFog.params.intensity : 0;
          if (newFog !== camera.state.fog) cameraState.fog = newFog;
        }
        if (props.cameraClipping !== undefined) {
          if (props.cameraClipping.far !== undefined && props.cameraClipping.far !== camera.state.clipFar) {
            cameraState.clipFar = props.cameraClipping.far;
          }
          if (props.cameraClipping.minNear !== undefined && props.cameraClipping.minNear !== camera.state.minNear) {
            cameraState.minNear = props.cameraClipping.minNear;
          }
          if (props.cameraClipping.radius !== undefined) {
            const radius = (getSceneRadius() / 100) * (100 - props.cameraClipping.radius);
            if (radius > 0 && radius !== cameraState.radius) {
              // if radius = 0, NaNs happen
              cameraState.radius = Math.max(radius, 0.01);
            }
          }
        }
        if (Object.keys(cameraState).length > 0) camera.setState(cameraState);

        if (props.camera?.helper) helper.camera.setProps(props.camera.helper);
        if (props.camera?.manualReset !== undefined) p.camera.manualReset = props.camera.manualReset;
        if (props.camera?.stereo !== undefined) {
          Object.assign(p.camera.stereo, props.camera.stereo);
          stereoCamera.setProps(p.camera.stereo.params);
        }
        if (props.cameraResetDurationMs !== undefined) p.cameraResetDurationMs = props.cameraResetDurationMs;
        if (props.cameraResetEasing !== undefined) p.cameraResetEasing = props.cameraResetEasing;
        if (props.cameraResetTrajectory !== undefined) p.cameraResetTrajectory = props.cameraResetTrajectory;
        if (props.transparentBackground !== undefined) p.transparentBackground = props.transparentBackground;
        if (props.checkeredTransparentBackground !== undefined)
          p.checkeredTransparentBackground = props.checkeredTransparentBackground;
        if (props.dpoitIterations !== undefined) p.dpoitIterations = props.dpoitIterations;
        if (props.volumeResolutionScale !== undefined) p.volumeResolutionScale = props.volumeResolutionScale;
        if (props.pickPadding !== undefined) {
          p.pickPadding = props.pickPadding;
          pickHelper.setPickPadding(p.pickPadding);
        }
        if (props.userInteractionReleaseMs !== undefined) p.userInteractionReleaseMs = props.userInteractionReleaseMs;
        if (props.viewport !== undefined) {
          const doNotUpdate =
            p.viewport === props.viewport ||
            (p.viewport.name === props.viewport.name && shallowEqual(p.viewport.params, props.viewport.params));

          if (!doNotUpdate) {
            p.viewport = props.viewport;
            updateViewport();
            syncViewport();
          }
        }

        if (props.postprocessing?.background) {
          Object.assign(p.postprocessing.background, props.postprocessing.background);
          passes.draw.postprocessing.background.update(camera, p.postprocessing.background, (changed) => {
            if (changed && !doNotRequestDraw) requestDraw();
          });
        }
        if (props.postprocessing) Object.assign(p.postprocessing, props.postprocessing);
        if (props.marking) Object.assign(p.marking, props.marking);
        if (props.illumination) Object.assign(p.illumination, props.illumination);
        if (props.multiSample) Object.assign(p.multiSample, props.multiSample);
        if (props.postprocessing || props.marking || props.multiSample) limitToSupported(p);
        if (props.hiZ) hiZ.setProps(props.hiZ);
        if (props.renderer || props.marking) {
          scene.setGlobals({
            dColorMarker: isMaterialColorMarker(
              {
                colorMarker: props.renderer?.colorMarker ?? renderer.props.colorMarker,
              },
              p.marking,
            ),
            dLightCount: Math.min(props.renderer?.light?.length ?? renderer.props.light.length, MaxLightCount),
          });
        }
        if (props.renderer) renderer.setProps(props.renderer);
        if (props.trackball) controls.setProps(props.trackball);
        if (props.interaction) interactionHelper.setProps(props.interaction);
        if (props.handle) helper.handle.setProps(props.handle);
        if (props.pointer) helper.pointer.setProps(props.pointer);
        if (props.xr) xrManager.setProps(props.xr);

        if (cameraState.mode === 'orthographic') {
          p.camera.stereo.name = 'off';
        }

        if (
          'transparentBackground' in props ||
          'checkeredTransparentBackground' in props ||
          (props.renderer && 'backgroundColor' in props.renderer)
        ) {
          syncCanvasBackground(canvas, p);
        }

        shaderManager.updateRequired(p);
        holdFramesUntilPrepared();
        scheduleTargetRelease();
        if (!doNotRequestDraw) {
          requestDraw();
        }
      },
      setAttribs: (attribs: PartialCanvas3DAttribs) => {
        busy();
        if (attribs.trackball) controls.setAttribs(attribs.trackball);
        if (attribs.xr) xrManager.setAttribs(attribs.xr);
        if (attribs.warmUp !== undefined) a.warmUp = attribs.warmUp;
      },
      getImagePass: (props: Partial<ImageProps> = {}) => {
        const pass = new ImagePass(graphics, assetManager, renderer, scene, camera, helper, props);
        imagePasses.add(new WeakRef(pass));
        return pass;
      },
      getRenderObjects(): GraphicsRenderObject[] {
        const renderObjects: GraphicsRenderObject[] = [];
        scene.forEach((_, ro) => renderObjects.push(ro));
        return renderObjects;
      },

      get props() {
        return getProps();
      },
      get attribs() {
        return {
          trackball: controls.attribs,
          xr: xrManager.attribs,
          warmUp: a.warmUp,
        };
      },
      get input() {
        return input;
      },
      get stats() {
        return renderer.stats;
      },
      get interaction() {
        return interactionHelper.events;
      },
      debugRegistry: helper.debug,
      xr,
      dispose: () => {
        canvasDisposed = true;
        warmUpCancel?.();
        warmUpCancel = undefined;
        passes.targets.cancelRelease();
        contextLostSub?.unsubscribe();
        ctxChangedSub?.unsubscribe();
        ctxContextChangedSub?.unsubscribe();
        for (const s of contextSubs) s.unsubscribe();

        for (const s of xrSubs) s.unsubscribe();
        xrSubs.length = 0;

        for (const s of interactionSubs) s.unsubscribe();
        interactionSubs.length = 0;

        cancelAnimationFrame(animationFrameHandle);
        animationFrameCB = undefined;

        markBuffer = [];

        scene.clear();
        helper.debug.clear();
        destroyPrecompileItems();
        controls.dispose();
        renderer.dispose();
        interactionHelper.dispose();
        hiZ.dispose();
        pickHelper.dispose();
        rayHelper.dispose();
        xrManager.dispose();

        reprCount.complete();
        interactionEvent.complete();
        didDraw.complete();
        wake.complete();
        resized.complete();
        commited.complete();
        commitQueueSize.complete();
        xr.isPresenting.complete();
        xr.isSupported.complete();
        xr.requestFailed.complete();
        contextChanged.complete();
        outOfMemory.complete();

        removeConsoleStatsProvider(consoleStats);
      },
    };

    function updateViewport() {
      const oldX = x,
        oldY = y,
        oldWidth = width,
        oldHeight = height;
      const drs = graphics.getDrawingBufferSize();

      if (p.viewport.name === 'canvas') {
        x = 0;
        y = 0;
        width = drs.width;
        height = drs.height;
      } else if (p.viewport.name === 'static-frame') {
        x = p.viewport.params.x * graphics.pixelRatio;
        height = p.viewport.params.height * graphics.pixelRatio;
        y = p.viewport.params.y * graphics.pixelRatio;
        width = p.viewport.params.width * graphics.pixelRatio;
      } else if (p.viewport.name === 'relative-frame') {
        x = Math.round(p.viewport.params.x * drs.width);
        height = Math.round(p.viewport.params.height * drs.height);
        y = Math.round(p.viewport.params.y * drs.height);
        width = Math.round(p.viewport.params.width * drs.width);
      }

      if (oldX !== x || oldY !== y || oldWidth !== width || oldHeight !== height) {
        forceNextRender = true;
      }
    }

    function syncViewport() {
      pickHelper.setViewport(x, y, width, height);
      renderer.setViewport(x, y, width, height);
      Viewport.set(camera.viewport, x, y, width, height);
      Viewport.set(controls.viewport, x, y, width, height);
      hiZ.setViewport(x, y, width, height);
    }
  }
}
