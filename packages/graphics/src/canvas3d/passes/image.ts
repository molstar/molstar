/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { type Renderer, RendererParams } from '@molstar/graphics/gl/renderer';
import { type Frame, createFrame } from '@molstar/graphics/gl/renderable';
import type { Scene } from '@molstar/graphics/gl/scene';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DrawPass } from './draw.js';
import { PostprocessingParams } from './postprocessing.js';
import { MultiSamplePass, MultiSampleParams, MultiSampleHelper } from './multi-sample.js';
import { createSupportedPropsLimiter } from './supported-props.js';
import { Camera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import { PixelData } from '@molstar/core/util/image';
import type { Helper } from '../helper/helper.js';
import { CameraHelper, CameraHelperParams } from '../helper/camera-helper.js';
import { MarkingParams } from './marking.js';
import type { AssetManager } from '@molstar/core/util/assets';
import { IlluminationParams, IlluminationPass } from './illumination.js';
import type { RuntimeContext } from '@molstar/core/task';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import { printTimerResults } from '@molstar/graphics/gl/timer';
import { ShaderManager } from '../helper/shader-manager.js';
import type { RenderContext } from '../util.js';
import { createPrepareBuilder, type PassBuilder } from './pass-builder.js';
import { PassTargets } from './pass-targets.js';

export const ImageParams = {
  transparentBackground: PD.Boolean(false),
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
  multiSample: PD.Group(MultiSampleParams),
  postprocessing: PD.Group(PostprocessingParams),
  marking: PD.Group(MarkingParams),
  illumination: PD.Group(IlluminationParams),

  cameraHelper: PD.Group(CameraHelperParams),
  renderer: PD.Group(RendererParams),
};
export type ImageProps = PD.Values<typeof ImageParams>;

/** Renders of `ImagePass.render` after the first one, while renders start pending work */
const MaxRerenders = 3;
/** context changes (`setContext`) a render or read starts over for, before it gives up */
const MaxContextChanges = 3;

export class ImagePass {
  private _width = 0;
  private _height = 0;
  private _camera = new Camera();
  /** bumped once per render() call, shared by every sub-render issued during that call */
  private _frame: Frame = createFrame();

  readonly props: ImageProps;

  private _colorTarget: RenderTarget;
  get colorTarget() {
    return this._colorTarget;
  }

  /** the targets of the passes, see `PassTargets` */
  private targets!: PassTargets;
  private drawPass!: DrawPass;
  private illuminationPass!: IlluminationPass;
  private multiSamplePass!: MultiSamplePass;
  private multiSampleHelper!: MultiSampleHelper;
  private helper!: Helper;
  private limitToSupported!: (props: ImageProps) => void;

  get width() {
    return this._width;
  }
  get height() {
    return this._height;
  }

  constructor(
    private graphics: GraphicsContext,
    private assetManager: AssetManager,
    private renderer: Renderer,
    private scene: Scene,
    private camera: Camera,
    helper: Helper,
    props: Partial<ImageProps>,
  ) {
    this.props = { ...PD.getDefaultValues(ImageParams), ...props };
    this.init(helper);
    this.setSize(1024, 768);
  }

  private init(helper: Helper) {
    const { graphics } = this;
    this.limitToSupported = createSupportedPropsLimiter(graphics);
    this.limitToSupported(this.props);

    this.targets = new PassTargets(graphics);
    this.drawPass = new DrawPass(this.targets, this.assetManager, 128, 128, this.scene.transparency);
    this.illuminationPass = new IlluminationPass(graphics, this.drawPass);
    this.multiSamplePass = new MultiSamplePass(graphics, this.drawPass);
    this.multiSampleHelper = new MultiSampleHelper(this.multiSamplePass);

    this.helper = {
      camera: new CameraHelper(graphics, this.props.cameraHelper),
      debug: helper.debug,
      handle: helper.handle,
      pointer: helper.pointer,
    };
  }

  /**
   * Renders with `graphics`, `renderer` and `helper` from now on, the passes are created again (the graphics context
   * changed, see `Canvas3D.contextChanged`)
   */
  /** advanced by `setContext`: a render or read awaited across it starts over on the new context */
  private contextVersion = 0;

  setContext(graphics: GraphicsContext, renderer: Renderer, helper: Helper) {
    this.contextVersion += 1;
    this.graphics = graphics;
    this.renderer = renderer;
    this.init(helper);
    const [width, height] = [this._width, this._height];
    this._width = 0;
    this._height = 0;
    this.setSize(width, height);
  }

  getByteCount() {
    return this.drawPass.getByteCount() + this.illuminationPass.getByteCount() + this.multiSamplePass.getByteCount();
  }

  /** Updates the background to the props and resolves when its images are loaded (pending work, see `render`) */
  async updateBackground() {
    this.drawPass.postprocessing.background.update(this.camera, this.props.postprocessing.background);
    await this.graphics.resources.pending.whenReady();
  }

  setSize(width: number, height: number) {
    if (width === this._width && height === this._height) return;

    this._width = width;
    this._height = height;

    this.drawPass.setSize(width, height);
    this.illuminationPass.setSize(width, height);
    this.multiSamplePass.syncSize();
  }

  setProps(props: Partial<ImageProps> = {}) {
    Object.assign(this.props, props);
    this.limitToSupported(this.props);
    if (props.cameraHelper) this.helper.camera.setProps(props.cameraHelper);
  }

  /**
   * Renders the image (`colorTarget`); the targets' storage is released a while after the last render (see
   * `PassTargets.scheduleRelease`), so that an image pass kept for later screenshots takes no memory meanwhile
   */
  async render(runtime: RuntimeContext) {
    this.targets.cancelRelease();
    try {
      await this.renderUnreleased(runtime);
    } finally {
      this.targets.scheduleRelease(() => this.targets.releaseAll());
    }
  }

  private async renderUnreleased(runtime: RuntimeContext) {
    // `prepare` makes what the render draws with ready, so that no draw is skipped while a WebGPU pipeline compiles;
    // should a render still start work (it may have settled during the render, e.g. between illumination
    // iterations), render again, so that the image is complete and doesn't depend on timing; the context changing
    // during the render (`setContext`: a backend switch or restore) starts it over on the new one
    let changes = 0;
    for (let i = 0; i <= MaxRerenders; ++i) {
      const version = this.contextVersion;
      await this.prepare();
      const { pending } = this.graphics.resources;
      const trackedCount = pending.trackedCount;
      if (version === this.contextVersion) await this.renderOnce(runtime, version);
      if (version !== this.contextVersion) {
        if (++changes > MaxContextChanges) {
          console.warn('ImagePass: the context kept changing during the render, the image may be incomplete');
          return;
        }
        i -= 1;
        continue;
      }
      if (pending.trackedCount === trackedCount) return;
      if (isDebugMode) {
        console.warn(
          `ImagePass: the render started work 'prepare' missed: ${pending.trackedSince(trackedCount).join(', ')}`,
        );
      }
    }
    if (isDebugMode) console.warn('ImagePass: rendering again started more work, the image may be incomplete');
  }

  /**
   * Makes the pipelines and resources `render` needs ready without drawing anything (see `PassBuilder`): the
   * programs of the scene's variants, the images of the background and lookups, and the pipelines of every pass the
   * props enable; resolves when they are
   */
  async prepare() {
    this.drawPass.setTransparency(this.scene.transparency);
    ShaderManager.ensureRequired(this.graphics, this.scene, this.props);
    // the background is updated before the render: the draw pass chooses its clear from the background's state
    // (`isEnabled`), which the postprocessing pass would only update during the render, and a loaded image enables
    // its draws
    await this.updateBackground();
    const b = createPrepareBuilder(this.graphics, this.renderer);
    this.plan(b, this.getRenderContext());
    await b.whenReady();
  }

  private plan(b: PassBuilder, ctx: RenderContext<Camera>) {
    if (this.illuminationPass.supported && this.props.illumination.enabled) {
      this.illuminationPass.plan(b, ctx, this.props, false);
    } else if (MultiSamplePass.isEnabled(this.props.multiSample)) {
      this.multiSampleHelper.plan(b, ctx, this.props, false);
    } else {
      this.drawPass.plan(b, ctx, this.props, false);
    }
  }

  private getRenderContext(): RenderContext<Camera> {
    Camera.copySnapshot(this._camera.state, this.camera.state);
    Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
    this._camera.update();
    this.renderer.setOcclusionTest(null);
    this._frame = createFrame();
    return {
      renderer: this.renderer,
      camera: this._camera,
      scene: this.scene,
      helper: this.helper,
      frame: this._frame,
    };
  }

  /** Renders on the context of `version`, returning early when it changed (`setContext`) */
  private async renderOnce(runtime: RuntimeContext, version: number) {
    const ctx = this.getRenderContext();
    if (this.illuminationPass.supported && this.props.illumination.enabled) {
      await runtime.update({
        message: 'Tracing...',
        current: 1,
        max: this.illuminationPass.getMaxIterations(this.props.illumination),
      });
      this.illuminationPass.restart(true);
      while (this.illuminationPass.shouldRender(this.props.illumination)) {
        if (isTimingMode) this.graphics.timer.group('ImagePass.render', { captureStats: true });
        this.illuminationPass.render(ctx, this.props, false);
        if (isTimingMode) this.graphics.timer.groupEnd('ImagePass.render');
        if (runtime.shouldUpdate) {
          await runtime.update({ current: this.illuminationPass.iteration });
          if (version !== this.contextVersion) return;
        }
        await this.graphics.waitForGpuCommandsComplete();
        if (version !== this.contextVersion) return;
      }
      this._colorTarget = this.illuminationPass.colorTarget;
    } else {
      if (isTimingMode) this.graphics.timer.group('ImagePass.render', { captureStats: true });
      if (MultiSamplePass.isEnabled(this.props.multiSample)) {
        this.multiSampleHelper.render(ctx, this.props, false);
        this._colorTarget = this.multiSamplePass.colorTarget;
      } else {
        this.drawPass.render(ctx, this.props, false);
        this._colorTarget = this.drawPass.getColorTarget(this.scene, this.props.postprocessing);
      }
      if (isTimingMode) this.graphics.timer.groupEnd('ImagePass.render');
    }

    if (isTimingMode) {
      const timerResults = this.graphics.timer.resolve();
      if (timerResults) {
        for (const result of timerResults) {
          printTimerResults([result]);
        }
      }
    }

    if (isDebugMode) {
      console.log(`image pass byte count ${(this.getByteCount() / 1024 / 1024).toFixed(3)} MiB`);
    }
  }

  async getImageData(runtime: RuntimeContext, width: number, height: number, viewport?: Viewport) {
    const w = viewport?.width ?? width,
      h = viewport?.height ?? height;
    const x = viewport?.x ?? 0;
    const y = viewport?.y ?? 0;
    // a read in flight while the context changes (it may fail with the previous one) renders again on the new one
    for (let changes = 0; ; ++changes) {
      const version = this.contextVersion;
      this.setSize(width, height);
      await this.render(runtime);
      let array: Uint8Array;
      try {
        array = await this.graphics.readRenderTarget(this.colorTarget, x, y, w, h);
      } catch (e) {
        if (version === this.contextVersion || changes >= MaxContextChanges) throw e;
        continue;
      }
      if (version !== this.contextVersion && changes < MaxContextChanges) continue;
      PixelData.divideByAlpha(PixelData.create(array, w, h));
      return new ImageData(new Uint8ClampedArray(array), w, h);
    }
  }
}
