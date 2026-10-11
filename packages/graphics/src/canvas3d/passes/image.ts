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

export const ImageParams = {
  transparentBackground: PD.Boolean(false),
  dpoitIterations: PD.Numeric(2, { min: 1, max: 10, step: 1 }),
  multiSample: PD.Group(MultiSampleParams),
  postprocessing: PD.Group(PostprocessingParams),
  marking: PD.Group(MarkingParams),
  illumination: PD.Group(IlluminationParams),

  cameraHelper: PD.Group(CameraHelperParams),
  renderer: PD.Group(RendererParams),
};
export type ImageProps = PD.Values<typeof ImageParams>;

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

  private readonly drawPass: DrawPass;
  private readonly illuminationPass: IlluminationPass;
  private readonly multiSamplePass: MultiSamplePass;
  private readonly multiSampleHelper: MultiSampleHelper;
  private readonly helper: Helper;
  private readonly limitToSupported: (props: ImageProps) => void;

  get width() {
    return this._width;
  }
  get height() {
    return this._height;
  }

  constructor(
    private graphics: GraphicsContext,
    assetManager: AssetManager,
    private renderer: Renderer,
    private scene: Scene,
    private camera: Camera,
    helper: Helper,
    props: Partial<ImageProps>,
  ) {
    this.props = { ...PD.getDefaultValues(ImageParams), ...props };
    this.limitToSupported = createSupportedPropsLimiter(graphics);
    this.limitToSupported(this.props);

    this.drawPass = new DrawPass(graphics, assetManager, 128, 128, scene.transparency);
    this.illuminationPass = new IlluminationPass(graphics, this.drawPass);
    this.multiSamplePass = new MultiSamplePass(graphics, this.drawPass);
    this.multiSampleHelper = new MultiSampleHelper(this.multiSamplePass);

    this.helper = {
      camera: new CameraHelper(graphics, this.props.cameraHelper),
      debug: helper.debug,
      handle: helper.handle,
      pointer: helper.pointer,
    };

    this.setSize(1024, 768);
  }

  getByteCount() {
    return this.drawPass.getByteCount() + this.illuminationPass.getByteCount() + this.multiSamplePass.getByteCount();
  }

  updateBackground() {
    return new Promise<void>((resolve) => {
      this.drawPass.postprocessing.background.update(this.camera, this.props.postprocessing.background, () => {
        resolve();
      });
    });
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

  async render(runtime: RuntimeContext) {
    // draws whose programs compile in the background (WebGPU pipelines) are skipped: render again if a render started
    // such work (it may have settled during the render, e.g. between illumination iterations), so the image is complete
    // and doesn't depend on timing. A render reaches all draws, so one more render suffices.
    const { pending } = this.graphics.resources;
    const trackedCount = pending.trackedCount;
    await this.renderOnce(runtime);
    if (pending.trackedCount === trackedCount) return;
    const retrackedCount = pending.trackedCount;
    await this.renderOnce(runtime);
    if (isDebugMode && pending.trackedCount !== retrackedCount) {
      console.warn('ImagePass: rendering again started more work, the image may be incomplete');
    }
  }

  private async renderOnce(runtime: RuntimeContext) {
    this.drawPass.setTransparency(this.scene.transparency);
    ShaderManager.ensureRequired(this.graphics, this.scene, this.props);
    // e.g. images still loading into textures (SMAA lookups), so that the image doesn't depend on timing
    await this.graphics.resources.pending.whenReady();
    Camera.copySnapshot(this._camera.state, this.camera.state);
    Viewport.set(this._camera.viewport, 0, 0, this._width, this._height);
    this._camera.update();
    this.renderer.setOcclusionTest(null);
    this._frame = createFrame();

    const ctx = {
      renderer: this.renderer,
      camera: this._camera,
      scene: this.scene,
      helper: this.helper,
      frame: this._frame,
    };
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
        }
        await this.graphics.waitForGpuCommandsComplete();
      }
      this._colorTarget = this.illuminationPass.colorTarget;
    } else {
      if (isTimingMode) this.graphics.timer.group('ImagePass.render', { captureStats: true });
      if (MultiSamplePass.isEnabled(this.props.multiSample)) {
        this.multiSampleHelper.render(ctx, this.props, false);
        this._colorTarget = this.multiSamplePass.colorTarget;
      } else {
        this.drawPass.render(ctx, this.props, false);
        this._colorTarget = this.drawPass.getColorTarget(this.props.postprocessing);
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
    this.setSize(width, height);
    await this.render(runtime);

    const w = viewport?.width ?? width,
      h = viewport?.height ?? height;
    // `viewport` has a top-left origin, the read rectangle a bottom-left one
    const x = viewport?.x ?? 0;
    const y = viewport ? height - viewport.y - viewport.height : 0;
    const array = await this.graphics.readRenderTarget(this.colorTarget, x, y, w, h);
    PixelData.divideByAlpha(PixelData.create(array, w, h));
    return new ImageData(new Uint8ClampedArray(array), w, h);
  }
}
