/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import { ValueCell } from '@molstar/core/util';
import { Vec2 } from '@molstar/core/math/linear-algebra';

import { StereoCamera } from '../camera/stereo.js';
import { WboitPass } from './wboit.js';
import { DpoitPass } from './dpoit.js';
import { AntialiasingPass, PostprocessingPass, type PostprocessingProps } from './postprocessing.js';
import { MarkingPass, type MarkingProps, type MarkingShading, SingleSample } from './marking.js';
import { type CopyRenderable, createCopyRenderable } from '@molstar/graphics/gl/compute/util';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import type { AssetManager } from '@molstar/core/util/assets';
import { DofPass } from './dof.js';
import { BloomPass } from './bloom.js';
import { SsaoPass } from './ssao.js';
import { VolumePass } from './volume.js';
import type { RenderContext } from '../util.js';
import { createPrepareBuilder, type PassBuilder } from './pass-builder.js';
import type { PassTargets } from './pass-targets.js';

type Props = {
  postprocessing: PostprocessingProps;
  marking: MarkingProps;
  transparentBackground: boolean;
  dpoitIterations: number;
  /** see `Canvas3DParams.volumeResolutionScale` */
  volumeResolutionScale: number;
};

export type TransparencyMode = 'wboit' | 'dpoit' | 'blended';

export class DrawPass {
  private readonly graphics: GraphicsContext;
  private readonly drawTarget: RenderTarget;

  readonly colorTarget: RenderTarget;
  readonly transparentColorTarget: RenderTarget;
  readonly depthTextureTransparent: Texture;
  readonly depthTextureOpaque: Texture;

  readonly packedDepth: boolean;

  readonly depthTargetTransparent: RenderTarget;
  private depthTargetOpaque: RenderTarget | null;

  private copyFbo: CopyRenderable;
  private readonly renderBuilder: PassBuilder;

  readonly wboit: WboitPass;
  readonly dpoit: DpoitPass;
  readonly marking: MarkingPass;
  readonly postprocessing: PostprocessingPass;
  private readonly ssaoShading: MarkingShading;
  readonly antialiasing: AntialiasingPass;
  readonly dof: DofPass;
  readonly volume: VolumePass;
  /** `Props.volumeResolutionScale` of the eye planned */
  private volumeScale = 1;

  private transparencyMode: TransparencyMode = 'blended';
  /** `transparency` as this context draws it: blended where the mode is not supported */
  resolveTransparency(transparency: TransparencyMode): TransparencyMode {
    if (transparency === 'wboit' && !this.wboit.supported) return 'blended';
    if (transparency === 'dpoit' && !this.dpoit.supported) return 'blended';
    return transparency;
  }
  setTransparency(transparency: TransparencyMode) {
    this.transparencyMode = this.resolveTransparency(transparency);
    if (isDebugMode && this.transparencyMode !== transparency) {
      console.log(`Missing "${transparency}" support, falling back to "blended".`);
    }
    this.resetTransparencyAttachments();
  }
  private resetTransparencyAttachments() {
    this.depthTextureOpaque.detachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
    this.marking.invalidate();
  }
  /**
   * Runs `plan` with the passes planning for `transparency` (as this context draws it) instead of the current mode,
   * e.g. to prepare for props not set yet; the current mode is restored after
   */
  withTransparency<T>(transparency: TransparencyMode, plan: () => T): T {
    const mode = this.resolveTransparency(transparency);
    const current = this.transparencyMode;
    if (mode === current) return plan();
    this.transparencyMode = mode;
    this.resetTransparencyAttachments();
    try {
      return plan();
    } finally {
      this.transparencyMode = current;
      this.resetTransparencyAttachments();
    }
  }
  get transparency() {
    return this.transparencyMode;
  }

  constructor(
    /** the targets of this pass and of the passes drawn with it, see `PassTargets` */
    readonly targets: PassTargets,
    assetManager: AssetManager,
    width: number,
    height: number,
    transparency: 'wboit' | 'dpoit' | 'blended',
  ) {
    const graphics = (this.graphics = targets.graphics);
    const { capabilities } = graphics;
    this.drawTarget = graphics.createDrawTarget();
    this.colorTarget = targets.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'linear');
    this.transparentColorTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');

    this.packedDepth = !capabilities.depthTexture;

    this.depthTargetTransparent = targets.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'nearest');
    this.depthTextureTransparent = this.depthTargetTransparent.texture;

    this.depthTargetOpaque = this.packedDepth ? targets.createRenderTarget(width, height, 'depth-stencil') : null;

    this.depthTextureOpaque = this.depthTargetOpaque
      ? this.depthTargetOpaque.texture
      : targets.texture('image-depth', 'depth-stencil', 'float-stencil', 'nearest');
    if (!this.packedDepth) {
      this.depthTextureOpaque.define(width, height);
    }

    // the transparency and marking passes draw with the depth-stencil buffer of the transparent depth target, used
    // before them in a frame: for the caps' stencil and the depth impostors write, not read across the passes
    const transparencyDepthStencil = this.depthTargetTransparent.depthRenderbuffer!;
    this.wboit = new WboitPass(targets, width, height, transparencyDepthStencil);
    this.dpoit = new DpoitPass(targets, width, height, transparencyDepthStencil);
    this.marking = new MarkingPass(targets, width, height, transparencyDepthStencil);
    this.postprocessing = new PostprocessingPass(graphics, assetManager, this);
    this.ssaoShading = {
      name: 'ssao',
      ssao: this.postprocessing.ssao.ssaoDepthTexture,
      depth: this.depthTextureOpaque,
      packedDepth: this.packedDepth,
    };
    this.antialiasing = new AntialiasingPass(targets, width, height);
    this.dof = new DofPass(targets, width, height);
    this.volume = new VolumePass(targets);

    this.copyFbo = createCopyRenderable(graphics, this.colorTarget.texture);
    this.renderBuilder = targets.createRenderBuilder();

    this.setTransparency(transparency);
  }

  getByteCount() {
    return (
      this.drawTarget.getByteCount() +
      this.colorTarget.getByteCount() +
      this.transparentColorTarget.getByteCount() +
      this.depthTargetTransparent.getByteCount() +
      (this.depthTargetOpaque ? this.depthTargetOpaque.getByteCount() : this.depthTextureOpaque.getByteCount()) +
      this.wboit.getByteCount() +
      this.dpoit.getByteCount() +
      this.marking.getByteCount() +
      this.postprocessing.getByteCount() +
      this.antialiasing.getByteCount() +
      this.dof.getByteCount() +
      this.volume.getByteCount()
    );
  }

  setSize(width: number, height: number) {
    const w = this.colorTarget.getWidth();
    const h = this.colorTarget.getHeight();

    if (width !== w || height !== h) {
      this.colorTarget.setSize(width, height);
      this.depthTargetTransparent.setSize(width, height);
      this.transparentColorTarget.setSize(width, height);

      if (this.depthTargetOpaque) {
        this.depthTargetOpaque.setSize(width, height);
      } else {
        this.depthTextureOpaque.define(width, height);
      }

      ValueCell.update(this.copyFbo.values.uTexSize, Vec2.set(this.copyFbo.values.uTexSize.ref.value, width, height));
    }

    if (this.wboit.supported) {
      this.wboit.setSize(width, height);
    }

    if (this.dpoit.supported) {
      this.dpoit.setSize(width, height);
    }

    this.marking.setSize(width, height);
    this.postprocessing.setSize(width, height);
    this.antialiasing.setSize(width, height);
    this.dof.setSize(width, height);
  }

  renderEmissiveBloom(b: PassBuilder, renderer: Renderer, camera: ICamera, scene: Scene, transparency: boolean): void {
    const bloom = this.postprocessing.bloom;
    b.begin(bloom.emissiveTarget, { label: 'emissive' });
    // (0,0,0,0) clear so the MAX blend in renderEmissiveTransparent isn't polluted by a white clear
    if (b.drawing) renderer.clear(false, true, true);
    // occlude emitters against real opaque depth so glow can't bleed through opaque foreground; packed depth builds its own
    const occludeWithOpaqueDepth = !this.packedDepth;
    if (occludeWithOpaqueDepth) {
      this.depthTextureOpaque.attachFramebuffer(bloom.emissiveTarget.framebuffer, 'depth');
      b.begin(bloom.emissiveTarget, { label: 'emissive' });
    }
    b.renderer({ name: 'emissiveOpaque', occludeWithOpaqueDepth }, scene.primitives, () =>
      renderer.renderEmissiveOpaque(scene.primitives, camera, this.depthTextureTransparent, occludeWithOpaqueDepth),
    );
    if (transparency && scene.opacityAverage < 1) {
      b.renderer({ name: 'emissiveTransparent' }, scene.primitives, () =>
        renderer.renderEmissiveTransparent(scene.primitives, camera, this.depthTextureTransparent),
      );
    }
    if (occludeWithOpaqueDepth) {
      this.depthTextureOpaque.detachFramebuffer(bloom.emissiveTarget.framebuffer, 'depth');
      bloom.emissiveTarget.depthRenderbuffer?.attachFramebuffer(bloom.emissiveTarget.framebuffer);
    }
  }

  private _renderBloom(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    postprocessingProps: PostprocessingProps,
  ): boolean {
    if (!BloomPass.isActive(scene, postprocessingProps) || postprocessingProps.bloom.name !== 'on') return false;
    const bloom = this.postprocessing.bloom;
    const params = postprocessingProps.bloom.params;
    const emissiveBloom = params.mode === 'emissive';

    if (emissiveBloom) {
      this.renderEmissiveBloom(b, renderer, camera, scene, params.transparency);
    }

    bloom.update(
      this.colorTarget.texture,
      this.transparentColorTarget.texture,
      bloom.emissiveTarget.texture,
      this.depthTextureOpaque,
      this.depthTextureTransparent,
      params,
      camera,
      true,
      // opaque only: the transparent targets are neither drawn nor read
      scene.opacityAverage < 1,
    );
    bloom.render(b, camera.viewport);
    return true;
  }

  private renderPostprocessing(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    transparentBackground: boolean,
    postprocessingProps: PostprocessingProps,
  ) {
    this.postprocessing.render(
      b,
      camera,
      scene,
      false,
      transparentBackground,
      renderer.props.backgroundColor,
      postprocessingProps,
      renderer.light,
      renderer.ambientColor,
      this._renderBloom(b, renderer, camera, scene, postprocessingProps),
    );
  }

  /**
   * Draws the volumes at the reduced resolution (`VolumePass`) onto `target`; the opaque depth, which the pass reads, is
   * detached from `target` meanwhile where `depthAttached`
   */
  private renderVolumesReduced(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    target: RenderTarget,
    depthAttached: boolean,
  ) {
    if (depthAttached) this.depthTextureOpaque.detachFramebuffer(target.framebuffer, 'depth');
    this.volume.render(
      b,
      renderer,
      camera,
      scene,
      this.depthTextureOpaque,
      target,
      this.volumeScale,
      this.transparencyMode,
    );
    if (depthAttached) {
      this.depthTextureOpaque.attachFramebuffer(target.framebuffer, 'depth');
      b.begin(target, { label: 'volume composite' });
    }
  }

  /**
   * Draws the volumes onto `target` (at the reduced resolution, or as the transparency mode draws them at the full
   * size); `depthAttached`: the opaque depth is attached to `target` (see `renderVolumesReduced`)
   */
  private renderVolumes(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    target: RenderTarget,
    depthAttached: boolean,
  ) {
    if (VolumePass.isEnabled(this.volumeScale, this.packedDepth)) {
      this.renderVolumesReduced(b, renderer, camera, scene, target, depthAttached);
    } else if (this.transparencyMode === 'wboit' && this.wboit.supported) {
      // a weighted blended layer of their own
      this.wboit.bind(b);
      b.renderer({ name: 'wboitTransparent' }, scene.volumes, () =>
        renderer.renderWboitTransparent(scene.volumes, camera, this.depthTextureOpaque),
      );
      b.begin(target, { label: 'wboit evaluate' });
      this.wboit.render(b);
    } else {
      b.begin(target, { label: 'volume' });
      b.renderer({ name: 'volume' }, scene.volumes, () =>
        renderer.renderVolume(scene.volumes, camera, this.depthTextureOpaque),
      );
    }
  }

  private _renderDpoit(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    iterations: number,
    transparentBackground: boolean,
    postprocessingProps: PostprocessingProps,
  ) {
    if (!this.dpoit.supported) throw new Error('expected dpoit to be supported');

    this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
    b.begin(this.colorTarget, { label: 'opaque' });
    if (b.drawing) renderer.clear(true);

    // render opaque primitives
    if (scene.hasOpaque) {
      b.renderer({ name: 'opaque' }, scene.primitives, () => renderer.renderOpaque(scene.primitives, camera));
    }

    this.depthTextureOpaque.detachFramebuffer(this.colorTarget.framebuffer, 'depth');

    if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
      b.begin(this.depthTargetTransparent, { label: 'transparent depth' });
      if (b.drawing) renderer.clearDepth(true);
      if (scene.opacityAverage < 1) {
        b.renderer({ name: 'depthTransparent' }, scene.primitives, () =>
          renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque),
        );
      }
    }

    // render transparent primitives
    const isPostprocessingEnabled = PostprocessingPass.isEnabled(scene, postprocessingProps);
    const hasVolumes = scene.volumes.renderables.length > 0;
    const hasTransparent = scene.opacityAverage < 1;
    if (hasTransparent) {
      const target = isPostprocessingEnabled ? this.transparentColorTarget : this.colorTarget;
      if (isPostprocessingEnabled && b.drawing) {
        b.begin(target, { label: 'transparent clear' });
        renderer.clear(false, false, true);
      }

      // volumes under the transparent primitives (one layer, not split at their surfaces)
      if (hasVolumes) this.renderVolumes(b, renderer, camera, scene, target, false);

      const dpoitTextures = this.dpoit.bind(b);
      b.renderer({ name: 'dpoitTransparent' }, scene.primitives, () =>
        renderer.renderDpoitTransparent(scene.primitives, camera, this.depthTextureOpaque, dpoitTextures),
      );

      // preparing: one peeling layer draws as all of them
      const layers = b.drawing ? iterations : Math.min(iterations, 1);
      for (let i = 0; i < layers; i++) {
        if (isTimingMode && b.drawing) this.graphics.timer.group('DpoitPass.layer');
        const dpoitTextures = this.dpoit.bindDualDepthPeeling(b);
        b.renderer({ name: 'dpoitTransparent' }, scene.primitives, () =>
          renderer.renderDpoitTransparent(scene.primitives, camera, this.depthTextureOpaque, dpoitTextures),
        );

        if (iterations > 1) {
          b.begin(target, { label: 'dpoit blend back' });
          this.dpoit.renderBlendBack(b);
        }
        if (isTimingMode && b.drawing) this.graphics.timer.groupEnd('DpoitPass.layer');
      }

      // evaluate dpoit
      b.begin(target, { label: 'dpoit evaluate' });
      this.dpoit.render(b);
    }

    if (PostprocessingPass.isEnabled(scene, postprocessingProps)) {
      this.renderPostprocessing(b, renderer, camera, scene, transparentBackground, postprocessingProps);
    }

    // without transparent primitives, volumes over the postprocessed image
    if (hasVolumes && !hasTransparent) {
      const target = isPostprocessingEnabled ? this.postprocessing.target : this.colorTarget;
      this.renderVolumes(b, renderer, camera, scene, target, target === this.colorTarget);
    }
  }

  private _renderWboit(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    transparentBackground: boolean,
    postprocessingProps: PostprocessingProps,
  ) {
    if (!this.wboit.supported) throw new Error('expected wboit to be supported');

    this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
    b.begin(this.colorTarget, { label: 'opaque' });
    if (b.drawing) renderer.clear(true);

    // render opaque primitives
    if (scene.hasOpaque) {
      b.renderer({ name: 'opaque' }, scene.primitives, () => renderer.renderOpaque(scene.primitives, camera));
    }

    if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
      b.begin(this.depthTargetTransparent, { label: 'transparent depth' });
      if (b.drawing) renderer.clearDepth(true);
      if (scene.opacityAverage < 1) {
        b.renderer({ name: 'depthTransparent' }, scene.primitives, () =>
          renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque),
        );
      }
    }

    // render transparent primitives
    const isPostprocessingEnabled = PostprocessingPass.isEnabled(scene, postprocessingProps);
    const hasVolumes = scene.volumes.renderables.length > 0;
    const hasTransparent = scene.opacityAverage < 1;
    if (hasTransparent) {
      const target = isPostprocessingEnabled ? this.transparentColorTarget : this.colorTarget;
      if (isPostprocessingEnabled && b.drawing) {
        b.begin(target, { label: 'transparent clear' });
        renderer.clear(false, false, true);
      }

      // volumes under the transparent primitives (one layer, not split at their surfaces)
      if (hasVolumes) this.renderVolumes(b, renderer, camera, scene, target, target === this.colorTarget);

      this.wboit.bind(b);
      b.renderer({ name: 'wboitTransparent' }, scene.primitives, () =>
        renderer.renderWboitTransparent(scene.primitives, camera, this.depthTextureOpaque),
      );

      // evaluate wboit
      b.begin(target, { label: 'wboit evaluate' });
      this.wboit.render(b);
    }

    if (PostprocessingPass.isEnabled(scene, postprocessingProps)) {
      this.renderPostprocessing(b, renderer, camera, scene, transparentBackground, postprocessingProps);
    }

    // without transparent primitives, volumes over the postprocessed image
    if (hasVolumes && !hasTransparent) {
      const target = isPostprocessingEnabled ? this.postprocessing.target : this.colorTarget;
      this.renderVolumes(b, renderer, camera, scene, target, target === this.colorTarget);
    }
  }

  /** Draws the volumes onto `target` of the blended mode, its opaque depth (texture or renderbuffer) detached meanwhile */
  private renderBlendedVolumes(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    target: RenderTarget,
  ) {
    if (!this.packedDepth) {
      this.depthTextureOpaque.detachFramebuffer(target.framebuffer, 'depth');
    } else {
      this.colorTarget.depthRenderbuffer?.detachFramebuffer(target.framebuffer);
    }
    this.renderVolumes(b, renderer, camera, scene, target, false);
    if (!this.packedDepth) {
      this.depthTextureOpaque.attachFramebuffer(target.framebuffer, 'depth');
    } else {
      this.colorTarget.depthRenderbuffer?.attachFramebuffer(target.framebuffer);
    }
    b.begin(target, { label: 'volume' });
  }

  private _renderBlended(
    b: PassBuilder,
    renderer: Renderer,
    camera: ICamera,
    scene: Scene,
    toDrawingBuffer: boolean,
    transparentBackground: boolean,
    postprocessingProps: PostprocessingProps,
  ) {
    if (toDrawingBuffer) {
      b.begin(this.drawTarget, { label: 'opaque' });
    } else {
      if (!this.packedDepth) {
        this.depthTextureOpaque.attachFramebuffer(this.colorTarget.framebuffer, 'depth');
      }
      b.begin(this.colorTarget, { label: 'opaque' });
    }

    if (b.drawing) renderer.clear(true);
    if (scene.hasOpaque) {
      b.renderer({ name: 'opaque' }, scene.primitives, () => renderer.renderOpaque(scene.primitives, camera));
    }

    if (!toDrawingBuffer) {
      // do a depth pass if not rendering to drawing buffer and
      // depthTexture is unsupported (i.e. depthTarget is set)
      if (this.depthTargetOpaque) {
        b.begin(this.depthTargetOpaque, { label: 'opaque depth' });
        if (b.drawing) renderer.clearDepth(true);
        b.renderer({ name: 'depthOpaque' }, scene.primitives, () =>
          renderer.renderDepthOpaque(scene.primitives, camera),
        );
        b.begin(this.colorTarget, { label: 'transparent' });
      }

      if (PostprocessingPass.isTransparentDepthRequired(scene, postprocessingProps)) {
        b.begin(this.depthTargetTransparent, { label: 'transparent depth' });
        if (b.drawing) renderer.clearDepth(true);
        if (scene.opacityAverage < 1) {
          b.renderer({ name: 'depthTransparent' }, scene.primitives, () =>
            renderer.renderDepthTransparent(scene.primitives, camera, this.depthTextureOpaque),
          );
        }
      }

      // render transparent primitives
      const isPostprocessingEnabled = PostprocessingPass.isEnabled(scene, postprocessingProps);
      const hasVolumes = scene.volumes.renderables.length > 0;
      const hasTransparent = scene.opacityAverage < 1;
      if (hasTransparent) {
        // volumes under the transparent primitives (one layer, not split at their surfaces)
        if (isPostprocessingEnabled) {
          b.begin(this.transparentColorTarget, { label: 'transparent' });
          if (b.drawing) renderer.clear(false, false, true);
          if (hasVolumes) this.renderVolumes(b, renderer, camera, scene, this.transparentColorTarget, false);

          if (!this.packedDepth) {
            this.depthTextureOpaque.attachFramebuffer(this.transparentColorTarget.framebuffer, 'depth');
          } else {
            this.colorTarget.depthRenderbuffer?.attachFramebuffer(this.transparentColorTarget.framebuffer);
          }
          b.begin(this.transparentColorTarget, { label: 'transparent' });
        } else if (hasVolumes) {
          this.renderBlendedVolumes(b, renderer, camera, scene, this.colorTarget);
          b.begin(this.colorTarget, { label: 'transparent' });
        }

        b.renderer({ name: 'blendedTransparent' }, scene.primitives, () =>
          renderer.renderBlendedTransparent(scene.primitives, camera),
        );

        if (isPostprocessingEnabled) {
          if (!this.packedDepth) {
            this.depthTextureOpaque.detachFramebuffer(this.transparentColorTarget.framebuffer, 'depth');
          } else {
            this.colorTarget.depthRenderbuffer?.detachFramebuffer(this.transparentColorTarget.framebuffer);
          }
        }
      }

      if (isPostprocessingEnabled) {
        if (!this.packedDepth) {
          this.depthTextureOpaque.detachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
        } else {
          this.colorTarget.depthRenderbuffer?.detachFramebuffer(this.postprocessing.target.framebuffer);
        }

        this.renderPostprocessing(b, renderer, camera, scene, transparentBackground, postprocessingProps);

        if (!this.packedDepth) {
          this.depthTextureOpaque.attachFramebuffer(this.postprocessing.target.framebuffer, 'depth');
        } else {
          this.colorTarget.depthRenderbuffer?.attachFramebuffer(this.postprocessing.target.framebuffer);
        }
      }

      // without transparent primitives, volumes over the postprocessed image
      if (hasVolumes && !hasTransparent) {
        this.renderBlendedVolumes(
          b,
          renderer,
          camera,
          scene,
          isPostprocessingEnabled ? this.postprocessing.target : this.colorTarget,
        );
      }
    } else if (scene.opacityAverage < 1) {
      b.renderer({ name: 'blendedTransparent' }, scene.primitives, () =>
        renderer.renderBlendedTransparent(scene.primitives, camera),
      );
    }
  }

  private planEye(
    b: PassBuilder,
    ctx: RenderContext<ICamera>,
    toDrawingBuffer: boolean,
    transparentBackground: boolean,
    props: Props,
    skipMarking: boolean,
    asIfMarked: boolean,
  ) {
    const { renderer, camera, scene, frame } = ctx;
    if (camera.disabled) return;
    this.volumeScale = props.volumeResolutionScale;

    const volumeRendering = scene.volumes.renderables.length > 0;
    const postprocessingEnabled = PostprocessingPass.isEnabled(scene, props.postprocessing);
    const antialiasingEnabled = AntialiasingPass.isEnabled(props.postprocessing);
    const dofEnabled = DofPass.isEnabled(props.postprocessing);
    const hasMarking = !skipMarking && (asIfMarked || MarkingPass.hasMarking(scene, props));
    // marking is composed onto the finished image, which therefore has to be readable
    const offscreen = !toDrawingBuffer || hasMarking;

    if (b.drawing) {
      const { x, y, width, height } = camera.viewport;
      renderer.setViewport(x, y, width, height);
      renderer.update(camera, scene, frame);
      // no pass open yet: GPU culling for the frame's passes with this camera (others draw the CPU culling's lists)
      renderer.prepareCull(scene.primitives);
    }

    if (transparentBackground && !antialiasingEnabled && !dofEnabled && !offscreen && !postprocessingEnabled) {
      b.begin(this.drawTarget, { label: 'background clear' });
      if (b.drawing) renderer.clear(false);
    }

    let oitEnabled = false;
    if (this.transparencyMode === 'wboit' && this.wboit.supported) {
      this._renderWboit(b, renderer, camera, scene, transparentBackground, props.postprocessing);
      oitEnabled = true;
    } else if (this.transparencyMode === 'dpoit' && this.dpoit.supported) {
      this._renderDpoit(b, renderer, camera, scene, props.dpoitIterations, transparentBackground, props.postprocessing);
      oitEnabled = true;
    } else {
      this._renderBlended(
        b,
        renderer,
        camera,
        scene,
        !volumeRendering && !postprocessingEnabled && !antialiasingEnabled && !dofEnabled && !offscreen,
        transparentBackground,
        props.postprocessing,
      );
    }

    const target = postprocessingEnabled
      ? this.postprocessing.target
      : offscreen || volumeRendering || oitEnabled || dofEnabled
        ? this.colorTarget
        : this.drawTarget;

    b.begin(target, { label: 'helpers' });
    this._renderHelpers(b, ctx, target);

    const output = this._renderOutput(b, camera, scene, props, {
      offscreen,
      postprocessingEnabled,
      antialiasingEnabled,
      dofEnabled,
      volumeRendering,
      oitEnabled,
    });

    if (!output) {
      if (b.drawing) this.marking.invalidate();
    } else if (!skipMarking) {
      this.marking.present(b, ctx, props, {
        base: output,
        toDrawingBuffer,
        offsets: SingleSample,
        restart: true,
        samples: 1,
        shading: this.getMarkingShading(props.postprocessing),
        asIfMarked,
      });
    }

    if (b.drawing) this.graphics.flush();
  }

  private _renderHelpers(b: PassBuilder, ctx: RenderContext<ICamera>, target: RenderTarget) {
    const { renderer, camera, helper, frame } = ctx;
    if (helper.debug.isEnabled || helper.pointer.isEnabled) {
      if (!this.packedDepth) {
        this.depthTextureOpaque.attachFramebuffer(target.framebuffer, 'depth');
        b.begin(target, { label: 'helpers' });
      }
      if (helper.debug.isEnabled) {
        if (b.drawing) helper.debug.syncVisibility();
        for (const scene of helper.debug.scenes) {
          this.renderBlended(b, renderer, scene, camera);
        }
      }
      if (helper.pointer.isEnabled) {
        if (b.drawing) {
          helper.pointer.setCamera(camera);
          renderer.update(helper.pointer.camera, helper.pointer.scene, frame);
        }
        this.renderBlended(b, renderer, helper.pointer.scene, helper.pointer.camera);
      }
      if (!this.packedDepth) {
        this.depthTextureOpaque.detachFramebuffer(target.framebuffer, 'depth');
        b.begin(target, { label: 'helpers' });
      }
    }
    if (helper.handle.isEnabled) {
      this.renderBlended(b, renderer, helper.handle.scene, camera);
    }
    if (helper.camera.isEnabled) {
      if (b.drawing) {
        helper.camera.update(camera);
        renderer.update(helper.camera.camera, helper.camera.scene, frame);
      }
      this.renderBlended(b, renderer, helper.camera.scene, helper.camera.camera);
    }
  }

  /** `renderer.renderBlended` as renderer passes */
  renderBlended(b: PassBuilder, renderer: Renderer, scene: Scene, camera: ICamera) {
    if (b.drawing) {
      renderer.renderBlended(scene, camera);
      return;
    }
    if (scene.hasOpaque) b.renderer({ name: 'opaque' }, scene, () => {});
    if (scene.opacityAverage < 1) b.renderer({ name: 'blendedTransparent' }, scene, () => {});
  }

  /** Returns the target the finished image ended up in, or undefined when it went to the drawing buffer. */
  private _renderOutput(
    b: PassBuilder,
    camera: ICamera,
    scene: Scene,
    props: Props,
    flags: {
      offscreen: boolean;
      postprocessingEnabled: boolean;
      antialiasingEnabled: boolean;
      dofEnabled: boolean;
      volumeRendering: boolean;
      oitEnabled: boolean;
    },
  ): RenderTarget | undefined {
    const { offscreen, postprocessingEnabled, antialiasingEnabled, dofEnabled, volumeRendering, oitEnabled } = flags;
    const toBuffer = !offscreen;

    let needsTargetCopy = false;

    if (antialiasingEnabled) {
      const input = postprocessingEnabled ? this.postprocessing.target.texture : this.colorTarget.texture;
      this.antialiasing.render(b, camera, input, toBuffer && !dofEnabled, props.postprocessing);
    } else if (toBuffer && !dofEnabled) {
      needsTargetCopy = true;
    }

    if (dofEnabled && props.postprocessing.dof.name === 'on') {
      const input = antialiasingEnabled
        ? this.antialiasing.target.texture
        : postprocessingEnabled
          ? this.postprocessing.target.texture
          : this.colorTarget.texture;
      this.dof.update(
        camera,
        input,
        this.depthTargetOpaque?.texture || this.depthTextureOpaque,
        this.depthTextureTransparent,
        props.postprocessing.dof.params,
        scene.boundingSphereVisible,
      );
      this.dof.render(b, camera.viewport, toBuffer ? undefined : this.getColorTarget(scene, props.postprocessing));
    } else if (toBuffer && !antialiasingEnabled) {
      needsTargetCopy = true;
    }

    if (needsTargetCopy) {
      if (postprocessingEnabled) {
        this.copyToDrawingBuffer(b, this.postprocessing.target, camera.viewport);
      } else if (volumeRendering || oitEnabled) {
        this.copyToDrawingBuffer(b, this.colorTarget, camera.viewport);
      }
    }

    return offscreen ? this.getColorTarget(scene, props.postprocessing) : undefined;
  }

  private copyToDrawingBuffer(b: PassBuilder, src: RenderTarget, viewport: Viewport) {
    if (this.copyFbo.values.tColor.ref.value !== src.texture) {
      ValueCell.update(this.copyFbo.values.tColor, src.texture);
      this.copyFbo.update();
    }

    b.begin(this.drawTarget, { label: 'copy' });
    if (b.drawing) {
      const { x, y, width, height } = viewport;
      this.graphics.setViewport(x, y, width, height);
      this.graphics.setScissor(x, y, width, height);
    }
    b.quad(this.copyFbo, RenderState.Default);
  }

  render(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, skipMarking = false) {
    if (isTimingMode) this.graphics.timer.group('DrawPass.render');
    this.plan(this.renderBuilder, ctx, props, toDrawingBuffer, skipMarking);
    if (isTimingMode) this.graphics.timer.groupEnd('DrawPass.render');
  }

  /**
   * Makes the pipelines and resources of `render` with these arguments ready without drawing (see `PassBuilder`);
   * resolves when they are
   */
  prepare(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, skipMarking = false): Promise<void> {
    const b = createPrepareBuilder(this.graphics, ctx.renderer);
    this.plan(b, ctx, props, toDrawingBuffer, skipMarking);
    return b.whenReady();
  }

  /**
   * What `render` draws, see `PassBuilder`; `asIfMarked` draws marking as if objects were marked (preparing marking,
   * which draws only marked objects, before the first hover)
   */
  plan(
    b: PassBuilder,
    ctx: RenderContext,
    props: Props,
    toDrawingBuffer: boolean,
    skipMarking = false,
    asIfMarked = false,
  ) {
    const { renderer, camera, scene } = ctx;

    this.postprocessing.setTransparentBackground(props.transparentBackground);
    const transparentBackground = this.isTransparentBackground(scene, props);

    renderer.setTransparentBackground(transparentBackground);
    if (b.drawing) {
      renderer.setDrawingBufferSize(this.colorTarget.getWidth(), this.colorTarget.getHeight());
      renderer.setPixelRatio(this.graphics.pixelRatio);
    }

    if (StereoCamera.is(camera)) {
      const timing = isTimingMode && b.drawing;
      if (timing) this.graphics.timer.group('StereoCamera.left');
      this.planEye(
        b,
        { ...ctx, camera: camera.left },
        toDrawingBuffer,
        transparentBackground,
        props,
        skipMarking,
        asIfMarked,
      );
      if (timing) this.graphics.timer.groupEnd('StereoCamera.left');
      if (!b.drawing) return; // the right eye draws as the left one
      if (timing) this.graphics.timer.group('StereoCamera.right');
      this.planEye(
        b,
        { ...ctx, camera: camera.right },
        toDrawingBuffer,
        transparentBackground,
        props,
        skipMarking,
        asIfMarked,
      );
      if (timing) this.graphics.timer.groupEnd('StereoCamera.right');
      // only the second eye is left in the target, which is not a valid recompose source
      this.marking.invalidate();
    } else {
      this.planEye(b, { ...ctx, camera }, toDrawingBuffer, transparentBackground, props, skipMarking, asIfMarked);
    }
  }

  private isTransparentBackground(scene: Scene, props: Props) {
    const pp = props.postprocessing;
    const backgroundEnabled = this.postprocessing.background.isEnabled(pp);
    // premultiplied scene so bloom can composite the solid background last
    const bloomCompositesBackground =
      !props.transparentBackground && !backgroundEnabled && BloomPass.isActive(scene, pp);
    return props.transparentBackground || backgroundEnabled || bloomCompositesBackground;
  }

  /** shading of the last render that marking tint and dim keep, or null if there is none */
  getMarkingShading(postprocessingProps: PostprocessingProps): MarkingShading | null {
    return SsaoPass.isEnabled(postprocessingProps) ? this.ssaoShading : null;
  }

  getColorTarget(scene: Scene, postprocessingProps: PostprocessingProps): RenderTarget {
    if (DofPass.isEnabled(postprocessingProps)) {
      return this.dof.target;
    } else if (AntialiasingPass.isEnabled(postprocessingProps)) {
      return this.antialiasing.target;
    } else if (PostprocessingPass.isEnabled(scene, postprocessingProps)) {
      return this.postprocessing.target;
    }
    return this.colorTarget;
  }
}
