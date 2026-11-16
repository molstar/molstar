/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { ICamera } from '../camera.js';
import type { RenderContext } from '../util.js';
import type { PassBuilder } from './pass-builder.js';
import type { PassTargets } from './pass-targets.js';
import { createRenderPassAttachments, type RenderPassAttachments } from '@molstar/graphics/gl/render-pass-attachments';

/**
 * Draws the overlay layer (`OverlayHelper`) and then the screen layer (`ScreenHelper`) onto a finished image: after the
 * scene's post-processing (and illumination's composite), before antialiasing, with a depth of their own (cleared
 * before each layer, so a gizmo occludes itself). The overlay with its camera (near and far planes fitting its bounds,
 * no fog), the scene's depth applying per object by its occluded alpha (`dOverlay`); each screen object with its own
 * (`ScreenCamera`).
 */
export class OverlayPass {
  /** the composite's target: the caller's color (`colorTexture`) with the pass's own depth, where it has one */
  private readonly composite: RenderPassAttachments;
  private readonly depthRenderbuffer: Renderbuffer;
  /** the color of `composite`, the target's of the last draw */
  private colorTexture: Texture | undefined;

  /** `transparentDepth`: the scene's front transparent depth, written while an overlay is over transparent objects */
  constructor(
    targets: PassTargets,
    width: number,
    height: number,
    private readonly transparentDepth: Texture,
  ) {
    this.depthRenderbuffer = targets.renderbuffer('depth32f', 'depth', width, height);
    this.composite = createRenderPassAttachments(targets.graphics.resources, [], this.depthRenderbuffer);
  }

  getByteCount() {
    return this.depthRenderbuffer.getByteCount();
  }

  setSize(width: number, height: number) {
    this.depthRenderbuffer.setSize(width, height);
  }

  /**
   * Draws the overlay and the screen layer onto `target`, the overlay over the scene's opaque depth `depth` (drawn
   * with `ctx.camera`) and transparent depth, which its occlusion tests against (`Renderer.renderOverlay`) and its picking after
   * (`OverlayHelper.occluderDepth`); `targetDepth`: with the target's own depth buffer instead of the pass's (the
   * drawing buffer, which takes no other attachment and whose depth nothing reads after; only without an overlay,
   * whose occlusion needs the opaque depth that drawing straight to the drawing buffer doesn't write)
   */
  render(b: PassBuilder, ctx: RenderContext<ICamera>, target: RenderTarget, depth: Texture, targetDepth = false) {
    const { overlay, screen } = ctx.helper;
    const overlayEnabled = overlay.isEnabled;
    const screenEnabled = screen.isEnabled;
    if (!overlayEnabled && !screenEnabled) return;

    if (targetDepth) {
      b.begin(target, { label: 'overlay' });
    } else {
      if (this.colorTexture !== target.texture) {
        this.composite.color = [target.texture];
        this.colorTexture = target.texture;
      }
      b.begin(this.composite, { label: 'overlay' });
    }

    if (overlayEnabled) this.renderOverlay(b, ctx, depth);
    if (screenEnabled) this.renderScreen(b, ctx);
  }

  private renderOverlay(b: PassBuilder, ctx: RenderContext<ICamera>, depth: Texture) {
    const { renderer, camera, helper, frame } = ctx;
    const { overlay } = helper;
    const { scene } = overlay;
    if (b.drawing) {
      const depthAlpha = ctx.scene.opacityAverage < 1 ? this.transparentDepth : undefined;
      overlay.occluderDepth = depth;
      overlay.occluderDepthAlpha = depthAlpha;
      renderer.clearDepth();
      const overlayCamera = overlay.getCamera(camera);
      renderer.update(overlayCamera, scene, frame);
      // what the scene's depth hides is drawn
      const occlusionTest = renderer.occlusionTest;
      renderer.setOcclusionTest(null);
      renderer.renderOverlay(scene, overlayCamera, { depth, depthAlpha, near: camera.near, far: camera.far });
      renderer.setOcclusionTest(occlusionTest);
    } else {
      b.renderer({ name: 'opaque' }, scene, () => {});
      // also the occluded parts of opaque objects
      b.renderer({ name: 'blendedTransparent' }, scene, () => {});
    }
  }

  private renderScreen(b: PassBuilder, ctx: RenderContext<ICamera>) {
    const { renderer, camera, helper, frame } = ctx;
    if (b.drawing) renderer.clearDepth();
    const { screen } = helper;
    // a layer drawn larger than the canvas (`ScreenHelper.pixelScale`): its line widths and point sizes by its pixels
    const scaled = b.drawing && screen.pixelScale !== 1;
    if (scaled) renderer.setPixelRatio(screen.pixelRatio);
    screen.forEach(camera, (scene, screenCamera) => {
      if (b.drawing) {
        renderer.update(screenCamera, scene, frame);
        // the occlusion test is of the scene's view
        const occlusionTest = renderer.occlusionTest;
        renderer.setOcclusionTest(null);
        renderer.renderBlended(scene, screenCamera);
        renderer.setOcclusionTest(occlusionTest);
      } else {
        if (scene.hasOpaque) b.renderer({ name: 'opaque' }, scene, () => {});
        if (scene.opacityAverage < 1) b.renderer({ name: 'blendedTransparent' }, scene, () => {});
      }
    });
    if (scaled) renderer.setPixelRatio(screen.pixelRatio / screen.pixelScale);
  }
}
