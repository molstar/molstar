/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isTimingMode } from '@molstar/core/util/debug';
import { type OccluderGrid, getOccluderGrid } from '@molstar/graphics/gl/occluder-grid';
import type { Renderer } from '@molstar/graphics/gl/renderer';
import type { Camera, ICamera } from './camera.js';
import type { StereoCamera } from './camera/stereo.js';
import type { Canvas3DParts } from './parts.js';
import { getJitterOffsets, getTemporalSamplesPerFrame } from './passes/jitter.js';
import { isMaterialColorMarker, SingleSample } from './passes/marking.js';
import { MultiSamplePass } from './passes/multi-sample.js';
import type { PassBuilder } from './passes/pass-builder.js';
import { isStereoPyramidOf } from './passes/hi-z.js';

/** What the frame render reads besides the canvas's parts */
export interface Canvas3DFrameRenderHost {
  readonly parts: Canvas3DParts;
  /** the GPU occlusion pyramid of the last frame, for the renderer's culling */
  readonly getOcclusionPyramid: NonNullable<Parameters<Renderer['setOcclusionPyramid']>[0]>;
  /** the user is interacting (illumination holds its iterations meanwhile) */
  isActivelyInteracting(): boolean;
  /** the controls animate the camera (no flicker reduction meanwhile) */
  isControlsAnimating(): boolean;
  /** the marks changed a helper's (no marking-only redraw over the last image) */
  isMarkingHelperUpdated(): boolean;
  /** another frame drawn at the next tick, also where nothing changed */
  requestFrame(): void;
}

/**
 * Draws a frame of the canvas (`Canvas3D` render loop): a marking-only redraw over the last image where that is all
 * that changed, else illumination or the draw pass (multi-sampled where on), the Hi-Z levels of its depth, and more
 * marking samples while the scene rests; returns whether it drew. When to draw is the canvas's. What changed, as the
 * render loop found it (arguments, not an object: a frame each tick): `shouldRender` the scene is drawn again (forced,
 * the camera, a resize, an XR view or an animation changed), `xrChanged` the XR views (drawn with the stereo camera),
 * `markingUpdated` the marks changed what marking draws
 */
export function createCanvas3DFrameRender(host: Canvas3DFrameRenderHost) {
  const { parts, getOcclusionPyramid } = host;
  /**
   * the pyramid of the frame's views: of the main camera alone, or of stereo's eyes both (`OcclusionPyramid.second`);
   * after a switch between them none, as the next frame builds one
   */
  const getMonoPyramid = () => {
    const pyramid = getOcclusionPyramid();
    return pyramid?.second ? undefined : pyramid;
  };
  const getStereoPyramid = () => {
    const pyramid = getOcclusionPyramid();
    return pyramid?.second ? pyramid : undefined;
  };
  /**
   * whether stereo's frame is culled against levels of other views (the eyes moved): without a second phase what came
   * into view is missing until a frame against levels of these views draws it (`requestFrame`)
   */
  const isStereoPyramidStale = (stereo: StereoCamera) => {
    const pyramid = getStereoPyramid();
    return !!pyramid && !isStereoPyramidOf(pyramid, stereo);
  };
  /** the GPU culling's second phase within the frame (`hiZ.twoPhase`), made once: frames draw every tick */
  const planOcclusion = (b: PassBuilder, c: ICamera, complete: boolean) => parts.hiZ.planFrame(b, c, complete);

  /** Marking drawn again over the last image (`restart`: from its first sample), as the last full render culled */
  function redrawMarking(
    illuminationEnabled: boolean,
    markingOffsets: ReturnType<typeof getJitterOffsets>,
    restart: boolean,
    samples: number,
  ) {
    const { graphics, renderer, passes, hiZ, scene, helper, camera, frame, props } = parts;
    // other passes (e.g. pick) reset this, it has to match the last full render
    renderer.setOcclusionTest(illuminationEnabled ? null : hiZ.occlusionTest);
    renderer.setOcclusionPyramid(illuminationEnabled ? null : getMonoPyramid);
    if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
    const redrawn = passes.draw.marking.redraw(
      { renderer, camera, scene, helper, frame },
      props,
      markingOffsets,
      restart,
      samples,
    );
    if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');
    return redrawn;
  }

  return function renderFrame(
    shouldRender: boolean,
    cameraChanged: boolean,
    xrChanged: boolean,
    xrFrame: XRFrame | undefined,
    markingUpdated: boolean,
  ): boolean {
    const { graphics, renderer, passes, hiZ, multiSampleHelper, pickHelper, scene, helper, camera, frame } = parts;
    // what a visible object's content hides (`OccluderGrid`): the depth pyramid primed from it
    let occluder: OccluderGrid | undefined;
    for (const r of scene.renderables) {
      if (!r.state.visible) continue;
      occluder = getOccluderGrid(r.state);
      if (occluder) break;
    }
    hiZ.setOccluder(occluder);
    const { stereoCamera, props: p } = parts;
    let didRender = false;

    const illuminationEnabled = passes.illumination.supported && p.illumination.enabled && !xrFrame;
    const multiSampleEnabled = MultiSamplePass.isEnabled(p.multiSample) && !xrFrame;

    const markingOffsets = (illuminationEnabled ? p.multiSample.mode === 'on' : multiSampleEnabled)
      ? getJitterOffsets(p.multiSample.sampleLevel)
      : SingleSample;

    // marking-only changes are redrawn over the image of the last render, avoiding a full
    // re-render; the per-object color marker is part of that image, so it needs a re-render
    const canRedrawMarking =
      !host.isMarkingHelperUpdated() &&
      !shouldRender &&
      !isMaterialColorMarker(renderer.props, p.marking) &&
      !xrFrame &&
      p.camera.stereo.name !== 'on' &&
      (!illuminationEnabled || passes.illumination.iteration > 0);

    if (canRedrawMarking && markingUpdated) {
      // illumination multi-samples are expensive (full depth pass each), so add them over frames
      const progressive = illuminationEnabled || p.multiSample.mode === 'temporal';
      if (redrawMarking(illuminationEnabled, markingOffsets, true, progressive ? 1 : markingOffsets.length))
        return true;
    }

    if (illuminationEnabled) {
      if (shouldRender || markingUpdated) passes.illumination.restart();

      if (
        passes.illumination.shouldRender(p.illumination) &&
        ((!host.isActivelyInteracting() && scene.count > 0) ||
          passes.illumination.iteration === 0 ||
          p.userInteractionReleaseMs === 0)
      ) {
        // the inputs (tracing, back depth) culled by occlusion only on the GPU with its second phase: the image
        // accumulates them, what an earlier frame's pyramid hides and the view reveals would stay missing; so no CPU
        // test, which lags behind (WebGL: none). Set per frame, other passes (e.g. pick) reset them
        const gpuOcclusion = hiZ.props.gpu && hiZ.props.twoPhase;
        renderer.setOcclusionTest(null);
        renderer.setOcclusionPyramid(gpuOcclusion ? getMonoPyramid : null);
        const occlusion = gpuOcclusion ? planOcclusion : undefined;
        if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
        const ctx = { renderer, camera, scene, helper, frame, occlusion };
        // the first iteration draws the inputs: the tracing input's depth is the opaque depth the levels are built from
        const inputDrawn = passes.illumination.iteration === 0;
        const skippedDraws = graphics.stats.skippedDraws;
        passes.illumination.render(ctx, p, true);
        if (inputDrawn && gpuOcclusion && renderer.gpuCulledCount > 0) {
          // for the next frame's first phase, not read back (no CPU test)
          hiZ.setGpuCulling(true);
          hiZ.render(camera, graphics.stats.skippedDraws === skippedDraws, false);
        }
        if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');

        // if only marking has updated, do not set the flag to dirty
        pickHelper.dirty = pickHelper.dirty || shouldRender;
        didRender = true;
      }
    } else {
      const multiSampleChanged = multiSampleHelper.update(markingUpdated || shouldRender, p.multiSample) && !xrFrame;

      if (shouldRender || multiSampleChanged || markingUpdated) {
        renderer.setOcclusionTest(hiZ.occlusionTest);
        renderer.setOcclusionPyramid(getMonoPyramid);

        let cam: Camera | StereoCamera = camera;
        if (p.camera.stereo.name === 'on' || xrChanged) {
          if (!xrChanged) stereoCamera.update();
          cam = stereoCamera;
          // the CPU test reads the main camera's levels, not an eye's; the GPU's tests both eyes' (`hiZ.render`): an
          // instance occluded only where occluded in both, against the last frame's (no second phase: what comes into
          // view shows a frame late)
          renderer.setOcclusionTest(null);
          renderer.setOcclusionPyramid(getStereoPyramid);
          // before the frame: its levels built after it reuse the pyramid's matrices
          if (isStereoPyramidStale(stereoCamera)) host.requestFrame();
        }

        if (isTimingMode) graphics.timer.group('Canvas3D.render', { captureStats: true });
        // GPU culling catching up within the frame (`hiZ.twoPhase`)
        const occlusion = cam === camera && !xrFrame && hiZ.props.gpu && hiZ.props.twoPhase ? planOcclusion : undefined;
        const ctx = { renderer, camera: cam, scene, helper, frame, occlusion };
        // the levels need the frame's depth complete, not every pending work done (a pick's readback would skip them)
        const skippedDraws = graphics.stats.skippedDraws;
        if (multiSampleEnabled) {
          const forceOn =
            p.multiSample.reduceFlicker &&
            isMaterialColorMarker(renderer.props, p.marking) &&
            !cameraChanged &&
            markingUpdated &&
            !host.isControlsAnimating();
          multiSampleHelper.render(ctx, p, true, forceOn);
        } else {
          passes.draw.render(ctx, p, true);
        }
        // the GPU's occlusion test of the objects culled there needs the levels (also with Hi-Z off)
        hiZ.setGpuCulling(renderer.gpuCulledCount > 0);
        // stereo's eyes: the levels of both (side by side; XR's views where laid out so); none of a frame drawn straight
        // to the drawing buffer: the texture the levels read holds an earlier frame's
        if (cam === camera && !xrFrame && passes.draw.opaqueDepthDrawn) {
          hiZ.render(camera, graphics.stats.skippedDraws === skippedDraws);
        } else if (cam === stereoCamera && passes.draw.opaqueDepthDrawn) {
          hiZ.render(stereoCamera, graphics.stats.skippedDraws === skippedDraws, false);
        } else hiZ.clear();
        if (isTimingMode) graphics.timer.groupEnd('Canvas3D.render');

        // if only marking has updated, do not set the flag to dirty
        pickHelper.dirty = pickHelper.dirty || shouldRender;
        didRender = true;
      }
    }

    // keep accumulating marking samples after the scene has settled
    if (!didRender && !markingUpdated && canRedrawMarking && passes.draw.marking.needsMoreSamples) {
      didRender = redrawMarking(
        illuminationEnabled,
        markingOffsets,
        false,
        getTemporalSamplesPerFrame(p.multiSample.sampleLevel),
      );
    }
    return didRender;
  };
}
