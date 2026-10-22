/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MediabunnyBackend } from './backends/mediabunny.js';
import { resolveVideoParams, type VideoRendererBackend, type VideoResult } from './backend.js';
import { Viewport } from '@molstar/graphics/canvas3d/camera/util';
import { PluginStateAnimation } from '@molstar/plugin/state/animation/model';
import { PluginContext } from '@molstar/plugin/context';
import { RuntimeContext } from '@molstar/core/task';
import { Color } from '@molstar/core/util/color';

export interface VideoEncoderParams<A extends PluginStateAnimation = PluginStateAnimation> {
  pass: {
    updateBackground(): Promise<void>;
    getImageData(
      ctx: RuntimeContext,
      width: number,
      height: number,
      viewport?: Viewport,
    ): Promise<{ data: Uint8Array | Uint8ClampedArray }>;
  };
  customBackground?: Color;
  animation: PluginStateAnimation.Instance<A>;
  width: number;
  height: number;
  viewport: Viewport;
  /** default is 30 */
  fps?: number;
  /** Defaults to Mediabunny; explicitly import a backend for environments without WebCodecs. */
  backend?: VideoRendererBackend;
  backendParams?: Record<string, any>;
}

export async function encodeVideoAnimation<A extends PluginStateAnimation>(
  plugin: PluginContext,
  ctx: RuntimeContext,
  params: VideoEncoderParams<A>,
): Promise<VideoResult> {
  await ctx.update({ message: 'Initializing...', isIndeterminate: true });

  validateViewport(params);
  const durationMs = PluginStateAnimation.getDuration(plugin, params.animation);
  if (durationMs === void 0) {
    throw new Error('The animation does not have the duration specified.');
  }

  const backend = params.backend ?? MediabunnyBackend;
  const fps = params.fps ?? 30;
  if (!Number.isFinite(fps) || fps <= 0) throw new Error('Frame rate must be positive.');
  if (!Number.isFinite(durationMs) || durationMs <= 0) throw new Error('Animation duration must be positive.');
  const { width, height } = params;
  const alignment = backend.dimensionAlignment;
  const vw = Math.floor(params.viewport.width / alignment) * alignment;
  const vh = Math.floor(params.viewport.height / alignment) * alignment;
  if (vw <= 0 || vh <= 0) throw new Error('Video dimensions must be positive.');
  const normalizedViewport: Viewport = { ...params.viewport, width: vw, height: vh };
  const config = { width: vw, height: vh, fps };
  const values = await resolveVideoParams(backend, config, params.backendParams);
  const encoder = await backend.create(config, values);

  const loop = plugin.animationLoop;
  const canvasProps = plugin.canvas3d?.props;
  const wasAnimating = loop.isAnimating;
  let stoppedAnimation = true;

  try {
    loop.stop();
    loop.resetTime(0);

    if (params.customBackground !== void 0) {
      plugin.canvas3d?.setProps(
        { renderer: { backgroundColor: params.customBackground }, transparentBackground: false },
        true,
      );
    }

    const N = Math.ceil((durationMs / 1000) * fps);
    const dt = 1000 / fps;

    await ctx.update({ message: 'Rendering...', isIndeterminate: false, current: 0, max: N });
    await params.pass.updateBackground();

    stoppedAnimation = false;
    await plugin.managers.animation.play(params.animation.definition, params.animation.params);
    for (let i = 0; i < N; i++) {
      await loop.tick(i * dt, {
        isSynchronous: true,
        animation: { currentFrame: i, frameCount: N - 1 },
        manualDraw: true,
        updateControls: true,
      });
      const image = await params.pass.getImageData(ctx, width, height, normalizedViewport);
      await encoder.addFrame(image.data, i / fps, Math.min(dt, durationMs - i * dt) / 1000);

      if (ctx.shouldUpdate) {
        await ctx.update({ current: i + 1 });
      }
    }
    await ctx.update({ message: 'Applying finishing touches...', isIndeterminate: true });
    await plugin.managers.animation.stop();
    stoppedAnimation = true;
    return await encoder.finalize();
  } finally {
    try {
      await encoder.dispose();
    } finally {
      if (params.customBackground !== void 0) {
        plugin.canvas3d?.setProps({
          renderer: { backgroundColor: canvasProps?.renderer!.backgroundColor },
          transparentBackground: canvasProps?.transparentBackground,
        });
      }
      try {
        if (!stoppedAnimation) await plugin.managers.animation.stop();
      } finally {
        if (wasAnimating) loop.start();
      }
    }
  }
}

function validateViewport(params: VideoEncoderParams) {
  if (
    [params.width, params.height, params.viewport.width, params.viewport.height].some(
      (v) => !Number.isInteger(v) || v <= 0,
    ) ||
    [params.viewport.x, params.viewport.y].some((v) => !Number.isInteger(v) || v < 0)
  ) {
    throw new Error('Invalid video dimensions or viewport.');
  }
  if (
    params.viewport.x + params.viewport.width > params.width ||
    params.viewport.y + params.viewport.height > params.height
  ) {
    throw new Error('Viewport exceeds the canvas dimensions.');
  }
}
