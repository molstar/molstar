/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { writeFile } from 'node:fs/promises';
import type { VideoRendererBackend } from './backend.js';
import { VideoExport } from './index.js';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { encodeVideoAnimation } from '@molstar/video-export-extension/encoder';
import type { PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { AnimateStateSnapshots } from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { RuntimeContext, Task } from '@molstar/core/task';

export interface VideoAnimationOptions {
  backend?: VideoRendererBackend;
  backendParams?: Record<string, any>;
  size?: { width: number; height: number };
  fps?: number;
  postprocessing?: Partial<PostprocessingProps>;
}

/** Explicit video integration for headless embedding and the rendering CLI. */
export class VideoHeadlessPluginContext extends HeadlessPluginContext {
  /** Render plugin state snapshots animation and return as encoded video and its MIME type/extension */
  async getAnimation(options?: VideoAnimationOptions) {
    if (!this.state.hasBehavior(VideoExport)) {
      throw new Error('PluginContext must have VideoExport extension registered in order to save animation.');
    }

    const task = Task.create('Export Animation', async (ctx) => {
      const { width, height } = options?.size ?? this.renderer.canvasSize;
      const movie = await encodeVideoAnimation(this, ctx, {
        animation: { definition: AnimateStateSnapshots, params: {} },
        width,
        height,
        viewport: { x: 0, y: 0, width, height },
        backend: options?.backend,
        backendParams: options?.backendParams,
        fps: options?.fps,
        pass: {
          getImageData: (
            runtime: RuntimeContext,
            _width: number,
            _height: number,
            viewport?: { width: number; height: number },
          ) => this.renderer.getImageRaw(runtime, viewport ?? { width, height }, options?.postprocessing),
          updateBackground: () => this.renderer.imagePass.updateBackground(),
        },
      });
      return movie;
    });
    return this.runTask(task, { useOverlay: true });
  }

  /** Render plugin state snapshots animation and save to a video file. */
  async saveAnimation(outPath: string, options?: VideoAnimationOptions) {
    const movie = await this.getAnimation(options);
    if (!outPath.toLowerCase().endsWith(`.${movie.extension}`)) {
      throw new Error(`Output filename must end with .${movie.extension}.`);
    }
    await writeFile(outPath, movie.data);
  }
}
