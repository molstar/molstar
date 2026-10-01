import fs from 'fs';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { encodeMp4Animation } from '@molstar/mp4-export-extension/encoder';
import { ImagePass } from '@molstar/graphics/canvas3d/passes/image';
import type { PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { AnimateStateSnapshots } from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { RuntimeContext, Task } from '@molstar/core/task';

/** Explicit MP4 integration for headless embedding and the rendering CLI. */
export class Mp4HeadlessPluginContext extends HeadlessPluginContext {
    /** Render plugin state snapshots animation and return as raw MP4 data */
    async getAnimation(options?: { quantization?: number, size?: { width: number, height: number }, fps?: number, postprocessing?: Partial<PostprocessingProps> }) {
        if (!this.state.hasBehavior('extension-mp4-export')) {
            throw new Error('PluginContext must have Mp4Export extension registered in order to save animation.');
        }

        const task = Task.create('Export Animation', async ctx => {
            const { width, height } = options?.size ?? this.renderer.canvasSize;
            const movie = await encodeMp4Animation(this, ctx, {
                animation: { definition: AnimateStateSnapshots, params: {} },
                width,
                height,
                viewport: { x: 0, y: 0, width, height },
                quantizationParameter: options?.quantization ?? 18,
                fps: options?.fps,
                pass: {
                    getImageData: (runtime: RuntimeContext, width: number, height: number) => this.renderer.getImageRaw(runtime, { width, height }, options?.postprocessing),
                    updateBackground: () => this.renderer.imagePass.updateBackground(),
                } as ImagePass,
            });
            return movie;
        });
        return this.runTask(task, { useOverlay: true });
    }

    /** Render plugin state snapshots animation and save to a MP4 file */
    async saveAnimation(outPath: string, options?: { quantization?: number, size?: { width: number, height: number }, fps?: number, postprocessing?: Partial<PostprocessingProps> }) {
        const movie = await this.getAnimation(options);
        await new Promise<void>(resolve => {
            fs.writeFile(outPath, movie, () => resolve());
        });
    }
}
