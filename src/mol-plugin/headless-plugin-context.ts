/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import fs from 'fs';
import { type BufferRet as JpegBufferRet } from 'jpeg-js'; // Only import type here, the actual import must be provided by the caller
import { type PNG } from 'pngjs'; // Only import type here, the actual import must be provided by the caller

import { encodeMp4Animation } from '../extensions/mp4-export/encoder';
import { ImagePass } from '../mol-canvas3d/passes/image';
import { Viewport } from '../mol-canvas3d/camera/util';
import { PostprocessingProps } from '../mol-canvas3d/passes/postprocessing';
import { AnimateStateSnapshots } from '../mol-plugin-state/animation/built-in/state-snapshots';
import { RuntimeContext, Task } from '../mol-task';
import { PluginContext } from './context';
import { PluginSpec } from './spec';
import { ExternalModules, HeadlessScreenshotHelper, HeadlessScreenshotHelperOptions, RawImageData } from './util/headless-screenshot';


/** PluginContext that can be used in Node.js (without DOM) */
export class HeadlessPluginContext extends PluginContext {
    renderer: HeadlessScreenshotHelper;

    /** Legacy synchronous constructor accepts `gl`. Use create() with a caller-provided WebGPU provider for native rendering. */
    constructor(externalModules: ExternalModules, spec: PluginSpec, canvasSize: { width: number, height: number } = { width: 640, height: 480 }, rendererOptions?: HeadlessScreenshotHelperOptions, renderer?: HeadlessScreenshotHelper) {
        super(spec);
        this.renderer = renderer ?? new HeadlessScreenshotHelper(externalModules, canvasSize, undefined, rendererOptions);
        this.setCanvas3D(this.renderer.canvas3d, this.renderer.context);
        const create = async () => {
            const previous = this.renderer;
            this.renderer = await HeadlessScreenshotHelper.create(externalModules, canvasSize, { ...rendererOptions, canvas: previous.canvas3d.props, imagePass: previous.imagePass.props });
            this.renderer.context?.setProps({ transparency: previous.context?.props.transparency });
            this.setCanvas3D(this.renderer.canvas3d, this.renderer.context);
            await previous.dispose();
            this.configureWebGPURecovery(create);
        };
        if (this.renderer.context?.webgpu) this.configureWebGPURecovery(create);
    }

    /** Creates a headless plugin with native WebGPU rendering by default. Call init() afterwards. */
    static async create(externalModules: ExternalModules, spec: PluginSpec, canvasSize: { width: number, height: number } = { width: 640, height: 480 }, rendererOptions?: HeadlessScreenshotHelperOptions) {
        const renderer = await HeadlessScreenshotHelper.create(externalModules, canvasSize, rendererOptions);
        try { return new HeadlessPluginContext(externalModules, spec, canvasSize, rendererOptions, renderer); } catch (error) { await renderer.dispose(); throw error; }
    }

    /** Render the current plugin state and save to a PNG or JPEG file */
    async saveImage(outPath: string, imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>, format?: 'png' | 'jpeg', jpegQuality = 90) {
        await this.webgpuRecovery;
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.saveImage(ctx, outPath, imageSize, props, format, jpegQuality);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as raw image data */
    async getImageRaw(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>): Promise<RawImageData> {
        await this.webgpuRecovery;
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.getImageRaw(ctx, imageSize, props);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as a PNG object */
    async getImagePng(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>): Promise<PNG> {
        await this.webgpuRecovery;
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.getImagePng(ctx, imageSize, props);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as a JPEG object */
    async getImageJpeg(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>, jpegQuality: number = 90): Promise<JpegBufferRet> {
        await this.webgpuRecovery;
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.getImageJpeg(ctx, imageSize, props);
        });
        return this.runTask(task);
    }

    /** Get the current plugin state */
    async getStateSnapshot() {
        this.canvas3d!.commit(true);
        return await this.managers.snapshot.getStateSnapshot({ params: {} });
    }

    /** Save the current plugin state to a MOLJ file */
    async saveStateSnapshot(outPath: string) {
        const snapshot = await this.getStateSnapshot();
        const snapshot_json = JSON.stringify(snapshot, null, 2);
        await new Promise<void>(resolve => {
            fs.writeFile(outPath, snapshot_json, () => resolve());
        });
    }

    /** Render plugin state snapshots animation and return as raw MP4 data */
    async getAnimation(options?: { quantization?: number, size?: { width: number, height: number }, fps?: number, postprocessing?: Partial<PostprocessingProps> }) {
        const { Mp4Export } = await import('../extensions/mp4-export');
        if (!this.state.hasBehavior(Mp4Export)) {
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
                    getImageData: (runtime: RuntimeContext, width: number, height: number, viewport?: Viewport) => this.renderer.getImageRaw(runtime, { width, height }, options?.postprocessing, viewport),
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
