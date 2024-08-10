/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import fs from 'fs';
import { type PNG } from 'pngjs'; // Only import type here, the actual import must be provided by the caller
import { type BufferRet as JpegBufferRet } from 'jpeg-js'; // Only import type here, the actual import must be provided by the caller

import { Canvas3D } from '../mol-canvas3d/canvas3d';
import { PostprocessingProps } from '../mol-canvas3d/passes/postprocessing';
import { PluginContext } from './context';
import { PluginSpec } from './spec';
import { HeadlessScreenshotHelper, HeadlessScreenshotHelperOptions, ExternalModules, RawImageData } from './util/headless-screenshot';
import { Task } from '../mol-task';


/** PluginContext that can be used in Node.js (without DOM) */
export class HeadlessPluginContext extends PluginContext {
    renderer: HeadlessScreenshotHelper;

    /** External modules (`gl` and optionally `pngjs` and `jpeg-js`) must be provided to the constructor (this is to avoid Mol* being dependent on these packages which are only used here) */
    constructor(externalModules: ExternalModules, spec: PluginSpec, canvasSize: { width: number, height: number } = { width: 640, height: 480 }, rendererOptions?: HeadlessScreenshotHelperOptions) {
        super(spec);
        this.renderer = new HeadlessScreenshotHelper(externalModules, canvasSize, undefined, rendererOptions);
        (this.canvas3d as Canvas3D) = this.renderer.canvas3d;
    }

    /** Render the current plugin state and save to a PNG or JPEG file */
    async saveImage(outPath: string, imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>, format?: 'png' | 'jpeg', jpegQuality = 90) {
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.saveImage(ctx, outPath, imageSize, props, format, jpegQuality);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as raw image data */
    async getImageRaw(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>): Promise<RawImageData> {
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.getImageRaw(ctx, imageSize, props);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as a PNG object */
    async getImagePng(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>): Promise<PNG> {
        const task = Task.create('Render Screenshot', async ctx => {
            this.canvas3d!.commit(true);
            return await this.renderer.getImagePng(ctx, imageSize, props);
        });
        return this.runTask(task);
    }

    /** Render the current plugin state and return as a JPEG object */
    async getImageJpeg(imageSize?: { width: number, height: number }, props?: Partial<PostprocessingProps>, jpegQuality: number = 90): Promise<JpegBufferRet> {
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
}
