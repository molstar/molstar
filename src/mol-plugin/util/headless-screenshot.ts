/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Jesse Liang <jesse.liang@rcsb.org>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Ke Ma <mark.ma@rcsb.org>
 * @author Adam Midlik <midlik@gmail.com>
 */

import fs from 'fs';
import path from 'path';
import { type BufferRet as JpegBufferRet } from 'jpeg-js'; // Only import type here, the actual import must be provided by the caller
import { type PNG } from 'pngjs'; // Only import type here, the actual import must be provided by the caller

import { Canvas3D, Canvas3DContext, Canvas3DProps, DefaultCanvas3DParams } from '../../mol-canvas3d/canvas3d';
import { ImagePass, ImageProps } from '../../mol-canvas3d/passes/image';
import { Passes } from '../../mol-canvas3d/passes/passes';
import { PostprocessingParams, PostprocessingProps } from '../../mol-canvas3d/passes/postprocessing';
import { createContext } from '../../mol-gl/webgl/context';
import { AssetManager } from '../../mol-util/assets';
import { ColorNames } from '../../mol-util/color/names';
import { PixelData } from '../../mol-util/image';
import { InputObserver } from '../../mol-util/input/input-observer';
import { ParamDefinition } from '../../mol-util/param-definition';
import { RuntimeContext } from '../../mol-task';

export interface ExternalModules {
    'gl': typeof import('gl'),
    'jpeg-js'?: typeof import('jpeg-js'),
    'pngjs'?: typeof import('pngjs'),
}

export type HeadlessScreenshotHelperOptions = {
    webgl?: WebGLContextAttributes,
    canvas?: Partial<Canvas3DProps>,
    imagePass?: Partial<ImageProps>,
}

export type RawImageData = {
    data: Uint8ClampedArray,
    width: number,
    height: number,
}

/** To render Canvas3D when running in Node.js (without DOM) */
export class HeadlessScreenshotHelper {
    readonly canvas3d: Canvas3D;
    readonly imagePass: ImagePass;

    constructor(readonly externalModules: ExternalModules, readonly canvasSize: { width: number, height: number }, canvas3d?: Canvas3D, options?: HeadlessScreenshotHelperOptions) {
        if (canvas3d) {
            this.canvas3d = canvas3d;
        } else {
            const glContext = this.externalModules.gl(this.canvasSize.width, this.canvasSize.height, options?.webgl ?? defaultWebGLAttributes());
            const webgl = createContext(glContext);
            const input = InputObserver.create();
            const attribs = { ...Canvas3DContext.DefaultAttribs };
            const props = { ...Canvas3DContext.DefaultProps };
            const assetManager = new AssetManager();
            const passes = new Passes(webgl, assetManager, props);
            const pixelScale = 1;
            const syncPixelScale = () => {};
            const setProps = () => {};
            const dispose = () => {
                input.dispose();
                webgl.destroy();
            };
            this.canvas3d = Canvas3D.create({ webgl, input, passes, attribs, props, assetManager, pixelScale, syncPixelScale, setProps, dispose }, options?.canvas ?? defaultCanvas3DParams());
        }

        this.imagePass = this.canvas3d.getImagePass(options?.imagePass ?? defaultImagePassParams());
        this.imagePass.setSize(this.canvasSize.width, this.canvasSize.height);
    }

    private async getImageData(runtime: RuntimeContext, width: number, height: number): Promise<RawImageData> {
        this.imagePass.setSize(width, height);
        await this.imagePass.render(runtime);
        this.imagePass.colorTarget.bind();

        const array = new Uint8Array(width * height * 4);
        this.canvas3d.webgl.readPixels(0, 0, width, height, array);
        const pixelData = PixelData.create(array, width, height);
        PixelData.flipY(pixelData);
        PixelData.divideByAlpha(pixelData);
        // ImageData is not defined in Node.js
        return { data: new Uint8ClampedArray(array), width, height };
    }

    async getImageRaw(runtime: RuntimeContext, imageSize?: { width: number, height: number }, postprocessing?: Partial<PostprocessingProps>): Promise<RawImageData> {
        const width = imageSize?.width ?? this.canvasSize.width;
        const height = imageSize?.height ?? this.canvasSize.height;
        this.canvas3d.commit(true);
        this.imagePass.setProps({
            postprocessing: ParamDefinition.merge(PostprocessingParams, this.canvas3d.props.postprocessing, postprocessing),
        });
        return this.getImageData(runtime, width, height);
    }

    async getImagePng(runtime: RuntimeContext, imageSize?: { width: number, height: number }, postprocessing?: Partial<PostprocessingProps>): Promise<PNG> {
        const imageData = await this.getImageRaw(runtime, imageSize, postprocessing);
        if (!this.externalModules.pngjs) {
            throw new Error("External module 'pngjs' was not provided. If you want to use getImagePng, you must import 'pngjs' and provide it to the HeadlessPluginContext/HeadlessScreenshotHelper constructor.");
        }
        const generatedPng = new this.externalModules.pngjs.PNG({ width: imageData.width, height: imageData.height });
        generatedPng.data = Buffer.from(imageData.data.buffer);
        return generatedPng;
    }

    async getImageJpeg(runtime: RuntimeContext, imageSize?: { width: number, height: number }, postprocessing?: Partial<PostprocessingProps>, jpegQuality: number = 90): Promise<JpegBufferRet> {
        const imageData = await this.getImageRaw(runtime, imageSize, postprocessing);
        if (!this.externalModules['jpeg-js']) {
            throw new Error("External module 'jpeg-js' was not provided. If you want to use getImageJpeg, you must import 'jpeg-js' and provide it to the HeadlessPluginContext/HeadlessScreenshotHelper constructor.");
        }
        const generatedJpeg = this.externalModules['jpeg-js'].encode(imageData, jpegQuality);
        return generatedJpeg;
    }

    async saveImage(runtime: RuntimeContext, outPath: string, imageSize?: { width: number, height: number }, postprocessing?: Partial<PostprocessingProps>, format?: 'png' | 'jpeg', jpegQuality = 90) {
        if (!format) {
            const extension = path.extname(outPath).toLowerCase();
            if (extension === '.png') format = 'png';
            else if (extension === '.jpg' || extension === '.jpeg') format = 'jpeg';
            else throw new Error(`Cannot guess image format from file path '${outPath}'. Specify format explicitly or use path with one of these extensions: .png, .jpg, .jpeg`);
        }
        if (format === 'png') {
            const generatedPng = await this.getImagePng(runtime, imageSize, postprocessing);
            await writePngFile(generatedPng, outPath);
        } else if (format === 'jpeg') {
            const generatedJpeg = await this.getImageJpeg(runtime, imageSize, postprocessing, jpegQuality);
            await writeJpegFile(generatedJpeg, outPath);
        } else {
            throw new Error(`Invalid format: ${format}`);
        }
    }
}

async function writePngFile(png: PNG, outPath: string) {
    await new Promise<void>(resolve => {
        png.pack().pipe(fs.createWriteStream(outPath)).on('finish', resolve);
    });
}
async function writeJpegFile(jpeg: JpegBufferRet, outPath: string) {
    await new Promise<void>(resolve => {
        fs.writeFile(outPath, jpeg.data, () => resolve());
    });
}

export function defaultCanvas3DParams(): Partial<Canvas3DProps> {
    return {
        camera: {
            mode: 'orthographic',
            helper: {
                axes: { name: 'off', params: {} }
            },
            stereo: {
                name: 'off', params: {}
            },
            fov: 90,
            manualReset: false,
            scale: 1,
        },
        cameraResetDurationMs: 0,
        cameraFog: {
            name: 'on',
            params: {
                intensity: 50
            }
        },
        renderer: {
            ...DefaultCanvas3DParams.renderer,
            backgroundColor: ColorNames.white,
        },
        postprocessing: {
            ...DefaultCanvas3DParams.postprocessing,
            occlusion: {
                name: 'off', params: {}
            },
            outline: {
                name: 'off', params: {}
            },
            antialiasing: {
                name: 'fxaa',
                params: {
                    edgeThresholdMin: 0.0312,
                    edgeThresholdMax: 0.063,
                    iterations: 12,
                    subpixelQuality: 0.3
                }
            },
            background: { variant: { name: 'off', params: {} } },
            shadow: { name: 'off', params: {} },
        }
    };
}

export function defaultWebGLAttributes(): WebGLContextAttributes {
    return {
        antialias: true,
        preserveDrawingBuffer: true,
        alpha: true, // the renderer requires an alpha channel
        depth: true, // the renderer requires a depth buffer
        premultipliedAlpha: true, // the renderer outputs PMA
    };
}

export function defaultImagePassParams(): Partial<ImageProps> {
    return {
        cameraHelper: {
            axes: { name: 'off', params: {} },
        },
        multiSample: {
            ...DefaultCanvas3DParams.multiSample,
            mode: 'on',
            sampleLevel: 4,
            reuseOcclusion: false,
        }
    };
}

export const STYLIZED_POSTPROCESSING: Partial<PostprocessingProps> = {
    occlusion: {
        name: 'on' as const, params: {
            samples: 32,
            multiScale: { name: 'off', params: {} },
            radius: 5,
            bias: 0.8,
            blurKernelSize: 15,
            blurDepthBias: 0.5,
            resolutionScale: 1,
            color: ColorNames.black,
            transparentThreshold: 0.4,
        }
    }, outline: {
        name: 'on' as const, params: {
            scale: 1,
            threshold: 0.95,
            color: ColorNames.black,
            includeTransparent: true,
        }
    }
};
