/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Viewport } from '../../mol-canvas3d/camera/util';
import { CameraHelperParams } from '../../mol-canvas3d/helper/camera-helper';
import { IlluminationProps } from '../../mol-canvas3d/passes/illumination';
import { ImagePass } from '../../mol-canvas3d/passes/image';
import { PostprocessingProps } from '../../mol-canvas3d/passes/postprocessing';
import { canvasToBlob } from '../../mol-canvas3d/util';
import { equalEps } from '../../mol-math/linear-algebra/3d/common';
import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateSelection } from '../../mol-state';
import { RuntimeContext, Task } from '../../mol-task';
import { Color } from '../../mol-util/color';
import { download } from '../../mol-util/download';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { SetUtils } from '../../mol-util/set';
import { PluginContext } from '../context';

export { ViewportScreenshotHelper, ViewportScreenshotHelperParams };

namespace ViewportScreenshotHelper {
    export type ResolutionSettings = PD.Values<ReturnType<ViewportScreenshotHelper['createParams']>>['resolution']
    export type ResolutionTypes = ResolutionSettings['name']
}

type ViewportScreenshotHelperParams = PD.Values<ReturnType<ViewportScreenshotHelper['createParams']>>

class ViewportScreenshotHelper extends PluginComponent {
    private createParams() {
        let max = 8192;
        if (this.plugin.canvas3d) {
            const { webgl } = this.plugin.canvas3d;
            max = Math.floor(Math.min(webgl.maxRenderbufferSize, webgl.maxTextureSize) / 2);
        }
        return {
            resolution: PD.MappedStatic('viewport', {
                viewport: PD.Group({}),
                hd: PD.Group({}),
                'full-hd': PD.Group({}),
                'ultra-hd': PD.Group({}),
                custom: PD.Group({
                    width: PD.Numeric(1920, { min: 128, max, step: 1 }),
                    height: PD.Numeric(1080, { min: 128, max, step: 1 }),
                }, { isFlat: true })
            }, {
                options: [
                    ['viewport', 'Viewport'],
                    ['hd', 'HD (1280 x 720)'],
                    ['full-hd', 'Full HD (1920 x 1080)'],
                    ['ultra-hd', 'Ultra HD (3840 x 2160)'],
                    ['custom', 'Custom']
                ]
            }),
            transparent: PD.Boolean(false),
            axes: CameraHelperParams.axes,
            illumination: PD.Group({
                extraIterations: PD.Numeric(1, { min: 0, max: 5, step: 1 }),
                targetIterationTimeMs: PD.Numeric(300, { min: 100, max: 3000, step: 10 }),
            }),
        };
    }
    private _params: ReturnType<ViewportScreenshotHelper['createParams']> = void 0 as any;
    get params() {
        if (this._params) return this._params;
        return this._params = this.createParams();
    }

    readonly behaviors = {
        values: this.ev.behavior<ViewportScreenshotHelperParams>({
            transparent: this.params.transparent.defaultValue,
            axes: { name: 'off', params: {} },
            resolution: this.params.resolution.defaultValue,
            illumination: this.params.illumination.defaultValue,
        }),
        cropParams: this.ev.behavior<{ auto: boolean, relativePadding: number }>({ auto: true, relativePadding: 0.1 }),
        relativeCrop: this.ev.behavior<Viewport>({ x: 0, y: 0, width: 1, height: 1 }),
    };

    readonly events = {
        previewed: this.ev<any>()
    };

    get values() {
        return this.behaviors.values.value;
    }

    get cropParams() {
        return this.behaviors.cropParams.value;
    }

    get relativeCrop() {
        return this.behaviors.relativeCrop.value;
    }

    private getCanvasSize() {
        return {
            width: this.plugin.canvas3d?.webgl.gl.drawingBufferWidth || 0,
            height: this.plugin.canvas3d?.webgl.gl.drawingBufferHeight || 0
        };
    }

    private getSize() {
        const values = this.values;
        switch (values.resolution.name) {
            case 'viewport': return this.getCanvasSize();
            case 'hd': return { width: 1280, height: 720 };
            case 'full-hd': return { width: 1920, height: 1080 };
            case 'ultra-hd': return { width: 3840, height: 2160 };
            default: return { width: values.resolution.params.width, height: values.resolution.params.height };
        }
    }

    private getPostprocessingProps() {
        const c = this.plugin.canvas3d!;
        const aoProps = c.props.postprocessing.occlusion;
        return {
            ...c.props.postprocessing,
            occlusion: aoProps.name === 'on'
                ? { name: 'on', params: { ...aoProps.params, samples: 128, resolutionScale: c.webgl.pixelRatio, transparentThreshold: 1 } }
                : aoProps
        } as PostprocessingProps;
    }

    private getIlluminationProps(isPreview: boolean) {
        const c = this.plugin.canvas3d!;
        const giProps = c.props.illumination;
        const { extraIterations, targetIterationTimeMs } = this.values.illumination;
        return {
            ...giProps,
            enabled: isPreview ? false : giProps.enabled,
            maxIterations: Math.ceil(Math.log2(Math.pow(2, giProps.maxIterations + extraIterations) * giProps.rendersPerFrame[1])),
            targetFps: 1000 / targetIterationTimeMs,
            denoiseThreshold: [giProps.denoiseThreshold[0], giProps.denoiseThreshold[0]],
            rendersPerFrame: [1, 1],
        } as IlluminationProps;
    }

    private createPass(isPreview: boolean) {
        const c = this.plugin.canvas3d!;
        const { colorBufferFloat, textureFloat } = c.webgl.extensions;
        return c.getImagePass({
            transparentBackground: this.values.transparent,
            cameraHelper: { axes: this.values.axes },
            multiSample: {
                ...c.props.multiSample,
                mode: isPreview ? 'off' : 'on',
                sampleLevel: colorBufferFloat && textureFloat ? 4 : 2,
                reuseOcclusion: false,
            },
            postprocessing: this.getPostprocessingProps(),
            marking: { ...c.props.marking },
            illumination: this.getIlluminationProps(isPreview),
        });
    }

    private _previewPass: ImagePass;
    private get previewPass() {
        return this._previewPass || (this._previewPass = this.createPass(true));
    }

    private _imagePass: ImagePass;
    get imagePass() {
        if (this._imagePass) {
            const c = this.plugin.canvas3d!;
            this._imagePass.setProps({
                cameraHelper: { axes: this.values.axes },
                transparentBackground: this.values.transparent,
                postprocessing: this.getPostprocessingProps(),
                marking: { ...c.props.marking },
                illumination: this.getIlluminationProps(false),
            });
            return this._imagePass;
        }
        return this._imagePass = this.createPass(false);
    }

    getFilename(extension = '.png') {
        const models = this.plugin.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model)).map(s => s.obj!.data);
        const uniqueIds = new Set<string>();
        models.forEach(m => uniqueIds.add(m.entryId.toUpperCase()));
        const idString = SetUtils.toArray(uniqueIds).join('-');
        return `${idString || 'molstar-image'}${extension}`;
    }

    private canvas = function () {
        const canvas = document.createElement('canvas');
        return canvas;
    }();

    private previewCanvas = function () {
        const canvas = document.createElement('canvas');
        return canvas;
    }();

    private previewData = {
        image: { data: new Uint8ClampedArray(1), width: 1, height: 0 } as ImageData,
        background: Color(0),
        transparent: false
    };

    resetCrop() {
        this.behaviors.relativeCrop.next({ x: 0, y: 0, width: 1, height: 1 });
    }

    toggleAutocrop() {
        if (this.cropParams.auto) {
            this.behaviors.cropParams.next({ ...this.cropParams, auto: false });
            this.resetCrop();
        } else {
            this.behaviors.cropParams.next({ ...this.cropParams, auto: true });
        }
    }

    get isFullFrame() {
        const crop = this.relativeCrop;
        return equalEps(crop.x, 0, 1e-5) && equalEps(crop.y, 0, 1e-5) && equalEps(crop.width, 1, 1e-5) && equalEps(crop.height, 1, 1e-5);
    }

    autocrop(relativePadding = this.cropParams.relativePadding) {
        const { data, width, height } = this.previewData.image;
        const isTransparent = this.previewData.transparent;
        const bgColor = isTransparent ? this.previewData.background : 0xff000000 | this.previewData.background;

        let l = width, r = 0, t = height, b = 0;

        for (let j = 0; j < height; j++) {
            const jj = j * width;
            for (let i = 0; i < width; i++) {
                const o = 4 * (jj + i);

                if (isTransparent) {
                    if (data[o + 3] === 0) continue;
                } else {
                    const c = (data[o] << 16) | (data[o + 1] << 8) | (data[o + 2]) | (data[o + 3] << 24);
                    if (c === bgColor) continue;
                }

                if (i < l) l = i;
                if (i > r) r = i;
                if (j < t) t = j;
                if (j > b) b = j;
            }
        }

        if (l > r) {
            const x = l;
            l = r;
            r = x;
        }

        if (t > b) {
            const x = t;
            t = b;
            b = x;
        }

        const tw = r - l + 1, th = b - t + 1;
        l -= relativePadding * tw;
        r += relativePadding * tw;
        t -= relativePadding * th;
        b += relativePadding * th;

        const crop: Viewport = {
            x: Math.max(0, l / width),
            y: Math.max(0, t / height),
            width: Math.min(1, (r - l + 1) / width),
            height: Math.min(1, (b - t + 1) / height)
        };

        this.behaviors.relativeCrop.next(crop);
    }

    async getPreview(ctx: RuntimeContext, maxDim = 320) {
        const { width, height } = this.getSize();
        if (width <= 0 || height <= 0) return;

        const f = width / height;

        let w = 0, h = 0;
        if (f > 1) {
            w = maxDim;
            h = Math.round(maxDim / f);
        } else {
            h = maxDim;
            w = Math.round(maxDim * f);
        }

        const canvasProps = this.plugin.canvas3d!.props;
        this.previewPass.setProps({
            cameraHelper: { axes: this.values.axes },
            transparentBackground: this.values.transparent,
            postprocessing: canvasProps.postprocessing,
            marking: canvasProps.marking,
        });
        const imageData = await this.previewPass.getImageData(ctx, w, h);
        const canvas = this.previewCanvas;
        canvas.width = imageData.width;
        canvas.height = imageData.height;

        this.previewData.image = imageData;
        this.previewData.background = canvasProps.renderer.backgroundColor;
        this.previewData.transparent = this.values.transparent;

        const canvasCtx = canvas.getContext('2d');
        if (!canvasCtx) throw new Error('Could not create canvas 2d context');
        canvasCtx.putImageData(imageData, 0, 0);
        if (this.cropParams.auto) this.autocrop();

        this.events.previewed.next(void 0);
        return { canvas, width: w, height: h };
    }

    getSizeAndViewport() {
        const { width, height } = this.getSize();
        const crop = this.relativeCrop;
        const viewport: Viewport = {
            x: Math.floor(crop.x * width),
            y: Math.floor(crop.y * height),
            width: Math.ceil(crop.width * width),
            height: Math.ceil(crop.height * height)
        };
        if (viewport.width + viewport.x > width) viewport.width = width - viewport.x;
        if (viewport.height + viewport.y > height) viewport.height = height - viewport.y;
        return { width, height, viewport };
    }

    private async draw(ctx: RuntimeContext) {
        const { width, height, viewport } = this.getSizeAndViewport();
        if (width <= 0 || height <= 0) return;

        this.plugin.canvas3d?.pause(true);
        try {
            await ctx.update('Rendering image...');
            const pass = this.imagePass;
            await pass.updateBackground();
            const imageData = await pass.getImageData(ctx, width, height, viewport);

            await ctx.update('Encoding image...');
            const canvas = this.canvas;
            canvas.width = imageData.width;
            canvas.height = imageData.height;
            const canvasCtx = canvas.getContext('2d');
            if (!canvasCtx) throw new Error('Could not create canvas 2d context');
            canvasCtx.putImageData(imageData, 0, 0);
        } finally {
            this.plugin.canvas3d?.animate();
        }
        return;
    }

    private copyToClipboardTask() {
        const cb = navigator.clipboard as any;

        if (!cb?.write) {
            this.plugin.log.error('clipboard.write not supported!');
            return;
        }

        return Task.create('Copy Image', async ctx => {
            await this.draw(ctx);
            await ctx.update('Converting image...');
            const blob = await canvasToBlob(this.canvas, 'png');
            const item = new ClipboardItem({ 'image/png': blob });
            await cb.write([item]);
            this.plugin.log.message('Image copied to clipboard.');
        });
    }

    getImageDataUri() {
        return this.plugin.runTask(Task.create('Generate Image', async ctx => {
            await this.draw(ctx);
            await ctx.update('Converting image...');
            return this.canvas.toDataURL('png');
        }));
    }

    copyToClipboard() {
        const task = this.copyToClipboardTask();
        if (!task) return;
        return this.plugin.runTask(task);
    }

    private downloadTask(filename?: string) {
        return Task.create('Download Image', async ctx => {
            await this.draw(ctx);
            await ctx.update('Downloading image...');
            const blob = await canvasToBlob(this.canvas, 'png');
            download(blob, filename ?? this.getFilename());
        });
    }

    download(filename?: string) {
        return this.plugin.runTask(this.downloadTask(filename), { useOverlay: true });
    }

    constructor(private plugin: PluginContext) {
        super();
    }
}

declare const ClipboardItem: any;