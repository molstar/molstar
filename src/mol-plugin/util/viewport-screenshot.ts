/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CameraHelperParams } from '../../mol-canvas3d/helper/camera-helper';
import { ImagePass } from '../../mol-canvas3d/passes/image';
import { canvasToBlob } from '../../mol-canvas3d/util';
import { PluginComponent } from '../../mol-plugin-state/component';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateSelection } from '../../mol-state';
import { RuntimeContext, Task } from '../../mol-task';
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
        const max = Math.min(this.plugin.canvas3d ? this.plugin.canvas3d.webgl.maxRenderbufferSize : 4096, 4096);
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
            resolution: this.params.resolution.defaultValue
        })
    };

    get values() {
        return this.behaviors.values.value;
    }

    private getCanvasSize() {
        return {
            width: this.plugin.canvas3d?.webgl.gl.drawingBufferWidth || 0,
            height: this.plugin.canvas3d?.webgl.gl.drawingBufferHeight || 0
        };
    }

    private getSize() {
        const values = this.values;
        switch (values.resolution.name ) {
            case 'viewport': return this.getCanvasSize();
            case 'hd': return { width: 1280, height: 720 };
            case 'full-hd': return { width: 1920, height: 1080 };
            case 'ultra-hd': return { width: 3840, height: 2160 };
            default: return { width: values.resolution.params.width, height: values.resolution.params.height };
        }
    }

    private createPass(mutlisample: boolean) {
        const c = this.plugin.canvas3d!;
        return this.plugin.canvas3d!.getImagePass({
            transparentBackground: this.values.transparent,
            cameraHelper: { axes: this.values.axes },
            multiSample: {
                mode: mutlisample ? 'on' : 'off',
                sampleLevel: c.webgl.extensions.colorBufferFloat ? 4 : 2
            },
            postprocessing: c.props.postprocessing
        });
    }

    private _previewPass: ImagePass;
    private get previewPass() {
        return this._previewPass || (this._previewPass = this.createPass(false));
    }

    private _imagePass: ImagePass;
    get imagePass() {
        return this._imagePass || (this._imagePass = this.createPass(true));
    }

    getFilename() {
        const models = this.plugin.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model)).map(s => s.obj!.data);
        const uniqueIds = new Set<string>();
        models.forEach(m => uniqueIds.add(m.entryId.toUpperCase()));
        const idString = SetUtils.toArray(uniqueIds).join('-');
        return `${idString || 'molstar-image'}.png`;
    }

    private canvas = function () {
        const canvas = document.createElement('canvas');
        return canvas;
    }();

    private previewCanvas = function () {
        const canvas = document.createElement('canvas');
        return canvas;
    }();

    getPreview(maxDim = 640) {
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

        this.previewPass.setProps({
            cameraHelper: { axes: this.values.axes },
            transparentBackground: this.values.transparent,
            // TODO: optimize because this creates a copy of a large object!
            postprocessing: this.plugin.canvas3d!.props.postprocessing
        });
        const imageData = this.previewPass.getImageData(w, h);
        const canvas = this.previewCanvas;
        canvas.width = imageData.width;
        canvas.height = imageData.height;
        const canvasCtx = canvas.getContext('2d');
        if (!canvasCtx) throw new Error('Could not create canvas 2d context');
        canvasCtx.putImageData(imageData, 0, 0);
        return { canvas, width: w, height: h };
    }

    private async draw(ctx: RuntimeContext) {
        const { width, height } = this.getSize();
        if (width <= 0 || height <= 0) return;

        await ctx.update('Rendering image...');
        this.imagePass.setProps({
            cameraHelper: { axes: this.values.axes },
            transparentBackground: this.values.transparent,
            // TODO: optimize because this creates a copy of a large object!
            postprocessing: this.plugin.canvas3d!.props.postprocessing
        });
        const imageData = this.imagePass.getImageData(width, height);

        await ctx.update('Encoding image...');
        const canvas = this.canvas;
        canvas.width = imageData.width;
        canvas.height = imageData.height;
        const canvasCtx = canvas.getContext('2d');
        if (!canvasCtx) throw new Error('Could not create canvas 2d context');
        canvasCtx.putImageData(imageData, 0, 0);
        return;
    }

    private copyToClipboardTask() {
        const cb = navigator.clipboard as any;

        if (!cb.write) {
            this.plugin.log.error('clipboard.write not supported!');
            return;
        }

        return Task.create('Copy Image', async ctx => {
            await this.draw(ctx);
            await ctx.update('Downloading image...');
            const blob = await canvasToBlob(this.canvas, 'png');
            const item = new ClipboardItem({ 'image/png': blob });
            cb.write([item]);
            this.plugin.log.message('Image copied to clipboard.');
        });
    }

    copyToClipboard() {
        const task = this.copyToClipboardTask();
        if (!task) return;
        return this.plugin.runTask(task);
    }

    private downloadTask() {
        return Task.create('Download Image', async ctx => {
            await this.draw(ctx);
            await ctx.update('Downloading image...');
            const blob = await canvasToBlob(this.canvas, 'png');
            download(blob, this.getFilename());
        });
    }

    download() {
        this.plugin.runTask(this.downloadTask());
    }

    constructor(private plugin: PluginContext) {
        super();
    }
}

declare const ClipboardItem: any;