/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginContext } from '../context';
import { ImagePass } from '../../mol-canvas3d/passes/image';
import { StateSelection } from '../../mol-state';
import { PluginStateObject } from '../state/objects';
import { Task, RuntimeContext } from '../../mol-task';
import { canvasToBlob } from '../../mol-canvas3d/util';
import { download } from '../../mol-util/download';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { SyncRuntimeContext } from '../../mol-task/execution/synchronous';

export class ViewportScreenshotWrapper {
    get params() {
        const max = Math.min(this.plugin.canvas3d ? this.plugin.canvas3d.webgl.maxRenderbufferSize : 4096, 4096)
        const { width, height } = this.size
        return {
            size: PD.MappedStatic(this.size.type, {
                canvas: PD.Group({}),
                custom: PD.Group({
                    width: PD.Numeric(width, { min: 128, max, step: 1 }),
                    height: PD.Numeric(height, { min: 128, max, step: 1 }),
                }, { isFlat: true })
            }, { options: [['canvas', 'Canvas'], ['custom', 'Custom']] })
        }
    }

    get values() {
        return this.size.type === 'canvas'
            ? { size: { name: 'canvas', params: {} } }
            : { size: { name: 'custom', params: { width: this.size.width, height: this.size.height } } }
    }

    private getCanvasSize() {
        return {
            width: this.plugin.canvas3d?.webgl.gl.drawingBufferWidth || 0,
            height: this.plugin.canvas3d?.webgl.gl.drawingBufferHeight || 0
        };
    }

    size = {
        type: 'custom' as 'canvas' | 'custom',
        width: 1920,
        height: 1080
    }

    getSize() {
        if (this.size.type === 'canvas') return this.getCanvasSize();
        return { width: this.size.width, height: this.size.height };
    }

    private _imagePass: ImagePass;

    get imagePass() {
        if (this._imagePass) return this._imagePass;

        this._imagePass = this.plugin.canvas3d!.getImagePass()
        this._imagePass.setProps({
            multiSample: { mode: 'on', sampleLevel: 2 },
            postprocessing: this.plugin.canvas3d!.props.postprocessing
        });
        return this._imagePass;
    }

    getFilename() {
        const models = this.plugin.state.dataState.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Model)).map(s => s.obj!.data)
        const uniqueIds = new Set<string>()
        models.forEach(m => uniqueIds.add(m.entryId.toUpperCase()))
        const idString = Array.from(uniqueIds).join('-')
        return `${idString || 'molstar-image'}.png`
    }

    private canvas = function () {
        const canvas = document.createElement('canvas');
        return canvas;
    }();

    private async draw(ctx: RuntimeContext) {
        const { width, height } = this.getSize();
        if (width <= 0 || height <= 0) return;

        await ctx.update('Rendering image...')
        const imageData = this.imagePass.getImageData(width, height);

        await ctx.update('Encoding image...')
        const canvas = this.canvas
        canvas.width = imageData.width
        canvas.height = imageData.height
        const canvasCtx = canvas.getContext('2d')
        if (!canvasCtx) throw new Error('Could not create canvas 2d context')
        canvasCtx.putImageData(imageData, 0, 0)
        return;
    }

    private downloadTask() {
        return Task.create('Download Image', async ctx => {
            this.draw(ctx);
            await ctx.update('Downloading image...')
            const blob = await canvasToBlob(this.canvas, 'png')
            download(blob, this.getFilename())
        })
    }

    async imageData() {
        await this.draw(SyncRuntimeContext)
        return this.canvas.toDataURL();
    }

    download() {
        this.plugin.runTask(this.downloadTask());
    }

    constructor(private plugin: PluginContext) {

    }
}