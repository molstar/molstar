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
import { Task } from '../../mol-task';
import { canvasToBlob } from '../../mol-canvas3d/util';
import { download } from '../../mol-util/download';

export class ViewportScreenshotWrapper {
    private getCanvasSize() {
        return {
            width: this.plugin.canvas3d?.webgl.gl.drawingBufferWidth || 0,
            height: this.plugin.canvas3d?.webgl.gl.drawingBufferHeight || 0
        };
    }

    size = this.getCanvasSize();

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

    private downloadTask() {
        return Task.create('Download Image', async ctx => {
            const { width, height } = this.size
            if (width <= 0 || height <= 0) return

            await ctx.update('Rendering image...')
            const imageData = this.imagePass.getImageData(width, height)

            await ctx.update('Encoding image...')
            const canvas = document.createElement('canvas')
            canvas.width = imageData.width
            canvas.height = imageData.height
            const canvasCtx = canvas.getContext('2d')
            if (!canvasCtx) throw new Error('Could not create canvas 2d context')
            canvasCtx.putImageData(imageData, 0, 0)

            await ctx.update('Downloading image...')
            const blob = await canvasToBlob(canvas, 'png')
            download(blob, this.getFilename())
        })
    }

    download() {
        this.plugin.runTask(this.downloadTask());
    }

    constructor(private plugin: PluginContext) {

    }
}