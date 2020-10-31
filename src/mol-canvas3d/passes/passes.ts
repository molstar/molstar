/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DrawPass } from './draw';
import { PickPass } from './pick';
import { PostprocessingPass } from './postprocessing';
import { MultiSamplePass } from './multi-sample';
import { WebGLContext } from '../../mol-gl/webgl/context';

export class Passes {
    readonly draw: DrawPass
    readonly pick: PickPass
    readonly postprocessing: PostprocessingPass
    readonly multiSample: MultiSamplePass

    constructor(private webgl: WebGLContext, attribs: Partial<{ pickScale: number }> = {}) {
        const { gl } = webgl;
        this.draw = new DrawPass(webgl, gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.pick = new PickPass(webgl, this.draw, attribs.pickScale || 0.25);
        this.postprocessing = new PostprocessingPass(webgl, this.draw);
        this.multiSample = new MultiSamplePass(webgl, this.draw, this.postprocessing);
    }

    updateSize() {
        const { gl } = this.webgl;
        this.draw.setSize(gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.pick.syncSize();
        this.postprocessing.syncSize();
        this.multiSample.syncSize();
    }
}