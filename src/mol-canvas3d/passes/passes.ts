/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DrawPass } from './draw';
import { PickPass } from './pick';
import { MultiSamplePass } from './multi-sample';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { AssetManager } from '../../mol-util/assets';
import { IlluminationPass } from './illumination';

export class Passes {
    readonly draw: DrawPass;
    readonly pick: PickPass;
    readonly multiSample: MultiSamplePass;
    readonly illumination: IlluminationPass;

    constructor(private webgl: WebGLContext, assetManager: AssetManager, attribs: Partial<{ pickScale: number, transparency: 'wboit' | 'dpoit' | 'blended' }> = {}) {
        const drs = this.webgl.getDrawingBufferSize();
        this.draw = new DrawPass(webgl, assetManager, drs.width, drs.height, attribs.transparency || 'blended');
        this.pick = new PickPass(webgl, drs.width, drs.height, attribs.pickScale || 0.25);
        this.multiSample = new MultiSamplePass(webgl, this.draw);
        this.illumination = new IlluminationPass(webgl, this.draw);
    }

    setPickScale(pickScale: number) {
        this.pick.setPickScale(pickScale);
    }

    setTransparency(transparency: 'wboit' | 'dpoit' | 'blended') {
        this.draw.setTransparency(transparency);
    }

    updateSize() {
        const drs = this.webgl.getDrawingBufferSize();
        // Avoid setting dimensions to 0x0 because it causes "empty textures are not allowed" error.
        const width = Math.max(drs.width, 2);
        const height = Math.max(drs.height, 2);
        this.draw.setSize(width, height);
        this.pick.setSize(width, height);
        this.multiSample.syncSize();
        this.illumination.setSize(width, height);
    }
}