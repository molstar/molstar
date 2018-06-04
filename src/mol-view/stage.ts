/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer'
import { StateContext } from './state/context';
import { Progress } from 'mol-task';
import { MmcifUrlToSpacefill } from './state/transform';
import { UrlEntity } from './state/entity';
import { SpacefillProps } from 'mol-geo/representation/structure/spacefill';

// export const ColorTheme = {
//     'atom-index': {},
//     'chain-id': {},
//     'element-symbol': {},
//     'instance-index': {},
//     'uniform': {}
// }
// export type ColorTheme = keyof typeof ColorTheme

const spacefillProps: SpacefillProps = {
    doubleSided: true,
    detail: 0,
    colorTheme: { name: 'atom-index' }
}

export class Stage {
    viewer: Viewer
    ctx = new StateContext(Progress.format)

    constructor() {

    }

    initRenderer (canvas: HTMLCanvasElement, container: HTMLDivElement) {
        this.viewer = Viewer.create(canvas, container)
        this.viewer.animate()
        this.ctx.viewer = this.viewer
        // this.loadPdbid('1crn')
        this.loadMmcifUrl(`../../examples/1cbs_full.bcif`)
    }

    loadMmcifUrl (url: string) {
        const urlEntity = UrlEntity.ofUrl(this.ctx, url)
        MmcifUrlToSpacefill.apply(this.ctx, urlEntity, spacefillProps)
    }

    loadPdbid (pdbid: string) {
        return this.loadMmcifUrl(`https://files.rcsb.org/download/${pdbid}.cif`)
    }

    dispose () {
        // TODO
    }
}