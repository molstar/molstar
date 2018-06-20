/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer'
import { StateContext } from './state/context';
import { Progress } from 'mol-task';
import { MmcifUrlToModel, ModelToStructure, StructureToSpacefill, StructureToBallAndStick } from './state/transform';
import { UrlEntity } from './state/entity';
import { SpacefillProps } from 'mol-geo/representation/structure/spacefill';
import { Context } from 'mol-app/context/context';
import { BallAndStickProps } from 'mol-geo/representation/structure/ball-and-stick';

const spacefillProps: SpacefillProps = {
    doubleSided: true,
    detail: 0,
    colorTheme: { name: 'atom-index' }
}

const ballAndStickProps: BallAndStickProps = {
    doubleSided: true,
    detail: 1,
    radialSegments: 8,
    colorTheme: { name: 'chain-id' },
    sizeTheme: { name: 'uniform', value: 0.25 },
}

export class Stage {
    viewer: Viewer
    ctx = new StateContext(Progress.format)

    constructor(public globalContext: Context) {

    }

    initRenderer (canvas: HTMLCanvasElement, container: HTMLDivElement) {
        this.viewer = Viewer.create(canvas, container)
        this.viewer.animate()
        this.ctx.viewer = this.viewer

        // this.loadPdbid('1jj2')
        this.loadPdbid('4umt') // ligand has bond with order 3
        // this.loadPdbid('1crn')
        // this.loadMmcifUrl(`../../examples/1cbs_full.bcif`)
    }

    async loadMmcifUrl (url: string) {
        const urlEntity = UrlEntity.ofUrl(this.ctx, url)
        const modelEntity = await MmcifUrlToModel.apply(this.ctx, urlEntity)
        const structureEntity = await ModelToStructure.apply(this.ctx, modelEntity)

        StructureToSpacefill.apply(this.ctx, structureEntity, { ...spacefillProps, visible: false })
        StructureToBallAndStick.apply(this.ctx, structureEntity, ballAndStickProps)

        this.globalContext.components.sequenceView.setState({ structure: structureEntity.value });
    }

    loadPdbid (pdbid: string) {
        return this.loadMmcifUrl(`http://www.ebi.ac.uk/pdbe/static/entry/${pdbid}_updated.cif`)
        // return this.loadMmcifUrl(`https://files.rcsb.org/download/${pdbid}.cif`)
    }

    dispose () {
        // TODO
    }
}