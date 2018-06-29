/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Viewer from 'mol-view/viewer'
import { StateContext } from './state/context';
import { Progress } from 'mol-task';
import { MmcifUrlToModel, ModelToStructure, StructureToSpacefill, StructureToBallAndStick, StructureToDistanceRestraint, StructureToCartoon, StructureToBackbone } from './state/transform';
import { UrlEntity } from './state/entity';
import { SpacefillProps } from 'mol-geo/representation/structure/spacefill';
import { Context } from 'mol-app/context/context';
import { BallAndStickProps } from 'mol-geo/representation/structure/ball-and-stick';
import { CartoonProps } from 'mol-geo/representation/structure/cartoon';
import { DistanceRestraintProps } from 'mol-geo/representation/structure/distance-restraint';
import { BackboneProps } from 'mol-geo/representation/structure/backbone';

const spacefillProps: SpacefillProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
}

const ballAndStickProps: BallAndStickProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    sizeTheme: { name: 'uniform', value: 0.25 },
    quality: 'auto',
    useFog: false
}

const distanceRestraintProps: DistanceRestraintProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    linkRadius: 0.5,
    quality: 'auto',
    useFog: false
}

const backboneProps: BackboneProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
}

const cartoonProps: CartoonProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
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
        // this.loadPdbid('4umt') // ligand has bond with order 3
        // this.loadPdbid('1crn') // small
        // this.loadPdbid('1rb8') // virus
        // this.loadPdbid('1blu') // metal coordination
        // this.loadPdbid('3pqr') // inter unit bonds
        // this.loadPdbid('4v5a') // ribosome
        // this.loadPdbid('3j3q') // ...
        this.loadPdbid('3sn6') // ...
        // this.loadMmcifUrl(`../../examples/1cbs_full.bcif`)

        // this.loadMmcifUrl(`../../../test/pdb-dev/PDBDEV_00000001.cif`)
    }

    async loadMmcifUrl (url: string) {
        const urlEntity = UrlEntity.ofUrl(this.ctx, url)
        const modelEntity = await MmcifUrlToModel.apply(this.ctx, urlEntity)
        const structureEntity = await ModelToStructure.apply(this.ctx, modelEntity)

        StructureToSpacefill.apply(this.ctx, structureEntity, { ...spacefillProps, visible: false })
        StructureToBallAndStick.apply(this.ctx, structureEntity, { ...ballAndStickProps, visible: false })
        StructureToDistanceRestraint.apply(this.ctx, structureEntity, { ...distanceRestraintProps, visible: false })
        StructureToBackbone.apply(this.ctx, structureEntity, { ...backboneProps, visible: true })
        StructureToCartoon.apply(this.ctx, structureEntity, { ...cartoonProps, visible: false })

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