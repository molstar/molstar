/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject } from 'rxjs';

// import { ValueCell } from 'mol-util/value-cell'

// import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import Viewer from 'mol-view/viewer'
// import { createColorTexture } from 'mol-gl/util';
// import Icosahedron from 'mol-geo/primitive/icosahedron'
// import Box from 'mol-geo/primitive/box'
import Spacefill, { SpacefillProps } from 'mol-geo/representation/structure/spacefill'
import Point, { PointProps } from 'mol-geo/representation/structure/point'

import { Run } from 'mol-task'
import { Symmetry, Structure, Model } from 'mol-model/structure'

// import mcubes from './utils/mcubes'
import { getModelFromPdbId, getModelFromFile, log, Volume, getVolumeFromEmdId } from './utils'
import { StructureRepresentation } from 'mol-geo/representation/structure';
import { Color } from 'mol-util/color';
import Surface, { SurfaceProps } from 'mol-geo/representation/volume/surface';
import { VolumeIsoValue } from 'mol-model/volume';
import { VolumeRepresentation } from 'mol-geo/representation/volume';
// import Cylinder from 'mol-geo/primitive/cylinder';


export const ColorTheme = {
    'atom-index': {},
    'chain-id': {},
    'element-symbol': {},
    'instance-index': {},
    'uniform': {}
}
export type ColorTheme = keyof typeof ColorTheme

export default class State {
    viewer: Viewer
    pdbId = ''
    // pdbId = '5ire'
    emdId = '8116'
    // pdbId = '6G1K'
    // emdId = '4339'
    // pdbId = '4cup'
    // emdId = ''
    model = new BehaviorSubject<Model | undefined>(undefined)
    volume = new BehaviorSubject<Volume | undefined>(undefined)
    initialized = new BehaviorSubject<boolean>(false)
    loading = new BehaviorSubject<boolean>(false)

    colorTheme = new BehaviorSubject<ColorTheme>('element-symbol')
    colorValue = new BehaviorSubject<Color>(0xFF4411)
    sphereDetail = new BehaviorSubject<number>(0)
    assembly = new BehaviorSubject<string>('')

    pointVisibility = new BehaviorSubject<boolean>(true)
    spacefillVisibility = new BehaviorSubject<boolean>(true)

    pointRepr: StructureRepresentation<PointProps>
    spacefillRepr: StructureRepresentation<SpacefillProps>
    surfaceRepr: VolumeRepresentation<SurfaceProps>

    constructor() {
        this.colorTheme.subscribe(() => this.update())
        this.colorValue.subscribe(() => this.update())
        this.sphereDetail.subscribe(() => this.update())
        this.assembly.subscribe(() => this.initStructure())

        this.pointVisibility.subscribe(() => this.updateVisibility())
        this.spacefillVisibility.subscribe(() => this.updateVisibility())
    }

    getSpacefillProps (): SpacefillProps {
        const colorThemeName = this.colorTheme.getValue()
        return {
            detail: this.sphereDetail.getValue(),
            colorTheme: colorThemeName === 'uniform' ?
                { name: colorThemeName, value: this.colorValue.getValue() } :
                { name: colorThemeName }
        }
    }

    getPointProps (): PointProps {
        const colorThemeName = this.colorTheme.getValue()
        return {
            sizeTheme: { name: 'uniform', value: 0.1 },
            colorTheme: colorThemeName === 'uniform' ?
                { name: colorThemeName, value: this.colorValue.getValue() } :
                { name: colorThemeName }
        }
    }

    async initRenderer (canvas: HTMLCanvasElement, container: HTMLDivElement) {
        this.viewer = Viewer.create(canvas, container)
        this.initialized.next(true)
        this.loadPdbId()
        this.loadEmdId()
        this.viewer.animate()
    }

    async getStructure () {
        const model = this.model.getValue()
        if (!model) return
        const assembly = this.assembly.getValue()
        let structure: Structure
        const assemblies = model.symmetry.assemblies
        if (assemblies.length) {
            structure = await Run(Symmetry.buildAssembly(Structure.ofModel(model), assembly || '1'), log, 500)
        } else {
            structure = Structure.ofModel(model)
        }
        return structure
    }

    async initStructure () {
        const { viewer } = this
        if (!viewer || !this.model.getValue()) return

        if (this.pointRepr) this.viewer.remove(this.pointRepr)
        if (this.spacefillRepr) this.viewer.remove(this.spacefillRepr)

        const structure = await this.getStructure()
        if (!structure) return

        this.pointRepr = StructureRepresentation(Point)
        await Run(this.pointRepr.create(structure, this.getPointProps()), log, 500)
        viewer.add(this.pointRepr)

        this.spacefillRepr = StructureRepresentation(Spacefill)
        await Run(this.spacefillRepr.create(structure, this.getSpacefillProps()), log, 500)
        viewer.add(this.spacefillRepr)

        this.updateVisibility()
        viewer.requestDraw()
        console.log(viewer.stats)
    }

    setModel(model: Model) {
        this.model.next(model)
        this.initStructure()
        this.loading.next(false)
    }

    async loadFile (file: File) {
        this.viewer.clear()
        this.loading.next(true)
        this.setModel((await getModelFromFile(file))[0])
    }

    async initVolume () {
        const { viewer } = this
        const v = this.volume.getValue()
        if (!viewer || !v) return

        if (this.surfaceRepr) this.viewer.remove(this.surfaceRepr)

        this.surfaceRepr = VolumeRepresentation(Surface)
        await Run(this.surfaceRepr.create(v.volume, {
            isoValue: VolumeIsoValue.relative(v.volume.dataStats, 3.0),
            alpha: 1.0
        }), log, 500)
        viewer.add(this.surfaceRepr)

        viewer.requestDraw()
        console.log(viewer.stats)
    }

    async loadPdbId () {
        if (this.pointRepr) this.viewer.remove(this.pointRepr)
        if (this.spacefillRepr) this.viewer.remove(this.spacefillRepr)
        if (this.pdbId.length !== 4) return
        this.loading.next(true)
        this.setModel((await getModelFromPdbId(this.pdbId))[0])
    }

    setVolume(volume: Volume) {
        this.volume.next(volume)
        this.initVolume()
        this.loading.next(false)
    }

    async loadEmdId () {
        if (this.surfaceRepr) this.viewer.remove(this.surfaceRepr)
        if (this.emdId.length !== 4) return
        this.loading.next(true)
        this.setVolume(await getVolumeFromEmdId(this.emdId))
    }

    async update () {
        if (!this.spacefillRepr) return
        await Run(this.spacefillRepr.update(this.getSpacefillProps()), log, 500)
        await Run(this.pointRepr.update(this.getPointProps()), log, 500)
        this.viewer.add(this.spacefillRepr)
        this.viewer.add(this.pointRepr)
        this.viewer.update()
        this.viewer.requestDraw()
        console.log(this.viewer.stats)
    }

    updateVisibility () {
        if (!this.viewer) return
        if (this.pointRepr) {
            if (this.pointVisibility.getValue()) {
                this.viewer.show(this.pointRepr)
            } else {
                this.viewer.hide(this.pointRepr)
            }
        }
        if (this.spacefillRepr) {
            if (this.spacefillVisibility.getValue()) {
                this.viewer.show(this.spacefillRepr)
            } else {
                this.viewer.hide(this.spacefillRepr)
            }
        }
        this.viewer.requestDraw()
    }
}