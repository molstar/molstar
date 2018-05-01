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
import { getModelFromPdbId, getModelFromFile, log } from './utils'
import { StructureRepresentation } from 'mol-geo/representation/structure';
import { Color } from 'mol-util/color';
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
    pdbId = '4cup'
    model = new BehaviorSubject<Model | undefined>(undefined)
    initialized = new BehaviorSubject<boolean>(false)
    loading = new BehaviorSubject<boolean>(false)

    colorTheme = new BehaviorSubject<ColorTheme>('element-symbol')
    colorValue = new BehaviorSubject<Color>(0xFF4411)
    detail = new BehaviorSubject<number>(0)
    assembly = new BehaviorSubject<string>('')

    pointVisibility = new BehaviorSubject<boolean>(true)
    spacefillVisibility = new BehaviorSubject<boolean>(true)

    pointRepr: StructureRepresentation<PointProps>
    spacefillRepr: StructureRepresentation<SpacefillProps>

    constructor() {
        this.colorTheme.subscribe(() => this.update())
        this.colorValue.subscribe(() => this.update())
        this.detail.subscribe(() => this.update())
        this.assembly.subscribe(() => this.initStructure())

        this.pointVisibility.subscribe(() => this.updateVisibility())
        this.spacefillVisibility.subscribe(() => this.updateVisibility())
    }

    getSpacefillProps (): SpacefillProps {
        const colorThemeName = this.colorTheme.getValue()
        return {
            detail: this.detail.getValue(),
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
        this.viewer.animate()
    }

    async getStructure () {
        const model = this.model.getValue()
        if (!model) return
        const assembly = this.assembly.getValue()
        let structure: Structure
        const assemblies = model.symmetry.assemblies
        if (assemblies.length) {
            structure = await Run(Symmetry.buildAssembly(Structure.ofModel(model), assembly || '1'), log, 100)
        } else {
            structure = Structure.ofModel(model)
        }
        return structure
    }

    async initStructure () {
        const { viewer, model } = this
        if (!viewer || !model) return

        viewer.clear()

        const structure = await this.getStructure()
        if (!structure) return

        this.pointRepr = StructureRepresentation(Point)
        await Run(this.pointRepr.create(structure, this.getPointProps()), log, 100)
        viewer.add(this.pointRepr)

        this.spacefillRepr = StructureRepresentation(Spacefill)
        await Run(this.spacefillRepr.create(structure, this.getSpacefillProps()), log, 100)
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

    async loadPdbId () {
        this.viewer.clear()
        if (this.pdbId.length !== 4) return
        this.loading.next(true)
        this.setModel((await getModelFromPdbId(this.pdbId))[0])
    }

    async update () {
        if (!this.spacefillRepr) return
        await Run(this.spacefillRepr.update(this.getSpacefillProps()), log, 100)
        await Run(this.pointRepr.update(this.getPointProps()), log, 100)
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