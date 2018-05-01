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
import { Symmetry, Structure } from 'mol-model/structure'

// import mcubes from './utils/mcubes'
import { getStructuresFromPdbId, getStructuresFromFile, log } from './utils'
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
    initialized = new BehaviorSubject<boolean>(false)
    loading = new BehaviorSubject<boolean>(false)

    colorTheme = new BehaviorSubject<ColorTheme>('element-symbol')
    colorValue = new BehaviorSubject<Color>(0xFF4411)
    detail = new BehaviorSubject<number>(0)

    pointVisibility = new BehaviorSubject<boolean>(true)
    spacefillVisibility = new BehaviorSubject<boolean>(true)

    pointRepr: StructureRepresentation<PointProps>
    spacefillRepr: StructureRepresentation<SpacefillProps>

    constructor() {
        this.colorTheme.subscribe(() => this.update())
        this.colorValue.subscribe(() => this.update())
        this.detail.subscribe(() => this.update())

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

    async initStructure (structure: Structure) {
        const { viewer, loading } = this
        viewer.clear()

        const struct = await Run(Symmetry.buildAssembly(structure, '1'), log, 100)

        this.pointRepr = StructureRepresentation(Point)
        await Run(this.pointRepr.create(struct, this.getPointProps()), log, 100)
        viewer.add(this.pointRepr)

        this.spacefillRepr = StructureRepresentation(Spacefill)
        await Run(this.spacefillRepr.create(struct, this.getSpacefillProps()), log, 100)
        viewer.add(this.spacefillRepr)

        this.updateVisibility()
        viewer.requestDraw()
        console.log(viewer.stats)

        loading.next(false)
    }

    async loadFile (file: File) {
        this.viewer.clear()
        this.loading.next(true)

        const structures = await getStructuresFromFile(file)
        this.initStructure(structures[0])
    }

    async loadPdbId () {
        this.viewer.clear()
        if (this.pdbId.length !== 4) return
        this.loading.next(true)

        const structures = await getStructuresFromPdbId(this.pdbId)
        this.initStructure(structures[0])
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



// async foo () {
//     const p1 = Vec3.create(0, 4, 0)
//     const p2 = Vec3.create(-3, 0, 0)

//     // const position = ValueCell.create(new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]))
//     // const normal = ValueCell.create(new Float32Array([0, 0, 0, 0, 0, 0, 0, 0, 0]))

//     const transformArray1 = ValueCell.create(new Float32Array(16))
//     const transformArray2 = ValueCell.create(new Float32Array(16 * 3))
//     const m4 = Mat4.identity()
//     Mat4.toArray(m4, transformArray1.ref.value, 0)
//     Mat4.toArray(m4, transformArray2.ref.value, 0)
//     Mat4.setTranslation(m4, p1)
//     Mat4.toArray(m4, transformArray2.ref.value, 16)
//     Mat4.setTranslation(m4, p2)
//     Mat4.toArray(m4, transformArray2.ref.value, 32)

//     const color = ValueCell.create(createColorTexture(3))
//     color.ref.value.set([
//         0, 0, 255,
//         0, 255, 0,
//         255, 0, 0
//     ])

//     // const points = createRenderObject('point', {
//     //     position,
//     //     transform: transformArray1
//     // })
//     // // renderer.add(points)

//     // const mesh = createRenderObject('mesh', {
//     //     position,
//     //     normal,
//     //     color,
//     //     transform: transformArray2
//     // })
//     // renderer.add(mesh)

//     // const cylinder = Cylinder({ height: 3, radiusBottom: 0.5, radiusTop: 0.5 })
//     // console.log(cylinder)
//     // const cylinderMesh = createRenderObject('mesh', {
//     //     position: ValueCell.create(cylinder.vertices),
//     //     normal: ValueCell.create(cylinder.normals),
//     //     color,
//     //     transform: transformArray2
//     // }, cylinder.indices)
//     // renderer.add(cylinderMesh)

//     // const sphere = Icosahedron()
//     // console.log(sphere)

//     // const box = Box()
//     // console.log(box)

//     // const points2 = createRenderObject('point', {
//     //     position: ValueCell.create(new Float32Array(box.vertices)),
//     //     transform: transformArray1
//     // })
//     // renderer.add(points2)

//     // let rr = 0.7;
//     // function cubesF(x: number, y: number, z: number) {
//     //     return x * x + y * y + z * z - rr * rr;
//     // }
//     // let cubes = await mcubes(cubesF);

//     // const makeCubesMesh = () => createRenderObject('mesh', {
//     //     position: cubes.surface.vertexBuffer,
//     //     normal: cubes.surface.normalBuffer,
//     //     color,
//     //     transform: transformArray2,
//     //     elements: cubes.surface.indexBuffer,

//     //     instanceCount: transformArray2.ref.value.length / 16,
//     //     elementCount: cubes.surface.triangleCount,
//     //     positionCount: cubes.surface.vertexCount
//     // }, {});
//     // const mesh2 = makeCubesMesh();
//     // renderer.add(mesh2)
// }