/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'

import Viewer from 'mol-view/viewer';
import CIF, { CifBlock } from 'mol-io/reader/cif'
// import { parse as parseObj } from 'mol-io/reader/obj/parser'
import { readUrlAs } from 'mol-util/read'
import { Model, Format, Structure, StructureSymmetry } from 'mol-model/structure';
import { CartoonRepresentation } from 'mol-geo/representation/structure/representation/cartoon';
import { BallAndStickRepresentation } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { EveryLoci } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/util/marker-data';
import { labelFirst } from 'mol-view/label';
// import { Queries as Q, StructureProperties as SP, StructureSelection, StructureQuery } from 'mol-model/structure';
import { MeshBuilder } from 'mol-geo/mesh/mesh-builder';
import { ShapeRepresentation } from 'mol-geo/representation/shape';
import { /*Vec3, Mat4,*/ Tensor } from 'mol-math/linear-algebra';
import { Shape } from 'mol-model/shape';
import { Color } from 'mol-util/color';
import { addSphere } from 'mol-geo/mesh/builder/sphere';
// import { Box } from 'mol-geo/primitive/box';
import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';
import { addCylinder } from 'mol-geo/mesh/builder/cylinder';
import { Table } from 'mol-data/db';

const container = document.getElementById('container')
if (!container) throw new Error('Can not find element with id "container".')

const canvas = document.getElementById('canvas') as HTMLCanvasElement
if (!canvas) throw new Error('Can not find element with id "canvas".')

const info = document.getElementById('info') as HTMLCanvasElement
if (!info) throw new Error('Can not find element with id "info".')

const viewer = Viewer.create(canvas, container)
viewer.animate()

viewer.input.resize.subscribe(() => {
    // do whatever appropriate
})

viewer.input.move.subscribe(({x, y, inside, buttons}) => {
    if (!inside || buttons) return
    const p = viewer.identify(x, y)
    const loci = viewer.getLoci(p)

    viewer.mark(EveryLoci, MarkerAction.RemoveHighlight)
    viewer.mark(loci, MarkerAction.Highlight)

    const label = labelFirst(loci)
    info.innerText = `${label}`
})


// async function getObjFromUrl(url: string) {
//     const data = await readUrlAs(url, false) as string
//     const comp = parseObj(data)
//     const parsed = await comp.run()
//     if (parsed.isError) throw parsed
//     return parsed.result
// }

async function getCifFromUrl(url: string) {
    const data = await readUrlAs(url, false)
    const comp = CIF.parse(data)
    const parsed = await comp.run()
    if (parsed.isError) throw parsed
    return parsed.result.blocks[0]
}

async function getModelFromMmcif(cif: CifBlock) {
    const models = await Model.create(Format.mmCIF(cif)).run()
    return models[0]
}

async function getStructureFromModel(model: Model, assembly = '1') {
    const assemblies = model.symmetry.assemblies
    if (assemblies.length) {
        return await StructureSymmetry.buildAssembly(Structure.ofModel(model), assembly).run()
    } else {
        return Structure.ofModel(model)
    }
}

function getAxesShape(featureId: number, assemblySymmetry: AssemblySymmetry) {
    const f = assemblySymmetry.db.rcsb_assembly_symmetry_feature
    const feature = Table.pickRow(f, i => f.id.value(i) === featureId)
    if (!feature) return

    const axes = assemblySymmetry.getAxes(featureId)
    if (!axes._rowCount) return

    const vectorSpace = AssemblySymmetry.Schema.rcsb_assembly_symmetry_axis.start.space;

    const colors: Color[] = []
    const labels: string[] = []

    const radius = 0.4
    const cylinderProps = { radiusTop: radius, radiusBottom: radius }
    const meshBuilder = MeshBuilder.create(256, 128)

    for (let i = 0, il = axes._rowCount; i < il; ++i) {
        const start = Tensor.toVec3(vectorSpace, axes.start.value(i))
        const end = Tensor.toVec3(vectorSpace, axes.end.value(i))
        meshBuilder.setGroup(i)
        addSphere(meshBuilder, start, radius, 2)
        addSphere(meshBuilder, end, radius, 2)
        addCylinder(meshBuilder, start, end, 1, cylinderProps)
        colors.push(Color(0xCCEE11))
        labels.push(`Axis ${i + 1} for ${feature.symmetry_value} ${feature.type.toLowerCase()} symmetry`)
    }
    const mesh = meshBuilder.getMesh()
    const shape = Shape.create('Axes', mesh, colors, labels)
    return shape
}

function getClusterColorTheme(featureId: number, assemblySymmetry: AssemblySymmetry) {
    const f = assemblySymmetry.db.rcsb_assembly_symmetry_feature
    const feature = Table.pickRow(f, i => f.id.value(i) === featureId)
    if (!feature) return

    const clusters = assemblySymmetry.getClusters(featureId)
    if (!clusters._rowCount) return

    for (let i = 0, il = clusters._rowCount; i < il; ++i) {
        console.log(clusters.members.value(i), clusters.avg_rmsd.value(i), feature.stoichiometry_value, feature.stoichiometry_description)
    }
}

async function init() {
    const assembly = '1'

    const cif = await getCifFromUrl('https://files.rcsb.org/download/4hhb.cif')
    const model = await getModelFromMmcif(cif)
    const structure = await getStructureFromModel(model, assembly)
    viewer.center(structure.boundary.sphere.center)

    // cartoon for whole structure
    const cartoonRepr = CartoonRepresentation()
    await cartoonRepr.create(structure, {
        colorTheme: { name: 'chain-id' },
        sizeTheme: { name: 'uniform', value: 0.2 },
        useFog: false // TODO fog not working properly
    }).run()
    viewer.add(cartoonRepr)

    // ball+stick for whole structure
    const ballStickRepr = BallAndStickRepresentation()
    await ballStickRepr.create(structure, {
        colorTheme: { name: 'element-symbol' },
        sizeTheme: { name: 'uniform', value: 0.1 },
        useFog: false // TODO fog not working properly
    }).run()
    viewer.add(ballStickRepr)

    // // create new structure via query
    // const q1 = Q.generators.atoms({
    //     residueTest: qtx => SP.residue.label_seq_id(qtx.element) < 7
    // });
    // const newStructure = StructureSelection.unionStructure(await StructureQuery.run(q1, structure));

    // // ball+stick for new structure
    // const newBallStickRepr = BallAndStickRepresentation()
    // await newBallStickRepr.create(newStructure, {
    //     colorTheme: { name: 'element-symbol' },
    //     sizeTheme: { name: 'uniform', value: 0.1 },
    //     useFog: false // TODO fog not working properly
    // }).run()
    // viewer.add(newBallStickRepr)

    // // create a mesh
    // const meshBuilder = MeshBuilder.create(256, 128)
    // const colors: Color[] = []
    // const labels: string[] = []
    // // red sphere
    // meshBuilder.setGroup(0)
    // colors[0] = Color(0xFF2233)
    // labels[0] = 'red sphere'
    // addSphere(meshBuilder, Vec3.create(0, 0, 0), 4, 2)
    // // green cube
    // meshBuilder.setGroup(1)
    // colors[1] = Color(0x2233FF)
    // labels[1] = 'blue cube'
    // const t = Mat4.identity()
    // Mat4.fromTranslation(t, Vec3.create(10, 0, 0))
    // Mat4.scale(t, t, Vec3.create(3, 3, 3))
    // meshBuilder.add(t, Box())
    // const mesh = meshBuilder.getMesh()
    // const mesh = getObjFromUrl('mesh.obj')

    // // create shape from mesh
    // const shape = Shape.create('myShape', mesh, colors, labels)

    // // add representation from shape
    // const customRepr = ShapeRepresentation()
    // await customRepr.create(shape, {
    //     colorTheme: { name: 'shape-group' },
    //     // colorTheme: { name: 'uniform', value: Color(0xFFCC22) },
    //     useFog: false // TODO fog not working properly
    // }).run()
    // viewer.add(customRepr)


    await AssemblySymmetry.attachFromCifOrAPI(model)
    const assemblySymmetry = AssemblySymmetry.get(model)
    console.log(assemblySymmetry)
    if (assemblySymmetry) {
        const features = assemblySymmetry.getFeatures(assembly)
        if (features._rowCount) {
            const axesShape = getAxesShape(features.id.value(1), assemblySymmetry)
            console.log(axesShape)
            if (axesShape) {
                const customRepr = ShapeRepresentation()
                await customRepr.create(axesShape, {
                    colorTheme: { name: 'shape-group' },
                    // colorTheme: { name: 'uniform', value: Color(0xFFCC22) },
                    useFog: false // TODO fog not working properly
                }).run()
                viewer.add(customRepr)
            }

            getClusterColorTheme(features.id.value(0), assemblySymmetry)
            getClusterColorTheme(features.id.value(1), assemblySymmetry)
        }
    }

    // ensure the added representations get rendered, i.e. without mouse input
    viewer.requestDraw()
}

init()