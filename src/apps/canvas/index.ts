/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'

import Viewer from 'mol-view/viewer';
import CIF, { CifBlock } from 'mol-io/reader/cif'
import { readUrlAs } from 'mol-util/read'
import { Model, Format, Structure, StructureSymmetry } from 'mol-model/structure';
import { CartoonRepresentation } from 'mol-geo/representation/structure/representation/cartoon';
import { BallAndStickRepresentation } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { EveryLoci } from 'mol-model/loci';
import { MarkerAction } from 'mol-geo/util/marker-data';
import { labelFirst } from 'mol-view/label';
import { Queries as Q, StructureProperties as SP, StructureSelection, StructureQuery } from 'mol-model/structure';
import { MeshBuilder } from 'mol-geo/shape/mesh-builder';
import { CustomRepresentation } from 'mol-geo/representation/custom';
import { Vec3 } from 'mol-math/linear-algebra';

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

async function init() {
    const cif = await getCifFromUrl('https://files.rcsb.org/download/1crn.cif')
    const model = await getModelFromMmcif(cif)
    const structure = await getStructureFromModel(model)

    viewer.center(structure.boundary.sphere.center)

    // cartoon for whole structure
    const cartoonRepr = CartoonRepresentation()
    await cartoonRepr.create(structure, {
        colorTheme: { name: 'chain-id' },
        sizeTheme: { name: 'uniform', value: 0.2 },
        useFog: false // TODO fog not working properly
    }).run()
    viewer.add(cartoonRepr)

    // create new structure via query
    const q1 = Q.generators.atoms({
        residueTest: qtx => SP.residue.label_seq_id(qtx.element) < 7
    });
    const newStructure = StructureSelection.unionStructure(await StructureQuery.run(q1, structure));

    // ball+stick for new structure
    const ballStickRepr = BallAndStickRepresentation()
    await ballStickRepr.create(newStructure, {
        colorTheme: { name: 'element-symbol' },
        sizeTheme: { name: 'uniform', value: 0.1 },
        useFog: false // TODO fog not working properly
    }).run()
    viewer.add(ballStickRepr)

    // create a custom mesh
    const meshBuilder = MeshBuilder.create(256, 128)
    meshBuilder.setId(0)
    meshBuilder.addSphere(Vec3.create(0, 0, 0), 4, 2)
    const mesh = meshBuilder.getMesh()
    // Mesh.computeNormalsImmediate(mesh)

    // add representation from custom mesh
    const customRepr = CustomRepresentation()
    await customRepr.create(mesh, {
        useFog: false // TODO fog not working properly
    }).run()
    viewer.add(customRepr)

    // ensure the added representations get rendered, i.e. without mouse input
    viewer.requestDraw()
}

init()