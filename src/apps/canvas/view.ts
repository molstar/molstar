/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, Structure } from "mol-model/structure";
import { CartoonRepresentation } from 'mol-geo/representation/structure/representation/cartoon';
import { BallAndStickRepresentation } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { getStructureFromModel } from "./util";

export interface DefaultView {
    readonly model: Model
    readonly structure: Structure
    readonly cartoon: CartoonRepresentation
    readonly ballAndStick: BallAndStickRepresentation
    
    setAssembly(assembly: string): void
}

export async function DefaultView(model: Model, assembly: string): Promise<DefaultView> {
    const cartoon = CartoonRepresentation()
    const ballAndStick = BallAndStickRepresentation()
    
    let structure: Structure

    async function setAssembly(assembly: string) {
        structure = await getStructureFromModel(model, assembly)
        await createRepr()
    }

    async function createRepr() {
        await cartoon.create(structure, {
            colorTheme: { name: 'chain-id' },
            sizeTheme: { name: 'uniform', value: 0.2 },
            useFog: false // TODO fog not working properly
        }).run()

        await ballAndStick.create(structure, {
            colorTheme: { name: 'element-symbol' },
            sizeTheme: { name: 'uniform', value: 0.1 },
            useFog: false // TODO fog not working properly
        }).run()
    }

    await setAssembly(assembly)

    return {
        model,
        get structure() { return structure },
        cartoon,
        ballAndStick,

        setAssembly
    }
}

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