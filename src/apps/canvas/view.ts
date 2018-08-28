/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, Structure } from 'mol-model/structure';
import { CartoonRepresentation } from 'mol-geo/representation/structure/representation/cartoon';
// import { BallAndStickRepresentation } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { getStructureFromModel } from './util';
import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';
import { ShapeRepresentation, ShapeProps } from 'mol-geo/representation/shape';
import { getAxesShape } from './assembly-symmetry';
import Viewer from 'mol-view/viewer';

export interface StructureView {
    readonly label: string
    readonly models: ReadonlyArray<Model>
    readonly structure: Structure | undefined
    readonly assemblySymmetry: AssemblySymmetry | undefined

    readonly cartoon: CartoonRepresentation
    // readonly ballAndStick: BallAndStickRepresentation
    readonly axes: ShapeRepresentation<ShapeProps>

    readonly modelId: number
    readonly assemblyId: string
    readonly symmetryFeatureId: number

    setAssembly(assembly: string): Promise<void>
    getAssemblyIds(): { id: string, label: string }[]
    setSymmetryFeature(symmetryFeature: number): Promise<void>
    getSymmetryFeatureIds(): { id: number, label: string }[]

    destroy: () => void
}

interface StructureViewProps {
    assembly?: string
    symmetryFeature?: number
}

export async function StructureView(viewer: Viewer, models: ReadonlyArray<Model>, props: StructureViewProps): Promise<StructureView> {
    const cartoon = CartoonRepresentation()
    // const ballAndStick = BallAndStickRepresentation()
    const axes = ShapeRepresentation()

    let label: string
    let model: Model | undefined
    let assemblySymmetry: AssemblySymmetry | undefined
    let structure: Structure | undefined

    let modelId: number
    let assemblyId: string
    let symmetryFeatureId: number

    async function setModel(newModelId: number, newAssemblyId?: string, newSymmetryFeatureId?: number) {
        modelId = newModelId

        model = models[modelId]
        await AssemblySymmetry.attachFromCifOrAPI(model)
        assemblySymmetry = AssemblySymmetry.get(model)

        await setAssembly(newAssemblyId, newSymmetryFeatureId)
    }

    async function setAssembly(newAssemblyId?: string, newSymmetryFeatureId?: number) {
        if (newAssemblyId !== undefined) {
            assemblyId = newAssemblyId
        } else if (model && model.symmetry.assemblies.length) {
            assemblyId = model.symmetry.assemblies[0].id
        } else if (model) {
            assemblyId = '0'
        } else {
            assemblyId = '-1'
        }
        await getStructure()
        await createStructureRepr()
        await setSymmetryFeature(newSymmetryFeatureId)
    }

    function getAssemblyIds() {
        const assemblyIds: { id: string, label: string }[] = [
            { id: '0', label: '0: model' }
        ]
        if (model) model.symmetry.assemblies.forEach(a => {
            assemblyIds.push({ id: a.id, label: `${a.id}: ${a.details}` })
        })
        return assemblyIds
    }

    async function setSymmetryFeature(newSymmetryFeatureId?: number) {
        if (newSymmetryFeatureId !== undefined) {
            symmetryFeatureId = newSymmetryFeatureId
        } else if (assemblySymmetry) {
            const f = assemblySymmetry.getFeatures(assemblyId)
            if (f._rowCount) {
                symmetryFeatureId = f.id.value(0)
            } else {
                symmetryFeatureId = -1
            }
        } else {
            symmetryFeatureId = -1
        }
        await createSymmetryRepr()
    }

    function getSymmetryFeatureIds() {
        const symmetryFeatureIds: { id: number, label: string }[] = []
        if (assemblySymmetry) {
            const symmetryFeatures = assemblySymmetry.getFeatures(assemblyId)
            for (let i = 0, il = symmetryFeatures._rowCount; i < il; ++i) {
                const id = symmetryFeatures.id.value(i)
                const symmetry = symmetryFeatures.symmetry_value.value(i)
                const type = symmetryFeatures.type.value(i)
                const stoichiometry = symmetryFeatures.stoichiometry_value.value(i)
                const label = `${id}: ${symmetry} ${type} ${stoichiometry}`
                symmetryFeatureIds.push({ id, label })
            }
        }
        return symmetryFeatureIds
    }

    async function getStructure() {
        if (model) structure = await getStructureFromModel(model, assemblyId)
        if (model && structure) {
            label = `${model.label} - Assembly ${assemblyId}`
        } else {
            label = ''
        }
        await createStructureRepr()
    }

    async function createStructureRepr() {
        if (structure) {
            await cartoon.create(structure, {
                colorTheme: { name: 'chain-id' },
                sizeTheme: { name: 'uniform', value: 0.2 },
                useFog: false // TODO fog not working properly
            }).run()

            // await ballAndStick.create(structure, {
            //     colorTheme: { name: 'element-symbol' },
            //     sizeTheme: { name: 'uniform', value: 0.1 },
            //     useFog: false // TODO fog not working properly
            // }).run()

            viewer.center(structure.boundary.sphere.center)
        } else {
            cartoon.destroy()
            // ballAndStick.destroy()
        }

        viewer.add(cartoon)
        // viewer.add(ballAndStick)
    }

    async function createSymmetryRepr() {
        if (assemblySymmetry) {
            const features = assemblySymmetry.getFeatures(assemblyId)
            if (features._rowCount) {
                const axesShape = getAxesShape(symmetryFeatureId, assemblySymmetry)
                if (axesShape) {
                    // getClusterColorTheme(symmetryFeatureId, assemblySymmetry)
                    await axes.create(axesShape, {
                        colorTheme: { name: 'shape-group' },
                        // colorTheme: { name: 'uniform', value: Color(0xFFCC22) },
                        useFog: false // TODO fog not working properly
                    }).run()
                } else {
                    axes.destroy()
                }
            } else {
                axes.destroy()
            }
        } else {
            axes.destroy()
        }
        viewer.add(axes)
        viewer.requestDraw()
    }

    await setModel(0, props.assembly, props.symmetryFeature)

    return {
        get label() { return label },
        models,
        get structure() { return structure },
        get assemblySymmetry() { return assemblySymmetry },

        cartoon,
        // ballAndStick,
        axes,

        get modelId() { return modelId },
        get assemblyId() { return assemblyId },
        get symmetryFeatureId() { return symmetryFeatureId },

        setAssembly,
        getAssemblyIds,
        setSymmetryFeature,
        getSymmetryFeatureIds,

        destroy: () => {
            viewer.remove(cartoon)
            // viewer.remove(ballAndStick)
            viewer.remove(axes)
            viewer.requestDraw()

            cartoon.destroy()
            // ballAndStick.destroy()
            axes.destroy()
        }
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