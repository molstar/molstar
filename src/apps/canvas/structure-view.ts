/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model, Structure } from 'mol-model/structure';
import { getStructureFromModel } from './util';
import { AssemblySymmetry } from 'mol-model-props/rcsb/symmetry';
import { getAxesShape } from './assembly-symmetry';
import Canvas3D from 'mol-canvas3d/canvas3d';
// import { MeshBuilder } from 'mol-geo/mesh/mesh-builder';
// import { addSphere } from 'mol-geo/mesh/builder/sphere';
// import { Shape } from 'mol-model/shape';
// import { Color } from 'mol-util/color';
// import { computeUnitBoundary } from 'mol-model/structure/structure/util/boundary';
// import { addBoundingBox } from 'mol-geo/mesh/builder/bounding-box';
import { BehaviorSubject } from 'rxjs';
import { App } from './app';
import { StructureRepresentation } from 'mol-repr/structure/representation';
import { ShapeRepresentation, ShapeProps } from 'mol-repr/shape/representation';

export interface StructureView {
    readonly app: App
    readonly canvas3d: Canvas3D

    readonly label: string
    readonly models: ReadonlyArray<Model>
    readonly structure: Structure | undefined
    readonly assemblySymmetry: AssemblySymmetry | undefined

    readonly active: { [k: string]: boolean }
    readonly structureRepresentations: { [k: string]: StructureRepresentation<any> }
    readonly updated: BehaviorSubject<null>
    readonly symmetryAxes: ShapeRepresentation<ShapeProps>

    setSymmetryAxes(value: boolean): void
    setStructureRepresentation(name: string, value: boolean): void

    readonly modelId: number
    readonly assemblyId: string
    readonly symmetryFeatureId: number

    setModel(modelId: number): Promise<void>
    getModelIds(): { id: number, label: string }[]
    setAssembly(assemblyId: string): Promise<void>
    getAssemblyIds(): { id: string, label: string }[]
    setSymmetryFeature(symmetryFeatureId: number): Promise<void>
    getSymmetryFeatureIds(): { id: number, label: string }[]

    destroy: () => void
}

interface StructureViewProps {
    assemblyId?: string
    symmetryFeatureId?: number
}

export async function StructureView(app: App, canvas3d: Canvas3D, models: ReadonlyArray<Model>, props: StructureViewProps = {}): Promise<StructureView> {
    const active: { [k: string]: boolean } = {
        'cartoon': true,
        'ball-and-stick': true,
        // point: false,
        // surface: false,
        // carbohydrate: false,
        // spacefill: false,
        // distanceRestraint: false,
        // symmetryAxes: true,
        // polymerSphere: false,
    }

    const structureRepresentations: { [k: string]: StructureRepresentation<any> } = {}

    const symmetryAxes = ShapeRepresentation()
    const polymerSphere = ShapeRepresentation()

    const updated: BehaviorSubject<null> = new BehaviorSubject<null>(null)

    let label: string
    let model: Model | undefined
    let assemblySymmetry: AssemblySymmetry | undefined
    let structure: Structure | undefined

    let modelId: number
    let assemblyId: string
    let symmetryFeatureId: number

    async function setSymmetryAxes(value: boolean) {
        if (!value) {
            assemblySymmetry = undefined
        } else {
            await app.runTask(AssemblySymmetry.attachFromCifOrAPI(models[modelId]), 'Load symmetry annotation')
            assemblySymmetry = AssemblySymmetry.get(models[modelId])
        }
        active.symmetryAxes = value
        await setSymmetryFeature()
    }

    async function setStructureRepresentation(k: string, value: boolean) {
        active[k] = value
        await createStructureRepr()
    }

    async function setModel(newModelId: number, newAssemblyId?: string, newSymmetryFeatureId?: number) {
        console.log('setModel', newModelId)
        modelId = newModelId
        model = models[modelId]
        if (active.symmetryAxes) {
            await AssemblySymmetry.attachFromCifOrAPI(model)
            assemblySymmetry = AssemblySymmetry.get(model)
        }
        await setAssembly(newAssemblyId, newSymmetryFeatureId)
    }

    function getModelIds() {
        const modelIds: { id: number, label: string }[] = []
        models.forEach((m, i) => {
            modelIds.push({ id: i, label: `${i}: ${m.label} #${m.modelNum}` })
        })
        return modelIds
    }

    async function setAssembly(newAssemblyId?: string, newSymmetryFeatureId?: number) {
        console.log('setAssembly', newAssemblyId)
        if (newAssemblyId !== undefined) {
            assemblyId = newAssemblyId
        } else if (model && model.symmetry.assemblies.length) {
            assemblyId = model.symmetry.assemblies[0].id
        } else if (model) {
            assemblyId = 'deposited'
        } else {
            assemblyId = '-1'
        }
        await getStructure()
        await setSymmetryFeature(newSymmetryFeatureId)
    }

    function getAssemblyIds() {
        const assemblyIds: { id: string, label: string }[] = [
            { id: 'deposited', label: 'deposited' }
        ]
        if (model) model.symmetry.assemblies.forEach(a => {
            assemblyIds.push({ id: a.id, label: `${a.id}: ${a.details}` })
        })
        return assemblyIds
    }

    async function setSymmetryFeature(newSymmetryFeatureId?: number) {
        console.log('setSymmetryFeature', newSymmetryFeatureId)
        if (newSymmetryFeatureId !== undefined) {
            symmetryFeatureId = newSymmetryFeatureId
        } else if (assemblySymmetry) {
            const s = assemblySymmetry.getSymmetries(assemblyId)
            if (s._rowCount) {
                symmetryFeatureId = s.id.value(0)
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
            const symmetries = assemblySymmetry.getSymmetries(assemblyId)
            for (let i = 0, il = symmetries._rowCount; i < il; ++i) {
                const id = symmetries.id.value(i)
                const kind = symmetries.kind.value(i)
                const type = symmetries.type.value(i)
                const stoichiometry = symmetries.stoichiometry.value(i)
                const label = `${id}: ${kind} ${type} ${stoichiometry}`
                symmetryFeatureIds.push({ id, label })
            }
        }
        return symmetryFeatureIds
    }

    async function getStructure() {
        if (model) structure = await app.runTask(getStructureFromModel(model, assemblyId), 'Build structure')
        if (model && structure) {
            label = `${model.label} - Assembly ${assemblyId}`
        } else {
            label = ''
        }
        await createStructureRepr()
    }

    async function createStructureRepr() {
        if (structure) {
            console.log('createStructureRepr')
            for (const k in active) {
                if (active[k]) {
                    let repr: StructureRepresentation
                    if (structureRepresentations[k]) {
                        repr = structureRepresentations[k]
                    } else {
                        repr = app.structureRepresentationRegistry.create(k, app.reprCtx, structure)
                        structureRepresentations[k] = repr
                    }
                    await app.runTask(repr.createOrUpdate(app.reprCtx, {}, {}, structure).run(
                        progress => app.log(progress)
                    ), 'Create/update representation')
                    canvas3d.add(repr)
                } else {
                    if (structureRepresentations[k]) {
                        canvas3d.remove(structureRepresentations[k])
                        structureRepresentations[k].destroy()
                        delete structureRepresentations[k]
                    }
                }
            }

            canvas3d.center(structure.boundary.sphere.center)

            // const mb = MeshBuilder.create()
            // mb.setGroup(0)
            // addSphere(mb, structure.boundary.sphere.center, structure.boundary.sphere.radius, 3)
            // addBoundingBox(mb, structure.boundary.box, 1, 2, 8)
            // for (let i = 0, il = structure.units.length; i < il; ++i) {
            //     mb.setGroup(1)
            //     const u = structure.units[i]
            //     const ci = u.model.atomicHierarchy.chainAtomSegments.index[u.elements[0]]
            //     const ek = u.model.atomicHierarchy.getEntityKey(ci)
            //     if (u.model.entities.data.type.value(ek) === 'water') continue
            //     const boundary = computeUnitBoundary(u)
            //     addSphere(mb, boundary.sphere.center, boundary.sphere.radius, 3)
            //     addBoundingBox(mb, boundary.box, 0.5, 2, 8)
            // }
            // const shape = Shape.create('boundary', mb.getMesh(), [Color(0xCC6633), Color(0x3366CC)], ['sphere boundary'])
            // await polymerSphere.createOrUpdate({
            //     alpha: 0.5,
            //     doubleSided: false,
            //     depthMask: false,
            //     useFog: false // TODO fog not working properly
            // }, shape).run()
        } else {
            for (const k in structureRepresentations) structureRepresentations[k].destroy()
            polymerSphere.destroy()
        }

        canvas3d.add(polymerSphere)

        updated.next(null)
        canvas3d.requestDraw(true)
        console.log('stats', canvas3d.stats)
    }

    async function createSymmetryRepr() {
        if (assemblySymmetry) {
            const symmetries = assemblySymmetry.getSymmetries(assemblyId)
            if (symmetries._rowCount) {
                const axesShape = getAxesShape(symmetryFeatureId, assemblySymmetry)
                if (axesShape) {
                    // const colorTheme = getClusterColorTheme(symmetryFeatureId, assemblySymmetry)
                    // await structureRepresentations['cartoon'].createOrUpdate({
                    //     colorTheme: 'custom',
                    //     colorFunction: colorTheme.color,
                    //     colorGranularity: colorTheme.granularity,
                    // }).run()
                    await symmetryAxes.createOrUpdate(app.reprCtx, {}, {}, axesShape).run()
                    canvas3d.add(symmetryAxes)
                } else {
                    canvas3d.remove(symmetryAxes)
                }
            } else {
                canvas3d.remove(symmetryAxes)
            }
        } else {
            canvas3d.remove(symmetryAxes)
        }
        updated.next(null)
        canvas3d.requestDraw(true)
    }

    await setModel(0, props.assemblyId, props.symmetryFeatureId)

    return {
        app,
        canvas3d,

        get label() { return label },
        models,
        get structure() { return structure },
        get assemblySymmetry() { return assemblySymmetry },

        active,
        structureRepresentations,
        updated,
        symmetryAxes,

        setSymmetryAxes,
        setStructureRepresentation,

        get modelId() { return modelId },
        get assemblyId() { return assemblyId },
        get symmetryFeatureId() { return symmetryFeatureId },

        setModel,
        getModelIds,
        setAssembly,
        getAssemblyIds,
        setSymmetryFeature,
        getSymmetryFeatureIds,

        destroy: () => {
            for (const k in structureRepresentations) {
                canvas3d.remove(structureRepresentations[k])
                structureRepresentations[k].destroy()
            }
            canvas3d.remove(polymerSphere)
            canvas3d.remove(symmetryAxes)
            canvas3d.requestDraw(true)

            polymerSphere.destroy()
            symmetryAxes.destroy()
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