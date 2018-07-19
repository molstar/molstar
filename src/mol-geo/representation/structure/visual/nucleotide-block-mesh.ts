/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from './util/common';
import { deepEqual } from 'mol-util';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getElementLoci, markElement } from './util/element';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Segmentation, SortedArray } from 'mol-data/int';
import { MoleculeType, isNucleic, isPurinBase, isPyrimidineBase } from 'mol-model/structure/model/types';
import { getElementIndexForAtomId } from 'mol-model/structure/util';
import { getElementIndexForResidueTypeAtomId } from './util/polymer';

async function createNucleotideBlockMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const builder = MeshBuilder.create(256, 128, mesh)

    const { elements, model } = unit
    const { chemicalComponentMap, modifiedResidues } = unit.model.properties
    const { chainAtomSegments, residueAtomSegments, residues } = unit.model.atomicHierarchy
    const { label_comp_id } = residues
    const pos = unit.conformation.invariantPosition

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)

    const p1 = Vec3.zero()
    const p2 = Vec3.zero()
    const p3 = Vec3.zero()
    const p4 = Vec3.zero()
    const p5 = Vec3.zero()
    const p6 = Vec3.zero()
    const v12 = Vec3.zero()
    const v34 = Vec3.zero()
    const vC = Vec3.zero()
    const center = Vec3.zero()
    const t = Mat4.identity()

    let i = 0
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();
            const cc = chemicalComponentMap.get(label_comp_id.value(residueIndex))
            const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown

            if (isNucleic(moleculeType)) {
                let compId = label_comp_id.value(residueIndex)
                const parentId = modifiedResidues.parentId.get(compId)
                if (parentId !== undefined) compId = parentId
                let idx1 = -1, idx2 = -1, idx3 = -1, idx4 = -1, idx5 = -1, idx6 = -1
                let width = 4.5, height = 4.5, depth = 0.5

                if (isPurinBase(compId)) {
                    height = 4.5
                    idx1 = getElementIndexForAtomId(model, residueIndex, 'N1')
                    idx2 = getElementIndexForAtomId(model, residueIndex, 'C4')
                    idx3 = getElementIndexForAtomId(model, residueIndex, 'C6')
                    idx4 = getElementIndexForAtomId(model, residueIndex, 'C2')
                    idx5 = getElementIndexForAtomId(model, residueIndex, 'N9')
                    idx6 = getElementIndexForResidueTypeAtomId(model, residueIndex, 'trace')
                } else if (isPyrimidineBase(compId)) {
                    height = 3.0
                    idx1 = getElementIndexForAtomId(model, residueIndex, 'N3')
                    idx2 = getElementIndexForAtomId(model, residueIndex, 'C6')
                    idx3 = getElementIndexForAtomId(model, residueIndex, 'C4')
                    idx4 = getElementIndexForAtomId(model, residueIndex, 'C2')
                    idx5 = getElementIndexForAtomId(model, residueIndex, 'N1')
                    idx6 = getElementIndexForResidueTypeAtomId(model, residueIndex, 'trace')
                }

                if (idx1 !== -1 && idx2 !== -1 && idx3 !== -1 && idx4 !== -1 && idx5 !== -1 && idx6 !== -1) {
                    pos(idx1, p1); pos(idx2, p2); pos(idx3, p3); pos(idx4, p4); pos(idx5, p5); pos(idx6, p6)
                    Vec3.normalize(v12, Vec3.sub(v12, p2, p1))
                    Vec3.normalize(v34, Vec3.sub(v34, p4, p3))
                    Vec3.normalize(vC, Vec3.cross(vC, v12, v34))
                    Mat4.targetTo(t, p1, p2, vC)
                    Vec3.scaleAndAdd(center, p1, v12, height / 2)
                    Mat4.setTranslation(t, center)
                    builder.setId(SortedArray.findPredecessorIndex(elements, idx6))
                    builder.addBox(t, { width: width, height: depth, depth: height })
                    builder.addCylinder(p5, p6, 1, { radiusTop: 0.2, radiusBottom: 0.2 })
                }
            }

            if (i % 10000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Gap mesh', current: i });
            }
            ++i
        }
    }

    return builder.getMesh()
}

export const DefaultNucleotideBlockProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type NucleotideBlockProps = Partial<typeof DefaultNucleotideBlockProps>

export function NucleotideBlockVisual(): UnitsVisual<NucleotideBlockProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultNucleotideBlockProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: NucleotideBlockProps = {}) {
            currentProps = Object.assign({}, DefaultNucleotideBlockProps, props)
            currentGroup = group

            const { colorTheme, unitKinds } = { ...DefaultNucleotideBlockProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createNucleotideBlockMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)
            // console.log(mesh)

            const transforms = createTransforms(group)
            const color = createColors(group, elementCount, colorTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts),
                aColor: ValueCell.create(new Float32Array(mesh.vertexCount * 3))
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: NucleotideBlockProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createNucleotideBlockMesh(ctx, unit, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                const elementCount = currentGroup.elements.length
                if (ctx.shouldUpdate) await ctx.update('Computing nucleotide block colors');
                createColors(currentGroup, elementCount, newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            return getElementLoci(renderObject.id, currentGroup, pickingId)
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(renderObject.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
