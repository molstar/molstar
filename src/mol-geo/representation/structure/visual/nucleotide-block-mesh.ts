/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { Mesh } from '../../../mesh/mesh';
import { MeshBuilder } from '../../../mesh/mesh-builder';
import { getElementLoci, markElement, StructureElementIterator } from './util/element';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Segmentation, SortedArray } from 'mol-data/int';
import { MoleculeType, isNucleic, isPurinBase, isPyrimidineBase } from 'mol-model/structure/model/types';
import { getElementIndexForAtomId, getElementIndexForAtomRole } from 'mol-model/structure/util';
import { DefaultUnitsMeshProps, UnitsMeshVisual } from '../units-visual';
import { addCylinder } from '../../../mesh/builder/cylinder';
import { Box } from '../../../primitive/box';

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
const sVec = Vec3.zero()
const box = Box()

// TODO define props, should be scalable
async function createNucleotideBlockMesh(ctx: RuntimeContext, unit: Unit, props: {}, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimate
    const builder = MeshBuilder.create(256, 128, mesh)

    const { elements, model } = unit
    const { chemicalComponentMap, modifiedResidues } = model.properties
    const { chainAtomSegments, residueAtomSegments, residues } = model.atomicHierarchy
    const { label_comp_id } = residues
    const pos = unit.conformation.invariantPosition

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)

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
                    idx6 = getElementIndexForAtomRole(model, residueIndex, 'trace')
                } else if (isPyrimidineBase(compId)) {
                    height = 3.0
                    idx1 = getElementIndexForAtomId(model, residueIndex, 'N3')
                    idx2 = getElementIndexForAtomId(model, residueIndex, 'C6')
                    idx3 = getElementIndexForAtomId(model, residueIndex, 'C4')
                    idx4 = getElementIndexForAtomId(model, residueIndex, 'C2')
                    idx5 = getElementIndexForAtomId(model, residueIndex, 'N1')
                    idx6 = getElementIndexForAtomRole(model, residueIndex, 'trace')
                }

                if (idx1 !== -1 && idx2 !== -1 && idx3 !== -1 && idx4 !== -1 && idx5 !== -1 && idx6 !== -1) {
                    pos(idx1, p1); pos(idx2, p2); pos(idx3, p3); pos(idx4, p4); pos(idx5, p5); pos(idx6, p6)
                    Vec3.normalize(v12, Vec3.sub(v12, p2, p1))
                    Vec3.normalize(v34, Vec3.sub(v34, p4, p3))
                    Vec3.normalize(vC, Vec3.cross(vC, v12, v34))
                    Mat4.targetTo(t, p1, p2, vC)
                    Vec3.scaleAndAdd(center, p1, v12, height / 2 - 0.2)
                    Mat4.scale(t, t, Vec3.set(sVec, width, depth, height))
                    Mat4.setTranslation(t, center)
                    builder.setGroup(SortedArray.findPredecessorIndex(elements, idx6))
                    builder.add(t, box)
                    addCylinder(builder, p5, p6, 1, { radiusTop: 0.2, radiusBottom: 0.2 })
                }
            }

            if (i % 10000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Nucleotide block mesh', current: i });
            }
            ++i
        }
    }

    return builder.getMesh()
}

export const DefaultNucleotideBlockProps = {
    ...DefaultUnitsMeshProps
}
export type NucleotideBlockProps = typeof DefaultNucleotideBlockProps

export function NucleotideBlockVisual(): UnitsVisual<NucleotideBlockProps> {
    return UnitsMeshVisual<NucleotideBlockProps>({
        defaultProps: DefaultNucleotideBlockProps,
        createMesh: createNucleotideBlockMesh,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: () => {}
    })
}