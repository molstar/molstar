/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../index';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { Segmentation } from 'mol-data/int';
import { MoleculeType, isNucleic, isPurinBase, isPyrimidineBase } from 'mol-model/structure/model/types';
import { getElementIndexForAtomRole } from 'mol-model/structure/util';
import { UnitsMeshVisual, UnitsMeshParams } from '../units-visual';
import { NucleotideLocationIterator, markNucleotideElement, getNucleotideElementLoci } from './util/nucleotide';
import { paramDefaultValues } from 'mol-util/parameter';
import { Box } from 'mol-geo/primitive/box';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { addCylinder } from 'mol-geo/geometry/mesh/builder/cylinder';
import { VisualContext } from 'mol-repr';
import { Theme } from 'mol-geo/geometry/geometry';

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
async function createNucleotideBlockMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: {}, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimate
    const builder = MeshBuilder.create(256, 128, mesh)

    const { elements, model } = unit
    const { chemicalComponentMap, modifiedResidues } = model.properties
    const { chainAtomSegments, residueAtomSegments, residues, index: atomicIndex } = model.atomicHierarchy
    const { label_comp_id } = residues
    const pos = unit.conformation.invariantPosition

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements)
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements)

    let i = 0
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();
            let compId = label_comp_id.value(residueIndex)
            const cc = chemicalComponentMap.get(compId)
            const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown

            if (isNucleic(moleculeType)) {
                const parentId = modifiedResidues.parentId.get(compId)
                if (parentId !== undefined) compId = parentId
                let idx1 = -1, idx2 = -1, idx3 = -1, idx4 = -1, idx5 = -1, idx6 = -1
                let width = 4.5, height = 4.5, depth = 0.5

                if (isPurinBase(compId)) {
                    height = 4.5
                    idx1 = atomicIndex.findAtomOnResidue(residueIndex, 'N1')
                    idx2 = atomicIndex.findAtomOnResidue(residueIndex, 'C4')
                    idx3 = atomicIndex.findAtomOnResidue(residueIndex, 'C6')
                    idx4 = atomicIndex.findAtomOnResidue(residueIndex, 'C2')
                    idx5 = atomicIndex.findAtomOnResidue(residueIndex, 'N9')
                    idx6 = getElementIndexForAtomRole(model, residueIndex, 'trace')
                } else if (isPyrimidineBase(compId)) {
                    height = 3.0
                    idx1 = atomicIndex.findAtomOnResidue(residueIndex, 'N3')
                    idx2 = atomicIndex.findAtomOnResidue(residueIndex, 'C6')
                    idx3 = atomicIndex.findAtomOnResidue(residueIndex, 'C4')
                    idx4 = atomicIndex.findAtomOnResidue(residueIndex, 'C2')
                    idx5 = atomicIndex.findAtomOnResidue(residueIndex, 'N1')
                    idx6 = getElementIndexForAtomRole(model, residueIndex, 'trace')
                }

                if (idx5 !== -1 && idx6 !== -1) {
                    pos(idx5, p5); pos(idx6, p6)
                    builder.setGroup(i)
                    addCylinder(builder, p5, p6, 1, { radiusTop: 0.2, radiusBottom: 0.2 })
                    if (idx1 !== -1 && idx2 !== -1 && idx3 !== -1 && idx4 !== -1) {
                        pos(idx1, p1); pos(idx2, p2); pos(idx3, p3); pos(idx4, p4);
                        Vec3.normalize(v12, Vec3.sub(v12, p2, p1))
                        Vec3.normalize(v34, Vec3.sub(v34, p4, p3))
                        Vec3.normalize(vC, Vec3.cross(vC, v12, v34))
                        Mat4.targetTo(t, p1, p2, vC)
                        Vec3.scaleAndAdd(center, p1, v12, height / 2 - 0.2)
                        Mat4.scale(t, t, Vec3.set(sVec, width, depth, height))
                        Mat4.setTranslation(t, center)
                        builder.add(t, box)
                    }
                }

                if (i % 10000 === 0 && ctx.runtime.shouldUpdate) {
                    await ctx.runtime.update({ message: 'Nucleotide block mesh', current: i });
                }
                ++i
            }
        }
    }

    return builder.getMesh()
}

export const NucleotideBlockParams = {
    ...UnitsMeshParams
}
export const DefaultNucleotideBlockProps = paramDefaultValues(NucleotideBlockParams)
export type NucleotideBlockProps = typeof DefaultNucleotideBlockProps

export function NucleotideBlockVisual(): UnitsVisual<NucleotideBlockProps> {
    return UnitsMeshVisual<NucleotideBlockProps>({
        defaultProps: DefaultNucleotideBlockProps,
        createGeometry: createNucleotideBlockMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        mark: markNucleotideElement,
        setUpdateState: () => {}
    })
}