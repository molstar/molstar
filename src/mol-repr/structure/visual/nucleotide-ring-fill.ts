/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { NumberArray } from '../../../mol-util/type-helpers';
import { VisualContext } from '../../visual';
import { Unit, Structure, ElementIndex } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../../mol-data/int';
import { isNucleic, isPurineBase, isPyrimidineBase } from '../../../mol-model/structure/model/types';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

// TODO support rings for multiple locations (including from microheterogeneity)

const pN1 = Vec3.zero();
const pC2 = Vec3.zero();
const pN3 = Vec3.zero();
const pC4 = Vec3.zero();
const pC5 = Vec3.zero();
const pC6 = Vec3.zero();
const pN7 = Vec3.zero();
const pC8 = Vec3.zero();
const pN9 = Vec3.zero();

const pC1_1 = Vec3.zero();
const pC2_1 = Vec3.zero();
const pC3_1 = Vec3.zero();
const pC4_1 = Vec3.zero();
const pO4_1 = Vec3.zero();

const mid = Vec3.zero();
const normal = Vec3.zero();
const shift = Vec3.zero();

export const NucleotideRingFillMeshParams = {
    nucleicRingThickness: PD.Numeric(0.5, { min: 0, max: 2, step: 0.01 }),
};
export const DefaultNucleotideRingFillMeshProps = PD.getDefaultValues(NucleotideRingFillMeshParams);
export type NucleotideRingFillProps = typeof DefaultNucleotideRingFillMeshProps

const positionsRing5_6 = new Float32Array(2 * 9 * 3);
const stripIndicesRing5_6 = new Uint32Array([0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 14, 15, 12, 13, 8, 9, 10, 11, 0, 1]);
const fanIndicesTopRing5_6 = new Uint32Array([8, 12, 14, 16, 6, 4, 2, 0, 10]);
const fanIndicesBottomRing5_6 = new Uint32Array([9, 11, 1, 3, 5, 7, 17, 15, 13]);

const positionsRing5 = new Float32Array(2 * 6 * 3);
const stripIndicesRing5 = new Uint32Array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3]);
const fanIndicesTopRing5 = new Uint32Array([0, 10, 8, 6, 4, 2, 10]);
const fanIndicesBottomRing5 = new Uint32Array([1, 3, 5, 7, 9, 11, 3]);

const positionsRing6 = new Float32Array(2 * 6 * 3);
const stripIndicesRing6 = new Uint32Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1]);
const fanIndicesTopRing6 = new Uint32Array([0, 10, 8, 6, 4, 2]);
const fanIndicesBottomRing6 = new Uint32Array([1, 3, 5, 7, 9, 11]);

const tmpShiftV = Vec3.zero();
function shiftPositions(out: NumberArray, dir: Vec3, ...positions: Vec3[]) {
    for (let i = 0, il = positions.length; i < il; ++i) {
        const v = positions[i];
        Vec3.toArray(Vec3.add(tmpShiftV, v, dir), out, (i * 2) * 3);
        Vec3.toArray(Vec3.sub(tmpShiftV, v, dir), out, (i * 2 + 1) * 3);
    }
}

function createNucleotideRingFillMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingFillProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { nucleicRingThickness } = props;

    const vertexCount = nucleotideElementCount * 25;
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments, atoms, index: atomicIndex } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;
    const { label_comp_id } = atoms;
    const pos = unit.conformation.invariantPosition;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);
    const thickness = nucleicRingThickness;

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const compId = label_comp_id.value(residueAtomSegments.offsets[residueIndex]);

                let idxN1: ElementIndex | -1 = -1, idxC2: ElementIndex | -1 = -1, idxN3: ElementIndex | -1 = -1, idxC4: ElementIndex | -1 = -1, idxC5: ElementIndex | -1 = -1, idxC6: ElementIndex | -1 = -1, idxN7: ElementIndex | -1 = -1, idxC8: ElementIndex | -1 = -1, idxN9: ElementIndex | -1 = -1,
                    idxC1_1: ElementIndex | -1 = -1, idxC2_1: ElementIndex | -1 = -1, idxC3_1: ElementIndex | -1 = -1, idxC4_1: ElementIndex | -1 = -1, idxO4_1: ElementIndex | -1 = -1;

                builderState.currentGroup = i;

                // sugar base
                idxC1_1 = atomicIndex.findAtomOnResidue(residueIndex, "C1'");
                idxC2_1 = atomicIndex.findAtomOnResidue(residueIndex, "C2'");
                idxC3_1 = atomicIndex.findAtomOnResidue(residueIndex, "C3'");
                idxC4_1 = atomicIndex.findAtomOnResidue(residueIndex, "C4'");
                idxO4_1 = atomicIndex.findAtomOnResidue(residueIndex, "O4'");
                if (idxC1_1 !== -1 && idxC2_1 !== -1 && idxC3_1 !== -1 && idxC4_1 !== -1 && idxO4_1 !== -1) {
                    pos(idxC1_1, pC1_1); pos(idxC2_1, pC2_1); pos(idxC3_1, pC3_1); pos(idxC4_1, pC4_1); pos(idxO4_1, pO4_1);

                    // sugar ring
                    Vec3.triangleNormal(normal, pC3_1, pC4_1, pC1_1);
                    Vec3.scale(mid, Vec3.add(mid, pO4_1, Vec3.add(mid, pC4_1, Vec3.add(mid, pC3_1, Vec3.add(mid, pC1_1, pC2_1)))), 0.2/* 1 / 5 */);

                    Vec3.scale(shift, normal, thickness);
                    shiftPositions(positionsRing5, shift, mid, pC3_1, pC4_1, pO4_1, pC1_1, pC2_1);

                    MeshBuilder.addTriangleStrip(builderState, positionsRing5, stripIndicesRing5);
                    MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5, fanIndicesTopRing5, normal);
                    Vec3.negate(normal, normal);
                    MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5, fanIndicesBottomRing5, normal);
                }

                let isPurine = isPurineBase(compId);
                let isPyrimidine = isPyrimidineBase(compId);

                if (!isPurine && !isPyrimidine) {
                    // detect Purine or Pyrimidin based on geometry
                    const idxC4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
                    const idxN9 = atomicIndex.findAtomOnResidue(residueIndex, 'N9');
                    if (idxC4 !== -1 && idxN9 !== -1 && Vec3.distance(pos(idxC4, pC4), pos(idxN9, pN9)) < 1.6) {
                        isPurine = true;
                    } else {
                        isPyrimidine = true;
                    }
                }

                if (isPurine) {
                    idxN1 = atomicIndex.findAtomOnResidue(residueIndex, 'N1');
                    idxC2 = atomicIndex.findAtomOnResidue(residueIndex, 'C2');
                    idxN3 = atomicIndex.findAtomOnResidue(residueIndex, 'N3');
                    idxC4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
                    idxC5 = atomicIndex.findAtomOnResidue(residueIndex, 'C5');
                    if (idxC5 === -1) {
                        // modified ring, e.g. DP
                        idxC5 = atomicIndex.findAtomOnResidue(residueIndex, 'N5');
                    }
                    idxC6 = atomicIndex.findAtomOnResidue(residueIndex, 'C6');
                    idxN7 = atomicIndex.findAtomOnResidue(residueIndex, 'N7');
                    if (idxN7 === -1) {
                        // modified ring, e.g. DP
                        idxN7 = atomicIndex.findAtomOnResidue(residueIndex, 'C7');
                    }
                    idxC8 = atomicIndex.findAtomOnResidue(residueIndex, 'C8');
                    idxN9 = atomicIndex.findAtomOnResidue(residueIndex, 'N9');

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1 && idxN7 !== -1 && idxC8 !== -1 && idxN9 !== -1) {
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6); pos(idxN7, pN7); pos(idxC8, pC8), pos(idxN9, pN9);

                        // base ring
                        Vec3.triangleNormal(normal, pN1, pC4, pC5);
                        Vec3.scale(shift, normal, thickness);
                        shiftPositions(positionsRing5_6, shift, pN1, pC2, pN3, pC4, pC5, pC6, pN7, pC8, pN9);

                        MeshBuilder.addTriangleStrip(builderState, positionsRing5_6, stripIndicesRing5_6);
                        MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5_6, fanIndicesTopRing5_6, normal);
                        Vec3.negate(normal, normal);
                        MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5_6, fanIndicesBottomRing5_6, normal);
                    }
                } else if (isPyrimidine) {
                    idxN1 = atomicIndex.findAtomOnResidue(residueIndex, 'N1');
                    if (idxN1 === -1) {
                        // modified ring, e.g. DZ
                        idxN1 = atomicIndex.findAtomOnResidue(residueIndex, 'C1');
                    }
                    idxC2 = atomicIndex.findAtomOnResidue(residueIndex, 'C2');
                    idxN3 = atomicIndex.findAtomOnResidue(residueIndex, 'N3');
                    idxC4 = atomicIndex.findAtomOnResidue(residueIndex, 'C4');
                    idxC5 = atomicIndex.findAtomOnResidue(residueIndex, 'C5');
                    idxC6 = atomicIndex.findAtomOnResidue(residueIndex, 'C6');

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1) {
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6);

                        // base ring
                        Vec3.triangleNormal(normal, pN1, pC4, pC5);
                        Vec3.scale(shift, normal, thickness);
                        shiftPositions(positionsRing6, shift, pN1, pC2, pN3, pC4, pC5, pC6);

                        MeshBuilder.addTriangleStrip(builderState, positionsRing6, stripIndicesRing6);
                        MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing6, fanIndicesTopRing6, normal);
                        Vec3.negate(normal, normal);
                        MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing6, fanIndicesBottomRing6, normal);
                    }
                }

                ++i;
            }
        }
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.nucleicRingThickness);
    m.setBoundingSphere(sphere);

    return m;
}

export const NucleotideRingFillParams = {
    ...UnitsMeshParams,
    ...NucleotideRingFillMeshParams
};
export type NucleotideRingFillParams = typeof NucleotideRingFillParams

export function NucleotideRingFillVisual(materialId: number): UnitsVisual<NucleotideRingFillParams> {
    return UnitsMeshVisual<NucleotideRingFillParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingFillParams),
        createGeometry: createNucleotideRingFillMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingFillParams>, currentProps: PD.Values<NucleotideRingFillParams>) => {
            state.createGeometry = (
                newProps.nucleicRingThickness !== currentProps.nucleicRingThickness
            );
        }
    }, materialId);
}
