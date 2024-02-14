/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { NumberArray } from '../../../mol-util/type-helpers';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../../mol-data/int';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement, getNucleotideBaseType, createNucleicIndices, setSugarIndices, hasSugarIndices, setPurinIndices, hasPyrimidineIndices, setPyrimidineIndices, hasPurinIndices } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { Sphere3D } from '../../../mol-math/geometry';

// TODO support ring-fills for multiple locations (including from microheterogeneity)

const pN1 = Vec3();
const pC2 = Vec3();
const pN3 = Vec3();
const pC4 = Vec3();
const pC5 = Vec3();
const pC6 = Vec3();
const pN7 = Vec3();
const pC8 = Vec3();
const pN9 = Vec3();

const pC1_1 = Vec3();
const pC2_1 = Vec3();
const pC3_1 = Vec3();
const pC4_1 = Vec3();
const pO4_1 = Vec3();

const mid = Vec3();
const normal = Vec3();
const shift = Vec3();

export const NucleotideAtomicRingFillMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    thicknessFactor: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
};
export const DefaultNucleotideAtomicRingFillMeshProps = PD.getDefaultValues(NucleotideAtomicRingFillMeshParams);
export type NucleotideAtomicRingFillProps = typeof DefaultNucleotideAtomicRingFillMeshProps

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

const tmpShiftV = Vec3();
function shiftPositions(out: NumberArray, dir: Vec3, ...positions: Vec3[]) {
    for (let i = 0, il = positions.length; i < il; ++i) {
        const v = positions[i];
        Vec3.toArray(Vec3.add(tmpShiftV, v, dir), out, (i * 2) * 3);
        Vec3.toArray(Vec3.sub(tmpShiftV, v, dir), out, (i * 2 + 1) * 3);
    }
}

function createNucleotideAtomicRingFillMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideAtomicRingFillProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, thicknessFactor } = props;

    const vertexCount = nucleotideElementCount * 25;
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const thickness = sizeFactor * thicknessFactor;

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const idx = createNucleicIndices();

                builderState.currentGroup = i;

                setSugarIndices(idx, unit, residueIndex);
                if (hasSugarIndices(idx)) {
                    c.invariantPosition(idx.C1_1, pC1_1); c.invariantPosition(idx.C2_1, pC2_1); c.invariantPosition(idx.C3_1, pC3_1); c.invariantPosition(idx.C4_1, pC4_1); c.invariantPosition(idx.O4_1, pO4_1);

                    // sugar ring
                    Vec3.triangleNormal(normal, pC3_1, pC4_1, pC1_1);
                    Vec3.scale(mid, Vec3.add(mid, pO4_1, Vec3.add(mid, pC4_1, Vec3.add(mid, pC3_1, Vec3.add(mid, pC1_1, pC2_1)))), 0.2 /* 1 / 5 */);

                    Vec3.scale(shift, normal, thickness);
                    shiftPositions(positionsRing5, shift, mid, pC3_1, pC4_1, pO4_1, pC1_1, pC2_1);

                    MeshBuilder.addTriangleStrip(builderState, positionsRing5, stripIndicesRing5);
                    MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5, fanIndicesTopRing5, normal);
                    Vec3.negate(normal, normal);
                    MeshBuilder.addTriangleFanWithNormal(builderState, positionsRing5, fanIndicesBottomRing5, normal);
                }

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8), c.invariantPosition(idx.N9, pN9);

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
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

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

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, thickness);
    m.setBoundingSphere(sphere);

    return m;
}

export const NucleotideAtomicRingFillParams = {
    ...UnitsMeshParams,
    ...NucleotideAtomicRingFillMeshParams
};
export type NucleotideAtomicRingFillParams = typeof NucleotideAtomicRingFillParams

export function NucleotideAtomicRingFillVisual(materialId: number): UnitsVisual<NucleotideAtomicRingFillParams> {
    return UnitsMeshVisual<NucleotideAtomicRingFillParams>({
        defaultProps: PD.getDefaultValues(NucleotideAtomicRingFillParams),
        createGeometry: createNucleotideAtomicRingFillMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideAtomicRingFillParams>, currentProps: PD.Values<NucleotideAtomicRingFillParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.thicknessFactor !== currentProps.thicknessFactor
            );
        }
    }, materialId);
}
