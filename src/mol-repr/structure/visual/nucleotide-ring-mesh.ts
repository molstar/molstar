/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
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
import { CylinderProps } from '../../../mol-geo/primitive/cylinder';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement, getNucleotideBaseType, createNucleicIndices, setPurinIndices, setPyrimidineIndices, hasPyrimidineIndices, hasPurinIndices } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';

// TODO support rings for multiple locations (including from microheterogeneity)

const pTrace = Vec3();
const pN1 = Vec3();
const pC2 = Vec3();
const pN3 = Vec3();
const pC4 = Vec3();
const pC5 = Vec3();
const pC6 = Vec3();
const pN7 = Vec3();
const pC8 = Vec3();
const pN9 = Vec3();
const normal = Vec3();

export const NucleotideRingMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    thicknessFactor: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
};
export const DefaultNucleotideRingMeshProps = PD.getDefaultValues(NucleotideRingMeshParams);
export type NucleotideRingProps = typeof DefaultNucleotideRingMeshProps

const positionsRing5_6 = new Float32Array(2 * 9 * 3);
const stripIndicesRing5_6 = new Uint32Array([0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 14, 15, 12, 13, 8, 9, 10, 11, 0, 1]);
const fanIndicesTopRing5_6 = new Uint32Array([8, 12, 14, 16, 6, 4, 2, 0, 10]);
const fanIndicesBottomRing5_6 = new Uint32Array([9, 11, 1, 3, 5, 7, 17, 15, 13]);

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

function createNucleotideRingMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, thicknessFactor, radialSegments, detail } = props;

    const vertexCount = nucleotideElementCount * (26 + radialSegments * 2);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const radius = 1 * sizeFactor;
    const thickness = thicknessFactor * sizeFactor;
    const cylinderProps: CylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments };

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const idx = createNucleicIndices();

                builderState.currentGroup = i;

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (idx.N9 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N9, pN9); c.invariantPosition(idx.trace, pTrace);
                        builderState.currentGroup = i;
                        addCylinder(builderState, pN9, pTrace, 1, cylinderProps);
                        addSphere(builderState, pN9, radius, detail);
                    }

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8);

                        Vec3.triangleNormal(normal, pN1, pC4, pC5);
                        Vec3.scale(normal, normal, thickness);
                        shiftPositions(positionsRing5_6, normal, pN1, pC2, pN3, pC4, pC5, pC6, pN7, pC8, pN9);

                        MeshBuilder.addTriangleStrip(builderState, positionsRing5_6, stripIndicesRing5_6);
                        MeshBuilder.addTriangleFan(builderState, positionsRing5_6, fanIndicesTopRing5_6);
                        MeshBuilder.addTriangleFan(builderState, positionsRing5_6, fanIndicesBottomRing5_6);
                    }
                } else if (isPyrimidine) {
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (idx.N1 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.trace, pTrace);
                        builderState.currentGroup = i;
                        addCylinder(builderState, pN1, pTrace, 1, cylinderProps);
                        addSphere(builderState, pN1, radius, detail);
                    }

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

                        Vec3.triangleNormal(normal, pN1, pC4, pC5);
                        Vec3.scale(normal, normal, thickness);
                        shiftPositions(positionsRing6, normal, pN1, pC2, pN3, pC4, pC5, pC6);

                        MeshBuilder.addTriangleStrip(builderState, positionsRing6, stripIndicesRing6);
                        MeshBuilder.addTriangleFan(builderState, positionsRing6, fanIndicesTopRing6);
                        MeshBuilder.addTriangleFan(builderState, positionsRing6, fanIndicesBottomRing6);
                    }
                }

                ++i;
            }
        }
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, radius);
    m.setBoundingSphere(sphere);

    return m;
}

export const NucleotideRingParams = {
    ...UnitsMeshParams,
    ...NucleotideRingMeshParams
};
export type NucleotideRingParams = typeof NucleotideRingParams

export function NucleotideRingVisual(materialId: number): UnitsVisual<NucleotideRingParams> {
    return UnitsMeshVisual<NucleotideRingParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingParams),
        createGeometry: createNucleotideRingMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingParams>, currentProps: PD.Values<NucleotideRingParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.thicknessFactor !== currentProps.thicknessFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        }
    }, materialId);
}