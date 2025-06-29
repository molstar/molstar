/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { VisualContext } from '../../visual';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../../mol-data/int';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsSpheresParams, UnitsSpheresVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement, getNucleotideBaseType, createNucleicIndices, setSugarIndices, hasSugarIndices, setPurinIndices, hasPurinIndices, setPyrimidineIndices, hasPyrimidineIndices } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Spheres } from '../../../mol-geo/geometry/spheres/spheres';
import { sphereVertexCount } from '../../../mol-geo/primitive/sphere';
import { SpheresBuilder } from '../../../mol-geo/geometry/spheres/spheres-builder';
import { checkSphereImpostorSupport, StructureGroup } from './util/common';

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

const pC1_1 = Vec3();
const pC2_1 = Vec3();
const pC3_1 = Vec3();
const pC4_1 = Vec3();
const pO4_1 = Vec3();

export const NucleotideAtomicElementParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true)
};
export type NucleotideAtomicElementParams = typeof NucleotideAtomicElementParams
interface NucleotideAtomicElementImpostorProps {
    sizeFactor: number,
}

export function NucleotideAtomicElementVisual(materialId: number, structure: Structure, props: PD.Values<NucleotideAtomicElementParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && checkSphereImpostorSupport(webgl)
        ? NucleotideAtomicElementImpostorVisual(materialId)
        : NucleotideAtomicElementMeshVisual(materialId);
}

function createNucleotideAtomicElementImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideAtomicElementImpostorProps, spheres?: Spheres) {
    if (!Unit.isAtomic(unit)) return Spheres.createEmpty(spheres);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Spheres.createEmpty(spheres);

    const spheresCountEstimate = nucleotideElementCount * 15; // 15 is the average purine (17) & pirimidine (13) bonds
    const builder = SpheresBuilder.create(spheresCountEstimate, spheresCountEstimate / 4, spheres);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const idx = createNucleicIndices();

                setSugarIndices(idx, unit, residueIndex);

                if (hasSugarIndices(idx)) {
                    c.invariantPosition(idx.C1_1, pC1_1); c.invariantPosition(idx.C2_1, pC2_1); c.invariantPosition(idx.C3_1, pC3_1); c.invariantPosition(idx.C4_1, pC4_1); c.invariantPosition(idx.O4_1, pO4_1);

                    // trace cylinder
                    c.invariantPosition(idx.trace, pTrace);
                    builder.add(pTrace[0], pTrace[1], pTrace[2], i);

                    // sugar ring
                    builder.add(pC3_1[0], pC3_1[1], pC3_1[2], i);
                    builder.add(pC4_1[0], pC4_1[1], pC4_1[2], i);
                    builder.add(pO4_1[0], pO4_1[1], pO4_1[2], i);
                    builder.add(pC1_1[0], pC1_1[1], pC1_1[2], i);
                    builder.add(pC2_1[0], pC2_1[1], pC2_1[2], i);
                }

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8); c.invariantPosition(idx.N9, pN9);

                        // base ring
                        builder.add(pN9[0], pN9[1], pN9[2], i);
                        builder.add(pC8[0], pC8[1], pC8[2], i);
                        builder.add(pN7[0], pN7[1], pN7[2], i);
                        builder.add(pC5[0], pC5[1], pC5[2], i);
                        builder.add(pC6[0], pC6[1], pC6[2], i);
                        builder.add(pN1[0], pN1[1], pN1[2], i);
                        builder.add(pC2[0], pC2[1], pC2[2], i);
                        builder.add(pN3[0], pN3[1], pN3[2], i);
                        builder.add(pC4[0], pC4[1], pC4[2], i);
                    }
                } else if (isPyrimidine) {
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

                        // base ring
                        builder.add(pN1[0], pN1[1], pN1[2], i);
                        builder.add(pC6[0], pC6[1], pC6[2], i);
                        builder.add(pC5[0], pC5[1], pC5[2], i);
                        builder.add(pC4[0], pC4[1], pC4[2], i);
                        builder.add(pN3[0], pN3[1], pN3[2], i);
                        builder.add(pC2[0], pC2[1], pC2[2], i);
                    }
                }

                ++i;
            }
        }
    }
    const s = builder.getSpheres();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    s.setBoundingSphere(sphere);

    return s;
}

export function NucleotideAtomicElementImpostorVisual(materialId: number): UnitsVisual<NucleotideAtomicElementParams> {
    return UnitsSpheresVisual<NucleotideAtomicElementParams>({
        defaultProps: PD.getDefaultValues(NucleotideAtomicElementParams),
        createGeometry: createNucleotideAtomicElementImpostor,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideAtomicElementParams>, currentProps: PD.Values<NucleotideAtomicElementParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideAtomicElementParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

interface NucleotideAtomicElementMeshProps {
    detail: number,
    sizeFactor: number,
}

function createNucleotideAtomicElementMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideAtomicElementMeshProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, detail } = props;

    const vertexCount = nucleotideElementCount * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const radius = 1 * sizeFactor;

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

                    // trace cylinder
                    c.invariantPosition(idx.trace, pTrace);
                    addSphere(builderState, pTrace, radius, detail);

                    // sugar ring
                    addSphere(builderState, pC4_1, radius, detail);
                    addSphere(builderState, pO4_1, radius, detail);
                    addSphere(builderState, pC1_1, radius, detail);
                    addSphere(builderState, pC2_1, radius, detail);
                    addSphere(builderState, pC3_1, radius, detail);
                }

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8); c.invariantPosition(idx.N9, pN9);

                        // base ring
                        addSphere(builderState, pC8, radius, detail);
                        addSphere(builderState, pN7, radius, detail);
                        addSphere(builderState, pC5, radius, detail);
                        addSphere(builderState, pC6, radius, detail);
                        addSphere(builderState, pN1, radius, detail);
                        addSphere(builderState, pC2, radius, detail);
                        addSphere(builderState, pN3, radius, detail);
                        addSphere(builderState, pC4, radius, detail);
                        addSphere(builderState, pC5, radius, detail);
                        addSphere(builderState, pN9, radius, detail);
                    }
                } else if (isPyrimidine) {
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

                        // base ring
                        addSphere(builderState, pC6, radius, detail);
                        addSphere(builderState, pC5, radius, detail);
                        addSphere(builderState, pC4, radius, detail);
                        addSphere(builderState, pN3, radius, detail);
                        addSphere(builderState, pC2, radius, detail);
                        addSphere(builderState, pN1, radius, detail);
                    }
                }

                ++i;
            }
        }
    }

    const m = MeshBuilder.getMesh(builderState);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}


export function NucleotideAtomicElementMeshVisual(materialId: number): UnitsVisual<NucleotideAtomicElementParams> {
    return UnitsMeshVisual<NucleotideAtomicElementParams>({
        defaultProps: PD.getDefaultValues(NucleotideAtomicElementParams),
        createGeometry: createNucleotideAtomicElementMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideAtomicElementParams>, currentProps: PD.Values<NucleotideAtomicElementParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideAtomicElementParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
