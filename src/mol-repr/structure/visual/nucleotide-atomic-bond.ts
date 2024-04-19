/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { CylinderProps } from '../../../mol-geo/primitive/cylinder';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsCylindersParams, UnitsCylindersVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement, getNucleotideBaseType, createNucleicIndices, setSugarIndices, hasSugarIndices, setPurinIndices, hasPurinIndices, setPyrimidineIndices, hasPyrimidineIndices } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';

import { WebGLContext } from '../../../mol-gl/webgl/context';

import { Cylinders } from '../../../mol-geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '../../../mol-geo/geometry/cylinders/cylinders-builder';
import { StructureGroup } from './util/common';

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

export const NucleotideAtomicBondParams = {
    ...UnitsMeshParams,
    ...UnitsCylindersParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true)
};
export type NucleotideAtomicBondParams = typeof NucleotideAtomicBondParams
interface NucleotideAtomicBondImpostorProps {
    sizeFactor: number,
}

export function NucleotideAtomicBondVisual(materialId: number, structure: Structure, props: PD.Values<NucleotideAtomicBondParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth
        ? NucleotideAtomicBondImpostorVisual(materialId)
        : NucleotideAtomicBondMeshVisual(materialId);
}

function createNucleotideAtomicBondImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideAtomicBondImpostorProps, cylinders?: Cylinders) {
    if (!Unit.isAtomic(unit)) return Cylinders.createEmpty(cylinders);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Cylinders.createEmpty(cylinders);

    const cylindersCountEstimate = nucleotideElementCount * 15; // 15 is the average purine (17) & pirimidine (13) bonds
    const builder = CylindersBuilder.create(cylindersCountEstimate, cylindersCountEstimate / 4, cylinders);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    let i = 0;
    const colorModeFlag = 2.0;
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
                    builder.add(pC3_1[0], pC3_1[1], pC3_1[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, colorModeFlag, i);

                    // sugar ring
                    builder.add(pC3_1[0], pC3_1[1], pC3_1[2], pC4_1[0], pC4_1[1], pC4_1[2], 1, true, true, colorModeFlag, i);
                    builder.add(pC4_1[0], pC4_1[1], pC4_1[2], pO4_1[0], pO4_1[1], pO4_1[2], 1, true, true, colorModeFlag, i);
                    builder.add(pO4_1[0], pO4_1[1], pO4_1[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, colorModeFlag, i);
                    builder.add(pC1_1[0], pC1_1[1], pC1_1[2], pC2_1[0], pC2_1[1], pC2_1[2], 1, true, true, colorModeFlag, i);
                    builder.add(pC2_1[0], pC2_1[1], pC2_1[2], pC3_1[0], pC3_1[1], pC3_1[2], 1, true, true, colorModeFlag, i);
                }

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (idx.C1_1 !== -1 && idx.N9 !== -1) {
                        c.invariantPosition(idx.C1_1, pC1_1); c.invariantPosition(idx.N9, pN9);
                        builder.add(pN9[0], pN9[1], pN9[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, colorModeFlag, i);
                    } else if (idx.N9 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N9, pN9); c.invariantPosition(idx.trace, pTrace);
                        builder.add(pN9[0], pN9[1], pN9[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, colorModeFlag, i);
                    }

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8); c.invariantPosition(idx.N9, pN9);

                        // base ring
                        builder.add(pN9[0], pN9[1], pN9[2], pC8[0], pC8[1], pC8[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC8[0], pC8[1], pC8[2], pN7[0], pN7[1], pN7[2], 1, true, true, colorModeFlag, i);
                        builder.add(pN7[0], pN7[1], pN7[2], pC5[0], pC5[1], pC5[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC5[0], pC5[1], pC5[2], pC6[0], pC6[1], pC6[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC6[0], pC6[1], pC6[2], pN1[0], pN1[1], pN1[2], 1, true, true, colorModeFlag, i);
                        builder.add(pN1[0], pN1[1], pN1[2], pC2[0], pC2[1], pC2[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC2[0], pC2[1], pC2[2], pN3[0], pN3[1], pN3[2], 1, true, true, colorModeFlag, i);
                        builder.add(pN3[0], pN3[1], pN3[2], pC4[0], pC4[1], pC4[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pC5[0], pC5[1], pC5[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pN9[0], pN9[1], pN9[2], 1, true, true, colorModeFlag, i);

                    }
                } else if (isPyrimidine) {
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (idx.C1_1 !== -1 && idx.N1 !== -1) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C1_1, pC1_1);
                        builder.add(pN1[0], pN1[1], pN1[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, colorModeFlag, i);
                    } else if (idx.N1 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.trace, pTrace);
                        builder.add(pN1[0], pN1[1], pN1[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, colorModeFlag, i);
                    }

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

                        // base ring
                        builder.add(pN1[0], pN1[1], pN1[2], pC6[0], pC6[1], pC6[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC6[0], pC6[1], pC6[2], pC5[0], pC5[1], pC5[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC5[0], pC5[1], pC5[2], pC4[0], pC4[1], pC4[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pN3[0], pN3[1], pN3[2], 1, true, true, colorModeFlag, i);
                        builder.add(pN3[0], pN3[1], pN3[2], pC2[0], pC2[1], pC2[2], 1, true, true, colorModeFlag, i);
                        builder.add(pC2[0], pC2[1], pC2[2], pN1[0], pN1[1], pN1[2], 1, true, true, colorModeFlag, i);
                    }
                }

                ++i;
            }
        }
    }
    const cy = builder.getCylinders();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    cy.setBoundingSphere(sphere);

    return cy;
}

export function NucleotideAtomicBondImpostorVisual(materialId: number): UnitsVisual<NucleotideAtomicBondParams> {
    return UnitsCylindersVisual<NucleotideAtomicBondParams>({
        defaultProps: PD.getDefaultValues(NucleotideAtomicBondParams),
        createGeometry: createNucleotideAtomicBondImpostor,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideAtomicBondParams>, currentProps: PD.Values<NucleotideAtomicBondParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideAtomicBondParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

interface NucleotideAtomicBondMeshProps {
    radialSegments: number,
    sizeFactor: number,
}

function createNucleotideAtomicBondMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideAtomicBondMeshProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, radialSegments } = props;

    const vertexCount = nucleotideElementCount * (radialSegments * 15); // 15 is the average purine (17) & pirimidine (13) bonds
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const cylinderProps: CylinderProps = { radiusTop: 1 * sizeFactor, radiusBottom: 1 * sizeFactor, radialSegments };

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
                    addCylinder(builderState, pC3_1, pTrace, 1, cylinderProps);

                    // sugar ring
                    addCylinder(builderState, pC3_1, pC4_1, 1, cylinderProps);
                    addCylinder(builderState, pC4_1, pO4_1, 1, cylinderProps);
                    addCylinder(builderState, pO4_1, pC1_1, 1, cylinderProps);
                    addCylinder(builderState, pC1_1, pC2_1, 1, cylinderProps);
                    addCylinder(builderState, pC2_1, pC3_1, 1, cylinderProps);
                }

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    setPurinIndices(idx, unit, residueIndex);

                    if (idx.C1_1 !== -1 && idx.N9 !== -1) {
                        c.invariantPosition(idx.C1_1, pC1_1); c.invariantPosition(idx.N9, pN9);
                        addCylinder(builderState, pN9, pC1_1, 1, cylinderProps);
                    } else if (idx.N9 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N9, pN9); c.invariantPosition(idx.trace, pTrace);
                        addCylinder(builderState, pN9, pTrace, 1, cylinderProps);
                    }

                    if (hasPurinIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6); c.invariantPosition(idx.N7, pN7); c.invariantPosition(idx.C8, pC8); c.invariantPosition(idx.N9, pN9);

                        // base ring
                        addCylinder(builderState, pN9, pC8, 1, cylinderProps);
                        addCylinder(builderState, pC8, pN7, 1, cylinderProps);
                        addCylinder(builderState, pN7, pC5, 1, cylinderProps);
                        addCylinder(builderState, pC5, pC6, 1, cylinderProps);
                        addCylinder(builderState, pC6, pN1, 1, cylinderProps);
                        addCylinder(builderState, pN1, pC2, 1, cylinderProps);
                        addCylinder(builderState, pC2, pN3, 1, cylinderProps);
                        addCylinder(builderState, pN3, pC4, 1, cylinderProps);
                        addCylinder(builderState, pC4, pC5, 1, cylinderProps);
                        addCylinder(builderState, pC4, pN9, 1, cylinderProps);
                    }
                } else if (isPyrimidine) {
                    setPyrimidineIndices(idx, unit, residueIndex);

                    if (idx.C1_1 !== -1 && idx.N1 !== -1) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C1_1, pC1_1);
                        addCylinder(builderState, pN1, pC1_1, 1, cylinderProps);
                    } else if (idx.N1 !== -1 && idx.trace !== -1) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.trace, pTrace);
                        addCylinder(builderState, pN1, pTrace, 1, cylinderProps);
                    }

                    if (hasPyrimidineIndices(idx)) {
                        c.invariantPosition(idx.N1, pN1); c.invariantPosition(idx.C2, pC2); c.invariantPosition(idx.N3, pN3); c.invariantPosition(idx.C4, pC4); c.invariantPosition(idx.C5, pC5); c.invariantPosition(idx.C6, pC6);

                        // base ring
                        addCylinder(builderState, pN1, pC6, 1, cylinderProps);
                        addCylinder(builderState, pC6, pC5, 1, cylinderProps);
                        addCylinder(builderState, pC5, pC4, 1, cylinderProps);
                        addCylinder(builderState, pC4, pN3, 1, cylinderProps);
                        addCylinder(builderState, pN3, pC2, 1, cylinderProps);
                        addCylinder(builderState, pC2, pN1, 1, cylinderProps);
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


export function NucleotideAtomicBondMeshVisual(materialId: number): UnitsVisual<NucleotideAtomicBondParams> {
    return UnitsMeshVisual<NucleotideAtomicBondParams>({
        defaultProps: PD.getDefaultValues(NucleotideAtomicBondParams),
        createGeometry: createNucleotideAtomicBondMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideAtomicBondParams>, currentProps: PD.Values<NucleotideAtomicBondParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideAtomicBondParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
