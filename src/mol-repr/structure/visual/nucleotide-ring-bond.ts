import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3 } from '../../mol-math/linear-algebra';
import { VisualContext } from '../visual';
import { Unit, Structure, ElementIndex } from '../../mol-model/structure';
import { Theme } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../mol-data/int';
import { CylinderProps } from '../../mol-geo/primitive/cylinder';
import { isNucleic, isPurineBase, isPyrimidineBase } from '../../mol-model/structure/model/types';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsCylindersParams, UnitsCylindersVisual } from '../structure/units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement } from '../structure/visual/util/nucleotide';
import { VisualUpdateState } from '../util';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Sphere3D } from '../../mol-math/geometry';

import { WebGLContext } from '../../mol-gl/webgl/context';

import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '../../mol-geo/geometry/cylinders/cylinders-builder';
import { StructureGroup } from './util/common';


// avoiding namespace lookup improved performance in Chrome (Aug 2020)

const pTrace = Vec3.zero();

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

export const NucleotideRingBondParams = {
    ...UnitsMeshParams,
    ...UnitsCylindersParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true)
};
export type NucleotideRingBondParams = typeof NucleotideRingBondParams
interface NucleotideRingBondImpostorProps {
    sizeFactor: number,
}

export function NucleotideRingBondVisual(materialId: number, structure: Structure, props: PD.Values<NucleotideRingBondParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth
        ? NucleotideRingBondImpostorVisual(materialId)
        : NucleotideRingBondMeshVisual(materialId);
}

function createNucleotideRingBondImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingBondImpostorProps, cylinders?: Cylinders) {
    if (!Unit.isAtomic(unit)) return Cylinders.createEmpty(cylinders);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Cylinders.createEmpty(cylinders);

    const cylindersCountEstimate = nucleotideElementCount * 15; // 15 is the average purine (17) & pirimidine (13) bonds
    const builder = CylindersBuilder.create(cylindersCountEstimate, cylindersCountEstimate / 4, cylinders);

    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments, atoms, index: atomicIndex } = model.atomicHierarchy;

    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const { label_comp_id } = atoms;
    const pos = unit.conformation.invariantPosition;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const compId = label_comp_id.value(residueAtomSegments.offsets[residueIndex]);

                let idxTrace: ElementIndex | -1 = -1, idxN1: ElementIndex | -1 = -1, idxC2: ElementIndex | -1 = -1, idxN3: ElementIndex | -1 = -1, idxC4: ElementIndex | -1 = -1, idxC5: ElementIndex | -1 = -1, idxC6: ElementIndex | -1 = -1, idxN7: ElementIndex | -1 = -1, idxC8: ElementIndex | -1 = -1, idxN9: ElementIndex | -1 = -1,
                    idxC1_1: ElementIndex | -1 = -1, idxC2_1: ElementIndex | -1 = -1, idxC3_1: ElementIndex | -1 = -1, idxC4_1: ElementIndex | -1 = -1, idxO4_1: ElementIndex | -1 = -1;

                idxTrace = traceElementIndex[residueIndex];

                // sugar base
                idxC1_1 = atomicIndex.findAtomOnResidue(residueIndex, "C1'");
                idxC2_1 = atomicIndex.findAtomOnResidue(residueIndex, "C2'");
                idxC3_1 = atomicIndex.findAtomOnResidue(residueIndex, "C3'");
                idxC4_1 = atomicIndex.findAtomOnResidue(residueIndex, "C4'");
                idxO4_1 = atomicIndex.findAtomOnResidue(residueIndex, "O4'");
                if (idxC1_1 !== -1 && idxC2_1 !== -1 && idxC3_1 !== -1 && idxC4_1 !== -1 && idxO4_1 !== -1) {
                    pos(idxC1_1, pC1_1); pos(idxC2_1, pC2_1); pos(idxC3_1, pC3_1); pos(idxC4_1, pC4_1); pos(idxO4_1, pO4_1);

                    // trace cylinder
                    if (idxTrace !== -1) {
                        pos(idxTrace, pTrace);
                        builder.add(pC3_1[0], pC3_1[1], pC3_1[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, i);
                    }

                    // sugar ring
                    builder.add(pC3_1[0], pC3_1[1], pC3_1[2], pC4_1[0], pC4_1[1], pC4_1[2], 1, true, true, i);
                    builder.add(pC4_1[0], pC4_1[1], pC4_1[2], pO4_1[0], pO4_1[1], pO4_1[2], 1, true, true, i);
                    builder.add(pO4_1[0], pO4_1[1], pO4_1[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, i);
                    builder.add(pC1_1[0], pC1_1[1], pC1_1[2], pC2_1[0], pC2_1[1], pC2_1[2], 1, true, true, i);
                    builder.add(pC2_1[0], pC2_1[1], pC2_1[2], pC3_1[0], pC3_1[1], pC3_1[2], 1, true, true, i);
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

                    if (idxC1_1 !== -1 && idxN9 !== -1) {
                        pos(idxC1_1, pC1_1); pos(idxN9, pN9);
                        builder.add(pN9[0], pN9[1], pN9[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, i);
                    } else if (idxN9 !== -1 && idxTrace !== -1) {
                        pos(idxN9, pN9); pos(idxTrace, pTrace);
                        builder.add(pN9[0], pN9[1], pN9[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, i);
                    }

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1 && idxN7 !== -1 && idxC8 !== -1 && idxN9 !== -1) {
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6); pos(idxN7, pN7); pos(idxC8, pC8);

                        // base ring
                        builder.add(pN9[0], pN9[1], pN9[2], pC8[0], pC8[1], pC8[2], 1, true, true, i);
                        builder.add(pC8[0], pC8[1], pC8[2], pN7[0], pN7[1], pN7[2], 1, true, true, i);
                        builder.add(pN7[0], pN7[1], pN7[2], pC5[0], pC5[1], pC5[2], 1, true, true, i);
                        builder.add(pC5[0], pC5[1], pC5[2], pC6[0], pC6[1], pC6[2], 1, true, true, i);
                        builder.add(pC6[0], pC6[1], pC6[2], pN1[0], pN1[1], pN1[2], 1, true, true, i);
                        builder.add(pN1[0], pN1[1], pN1[2], pC2[0], pC2[1], pC2[2], 1, true, true, i);
                        builder.add(pC2[0], pC2[1], pC2[2], pN3[0], pN3[1], pN3[2], 1, true, true, i);
                        builder.add(pN3[0], pN3[1], pN3[2], pC4[0], pC4[1], pC4[2], 1, true, true, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pC5[0], pC5[1], pC5[2], 1, true, true, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pN9[0], pN9[1], pN9[2], 1, true, true, i);

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

                    if (idxC1_1 !== -1 && idxN1 !== -1) {
                        pos(idxN1, pN1); pos(idxC1_1, pC1_1);
                        builder.add(pN1[0], pN1[1], pN1[2], pC1_1[0], pC1_1[1], pC1_1[2], 1, true, true, i);
                    } else if (idxN1 !== -1 && idxTrace !== -1) {
                        pos(idxN1, pN1); pos(idxTrace, pTrace);
                        builder.add(pN1[0], pN1[1], pN1[2], pTrace[0], pTrace[1], pTrace[2], 1, true, true, i);
                    }

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1) {
                        pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6);

                        // base ring
                        builder.add(pN1[0], pN1[1], pN1[2], pC6[0], pC6[1], pC6[2], 1, true, true, i);
                        builder.add(pC6[0], pC6[1], pC6[2], pC5[0], pC5[1], pC5[2], 1, true, true, i);
                        builder.add(pC5[0], pC5[1], pC5[2], pC4[0], pC4[1], pC4[2], 1, true, true, i);
                        builder.add(pC4[0], pC4[1], pC4[2], pN3[0], pN3[1], pN3[2], 1, true, true, i);
                        builder.add(pN3[0], pN3[1], pN3[2], pC2[0], pC2[1], pC2[2], 1, true, true, i);
                        builder.add(pC2[0], pC2[1], pC2[2], pN1[0], pN1[1], pN1[2], 1, true, true, i);
                    }
                }

                ++i;
            }

        }
    }
    const c = builder.getCylinders();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    c.setBoundingSphere(sphere);

    return c;
}

export function NucleotideRingBondImpostorVisual(materialId: number): UnitsVisual<NucleotideRingBondParams> {
    return UnitsCylindersVisual<NucleotideRingBondParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingBondParams),
        createGeometry: createNucleotideRingBondImpostor,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingBondParams>, currentProps: PD.Values<NucleotideRingBondParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideRingBondParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

interface NucleotideRingBondMeshProps {
    radialSegments: number,
    sizeFactor: number,
}

function createNucleotideRingBondMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingBondMeshProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, radialSegments } = props;

    const vertexCount = nucleotideElementCount * (radialSegments * 15); // 15 is the average purine (17) & pirimidine (13) bonds
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments, atoms, index: atomicIndex } = model.atomicHierarchy;
    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const { label_comp_id } = atoms;
    const pos = unit.conformation.invariantPosition;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const cylinderProps: CylinderProps = { radiusTop: 1 * sizeFactor, radiusBottom: 1 * sizeFactor, radialSegments };

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const compId = label_comp_id.value(residueAtomSegments.offsets[residueIndex]);

                let idxTrace: ElementIndex | -1 = -1, idxN1: ElementIndex | -1 = -1, idxC2: ElementIndex | -1 = -1, idxN3: ElementIndex | -1 = -1, idxC4: ElementIndex | -1 = -1, idxC5: ElementIndex | -1 = -1, idxC6: ElementIndex | -1 = -1, idxN7: ElementIndex | -1 = -1, idxC8: ElementIndex | -1 = -1, idxN9: ElementIndex | -1 = -1,
                    idxC1_1: ElementIndex | -1 = -1, idxC2_1: ElementIndex | -1 = -1, idxC3_1: ElementIndex | -1 = -1, idxC4_1: ElementIndex | -1 = -1, idxO4_1: ElementIndex | -1 = -1;

                builderState.currentGroup = i;

                idxTrace = traceElementIndex[residueIndex];

                // sugar base
                idxC1_1 = atomicIndex.findAtomOnResidue(residueIndex, "C1'");
                idxC2_1 = atomicIndex.findAtomOnResidue(residueIndex, "C2'");
                idxC3_1 = atomicIndex.findAtomOnResidue(residueIndex, "C3'");
                idxC4_1 = atomicIndex.findAtomOnResidue(residueIndex, "C4'");
                idxO4_1 = atomicIndex.findAtomOnResidue(residueIndex, "O4'");
                if (idxC1_1 !== -1 && idxC2_1 !== -1 && idxC3_1 !== -1 && idxC4_1 !== -1 && idxO4_1 !== -1) {
                    pos(idxC1_1, pC1_1); pos(idxC2_1, pC2_1); pos(idxC3_1, pC3_1); pos(idxC4_1, pC4_1); pos(idxO4_1, pO4_1);

                    // trace cylinder
                    if (idxTrace !== -1) {
                        pos(idxTrace, pTrace);
                        addCylinder(builderState, pC3_1, pTrace, 1, cylinderProps);
                    }

                    // sugar ring
                    addCylinder(builderState, pC3_1, pC4_1, 1, cylinderProps);
                    addCylinder(builderState, pC4_1, pO4_1, 1, cylinderProps);
                    addCylinder(builderState, pO4_1, pC1_1, 1, cylinderProps);
                    addCylinder(builderState, pC1_1, pC2_1, 1, cylinderProps);
                    addCylinder(builderState, pC2_1, pC3_1, 1, cylinderProps);
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

                    if (idxC1_1 !== -1 && idxN9 !== -1) {
                        pos(idxC1_1, pC1_1); pos(idxN9, pN9);
                        addCylinder(builderState, pN9, pC1_1, 1, cylinderProps);
                    } else if (idxN9 !== -1 && idxTrace !== -1) {
                        pos(idxN9, pN9); pos(idxTrace, pTrace);
                        addCylinder(builderState, pN9, pTrace, 1, cylinderProps);
                    }

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1 && idxN7 !== -1 && idxC8 !== -1 && idxN9 !== -1) {
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6); pos(idxN7, pN7); pos(idxC8, pC8);

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

                    if (idxC1_1 !== -1 && idxN1 !== -1) {
                        pos(idxN1, pN1); pos(idxC1_1, pC1_1);
                        addCylinder(builderState, pN1, pC1_1, 1, cylinderProps);
                    } else if (idxN1 !== -1 && idxTrace !== -1) {
                        pos(idxN1, pN1); pos(idxTrace, pTrace);
                        addCylinder(builderState, pN1, pTrace, 1, cylinderProps);
                    }

                    if (idxN1 !== -1 && idxC2 !== -1 && idxN3 !== -1 && idxC4 !== -1 && idxC5 !== -1 && idxC6 !== -1) {
                        pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6);

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


export function NucleotideRingBondMeshVisual(materialId: number): UnitsVisual<NucleotideRingBondParams> {
    return UnitsMeshVisual<NucleotideRingBondParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingBondParams),
        createGeometry: createNucleotideRingBondMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingBondParams>, currentProps: PD.Values<NucleotideRingBondParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideRingBondParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
