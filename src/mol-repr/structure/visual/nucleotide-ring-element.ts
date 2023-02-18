import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Vec3 } from '../../mol-math/linear-algebra';
import { VisualContext } from '../visual';
import { Unit, Structure, ElementIndex } from '../../mol-model/structure';
import { Theme } from '../../mol-theme/theme';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../mol-data/int';
import { isNucleic, isPurineBase, isPyrimidineBase } from '../../mol-model/structure/model/types';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsSpheresParams, UnitsSpheresVisual } from '../structure/units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement } from '../structure/visual/util/nucleotide';
import { VisualUpdateState } from '../util';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Sphere3D } from '../../mol-math/geometry';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { sphereVertexCount } from '../../mol-geo/primitive/sphere';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { StructureGroup } from '../structure/visual/util/common';

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

export const NucleotideRingElementParams = {
    ...UnitsMeshParams,
    ...UnitsSpheresParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    tryUseImpostor: PD.Boolean(true)
};
export type NucleotideRingElementParams = typeof NucleotideRingElementParams
interface NucleotideRingElementImpostorProps {
    sizeFactor: number,
}

export function NucleotideRingElementVisual(materialId: number, structure: Structure, props: PD.Values<NucleotideRingElementParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth
        ? NucleotideRingElementImpostorVisual(materialId)
        : NucleotideRingElementMeshVisual(materialId);
}

function createNucleotideRingElementImpostor(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingElementImpostorProps, spheres?: Spheres) {
    if (!Unit.isAtomic(unit)) return Spheres.createEmpty(spheres);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Spheres.createEmpty(spheres);

    const spheresCountEstimate = nucleotideElementCount * 15; // 15 is the average purine (17) & pirimidine (13) bonds
    const builder = SpheresBuilder.create(spheresCountEstimate, spheresCountEstimate / 4, spheres);

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
                        builder.add(pTrace[0], pTrace[1], pTrace[2], i);
                    }

                    // sugar ring
                    builder.add(pC3_1[0], pC3_1[1], pC3_1[2], i);
                    builder.add(pC4_1[0], pC4_1[1], pC4_1[2], i);
                    builder.add(pO4_1[0], pO4_1[1], pO4_1[2], i);
                    builder.add(pC1_1[0], pC1_1[1], pC1_1[2], i);
                    builder.add(pC2_1[0], pC2_1[1], pC2_1[2], i);
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
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6); pos(idxN7, pN7); pos(idxC8, pC8); pos(idxN9, pN9);

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
    const c = builder.getSpheres();

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * props.sizeFactor);
    c.setBoundingSphere(sphere);

    return c;
}

export function NucleotideRingElementImpostorVisual(materialId: number): UnitsVisual<NucleotideRingElementParams> {
    return UnitsSpheresVisual<NucleotideRingElementParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingElementParams),
        createGeometry: createNucleotideRingElementImpostor,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingElementParams>, currentProps: PD.Values<NucleotideRingElementParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideRingElementParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

interface NucleotideRingElementMeshProps {
    detail: number,
    sizeFactor: number,
}

function createNucleotideRingElementMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideRingElementMeshProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, detail } = props;

    const vertexCount = nucleotideElementCount * sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);

    const { elements, model } = unit;
    const { chainAtomSegments, residueAtomSegments, atoms, index: atomicIndex } = model.atomicHierarchy;
    const { moleculeType, traceElementIndex } = model.atomicHierarchy.derived.residue;
    const { label_comp_id } = atoms;
    const pos = unit.conformation.invariantPosition;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const radius = 1 * sizeFactor;

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
                        addSphere(builderState, pTrace, radius, detail);
                    }

                    // sugar ring
                    addSphere(builderState, pC4_1, radius, detail);
                    addSphere(builderState, pO4_1, radius, detail);
                    addSphere(builderState, pC1_1, radius, detail);
                    addSphere(builderState, pC2_1, radius, detail);
                    addSphere(builderState, pC3_1, radius, detail);
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
                        pos(idxN1, pN1); pos(idxC2, pC2); pos(idxN3, pN3); pos(idxC4, pC4); pos(idxC5, pC5); pos(idxC6, pC6); pos(idxN7, pN7); pos(idxC8, pC8); pos(idxN9, pN9);

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


export function NucleotideRingElementMeshVisual(materialId: number): UnitsVisual<NucleotideRingElementParams> {
    return UnitsMeshVisual<NucleotideRingElementParams>({
        defaultProps: PD.getDefaultValues(NucleotideRingElementParams),
        createGeometry: createNucleotideRingElementMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideRingElementParams>, currentProps: PD.Values<NucleotideRingElementParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<NucleotideRingElementParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
