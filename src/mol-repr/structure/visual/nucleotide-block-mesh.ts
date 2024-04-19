/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Vec3, Mat4 } from '../../../mol-math/linear-algebra';
import { Box } from '../../../mol-geo/primitive/box';
import { VisualContext } from '../../visual';
import { Unit, Structure, ElementIndex } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { Segmentation } from '../../../mol-data/int';
import { CylinderProps } from '../../../mol-geo/primitive/cylinder';
import { isNucleic } from '../../../mol-model/structure/model/types';
import { addCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual } from '../units-visual';
import { NucleotideLocationIterator, getNucleotideElementLoci, eachNucleotideElement, getNucleotideBaseType, createNucleicIndices, setPurinIndices, setPyrimidineIndices } from './util/nucleotide';
import { VisualUpdateState } from '../../util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Sphere3D } from '../../../mol-math/geometry';

// TODO support blocks for multiple locations (including from microheterogeneity)

const p1 = Vec3();
const p2 = Vec3();
const p3 = Vec3();
const p4 = Vec3();
const p5 = Vec3();
const pt = Vec3();
const v12 = Vec3();
const v34 = Vec3();
const vC = Vec3();
const center = Vec3();
const t = Mat4.identity();
const sVec = Vec3();
const box = Box();

export const NucleotideBlockMeshParams = {
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    thicknessFactor: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
};
export const DefaultNucleotideBlockMeshProps = PD.getDefaultValues(NucleotideBlockMeshParams);
export type NucleotideBlockMeshProps = typeof DefaultNucleotideBlockMeshProps

function createNucleotideBlockMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: NucleotideBlockMeshProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const nucleotideElementCount = unit.nucleotideElements.length;
    if (!nucleotideElementCount) return Mesh.createEmpty(mesh);

    const { sizeFactor, thicknessFactor, radialSegments } = props;

    const vertexCount = nucleotideElementCount * (box.vertices.length / 3 + radialSegments * 2);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 4, mesh);

    const { elements, model, conformation: c } = unit;
    const { chainAtomSegments, residueAtomSegments } = model.atomicHierarchy;
    const { moleculeType } = model.atomicHierarchy.derived.residue;

    const chainIt = Segmentation.transientSegments(chainAtomSegments, elements);
    const residueIt = Segmentation.transientSegments(residueAtomSegments, elements);

    const radius = 1 * sizeFactor;
    const width = 4.5;
    const depth = thicknessFactor * sizeFactor * 2;

    const cylinderProps: CylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments, bottomCap: true };

    let i = 0;
    while (chainIt.hasNext) {
        residueIt.setSegment(chainIt.move());

        while (residueIt.hasNext) {
            const { index: residueIndex } = residueIt.move();

            if (isNucleic(moleculeType[residueIndex])) {
                const idx = createNucleicIndices();
                let idx1: ElementIndex | -1 = -1, idx2: ElementIndex | -1 = -1, idx3: ElementIndex | -1 = -1, idx4: ElementIndex | -1 = -1, idx5: ElementIndex | -1 = -1;

                let height = 4.5;

                const { isPurine, isPyrimidine } = getNucleotideBaseType(unit, residueIndex);

                if (isPurine) {
                    height = 4.5;
                    setPurinIndices(idx, unit, residueIndex);
                    idx1 = idx.N1; idx2 = idx.C4; idx3 = idx.C6; idx4 = idx.C2; idx5 = idx.N9;
                } else if (isPyrimidine) {
                    height = 3.0;
                    setPyrimidineIndices(idx, unit, residueIndex);
                    idx1 = idx.N3; idx2 = idx.C6; idx3 = idx.C4; idx4 = idx.C2; idx5 = idx.N1;
                }

                if (idx5 !== -1 && idx.trace !== -1) {
                    c.invariantPosition(idx5, p5); c.invariantPosition(idx.trace, pt);
                    builderState.currentGroup = i;
                    addCylinder(builderState, p5, pt, 1, cylinderProps);
                    if (idx1 !== -1 && idx2 !== -1 && idx3 !== -1 && idx4 !== -1) {
                        c.invariantPosition(idx1, p1); c.invariantPosition(idx2, p2); c.invariantPosition(idx3, p3); c.invariantPosition(idx4, p4);
                        Vec3.normalize(v12, Vec3.sub(v12, p2, p1));
                        Vec3.normalize(v34, Vec3.sub(v34, p4, p3));
                        Vec3.normalize(vC, Vec3.cross(vC, v12, v34));
                        Mat4.targetTo(t, p1, p2, vC);
                        Vec3.scaleAndAdd(center, p1, v12, height / 2 - 0.2);
                        Mat4.scale(t, t, Vec3.set(sVec, width, depth, height));
                        Mat4.setTranslation(t, center);
                        MeshBuilder.addPrimitive(builderState, t, box);
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

export const NucleotideBlockParams = {
    ...UnitsMeshParams,
    ...NucleotideBlockMeshParams
};
export type NucleotideBlockParams = typeof NucleotideBlockParams

export function NucleotideBlockVisual(materialId: number): UnitsVisual<NucleotideBlockParams> {
    return UnitsMeshVisual<NucleotideBlockParams>({
        defaultProps: PD.getDefaultValues(NucleotideBlockParams),
        createGeometry: createNucleotideBlockMesh,
        createLocationIterator: NucleotideLocationIterator.fromGroup,
        getLoci: getNucleotideElementLoci,
        eachLocation: eachNucleotideElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NucleotideBlockParams>, currentProps: PD.Values<NucleotideBlockParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.thicknessFactor !== currentProps.thicknessFactor ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        }
    }, materialId);
}