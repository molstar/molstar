/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { Box, PerforatedBox } from '../../../mol-geo/primitive/box';
import { OctagonalPyramid, PerforatedOctagonalPyramid } from '../../../mol-geo/primitive/pyramid';
import { Star } from '../../../mol-geo/primitive/star';
import { Octahedron, PerforatedOctahedron } from '../../../mol-geo/primitive/octahedron';
import { DiamondPrism, PentagonalPrism, ShiftedHexagonalPrism, HexagonalPrism, HeptagonalPrism } from '../../../mol-geo/primitive/prism';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { getSaccharideShape, SaccharideShape } from '../../../mol-model/structure/structure/carbohydrates/constants';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { OrderedSet, Interval } from '../../../mol-data/int';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { getAltResidueLociFromId } from './util/common';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

const t = Mat4.identity();
const sVec = Vec3.zero();
const pd = Vec3.zero();

const SideFactor = 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)

const box = Box();
const perforatedBox = PerforatedBox();
const octagonalPyramid = OctagonalPyramid();
const perforatedOctagonalPyramid = PerforatedOctagonalPyramid();
const star = Star({ outerRadius: 1, innerRadius: 0.5, thickness: 0.5, pointCount: 5 });
const octahedron = Octahedron();
const perforatedOctahedron = PerforatedOctahedron();
const diamondPrism = DiamondPrism();
const pentagonalPrism = PentagonalPrism();
const hexagonalPrism = HexagonalPrism();
const shiftedHexagonalPrism = ShiftedHexagonalPrism();
const heptagonalPrism = HeptagonalPrism();

function createCarbohydrateSymbolMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CarbohydrateSymbolParams>, mesh?: Mesh) {
    const builderState = MeshBuilder.createState(256, 128, mesh);

    const { detail, sizeFactor } = props;

    const carbohydrates = structure.carbohydrates;
    const n = carbohydrates.elements.length;
    const l = StructureElement.Location.create(structure);

    for (let i = 0; i < n; ++i) {
        const c = carbohydrates.elements[i];
        const ring = c.unit.rings.all[c.ringIndex];
        const shapeType = getSaccharideShape(c.component.type, ring.length);

        l.unit = c.unit;
        l.element = c.unit.elements[ring[0]];
        const size = theme.size.size(l);
        const radius = size * sizeFactor;
        const side = size * sizeFactor * SideFactor;

        const { center, normal, direction } = c.geometry;
        Vec3.add(pd, center, direction);
        Mat4.targetTo(t, center, pd, normal);
        Mat4.setTranslation(t, center);

        builderState.currentGroup = i * 2;

        switch (shapeType) {
            case SaccharideShape.FilledSphere:
                addSphere(builderState, center, radius, detail);
                break;
            case SaccharideShape.FilledCube:
                Mat4.scaleUniformly(t, t, side);
                MeshBuilder.addPrimitive(builderState, t, box);
                break;
            case SaccharideShape.CrossedCube:
                Mat4.scaleUniformly(t, t, side);
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                Mat4.mul(t, t, Mat4.rotZ90X180);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                break;
            case SaccharideShape.FilledCone:
                Mat4.scaleUniformly(t, t, side * 1.2);
                MeshBuilder.addPrimitive(builderState, t, octagonalPyramid);
                break;
            case SaccharideShape.DevidedCone:
                Mat4.scaleUniformly(t, t, side * 1.2);
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid);
                Mat4.mul(t, t, Mat4.rotZ90);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid);
                break;
            case SaccharideShape.FlatBox:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, box);
                break;
            case SaccharideShape.FilledStar:
                Mat4.scaleUniformly(t, t, side);
                Mat4.mul(t, t, Mat4.rotZY90);
                MeshBuilder.addPrimitive(builderState, t, star);
                break;
            case SaccharideShape.FilledDiamond:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4));
                MeshBuilder.addPrimitive(builderState, t, octahedron);
                break;
            case SaccharideShape.DividedDiamond:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4));
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron);
                Mat4.mul(t, t, Mat4.rotY90);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron);
                break;
            case SaccharideShape.FlatDiamond:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side / 2, side / 2));
                MeshBuilder.addPrimitive(builderState, t, diamondPrism);
                break;
            case SaccharideShape.DiamondPrism:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, diamondPrism);
                break;
            case SaccharideShape.PentagonalPrism:
            case SaccharideShape.Pentagon:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, pentagonalPrism);
                break;
            case SaccharideShape.HexagonalPrism:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, hexagonalPrism);
                break;
            case SaccharideShape.HeptagonalPrism:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, heptagonalPrism);
                break;
            case SaccharideShape.FlatHexagon:
            default:
                Mat4.mul(t, t, Mat4.rotZYZ90);
                Mat4.scale(t, t, Vec3.set(sVec, side / 1.5, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, shiftedHexagonalPrism);
                break;
        }
    }

    return MeshBuilder.getMesh(builderState);
}

export const CarbohydrateSymbolParams = {
    ...ComplexMeshParams,
    detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    sizeFactor: PD.Numeric(1.75, { min: 0, max: 10, step: 0.01 }),
};
export type CarbohydrateSymbolParams = typeof CarbohydrateSymbolParams

export function CarbohydrateSymbolVisual(materialId: number): ComplexVisual<CarbohydrateSymbolParams> {
    return ComplexMeshVisual<CarbohydrateSymbolParams>({
        defaultProps: PD.getDefaultValues(CarbohydrateSymbolParams),
        createGeometry: createCarbohydrateSymbolMesh,
        createLocationIterator: CarbohydrateElementIterator,
        getLoci: getCarbohydrateLoci,
        eachLocation: eachCarbohydrate,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CarbohydrateSymbolParams>, currentProps: PD.Values<CarbohydrateSymbolParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.detail !== currentProps.detail
            );
        }
    }, materialId);
}

function CarbohydrateElementIterator(structure: Structure): LocationIterator {
    const carbElements = structure.carbohydrates.elements;
    const groupCount = carbElements.length * 2;
    const instanceCount = 1;
    const location = StructureElement.Location.create(structure);
    function getLocation (groupIndex: number, instanceIndex: number) {
        const carb = carbElements[Math.floor(groupIndex / 2)];
        const ring = carb.unit.rings.all[carb.ringIndex];
        location.unit = carb.unit;
        location.element = carb.unit.elements[ring[0]];
        return location;
    }
    function isSecondary (elementIndex: number, instanceIndex: number) {
        return (elementIndex % 2) === 1;
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true, isSecondary);
}

/** Return a Loci for the elements of the whole residue of a carbohydrate. */
function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const carb = structure.carbohydrates.elements[Math.floor(groupId / 2)];
        return getAltResidueLociFromId(structure, carb.unit, carb.residueIndex, carb.altId);
    }
    return EmptyLoci;
}

/** For each carbohydrate (usually a monosaccharide) when all its residue's elements are in a loci. */
function eachCarbohydrate(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getElementIndices } = structure.carbohydrates;
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;

        OrderedSet.forEach(indices, v => {
            // TODO avoid duplicate calls to apply
            const elementIndices = getElementIndices(unit, unit.elements[v]);
            for (let i = 0, il = elementIndices.length; i < il; ++i) {
                if (apply(Interval.ofSingleton(elementIndices[i] * 2))) changed = true;
            }
        });
    }
    return changed;
}