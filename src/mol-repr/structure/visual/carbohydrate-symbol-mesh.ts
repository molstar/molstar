/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Interval, OrderedSet } from '../../../mol-data/int';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { addSphere } from '../../../mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { PerforatedBox } from '../../../mol-geo/primitive/box';
import { PerforatedOctahedron } from '../../../mol-geo/primitive/octahedron';
import { PolygonalPrism, } from '../../../mol-geo/primitive/prism';
import { PerforatedOctagonalPyramid } from '../../../mol-geo/primitive/pyramid';
import { Star } from '../../../mol-geo/primitive/star';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { getSaccharideShape, SaccharideShape } from '../../../mol-model/structure/structure/carbohydrates/constants';
import { VisualContext } from '../../../mol-repr/visual';
import { Theme } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualUpdateState } from '../../util';
import { ComplexMeshParams, ComplexMeshVisual } from '../complex-visual';
import { ComplexVisual } from '../representation';
import { getAltResidueLociFromId } from './util/common';


const t = Mat4.identity();
const sVec = Vec3();
const pd = Vec3();

const SideFactor = 2 * 0.806; // 0.806 == Math.cos(Math.PI / 4)

const perforatedBox = PerforatedBox();
const perforatedOctagonalPyramid = PerforatedOctagonalPyramid();
const perforatedStar = Star({ outerRadius: 1, innerRadius: 0.5, thickness: 0.5, pointCount: 5, perforated: true });
const perforatedOctahedron = PerforatedOctahedron();
const diamondPrism = {
    caps: PolygonalPrism(4, { subset: 'caps' }),
    sides: PolygonalPrism(4, { subset: 'sides' }),
};
const pentagonalPrism = {
    caps: PolygonalPrism(5, { subset: 'caps' }),
    sides: PolygonalPrism(5, { subset: 'sides' }),
};
const hexagonalPrism = {
    caps: PolygonalPrism(6, { subset: 'caps' }),
    sides: PolygonalPrism(6, { subset: 'sides' }),
};
const shiftedHexagonalPrism = {
    caps: PolygonalPrism(6, { shifted: true, subset: 'caps' }),
    sides: PolygonalPrism(6, { shifted: true, subset: 'sides' }),
};
const heptagonalPrism = {
    caps: PolygonalPrism(7, { subset: 'caps' }),
    sides: PolygonalPrism(7, { subset: 'sides' }),
};

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
            case SaccharideShape.FilledSphere: // e.g. 3d11
                addSphere(builderState, center, radius, detail, { subset: 'ring' });
                builderState.currentGroup += 1;
                addSphere(builderState, center, radius, detail, { subset: 'caps' });
                break;
            case SaccharideShape.FilledCube: // e.g. 3d11
            case SaccharideShape.CrossedCube: // e.g. 5hwa
                Mat4.scaleUniformly(t, t, side);
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                Mat4.mul(t, t, Mat4.rotY90Z180);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                break;
            case SaccharideShape.FilledCone: // e.g. 9k1g
            case SaccharideShape.DevidedCone: // e.g. 4y9v
                Mat4.scaleUniformly(t, t, side * 1.2);
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid);
                Mat4.mul(t, t, Mat4.rotZ90);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedOctagonalPyramid);
                break;
            case SaccharideShape.FlatBox: // e.g. 1mfd
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.2, side * 0.6, side * 0.6));
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                Mat4.mul(t, t, Mat4.rotY90Z180);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedBox);
                break;
            case SaccharideShape.FilledStar: // e.g. 6rv7
                Mat4.scaleUniformly(t, t, side);
                Mat4.mul(t, t, Mat4.rotZY90);
                MeshBuilder.addPrimitive(builderState, t, perforatedStar);
                Mat4.mul(t, t, Mat4.rotX180);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedStar);
                break;
            case SaccharideShape.FilledDiamond: // e.g. 9q5e
            case SaccharideShape.DividedDiamond: // e.g. 1hv6
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side * 1.4, side * 1.4, side * 1.4));
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron);
                Mat4.mul(t, t, Mat4.rotY90);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, perforatedOctahedron);
                break;
            case SaccharideShape.FlatDiamond: // no CCD codes mapped
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side / 2, side / 2));
                MeshBuilder.addPrimitive(builderState, t, diamondPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, diamondPrism.sides);
                break;
            case SaccharideShape.DiamondPrism:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, diamondPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, diamondPrism.sides);
                break;
            case SaccharideShape.PentagonalPrism:
            case SaccharideShape.Pentagon: // e.g. 8jq5
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, pentagonalPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, pentagonalPrism.sides);
                break;
            case SaccharideShape.HexagonalPrism: // e.g. 9k1g
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, hexagonalPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, hexagonalPrism.sides);
                break;
            case SaccharideShape.HeptagonalPrism:
                Mat4.mul(t, t, Mat4.rotZY90);
                Mat4.scale(t, t, Vec3.set(sVec, side, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, heptagonalPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, heptagonalPrism.sides);
                break;
            case SaccharideShape.FlatHexagon: // e.g. 3k8d
            default:
                Mat4.mul(t, t, Mat4.rotZYZ90);
                Mat4.scale(t, t, Vec3.set(sVec, side / 1.5, side, side / 2));
                MeshBuilder.addPrimitive(builderState, t, shiftedHexagonalPrism.caps);
                builderState.currentGroup += 1;
                MeshBuilder.addPrimitive(builderState, t, shiftedHexagonalPrism.sides);
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
    function getLocation(groupIndex: number, instanceIndex: number) {
        const carb = carbElements[Math.floor(groupIndex / 2)];
        const ring = carb.unit.rings.all[carb.ringIndex];
        location.unit = carb.unit;
        location.element = carb.unit.elements[ring[0]];
        return location;
    }
    function isSecondary(elementIndex: number, instanceIndex: number) {
        return (elementIndex % 2) === 1;
    }
    return LocationIterator(groupCount, instanceCount, 1, getLocation, true, isSecondary);
}

/** Return a Loci for the elements of the whole residue of a carbohydrate. */
function getCarbohydrateLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        if (groupId === PickingId.Null) {
            return Structure.Loci(structure);
        } else {
            const carb = structure.carbohydrates.elements[Math.floor(groupId / 2)];
            return getAltResidueLociFromId(structure, carb.unit, carb.residueIndex, carb.altId);
        }
    }
    return EmptyLoci;
}

const __elementIndicesSet = new Set<number>();

/** For each carbohydrate (usually a monosaccharide) when all its residue's elements are in a loci. */
function eachCarbohydrate(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    const { getElementIndices } = structure.carbohydrates;
    let changed = false;
    if (!StructureElement.Loci.is(loci)) return false;
    if (!Structure.areEquivalent(loci.structure, structure)) return false;

    for (const { unit, indices } of loci.elements) {
        if (!Unit.isAtomic(unit)) continue;

        __elementIndicesSet.clear();
        OrderedSet.forEach(indices, v => {
            const elementIndices = getElementIndices(unit, unit.elements[v]);
            for (let i = 0, il = elementIndices.length; i < il; ++i) {
                if (!__elementIndicesSet.has(elementIndices[i])) {
                    __elementIndicesSet.add(elementIndices[i]);
                    const firstGroup = elementIndices[i] * 2;
                    const groupInterval = Interval.ofRange(firstGroup, firstGroup + 1); // each residue has up to 2 groups
                    if (apply(groupInterval)) changed = true;
                }
            }
        });
    }
    return changed;
}