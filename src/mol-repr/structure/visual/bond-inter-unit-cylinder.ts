/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Bond, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags, arrayEqual } from '../../../mol-util';
import { createLinkCylinderImpostors, createLinkCylinderMesh, LinkBuilderProps, LinkStyle } from './util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual, ComplexCylindersParams, ComplexCylindersVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondCylinderParams, BondIterator, getInterBondLoci, eachInterBond, makeInterBondIgnoreTest, hasStructureVisibleBonds } from './util/bond';
import { Sphere3D } from '../../../mol-math/geometry';
import { Cylinders } from '../../../mol-geo/geometry/cylinders/cylinders';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { SortedArray } from '../../../mol-data/int/sorted-array';
import { SizeTheme } from '../../../mol-theme/size';
import { EmptyLocationIterator } from '../../../mol-geo/util/location-iterator';
import { checkCylinderImpostorSupport } from './util/common';

const tmpRefPosBondIt = new Bond.ElementBondIterator();
function setRefPosition(pos: Vec3, structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    tmpRefPosBondIt.setElement(structure, unit, index);
    while (tmpRefPosBondIt.hasNext) {
        const bA = tmpRefPosBondIt.move();
        bA.otherUnit.conformation.position(bA.otherUnit.elements[bA.otherIndex], pos);
        return pos;
    }
    return null;
}

const tmpRef = Vec3();

function getInterUnitBondCylinderBuilderProps(structure: Structure, theme: Theme, props: PD.Values<InterUnitBondCylinderParams>): LinkBuilderProps {
    const locE = StructureElement.Location.create(structure);
    const locB = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);

    const bonds = structure.interUnitBonds;
    const { edgeCount, edges } = bonds;
    const { sizeFactor, sizeAspectRatio, adjustCylinderLength, aromaticBonds, multipleBonds } = props;

    const mbOff = multipleBonds === 'off';
    const mbSymmetric = multipleBonds === 'symmetric';

    const delta = Vec3();

    let stub: undefined | ((edgeIndex: number) => boolean);

    const { child } = structure;
    if (props.includeParent && child) {
        stub = (edgeIndex: number) => {
            const b = edges[edgeIndex];
            const childUnitA = child.unitMap.get(b.unitA);
            const childUnitB = child.unitMap.get(b.unitB);

            const unitA = structure.unitMap.get(b.unitA);
            const eA = unitA.elements[b.indexA];
            const unitB = structure.unitMap.get(b.unitB);
            const eB = unitB.elements[b.indexB];

            return (
                childUnitA && SortedArray.has(childUnitA.elements, eA) &&
                (!childUnitB || !SortedArray.has(childUnitB.elements, eB))
            );
        };
    }

    const radius = (edgeIndex: number) => {
        const b = edges[edgeIndex];
        locB.aUnit = structure.unitMap.get(b.unitA);
        locB.aIndex = b.indexA;
        locB.bUnit = structure.unitMap.get(b.unitB);
        locB.bIndex = b.indexB;
        return theme.size.size(locB) * sizeFactor;
    };

    const radiusA = (edgeIndex: number) => {
        const b = edges[edgeIndex];
        locE.unit = structure.unitMap.get(b.unitA);
        locE.element = locE.unit.elements[b.indexA];
        return theme.size.size(locE) * sizeFactor;
    };

    const radiusB = (edgeIndex: number) => {
        const b = edges[edgeIndex];
        locE.unit = structure.unitMap.get(b.unitB);
        locE.element = locE.unit.elements[b.indexB];
        return theme.size.size(locE) * sizeFactor;
    };

    return {
        linkCount: edgeCount,
        referencePosition: (edgeIndex: number) => {
            const b = edges[edgeIndex];
            let unitA: Unit.Atomic, unitB: Unit.Atomic;
            let indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex;
            if (b.unitA < b.unitB) {
                unitA = structure.unitMap.get(b.unitA) as Unit.Atomic;
                unitB = structure.unitMap.get(b.unitB) as Unit.Atomic;
                indexA = b.indexA;
                indexB = b.indexB;
            } else if (b.unitA > b.unitB) {
                unitA = structure.unitMap.get(b.unitB) as Unit.Atomic;
                unitB = structure.unitMap.get(b.unitA) as Unit.Atomic;
                indexA = b.indexB;
                indexB = b.indexA;
            } else {
                throw new Error('same units in createInterUnitBondCylinderMesh');
            }
            return setRefPosition(tmpRef, structure, unitA, indexA) || setRefPosition(tmpRef, structure, unitB, indexB);
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number, adjust: boolean) => {
            const b = edges[edgeIndex];
            const uA = structure.unitMap.get(b.unitA);
            const uB = structure.unitMap.get(b.unitB);

            uA.conformation.position(uA.elements[b.indexA], posA);
            uB.conformation.position(uB.elements[b.indexB], posB);

            if (adjust && adjustCylinderLength) {
                const rA = radiusA(edgeIndex), rB = radiusB(edgeIndex);
                const r = Math.min(rA, rB) * sizeAspectRatio;
                const oA = Math.sqrt(Math.max(0, rA * rA - r * r)) - 0.05;
                const oB = Math.sqrt(Math.max(0, rB * rB - r * r)) - 0.05;
                if (oA <= 0.01 && oB <= 0.01) return;

                Vec3.normalize(delta, Vec3.sub(delta, posB, posA));
                Vec3.scaleAndAdd(posA, posA, delta, oA);
                Vec3.scaleAndAdd(posB, posB, delta, -oB);
            }
        },
        style: (edgeIndex: number) => {
            const o = edges[edgeIndex].props.order;
            const f = BitFlags.create(edges[edgeIndex].props.flag);
            if (BondType.is(f, BondType.Flag.MetallicCoordination) || BondType.is(f, BondType.Flag.HydrogenBond)) {
                // show metallic coordinations and hydrogen bonds with dashed cylinders
                return LinkStyle.Dashed;
            } else if (o === 3) {
                return mbOff ? LinkStyle.Solid :
                    mbSymmetric ? LinkStyle.Triple :
                        LinkStyle.OffsetTriple;
            } else if (aromaticBonds && BondType.is(f, BondType.Flag.Aromatic)) {
                return LinkStyle.Aromatic;
            }

            return (o !== 2 || mbOff) ? LinkStyle.Solid :
                mbSymmetric ? LinkStyle.Double :
                    LinkStyle.OffsetDouble;
        },
        radius: (edgeIndex: number) => {
            return radius(edgeIndex) * sizeAspectRatio;
        },
        ignore: makeInterBondIgnoreTest(structure, props),
        stub
    };
}

function createInterUnitBondCylinderImpostors(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitBondCylinderParams>, cylinders?: Cylinders) {
    if (!hasStructureVisibleBonds(structure, props)) return Cylinders.createEmpty(cylinders);
    if (!structure.interUnitBonds.edgeCount) return Cylinders.createEmpty(cylinders);

    const builderProps = getInterUnitBondCylinderBuilderProps(structure, theme, props);
    const { cylinders: c, boundingSphere } = createLinkCylinderImpostors(ctx, builderProps, props, cylinders);

    if (boundingSphere) {
        c.setBoundingSphere(boundingSphere);
    } else if (c.cylinderCount > 0) {
        const { child } = structure;
        const sphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * props.sizeFactor);
        c.setBoundingSphere(sphere);
    }

    return c;
}

function createInterUnitBondCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitBondCylinderParams>, mesh?: Mesh) {
    if (!hasStructureVisibleBonds(structure, props)) return Mesh.createEmpty(mesh);
    if (!structure.interUnitBonds.edgeCount) return Mesh.createEmpty(mesh);

    const builderProps = getInterUnitBondCylinderBuilderProps(structure, theme, props);
    const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    if (boundingSphere) {
        m.setBoundingSphere(boundingSphere);
    } else if (m.triangleCount > 0) {
        const { child } = structure;
        const sphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * props.sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

export const InterUnitBondCylinderParams = {
    ...ComplexMeshParams,
    ...ComplexCylindersParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0, max: 3, step: 0.01 }),
    tryUseImpostor: PD.Boolean(true),
    includeParent: PD.Boolean(false),
};
export type InterUnitBondCylinderParams = typeof InterUnitBondCylinderParams

export function InterUnitBondCylinderVisual(materialId: number, structure: Structure, props: PD.Values<InterUnitBondCylinderParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && checkCylinderImpostorSupport(webgl)
        ? InterUnitBondCylinderImpostorVisual(materialId)
        : InterUnitBondCylinderMeshVisual(materialId);
}

export function InterUnitBondCylinderImpostorVisual(materialId: number): ComplexVisual<InterUnitBondCylinderParams> {
    return ComplexCylindersVisual<InterUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(InterUnitBondCylinderParams),
        createGeometry: createInterUnitBondCylinderImpostors,
        createLocationIterator: (structure: Structure, props: PD.Values<InterUnitBondCylinderParams>) => {
            return !hasStructureVisibleBonds(structure, props)
                ? EmptyLocationIterator
                : BondIterator.fromStructure(structure, { includeLocation2: props.colorMode === 'interpolate' });
        },
        getLoci: getInterBondLoci,
        eachLocation: eachInterBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitBondCylinderParams>, currentProps: PD.Values<InterUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.sizeAspectRatio !== currentProps.sizeAspectRatio ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.linkCap !== currentProps.linkCap ||
                newProps.aromaticScale !== currentProps.aromaticScale ||
                newProps.aromaticSpacing !== currentProps.aromaticSpacing ||
                newProps.aromaticDashCount !== currentProps.aromaticDashCount ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.dashScale !== currentProps.dashScale ||
                newProps.dashCap !== currentProps.dashCap ||
                newProps.stubCap !== currentProps.stubCap ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes) ||
                newProps.adjustCylinderLength !== currentProps.adjustCylinderLength ||
                newProps.multipleBonds !== currentProps.multipleBonds
            );

            if (newProps.colorMode !== currentProps.colorMode) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
            }

            if (hasStructureVisibleBonds(newStructure, newProps) && newStructure.interUnitBonds !== currentStructure.interUnitBonds) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
                state.updateSize = true;
            }
        },
        mustRecreate: (structure: Structure, props: PD.Values<InterUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function InterUnitBondCylinderMeshVisual(materialId: number): ComplexVisual<InterUnitBondCylinderParams> {
    return ComplexMeshVisual<InterUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(InterUnitBondCylinderParams),
        createGeometry: createInterUnitBondCylinderMesh,
        createLocationIterator: (structure: Structure, props: PD.Values<InterUnitBondCylinderParams>) => {
            return !hasStructureVisibleBonds(structure, props)
                ? EmptyLocationIterator
                : BondIterator.fromStructure(structure);
        },
        getLoci: getInterBondLoci,
        eachLocation: eachInterBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitBondCylinderParams>, currentProps: PD.Values<InterUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.sizeAspectRatio !== currentProps.sizeAspectRatio ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                newProps.linkCap !== currentProps.linkCap ||
                newProps.aromaticScale !== currentProps.aromaticScale ||
                newProps.aromaticSpacing !== currentProps.aromaticSpacing ||
                newProps.aromaticDashCount !== currentProps.aromaticDashCount ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.dashScale !== currentProps.dashScale ||
                newProps.dashCap !== currentProps.dashCap ||
                newProps.stubCap !== currentProps.stubCap ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes) ||
                newProps.adjustCylinderLength !== currentProps.adjustCylinderLength ||
                newProps.multipleBonds !== currentProps.multipleBonds ||
                newProps.adjustCylinderLength && !SizeTheme.areEqual(newTheme.size, currentTheme.size)
            );

            if (hasStructureVisibleBonds(newStructure, newProps) && newStructure.interUnitBonds !== currentStructure.interUnitBonds) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
                state.updateSize = true;
            }
        },
        mustRecreate: (structure: Structure, props: PD.Values<InterUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
