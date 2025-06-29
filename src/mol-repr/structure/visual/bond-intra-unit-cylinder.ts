/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Herman Bergwerf <post@hbergwerf.nl>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement, Bond } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { arrayEqual } from '../../../mol-util';
import { createLinkCylinderImpostors, createLinkCylinderMesh, EmptyLinkBuilderProps, LinkBuilderProps, LinkStyle } from './util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, UnitsCylindersParams, UnitsCylindersVisual } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondCylinderParams, BondIterator, eachIntraBond, eachStructureGroupsBond, getIntraBondLoci, getStructureGroupsBondLoci, hasStructureVisibleBonds, hasUnitVisibleBonds, ignoreBondType, makeIntraBondIgnoreTest } from './util/bond';
import { Sphere3D } from '../../../mol-math/geometry';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { Cylinders } from '../../../mol-geo/geometry/cylinders/cylinders';
import { SortedArray } from '../../../mol-data/int';
import { arrayIntersectionSize } from '../../../mol-util/array';
import { checkCylinderImpostorSupport, StructureGroup } from './util/common';
import { SizeTheme } from '../../../mol-theme/size';
import { ComplexCylindersParams, ComplexMeshParams, ComplexCylindersVisual, ComplexMeshVisual, ComplexVisual } from '../complex-visual';
import { EmptyLocationIterator } from '../../../mol-geo/util/location-iterator';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const isBondType = BondType.is;

function getIntraUnitBondCylinderBuilderProps(unit: Unit.Atomic, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondCylinderParams>): LinkBuilderProps {
    const elements = unit.elements;
    const bonds = unit.bonds;
    const { edgeCount, a, b, edgeProps, offset } = bonds;
    const { order: _order, flags: _flags } = edgeProps;
    const { sizeFactor, sizeAspectRatio, adjustCylinderLength, aromaticBonds, includeTypes, excludeTypes, multipleBonds } = props;

    const mbOff = multipleBonds === 'off';
    const mbSymmetric = multipleBonds === 'symmetric';

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);
    const ignoreComputedAromatic = ignoreBondType(include, exclude, BondType.Flag.Computed);

    const vRef = Vec3(), delta = Vec3();
    const c = unit.conformation;

    let stub: undefined | ((edgeIndex: number) => boolean);

    const locE = StructureElement.Location.create(structure, unit);
    const locB = Bond.Location(structure, unit, undefined, structure, unit, undefined);

    const { child } = structure;
    if (props.includeParent && child) {
        const childUnit = child.unitMap.get(unit.id);
        if (!childUnit) throw new Error('expected childUnit to exist');

        stub = (edgeIndex: number) => {
            const eA = elements[a[edgeIndex]];
            const eB = elements[b[edgeIndex]];
            return SortedArray.has(childUnit.elements, eA) && !SortedArray.has(childUnit.elements, eB);
        };
    }

    const radius = (edgeIndex: number) => {
        locB.aIndex = a[edgeIndex];
        locB.bIndex = b[edgeIndex];
        return theme.size.size(locB) * sizeFactor;
    };

    const radiusA = (edgeIndex: number) => {
        locE.element = elements[a[edgeIndex]];
        return theme.size.size(locE) * sizeFactor;
    };

    const radiusB = (edgeIndex: number) => {
        locE.element = elements[b[edgeIndex]];
        return theme.size.size(locE) * sizeFactor;
    };

    const { elementRingIndices, elementAromaticRingIndices } = unit.rings;
    const deloTriplets = aromaticBonds ? unit.resonance.delocalizedTriplets : undefined;

    return {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];

            const rI = deloTriplets?.getThirdElement(aI, bI);
            if (rI !== undefined) return c.invariantPosition(elements[rI], vRef);

            if (aI > bI) [aI, bI] = [bI, aI];
            if (offset[aI + 1] - offset[aI] === 1) [aI, bI] = [bI, aI];

            const aR = elementAromaticRingIndices.get(aI) || elementRingIndices.get(aI);
            let maxSize = 0;

            for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
                const _bI = b[i];
                if (_bI !== bI && _bI !== aI) {
                    if (aR) {
                        const _bR = elementAromaticRingIndices.get(_bI) || elementRingIndices.get(_bI);
                        if (!_bR) continue;

                        const size = arrayIntersectionSize(aR, _bR);
                        if (size > maxSize) {
                            maxSize = size;
                            c.invariantPosition(elements[_bI], vRef);
                        }
                    } else {
                        return c.invariantPosition(elements[_bI], vRef);
                    }
                }
            }
            return maxSize > 0 ? vRef : null;
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number, adjust: boolean) => {
            c.invariantPosition(elements[a[edgeIndex]], posA);
            c.invariantPosition(elements[b[edgeIndex]], posB);

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
            const o = _order[edgeIndex];
            const f = _flags[edgeIndex];
            if (isBondType(f, BondType.Flag.MetallicCoordination) || isBondType(f, BondType.Flag.HydrogenBond)) {
                // show metallic coordinations and hydrogen bonds with dashed cylinders
                return LinkStyle.Dashed;
            } else if (o === 3) {
                return mbOff ? LinkStyle.Solid :
                    mbSymmetric ? LinkStyle.Triple :
                        LinkStyle.OffsetTriple;
            } else if (aromaticBonds) {
                const aI = a[edgeIndex], bI = b[edgeIndex];
                const aR = elementAromaticRingIndices.get(aI);
                const bR = elementAromaticRingIndices.get(bI);
                const arCount = (aR && bR) ? arrayIntersectionSize(aR, bR) : 0;

                if (isBondType(f, BondType.Flag.Aromatic) || (arCount && !ignoreComputedAromatic)) {
                    if (arCount === 2) {
                        return LinkStyle.MirroredAromatic;
                    } else {
                        return LinkStyle.Aromatic;
                    }
                }
            }

            return (o !== 2 || mbOff) ? LinkStyle.Solid :
                mbSymmetric ? LinkStyle.Double :
                    LinkStyle.OffsetDouble;
        },
        radius: (edgeIndex: number) => {
            return radius(edgeIndex) * sizeAspectRatio;
        },
        ignore: makeIntraBondIgnoreTest(structure, unit, props),
        stub
    };
}

function createIntraUnitBondCylinderImpostors(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondCylinderParams>, cylinders?: Cylinders): Cylinders {
    if (!Unit.isAtomic(unit)) return Cylinders.createEmpty(cylinders);
    if (!hasUnitVisibleBonds(unit, props)) return Cylinders.createEmpty(cylinders);
    if (!unit.bonds.edgeCount) return Cylinders.createEmpty(cylinders);

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Cylinders.createEmpty(cylinders);

    const builderProps = getIntraUnitBondCylinderBuilderProps(unit, structure, theme, props);
    const { cylinders: c, boundingSphere } = createLinkCylinderImpostors(ctx, builderProps, props, cylinders);

    if (boundingSphere) {
        c.setBoundingSphere(boundingSphere);
    } else if (c.cylinderCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * props.sizeFactor);
        c.setBoundingSphere(sphere);
    }

    return c;
}

function createIntraUnitBondCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondCylinderParams>, mesh?: Mesh): Mesh {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);
    if (!hasUnitVisibleBonds(unit, props)) return Mesh.createEmpty(mesh);
    if (!unit.bonds.edgeCount) return Mesh.createEmpty(mesh);

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Mesh.createEmpty(mesh);

    const builderProps = getIntraUnitBondCylinderBuilderProps(unit, structure, theme, props);
    const { mesh: m, boundingSphere } = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    if (boundingSphere) {
        m.setBoundingSphere(boundingSphere);
    } else if (m.triangleCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * props.sizeFactor);
        m.setBoundingSphere(sphere);
    }

    return m;
}

export const IntraUnitBondCylinderParams = {
    ...UnitsMeshParams,
    ...UnitsCylindersParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0, max: 3, step: 0.01 }),
    tryUseImpostor: PD.Boolean(true),
    includeParent: PD.Boolean(false),
};
export type IntraUnitBondCylinderParams = typeof IntraUnitBondCylinderParams

export function IntraUnitBondCylinderVisual(materialId: number, structure: Structure, props: PD.Values<IntraUnitBondCylinderParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && checkCylinderImpostorSupport(webgl)
        ? IntraUnitBondCylinderImpostorVisual(materialId)
        : IntraUnitBondCylinderMeshVisual(materialId);
}

export function IntraUnitBondCylinderImpostorVisual(materialId: number): UnitsVisual<IntraUnitBondCylinderParams> {
    return UnitsCylindersVisual<IntraUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(IntraUnitBondCylinderParams),
        createGeometry: createIntraUnitBondCylinderImpostors,
        createLocationIterator: (structureGroup: StructureGroup, props) => BondIterator.fromGroup(structureGroup, { includeLocation2: props.colorMode === 'interpolate' }),
        getLoci: getIntraBondLoci,
        eachLocation: eachIntraBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitBondCylinderParams>, currentProps: PD.Values<IntraUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
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
                newProps.aromaticBonds !== currentProps.aromaticBonds ||
                newProps.multipleBonds !== currentProps.multipleBonds ||
                newProps.adjustCylinderLength && !SizeTheme.areEqual(newTheme.size, currentTheme.size)
            );

            if (newProps.colorMode !== currentProps.colorMode) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
            }

            const newUnit = newStructureGroup.group.units[0];
            const currentUnit = currentStructureGroup.group.units[0];
            if (Unit.isAtomic(newUnit) && Unit.isAtomic(currentUnit)) {
                if (!IntAdjacencyGraph.areEqual(newUnit.bonds, currentUnit.bonds)) {
                    state.createGeometry = true;
                    state.updateTransform = true;
                    state.updateColor = true;
                    state.updateSize = true;
                }
            }
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<IntraUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function IntraUnitBondCylinderMeshVisual(materialId: number): UnitsVisual<IntraUnitBondCylinderParams> {
    return UnitsMeshVisual<IntraUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(IntraUnitBondCylinderParams),
        createGeometry: createIntraUnitBondCylinderMesh,
        createLocationIterator: (structureGroup: StructureGroup) => BondIterator.fromGroup(structureGroup),
        getLoci: getIntraBondLoci,
        eachLocation: eachIntraBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitBondCylinderParams>, currentProps: PD.Values<IntraUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
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
                newProps.aromaticBonds !== currentProps.aromaticBonds ||
                newProps.multipleBonds !== currentProps.multipleBonds
            );

            const newUnit = newStructureGroup.group.units[0];
            const currentUnit = currentStructureGroup.group.units[0];
            if (Unit.isAtomic(newUnit) && Unit.isAtomic(currentUnit)) {
                if (!IntAdjacencyGraph.areEqual(newUnit.bonds, currentUnit.bonds)) {
                    state.createGeometry = true;
                    state.updateTransform = true;
                    state.updateColor = true;
                    state.updateSize = true;
                }
            }
        },
        mustRecreate: (structureGroup: StructureGroup, props: PD.Values<IntraUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}

//

function getStructureIntraUnitBondCylinderBuilderProps(structure: Structure, theme: Theme, props: PD.Values<StructureIntraUnitBondCylinderParams>): LinkBuilderProps {
    const intraUnitProps: { group: Unit.SymmetryGroup, props: LinkBuilderProps}[] = [];

    const { bondCount, unitIndex, unitEdgeIndex, unitGroupIndex } = structure.intraUnitBondMapping;
    const { child } = structure;

    for (const ug of structure.unitSymmetryGroups) {
        const unit = ug.units[0];
        const childUnit = child?.unitMap.get(unit.id);

        const p = Unit.isAtomic(unit) && !(child && !childUnit)
            ? getIntraUnitBondCylinderBuilderProps(unit, structure, theme, props)
            : EmptyLinkBuilderProps;
        intraUnitProps.push({ group: ug, props: p });
    }

    return {
        linkCount: bondCount,
        referencePosition: (edgeIndex: number) => {
            const { group, props } = intraUnitProps[unitIndex[edgeIndex]];
            if (!props.referencePosition) return null;

            const v = props.referencePosition(unitEdgeIndex[edgeIndex]);
            if (!v) return null;

            const u = group.units[unitGroupIndex[edgeIndex]];
            Vec3.transformMat4(v, v, u.conformation.operator.matrix);
            return v;
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number, adjust: boolean) => {
            const { group, props } = intraUnitProps[unitIndex[edgeIndex]];
            props.position(posA, posB, unitEdgeIndex[edgeIndex], adjust);
            const u = group.units[unitGroupIndex[edgeIndex]];
            Vec3.transformMat4(posA, posA, u.conformation.operator.matrix);
            Vec3.transformMat4(posB, posB, u.conformation.operator.matrix);
        },
        style: (edgeIndex: number) => {
            const { props } = intraUnitProps[unitIndex[edgeIndex]];
            return props.style ? props.style(unitEdgeIndex[edgeIndex]) : LinkStyle.Solid;
        },
        radius: (edgeIndex: number) => {
            const { props } = intraUnitProps[unitIndex[edgeIndex]];
            return props.radius(unitEdgeIndex[edgeIndex]);
        },
        ignore: (edgeIndex: number) => {
            const { props } = intraUnitProps[unitIndex[edgeIndex]];
            return props.ignore ? props.ignore(unitEdgeIndex[edgeIndex]) : false;
        },
        stub: (edgeIndex: number) => {
            const { props } = intraUnitProps[unitIndex[edgeIndex]];
            return props.stub ? props.stub(unitEdgeIndex[edgeIndex]) : false;
        }
    };
}

function createStructureIntraUnitBondCylinderImpostors(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureIntraUnitBondCylinderParams>, cylinders?: Cylinders) {
    if (!hasStructureVisibleBonds(structure, props)) return Cylinders.createEmpty(cylinders);
    if (!structure.intraUnitBondMapping.bondCount) return Cylinders.createEmpty(cylinders);

    const builderProps = getStructureIntraUnitBondCylinderBuilderProps(structure, theme, props);
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

function createStructureIntraUnitBondCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureIntraUnitBondCylinderParams>, mesh?: Mesh) {
    if (!hasStructureVisibleBonds(structure, props)) return Mesh.createEmpty(mesh);
    if (!structure.intraUnitBondMapping.bondCount) return Mesh.createEmpty(mesh);

    const builderProps = getStructureIntraUnitBondCylinderBuilderProps(structure, theme, props);
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

export const StructureIntraUnitBondCylinderParams = {
    ...ComplexMeshParams,
    ...ComplexCylindersParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0, max: 3, step: 0.01 }),
    tryUseImpostor: PD.Boolean(true),
    includeParent: PD.Boolean(false),
};
export type StructureIntraUnitBondCylinderParams = typeof StructureIntraUnitBondCylinderParams

export function StructureIntraUnitBondCylinderVisual(materialId: number, structure: Structure, props: PD.Values<StructureIntraUnitBondCylinderParams>, webgl?: WebGLContext) {
    return props.tryUseImpostor && webgl && webgl.extensions.fragDepth
        ? StructureIntraUnitBondCylinderImpostorVisual(materialId)
        : StructureIntraUnitBondCylinderMeshVisual(materialId);
}

export function StructureIntraUnitBondCylinderImpostorVisual(materialId: number): ComplexVisual<StructureIntraUnitBondCylinderParams> {
    return ComplexCylindersVisual<StructureIntraUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(StructureIntraUnitBondCylinderParams),
        createGeometry: createStructureIntraUnitBondCylinderImpostors,
        createLocationIterator: (structure: Structure, props: PD.Values<StructureIntraUnitBondCylinderParams>) => {
            return !hasStructureVisibleBonds(structure, props)
                ? EmptyLocationIterator
                : BondIterator.fromStructureGroups(structure, { includeLocation2: props.colorMode === 'interpolate' });
        },
        getLoci: getStructureGroupsBondLoci,
        eachLocation: eachStructureGroupsBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureIntraUnitBondCylinderParams>, currentProps: PD.Values<StructureIntraUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
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
        mustRecreate: (structure: Structure, props: PD.Values<StructureIntraUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return !props.tryUseImpostor || !webgl;
        }
    }, materialId);
}

export function StructureIntraUnitBondCylinderMeshVisual(materialId: number): ComplexVisual<StructureIntraUnitBondCylinderParams> {
    return ComplexMeshVisual<StructureIntraUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(StructureIntraUnitBondCylinderParams),
        createGeometry: createStructureIntraUnitBondCylinderMesh,
        createLocationIterator: (structure: Structure, props: PD.Values<StructureIntraUnitBondCylinderParams>) => {
            return !hasStructureVisibleBonds(structure, props)
                ? EmptyLocationIterator
                : BondIterator.fromStructureGroups(structure);
        },
        getLoci: getStructureGroupsBondLoci,
        eachLocation: eachStructureGroupsBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureIntraUnitBondCylinderParams>, currentProps: PD.Values<StructureIntraUnitBondCylinderParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
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
        mustRecreate: (structure: Structure, props: PD.Values<StructureIntraUnitBondCylinderParams>, webgl?: WebGLContext) => {
            return props.tryUseImpostor && !!webgl;
        }
    }, materialId);
}
