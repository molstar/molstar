/**
 * Copyright (c) 2020-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { arrayEqual } from '../../../mol-util';
import { LinkStyle, createLinkLines, LinkBuilderProps, EmptyLinkBuilderProps } from './util/link';
import { UnitsVisual, UnitsLinesParams, UnitsLinesVisual } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondIterator, BondLineParams, getIntraBondLoci, eachIntraBond, makeIntraBondIgnoreTest, ignoreBondType, hasUnitVisibleBonds, hasStructureVisibleBonds, getStructureGroupsBondLoci, eachStructureGroupsBond } from './util/bond';
import { Sphere3D } from '../../../mol-math/geometry';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { arrayIntersectionSize } from '../../../mol-util/array';
import { StructureGroup } from './util/common';
import { ComplexLinesParams, ComplexLinesVisual, ComplexVisual } from '../complex-visual';
import { EmptyLocationIterator } from '../../../mol-geo/util/location-iterator';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const isBondType = BondType.is;

function getIntraUnitBondLineBuilderProps(unit: Unit.Atomic, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondLineParams>): LinkBuilderProps {
    const location = StructureElement.Location.create(structure, unit);

    const elements = unit.elements;
    const bonds = unit.bonds;
    const { edgeCount, a, b, edgeProps, offset } = bonds;

    const { order: _order, flags: _flags } = edgeProps;
    const { sizeFactor, aromaticBonds, includeTypes, excludeTypes, multipleBonds } = props;

    const mbOff = multipleBonds === 'off';
    const mbSymmetric = multipleBonds === 'symmetric';

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);
    const ignoreComputedAromatic = ignoreBondType(include, exclude, BondType.Flag.Computed);

    const vRef = Vec3();
    const c = unit.conformation;

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
        position: (posA: Vec3, posB: Vec3, edgeIndex: number, _adjust: boolean) => {
            c.invariantPosition(elements[a[edgeIndex]], posA);
            c.invariantPosition(elements[b[edgeIndex]], posB);
        },
        style: (edgeIndex: number) => {
            const o = _order[edgeIndex];
            const f = _flags[edgeIndex];
            if (isBondType(f, BondType.Flag.MetallicCoordination) || isBondType(f, BondType.Flag.HydrogenBond)) {
                // show metallic coordinations and hydrogen bonds with dashed lines
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
            location.element = elements[a[edgeIndex]];
            const sizeA = theme.size.size(location);
            location.element = elements[b[edgeIndex]];
            const sizeB = theme.size.size(location);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: makeIntraBondIgnoreTest(structure, unit, props)
    };
}

function createIntraUnitBondLines(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondLineParams>, lines?: Lines) {
    if (!Unit.isAtomic(unit)) return Lines.createEmpty(lines);
    if (!hasUnitVisibleBonds(unit, props)) return Lines.createEmpty(lines);
    if (!unit.bonds.edgeCount) return Lines.createEmpty(lines);

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Lines.createEmpty(lines);

    const builderProps = getIntraUnitBondLineBuilderProps(unit, structure, theme, props);
    const { lines: l, boundingSphere } = createLinkLines(ctx, builderProps, props, lines);

    if (boundingSphere) {
        l.setBoundingSphere(boundingSphere);
    } else if (l.lineCount > 0) {
        const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * props.sizeFactor);
        l.setBoundingSphere(sphere);
    }

    return l;
}

export const IntraUnitBondLineParams = {
    ...UnitsLinesParams,
    ...BondLineParams,
    includeParent: PD.Boolean(false),
};
export type IntraUnitBondLineParams = typeof IntraUnitBondLineParams

export function IntraUnitBondLineVisual(materialId: number): UnitsVisual<IntraUnitBondLineParams> {
    return UnitsLinesVisual<IntraUnitBondLineParams>({
        defaultProps: PD.getDefaultValues(IntraUnitBondLineParams),
        createGeometry: createIntraUnitBondLines,
        createLocationIterator: (structureGroup: StructureGroup) => BondIterator.fromGroup(structureGroup),
        getLoci: getIntraBondLoci,
        eachLocation: eachIntraBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitBondLineParams>, currentProps: PD.Values<IntraUnitBondLineParams>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.aromaticDashCount !== currentProps.aromaticDashCount ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes) ||
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
        }
    }, materialId);
}

//

function getStructureIntraUnitBondLineBuilderProps(structure: Structure, theme: Theme, props: PD.Values<StructureIntraUnitBondLineParams>): LinkBuilderProps {
    const intraUnitProps: { group: Unit.SymmetryGroup, props: LinkBuilderProps}[] = [];

    const { bondCount, unitIndex, unitEdgeIndex, unitGroupIndex } = structure.intraUnitBondMapping;
    const { child } = structure;

    for (const ug of structure.unitSymmetryGroups) {
        const unit = ug.units[0];
        const childUnit = child?.unitMap.get(unit.id);
        const p = Unit.isAtomic(unit) && !(child && !childUnit)
            ? getIntraUnitBondLineBuilderProps(unit, structure, theme, props)
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

function createStructureIntraUnitBondLines(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<StructureIntraUnitBondLineParams>, lines?: Lines) {
    if (!hasStructureVisibleBonds(structure, props)) return Lines.createEmpty(lines);
    if (!structure.intraUnitBondMapping.bondCount) return Lines.createEmpty(lines);

    const builderProps = getStructureIntraUnitBondLineBuilderProps(structure, theme, props);
    const { lines: l, boundingSphere } = createLinkLines(ctx, builderProps, props, lines);

    if (boundingSphere) {
        l.setBoundingSphere(boundingSphere);
    } else if (l.lineCount > 0) {
        const { child } = structure;
        const sphere = Sphere3D.expand(Sphere3D(), (child ?? structure).boundary.sphere, 1 * props.sizeFactor);
        l.setBoundingSphere(sphere);
    }

    return l;
}

export const StructureIntraUnitBondLineParams = {
    ...ComplexLinesParams,
    ...BondLineParams,
    includeParent: PD.Boolean(false),
};
export type StructureIntraUnitBondLineParams = typeof StructureIntraUnitBondLineParams

export function StructureIntraUnitBondLineVisual(materialId: number): ComplexVisual<StructureIntraUnitBondLineParams> {
    return ComplexLinesVisual<StructureIntraUnitBondLineParams>({
        defaultProps: PD.getDefaultValues(StructureIntraUnitBondLineParams),
        createGeometry: createStructureIntraUnitBondLines,
        createLocationIterator: (structure: Structure, props: PD.Values<StructureIntraUnitBondLineParams>) => {
            return !hasStructureVisibleBonds(structure, props)
                ? EmptyLocationIterator
                : BondIterator.fromStructureGroups(structure);
        },
        getLoci: getStructureGroupsBondLoci,
        eachLocation: eachStructureGroupsBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<StructureIntraUnitBondLineParams>, currentProps: PD.Values<StructureIntraUnitBondLineParams>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.aromaticDashCount !== currentProps.aromaticDashCount ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes) ||
                newProps.multipleBonds !== currentProps.multipleBonds
            );

            if (hasStructureVisibleBonds(newStructure, newProps) && newStructure.interUnitBonds !== currentStructure.interUnitBonds) {
                state.createGeometry = true;
                state.updateTransform = true;
                state.updateColor = true;
                state.updateSize = true;
            }
        }
    }, materialId);
}
