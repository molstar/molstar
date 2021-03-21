/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { arrayEqual } from '../../../mol-util';
import { LinkStyle, createLinkLines } from './util/link';
import { UnitsVisual, UnitsLinesParams, UnitsLinesVisual, StructureGroup } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondIterator, BondLineParams, getIntraBondLoci, eachIntraBond, makeIntraBondIgnoreTest } from './util/bond';
import { Sphere3D } from '../../../mol-math/geometry';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { IntAdjacencyGraph } from '../../../mol-math/graph';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const isBondType = BondType.is;

function createIntraUnitBondLines(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondLineParams>, lines?: Lines) {
    if (!Unit.isAtomic(unit)) return Lines.createEmpty(lines);

    const child = Structure.WithChild.getChild(structure);
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) return Lines.createEmpty(lines);

    if (props.includeParent) {
        const child = Structure.WithChild.getChild(structure);
        if (!child) throw new Error('expected child to exist');
    }

    const location = StructureElement.Location.create(structure, unit);

    const elements = unit.elements;
    const bonds = unit.bonds;
    const { edgeCount, a, b, edgeProps, offset } = bonds;
    const { order: _order, flags: _flags } = edgeProps;
    const { sizeFactor } = props;

    if (!edgeCount) return Lines.createEmpty(lines);

    const vRef = Vec3();
    const pos = unit.conformation.invariantPosition;

    const builderProps = {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];

            if (aI > bI) [aI, bI] = [bI, aI];
            if (offset[aI + 1] - offset[aI] === 1) [aI, bI] = [bI, aI];
            // TODO prefer reference atoms within rings

            for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
                const _bI = b[i];
                if (_bI !== bI && _bI !== aI) return pos(elements[_bI], vRef);
            }
            for (let i = offset[bI], il = offset[bI + 1]; i < il; ++i) {
                const _aI = a[i];
                if (_aI !== aI && _aI !== bI) return pos(elements[_aI], vRef);
            }
            return null;
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            pos(elements[a[edgeIndex]], posA);
            pos(elements[b[edgeIndex]], posB);
        },
        style: (edgeIndex: number) => {
            const o = _order[edgeIndex];
            const f = _flags[edgeIndex];
            if (isBondType(f, BondType.Flag.MetallicCoordination) || isBondType(f, BondType.Flag.HydrogenBond)) {
                // show metall coordinations and hydrogen bonds with dashed cylinders
                return LinkStyle.Dashed;
            } else if (o === 2) {
                return LinkStyle.Double;
            } else if (o === 3) {
                return LinkStyle.Triple;
            } else {
                return LinkStyle.Solid;
            }
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

    const l = createLinkLines(ctx, builderProps, props, lines);

    const sphere = Sphere3D.expand(Sphere3D(), (childUnit ?? unit).boundary.sphere, 1 * sizeFactor);
    l.setBoundingSphere(sphere);

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
        createLocationIterator: BondIterator.fromGroup,
        getLoci: getIntraBondLoci,
        eachLocation: eachIntraBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitBondLineParams>, currentProps: PD.Values<IntraUnitBondLineParams>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes)
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