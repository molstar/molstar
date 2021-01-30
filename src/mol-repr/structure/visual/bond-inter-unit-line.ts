/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Bond, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags, arrayEqual } from '../../../mol-util';
import { LinkStyle, createLinkLines } from './util/link';
import { ComplexVisual, ComplexLinesVisual, ComplexLinesParams } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondIterator, getInterBondLoci, eachInterBond, BondLineParams, makeInterBondIgnoreTest } from './util/bond';
import { Lines } from '../../../mol-geo/geometry/lines/lines';
import { Sphere3D } from '../../../mol-math/geometry';

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
const tmpLoc = StructureElement.Location.create(void 0);

function createInterUnitBondLines(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitBondLineParams>, lines?: Lines) {
    const bonds = structure.interUnitBonds;
    const { edgeCount, edges } = bonds;
    const { sizeFactor } = props;

    if (!edgeCount) return Lines.createEmpty(lines);

    const builderProps = {
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
                throw new Error('same units in createInterUnitBondLines');
            }
            return setRefPosition(tmpRef, structure, unitA, indexA) || setRefPosition(tmpRef, structure, unitB, indexB);
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = edges[edgeIndex];
            const uA = structure.unitMap.get(b.unitA);
            const uB = structure.unitMap.get(b.unitB);
            uA.conformation.position(uA.elements[b.indexA], posA);
            uB.conformation.position(uB.elements[b.indexB], posB);
        },
        style: (edgeIndex: number) => {
            const o = edges[edgeIndex].props.order;
            const f = BitFlags.create(edges[edgeIndex].props.flag);
            if (BondType.is(f, BondType.Flag.MetallicCoordination) || BondType.is(f, BondType.Flag.HydrogenBond)) {
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
            const b = edges[edgeIndex];
            tmpLoc.structure = structure;
            tmpLoc.unit = structure.unitMap.get(b.unitA);
            tmpLoc.element = tmpLoc.unit.elements[b.indexA];
            const sizeA = theme.size.size(tmpLoc);
            tmpLoc.unit = structure.unitMap.get(b.unitB);
            tmpLoc.element = tmpLoc.unit.elements[b.indexB];
            const sizeB = theme.size.size(tmpLoc);
            return Math.min(sizeA, sizeB) * sizeFactor;
        },
        ignore: makeInterBondIgnoreTest(structure, props)
    };

    const l = createLinkLines(ctx, builderProps, props, lines);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, 1 * sizeFactor);
    l.setBoundingSphere(sphere);

    return l;
}

export const InterUnitBondLineParams = {
    ...ComplexLinesParams,
    ...BondLineParams,
};
export type InterUnitBondLineParams = typeof InterUnitBondLineParams

export function InterUnitBondLineVisual(materialId: number): ComplexVisual<InterUnitBondLineParams> {
    return ComplexLinesVisual<InterUnitBondLineParams>({
        defaultProps: PD.getDefaultValues(InterUnitBondLineParams),
        createGeometry: createInterUnitBondLines,
        createLocationIterator: BondIterator.fromStructure,
        getLoci: getInterBondLoci,
        eachLocation: eachInterBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitBondLineParams>, currentProps: PD.Values<InterUnitBondLineParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.dashCount !== currentProps.dashCount ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes)
            );
        }
    }, materialId);
}
