/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, StructureElement, Bond, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags, arrayEqual } from '../../../mol-util';
import { createLinkCylinderMesh, LinkStyle } from './util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { BondType } from '../../../mol-model/structure/model/types';
import { BondCylinderParams, BondIterator, getInterBondLoci, eachInterBond, makeInterBondIgnoreTest } from './util/bond';

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

function createInterUnitBondCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitBondCylinderParams>, mesh?: Mesh) {
    const bonds = structure.interUnitBonds;
    const { edgeCount, edges } = bonds;
    const { sizeFactor, sizeAspectRatio } = props;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const builderProps = {
        linkCount: edgeCount,
        referencePosition: (edgeIndex: number) => {
            const b = edges[edgeIndex];
            let unitA: Unit, unitB: Unit;
            let indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex;
            if (b.unitA.id < b.unitB.id) {
                unitA = b.unitA, unitB = b.unitB;
                indexA = b.indexA, indexB = b.indexB;
            } else if (b.unitA.id > b.unitB.id) {
                unitA = b.unitB, unitB = b.unitA;
                indexA = b.indexB, indexB = b.indexA;
            } else {
                throw new Error('same units in createInterUnitBondCylinderMesh');
            }
            return setRefPosition(tmpRef, structure, unitA, indexA) || setRefPosition(tmpRef, structure, unitB, indexB);
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = edges[edgeIndex];
            const uA = b.unitA, uB = b.unitB;
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
            tmpLoc.unit = b.unitA;
            tmpLoc.element = b.unitA.elements[b.indexA];
            const sizeA = theme.size.size(tmpLoc);
            tmpLoc.unit = b.unitB;
            tmpLoc.element = b.unitB.elements[b.indexB];
            const sizeB = theme.size.size(tmpLoc);
            return Math.min(sizeA, sizeB) * sizeFactor * sizeAspectRatio;
        },
        ignore: makeInterBondIgnoreTest(structure, props)
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const InterUnitBondCylinderParams = {
    ...ComplexMeshParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0, max: 3, step: 0.01 }),
};
export type InterUnitBondCylinderParams = typeof InterUnitBondCylinderParams

export function InterUnitBondCylinderVisual(materialId: number): ComplexVisual<InterUnitBondCylinderParams> {
    return ComplexMeshVisual<InterUnitBondCylinderParams>({
        defaultProps: PD.getDefaultValues(InterUnitBondCylinderParams),
        createGeometry: createInterUnitBondCylinderMesh,
        createLocationIterator: BondIterator.fromStructure,
        getLoci: getInterBondLoci,
        eachLocation: eachInterBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitBondCylinderParams>, currentProps: PD.Values<InterUnitBondCylinderParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.sizeAspectRatio !== currentProps.sizeAspectRatio ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
                newProps.linkCap !== currentProps.linkCap ||
                !arrayEqual(newProps.includeTypes, currentProps.includeTypes) ||
                !arrayEqual(newProps.excludeTypes, currentProps.excludeTypes)
            );
        }
    }, materialId);
}
