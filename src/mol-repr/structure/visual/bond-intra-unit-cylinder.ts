/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Unit, Structure, StructureElement, Bond } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { BitFlags, arrayEqual } from '../../../mol-util';
import { createLinkCylinderMesh, LinkCylinderStyle } from './util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../units-visual';
import { VisualUpdateState } from '../../util';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { isHydrogen } from './util/common';
import { BondType } from '../../../mol-model/structure/model/types';
import { ignoreBondType, BondCylinderParams, BondIterator } from './util/bond';
import { Sphere3D } from '../../../mol-math/geometry';

function createIntraUnitBondCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitBondParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const location = StructureElement.Location.create(structure, unit);

    const elements = unit.elements;
    const bonds = unit.bonds;
    const { edgeCount, a, b, edgeProps, offset } = bonds;
    const { order: _order, flags: _flags } = edgeProps;
    const { sizeFactor, sizeAspectRatio, ignoreHydrogens, includeTypes, excludeTypes } = props;

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);

    const ignoreHydrogen = ignoreHydrogens ? (edgeIndex: number) => {
        return isHydrogen(unit, elements[a[edgeIndex]]) || isHydrogen(unit, elements[b[edgeIndex]]);
    } : () => false;

    if (!edgeCount) return Mesh.createEmpty(mesh);

    const vRef = Vec3.zero();
    const pos = unit.conformation.invariantPosition;

    const builderProps = {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];

            if (aI > bI) [aI, bI] = [bI, aI];
            if (offset[aI + 1] - offset[aI] === 1) [aI, bI] = [bI, aI];
            // TODO prefer reference atoms in rings

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
            const f = BitFlags.create(_flags[edgeIndex]);
            if (BondType.is(f, BondType.Flag.MetallicCoordination) || BondType.is(f, BondType.Flag.HydrogenBond)) {
                // show metall coordinations and hydrogen bonds with dashed cylinders
                return LinkCylinderStyle.Dashed;
            } else if (o === 2) {
                return LinkCylinderStyle.Double;
            } else if (o === 3) {
                return LinkCylinderStyle.Triple;
            } else {
                return LinkCylinderStyle.Solid;
            }
        },
        radius: (edgeIndex: number) => {
            location.element = elements[a[edgeIndex]];
            const sizeA = theme.size.size(location);
            location.element = elements[b[edgeIndex]];
            const sizeB = theme.size.size(location);
            return Math.min(sizeA, sizeB) * sizeFactor * sizeAspectRatio;
        },
        ignore: (edgeIndex: number) => ignoreHydrogen(edgeIndex) || ignoreBondType(include, exclude, _flags[edgeIndex])
    };

    const m = createLinkCylinderMesh(ctx, builderProps, props, mesh);

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1 * sizeFactor);
    m.setBoundingSphere(sphere);

    return m;
}

export const IntraUnitBondParams = {
    ...UnitsMeshParams,
    ...BondCylinderParams,
    sizeFactor: PD.Numeric(0.3, { min: 0, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2 / 3, { min: 0, max: 3, step: 0.01 }),
    ignoreHydrogens: PD.Boolean(false),
};
export type IntraUnitBondParams = typeof IntraUnitBondParams

export function IntraUnitBondVisual(materialId: number): UnitsVisual<IntraUnitBondParams> {
    return UnitsMeshVisual<IntraUnitBondParams>({
        defaultProps: PD.getDefaultValues(IntraUnitBondParams),
        createGeometry: createIntraUnitBondCylinderMesh,
        createLocationIterator: BondIterator.fromGroup,
        getLoci: getBondLoci,
        eachLocation: eachBond,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitBondParams>, currentProps: PD.Values<IntraUnitBondParams>) => {
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

function getBondLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        if (Unit.isAtomic(unit)) {
            return Bond.Loci(structure, [
                Bond.Location(
                    structure, unit, unit.bonds.a[groupId] as StructureElement.UnitIndex,
                    structure, unit, unit.bonds.b[groupId] as StructureElement.UnitIndex
                ),
                Bond.Location(
                    structure, unit, unit.bonds.b[groupId] as StructureElement.UnitIndex,
                    structure, unit, unit.bonds.a[groupId] as StructureElement.UnitIndex
                )
            ]);
        }
    }
    return EmptyLoci;
}

function eachBond(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean, isMarking: boolean) {
    let changed = false;
    if (Bond.isLoci(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        const unit = group.units[0];
        if (!Unit.isAtomic(unit)) return false;
        const groupCount = unit.bonds.edgeCount * 2;
        for (const b of loci.bonds) {
            const unitIdx = group.unitIndexMap.get(b.aUnit.id);
            if (unitIdx !== undefined) {
                const idx = unit.bonds.getDirectedEdgeIndex(b.aIndex, b.bIndex);
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        const unit = group.units[0];
        if (!Unit.isAtomic(unit)) return false;
        const groupCount = unit.bonds.edgeCount * 2;
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id);
            if (unitIdx !== undefined) {
                const { offset, b } = unit.bonds;
                OrderedSet.forEach(e.indices, v => {
                    for (let t = offset[v], _t = offset[v + 1]; t < _t; t++) {
                        if (!isMarking || OrderedSet.has(e.indices, b[t])) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + t))) changed = true;
                        }
                    }
                });
            }
        }
    }
    return changed;
}