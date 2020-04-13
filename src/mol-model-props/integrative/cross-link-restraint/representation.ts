/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Theme } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval } from '../../../mol-data/int';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { VisualContext } from '../../../mol-repr/visual';
import { createLinkCylinderMesh, LinkCylinderParams } from '../../../mol-repr/structure/visual/util/link';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider } from '../../../mol-repr/structure/representation';
import { CustomProperty } from '../../common/custom-property';
import { CrossLinkRestraintProvider, CrossLinkRestraint } from './property';

function createCrossLinkRestraintCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<CrossLinkRestraintCylinderParams>, mesh?: Mesh) {

    const crossLinks = CrossLinkRestraintProvider.get(structure).value!;
    if (!crossLinks.count) return Mesh.createEmpty(mesh);
    const { sizeFactor } = props;

    const location = StructureElement.Location.create(structure);

    const builderProps = {
        linkCount: crossLinks.count,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex];
            const uA = b.unitA, uB = b.unitB;
            uA.conformation.position(uA.elements[b.indexA], posA);
            uB.conformation.position(uB.elements[b.indexB], posB);
        },
        radius: (edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex];
            location.unit = b.unitA;
            location.element = b.unitA.elements[b.indexA];
            return theme.size.size(location) * sizeFactor;
        },
    };

    return createLinkCylinderMesh(ctx, builderProps, props, mesh);
}

export const CrossLinkRestraintCylinderParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    sizeFactor: PD.Numeric(0.5, { min: 0, max: 10, step: 0.1 }),
};
export type CrossLinkRestraintCylinderParams = typeof CrossLinkRestraintCylinderParams

export function CrossLinkRestraintVisual(materialId: number): ComplexVisual<CrossLinkRestraintCylinderParams> {
    return ComplexMeshVisual<CrossLinkRestraintCylinderParams>({
        defaultProps: PD.getDefaultValues(CrossLinkRestraintCylinderParams),
        createGeometry: createCrossLinkRestraintCylinderMesh,
        createLocationIterator: createCrossLinkRestraintIterator,
        getLoci: getLinkLoci,
        eachLocation: eachCrossLink,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<CrossLinkRestraintCylinderParams>, currentProps: PD.Values<CrossLinkRestraintCylinderParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkCap !== currentProps.linkCap
            );
        }
    }, materialId);
}

function createCrossLinkRestraintIterator(structure: Structure): LocationIterator {
    const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
    const { pairs } = crossLinkRestraints;
    const groupCount = pairs.length;
    const instanceCount = 1;
    const location = CrossLinkRestraint.Location(crossLinkRestraints, structure);
    const getLocation = (groupIndex: number) => {
        location.element = groupIndex;
        return location;
    };
    return LocationIterator(groupCount, instanceCount, getLocation, true);
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
        const pair = crossLinkRestraints.pairs[groupId];
        if (pair) {
            return CrossLinkRestraint.Loci(structure, crossLinkRestraints, [groupId]);
        }
    }
    return EmptyLoci;
}

function eachCrossLink(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false;
    if (CrossLinkRestraint.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.data.structure, structure)) return false;
        const crossLinkRestraints = CrossLinkRestraintProvider.get(structure).value!;
        if (loci.data.crossLinkRestraints !== crossLinkRestraints) return false;

        for (const e of loci.elements) {
            if (apply(Interval.ofSingleton(e))) changed = true;
        }
    }
    return changed;
}

//

const CrossLinkRestraintVisuals = {
    'cross-link-restraint': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CrossLinkRestraintCylinderParams>) => ComplexRepresentation('Cross-link restraint', ctx, getParams, CrossLinkRestraintVisual),
};

export const CrossLinkRestraintParams = {
    ...CrossLinkRestraintCylinderParams,
};
export type CrossLinkRestraintParams = typeof CrossLinkRestraintParams
export function getCrossLinkRestraintParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(CrossLinkRestraintParams);
}

export type CrossLinkRestraintRepresentation = StructureRepresentation<CrossLinkRestraintParams>
export function CrossLinkRestraintRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CrossLinkRestraintParams>): CrossLinkRestraintRepresentation {
    return Representation.createMulti('CrossLinkRestraint', ctx, getParams, StructureRepresentationStateBuilder, CrossLinkRestraintVisuals as unknown as Representation.Def<Structure, CrossLinkRestraintParams>);
}

export const CrossLinkRestraintRepresentationProvider = StructureRepresentationProvider({
    name: CrossLinkRestraint.Tag.CrossLinkRestraint,
    label: 'Cross Link Restraint',
    description: 'Displays cross-link restraints.',
    factory: CrossLinkRestraintRepresentation,
    getParams: getCrossLinkRestraintParams,
    defaultValues: PD.getDefaultValues(CrossLinkRestraintParams),
    defaultColorTheme: { name: CrossLinkRestraint.Tag.CrossLinkRestraint },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => CrossLinkRestraint.isApplicable(structure),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, structure: Structure) => CrossLinkRestraintProvider.attach(ctx, structure, void 0, true),
        detach: (data) => CrossLinkRestraintProvider.ref(data, false)
    }
});