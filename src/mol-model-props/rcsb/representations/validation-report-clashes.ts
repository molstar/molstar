/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Unit, Structure, StructureElement } from '../../../mol-model/structure';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci, createDataLoci } from '../../../mol-model/loci';
import { Interval, OrderedSet } from '../../../mol-data/int';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../../../mol-repr/representation';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider, ComplexRepresentation } from '../../../mol-repr/structure/representation';
import { UnitKind, UnitKindOptions } from '../../../mol-repr/structure/visual/util/common';
import { VisualContext } from '../../../mol-repr/visual';
import { createLinkCylinderMesh, LinkCylinderParams, LinkCylinderStyle } from '../../../mol-repr/structure/visual/util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../../../mol-repr/structure/units-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { ClashesProvider } from '../validation-report';
import { CustomProperty } from '../../common/custom-property';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { Color } from '../../../mol-util/color';

//

function createIntraUnitClashCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitClashParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const clashes = ClashesProvider.get(structure).value!.intraUnit.get(unit.id)
    const { edgeCount, a, b, edgeProps } = clashes
    const { magnitude } = edgeProps
    const { sizeFactor } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const { elements } = unit
    const pos = unit.conformation.invariantPosition

    const builderProps = {
        linkCount: edgeCount * 2,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            pos(elements[a[edgeIndex]], posA)
            pos(elements[b[edgeIndex]], posB)
        },
        style: (edgeIndex: number) => LinkCylinderStyle.Disk,
        radius: (edgeIndex: number) => magnitude[edgeIndex] * sizeFactor,
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const IntraUnitClashParams = {
    ...UnitsMeshParams,
    ...LinkCylinderParams,
    linkCap: PD.Boolean(true),
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.01 }),
}
export type IntraUnitClashParams = typeof IntraUnitClashParams

export function IntraUnitClashVisual(materialId: number): UnitsVisual<IntraUnitClashParams> {
    return UnitsMeshVisual<IntraUnitClashParams>({
        defaultProps: PD.getDefaultValues(IntraUnitClashParams),
        createGeometry: createIntraUnitClashCylinderMesh,
        createLocationIterator: createIntraClashIterator,
        getLoci: getIntraClashLoci,
        eachLocation: eachIntraClash,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<IntraUnitClashParams>, currentProps: PD.Values<IntraUnitClashParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.linkCap !== currentProps.linkCap
            )
        }
    }, materialId)
}

function getIntraClashLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId
    if (id === objectId) {
        const { structure, group } = structureGroup
        const unit = group.units[instanceId]
        if (Unit.isAtomic(unit)) {
            const clashes = ClashesProvider.get(structure).value!.intraUnit.get(unit.id)
            return createDataLoci(clashes, 'clashes', OrderedSet.ofSingleton(groupId))
        }
    }
    return EmptyLoci
}

function eachIntraClash(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    let changed = false
    // TODO
    return changed
}

function createIntraClashIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup
    const unit = group.units[0]
    const clashes = ClashesProvider.get(structure).value!.intraUnit.get(unit.id)
    const { a } = clashes
    const groupCount = clashes.edgeCount * 2
    const instanceCount = group.units.length
    const location = StructureElement.Location.create()
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        const unit = group.units[instanceIndex]
        location.unit = unit
        location.element = unit.elements[a[groupIndex]]
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation)
}

//

function createInterUnitClashCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitClashParams>, mesh?: Mesh) {
    const clashes = ClashesProvider.get(structure).value!.interUnit
    const { edges, edgeCount } = clashes
    const { sizeFactor } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const builderProps = {
        linkCount: edgeCount,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = edges[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        style: (edgeIndex: number) => LinkCylinderStyle.Disk,
        radius: (edgeIndex: number) => edges[edgeIndex].props.magnitude * sizeFactor
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const InterUnitClashParams = {
    ...ComplexMeshParams,
    ...LinkCylinderParams,
    linkCap: PD.Boolean(true),
    sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.01 }),
}
export type InterUnitClashParams = typeof InterUnitClashParams

export function InterUnitClashVisual(materialId: number): ComplexVisual<InterUnitClashParams> {
    return ComplexMeshVisual<InterUnitClashParams>({
        defaultProps: PD.getDefaultValues(InterUnitClashParams),
        createGeometry: createInterUnitClashCylinderMesh,
        createLocationIterator: createInterClashIterator,
        getLoci: getInterClashLoci,
        eachLocation: eachInterClash,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<InterUnitClashParams>, currentProps: PD.Values<InterUnitClashParams>) => {
            state.createGeometry = (
                newProps.sizeFactor !== currentProps.sizeFactor ||
                newProps.radialSegments !== currentProps.radialSegments ||
                newProps.linkScale !== currentProps.linkScale ||
                newProps.linkSpacing !== currentProps.linkSpacing ||
                newProps.linkCap !== currentProps.linkCap
            )
        }
    }, materialId)
}

function getInterClashLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId
    if (id === objectId) {
        structure
        groupId
        // TODO
    }
    return EmptyLoci
}

function eachInterClash(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
    let changed = false
    // TODO
    return changed
}

function createInterClashIterator(structure: Structure): LocationIterator {
    const clashes = ClashesProvider.get(structure).value!.interUnit
    const groupCount = clashes.edgeCount
    const instanceCount = 1
    const location = StructureElement.Location.create()
    const getLocation = (groupIndex: number) => {
        const clash = clashes.edges[groupIndex]
        location.unit = clash.unitA
        location.element = clash.unitA.elements[clash.indexA]
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation, true)
}

//

const ClashesVisuals = {
    'intra-clash': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, IntraUnitClashParams>) => UnitsRepresentation('Intra-unit clash cylinder', ctx, getParams, IntraUnitClashVisual),
    'inter-clash': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, InterUnitClashParams>) => ComplexRepresentation('Inter-unit clash cylinder', ctx, getParams, InterUnitClashVisual),
}

export const ClashesParams = {
    ...IntraUnitClashParams,
    ...InterUnitClashParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic'], UnitKindOptions),
    visuals: PD.MultiSelect(['intra-clash', 'inter-clash'], PD.objectToOptions(ClashesVisuals))
}
export type ClashesParams = typeof ClashesParams
export function getClashesParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(ClashesParams)
}

export type ClashesRepresentation = StructureRepresentation<ClashesParams>
export function ClashesRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ClashesParams>): ClashesRepresentation {
    return Representation.createMulti('Clashes', ctx, getParams, StructureRepresentationStateBuilder, ClashesVisuals as unknown as Representation.Def<Structure, ClashesParams>)
}

export const ClashesRepresentationProvider: StructureRepresentationProvider<ClashesParams> = {
    label: 'RCSB Clashes',
    description: 'Displays clashes between atoms as disks.',
    factory: ClashesRepresentation,
    getParams: getClashesParams,
    defaultValues: PD.getDefaultValues(ClashesParams),
    defaultColorTheme: { name: 'uniform', props: { value: Color(0xFA28FF) } },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    ensureCustomProperties: (ctx: CustomProperty.Context, structure: Structure) => {
        return ClashesProvider.attach(ctx, structure)
    }
}