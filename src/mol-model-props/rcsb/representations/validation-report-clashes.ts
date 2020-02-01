/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Unit, Structure } from '../../../mol-model/structure';
import { Theme, ThemeRegistryContext, ThemeDataContext } from '../../../mol-theme/theme';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Interval, SortedArray } from '../../../mol-data/int';
import { RepresentationContext, RepresentationParamsGetter, Representation } from '../../../mol-repr/representation';
import { UnitsRepresentation, StructureRepresentation, StructureRepresentationStateBuilder, StructureRepresentationProvider, ComplexRepresentation } from '../../../mol-repr/structure/representation';
import { UnitKind, UnitKindOptions } from '../../../mol-repr/structure/visual/util/common';
import { VisualContext } from '../../../mol-repr/visual';
import { createLinkCylinderMesh, LinkCylinderParams, LinkCylinderStyle } from '../../../mol-repr/structure/visual/util/link';
import { UnitsMeshParams, UnitsVisual, UnitsMeshVisual, StructureGroup } from '../../../mol-repr/structure/units-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { DataLocation } from '../../../mol-model/location';
import { ValidationReportProvider, ValidationReport } from '../validation-report';
import { CustomProperty } from '../../common/custom-property';
import { ComplexMeshParams, ComplexVisual, ComplexMeshVisual } from '../../../mol-repr/structure/complex-visual';
import { arrayMax } from '../../../mol-util/array';
import { UnitIndex } from '../../../mol-model/structure/structure/element/element';
import { InterUnitGraph } from '../../../mol-math/graph/inter-unit-graph';
import { ColorTheme } from '../../../mol-theme/color';
import { ColorNames } from '../../../mol-util/color/names';

//

function createIntraUnitClashCylinderMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<IntraUnitClashParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const validationReport = ValidationReportProvider.get(unit.model).value!
    const { clashes } = validationReport

    const { edgeCount, a, b, edgeProps } = clashes
    const { magnitude } = edgeProps
    const { sizeFactor } = props

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const pos = unit.conformation.invariantPosition

    const builderProps = {
        linkCount: edgeCount * 2,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            pos(a[edgeIndex], posA)
            pos(b[edgeIndex], posB)
        },
        style: (edgeIndex: number) => LinkCylinderStyle.Disk,
        radius: (edgeIndex: number) => magnitude[edgeIndex] * sizeFactor,
        ignore: (edgeIndex: number) => {
            return (
                // TODO create lookup
                !SortedArray.has(unit.elements, a[edgeIndex]) ||
                !SortedArray.has(unit.elements, b[edgeIndex])
            )
        }
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
            structure
            groupId
            // TODO
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
    const { group } = structureGroup
    const unit = group.units[0]
    const validationReport = ValidationReportProvider.get(unit.model).value!
    const { clashes } = validationReport
    const groupCount = clashes.edgeCount * 2
    const instanceCount = group.units.length
    const location = DataLocation(validationReport, 'clashes')
    const getLocation = (groupIndex: number, instanceIndex: number) => {
        location.index = groupIndex + instanceIndex * groupCount
        return location
    }
    return LocationIterator(groupCount, instanceCount, getLocation)
}

//

type InterUnitClashesProps = { id: number, magnitude: number, distance: number }

function createInterUnitClashes(structure: Structure, clashes: ValidationReport['clashes']) {
    const builder = new InterUnitGraph.Builder<Unit.Atomic, UnitIndex, InterUnitClashesProps>()
    const { a, b, edgeProps: { id, magnitude, distance } } = clashes

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        const elementsA = unitA.elements
        const elementsB = unitB.elements

        builder.startUnitPair(unitA as Unit.Atomic, unitB as Unit.Atomic)

        for (let i = 0, il = clashes.edgeCount * 2; i < il; ++i) {
            // TODO create lookup
            let indexA = SortedArray.indexOf(elementsA, a[i])
            let indexB = SortedArray.indexOf(elementsB, b[i])

            if (indexA !== -1 && indexB !== -1) {
                builder.add(indexA as UnitIndex, indexB as UnitIndex, {
                    id: id[i],
                    magnitude: magnitude[i],
                    distance: distance[i]
                })
            }
        }

        builder.finishUnitPair()
    }, {
        maxRadius: arrayMax(clashes.edgeProps.distance),
        validUnit: (unit: Unit) => Unit.isAtomic(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => unitA.model === unitB.model
    })

    return new InterUnitGraph(builder.getMap())
}

function createInterUnitClashCylinderMesh(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<InterUnitClashParams>, mesh?: Mesh) {
    const validationReport = ValidationReportProvider.get(structure.models[0]).value!
    const clashes = createInterUnitClashes(structure, validationReport.clashes)

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
    const location = DataLocation({}, 'clashes')
    return LocationIterator(1, 1, () => location)
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
    defaultColorTheme: ValidationReport.Tag.Clashes,
    defaultSizeTheme: 'physical',
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    ensureCustomProperties: (ctx: CustomProperty.Context, structure: Structure) => {
        return ValidationReportProvider.attach(ctx, structure.models[0])
    }
}

//

function ClashesColorTheme(ctx: ThemeDataContext, props: {}): ColorTheme<{}> {
    return {
        factory: ClashesColorTheme,
        granularity: 'uniform',
        color: () => ColorNames.hotpink,
        props,
        description: 'Uniform color for clashes',
    }
}

export const ClashesColorThemeProvider: ColorTheme.Provider<{}> = {
    label: 'RCSB Clashes',
    factory: ClashesColorTheme,
    getParams: () => ({}),
    defaultValues: PD.getDefaultValues({}),
    isApplicable: (ctx: ThemeDataContext) => false,
}