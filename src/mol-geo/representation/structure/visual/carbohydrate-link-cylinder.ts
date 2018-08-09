/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Structure, Link, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, StructureVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createIdentityTransform, createColors } from './util/common';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { createMarkers, MarkerAction, MarkerData, applyMarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { Vec3 } from 'mol-math/linear-algebra';
import { deepEqual } from 'mol-util';
import { LocationIterator } from './util/location-iterator';
import { createLinkCylinderMesh, DefaultLinkCylinderProps, LinkCylinderProps } from './util/link';
import { OrderedSet } from 'mol-data/int';

// TODO create seperate visual
// for (let i = 0, il = carbohydrates.terminalLinks.length; i < il; ++i) {
//     const tl = carbohydrates.terminalLinks[i]
//     const center = carbohydrates.elements[tl.carbohydrateIndex].geometry.center
//     tl.elementUnit.conformation.position(tl.elementUnit.elements[tl.elementIndex], p)
//     if (tl.fromCarbohydrate) {
//         builder.addCylinder(center, p, 0.5, linkParams)
//     } else {
//         builder.addCylinder(p, center, 0.5, linkParams)
//     }
// }

async function createCarbohydrateLinkCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {
    const { links, elements } = structure.carbohydrates

    const builderProps = {
        linkCount: links.length,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const l = links[edgeIndex]
            Vec3.copy(posA, elements[l.carbohydrateIndexA].geometry.center)
            Vec3.copy(posB, elements[l.carbohydrateIndexB].geometry.center)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => 0
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultCarbohydrateLinkProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type CarbohydrateLinkProps = Partial<typeof DefaultCarbohydrateLinkProps>

export function CarbohydrateLinkVisual(): StructureVisual<CarbohydrateLinkProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultCarbohydrateLinkProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, structure: Structure, props: CarbohydrateLinkProps = {}) {
            currentProps = Object.assign({}, DefaultCarbohydrateLinkProps, props)
            currentStructure = structure

            const { colorTheme } = { ...DefaultCarbohydrateLinkProps, ...props }
            const instanceCount = 1
            const elementCount = currentStructure.elementCount

            mesh = await createCarbohydrateLinkCylinderMesh(ctx, currentStructure, currentProps, mesh)
            // console.log(mesh)

            const transforms = createIdentityTransform()
            const color = createColors(createCarbohydrateLinkIterator(structure), colorTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts)
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: CarbohydrateLinkProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                createColors(createCarbohydrateLinkIterator(currentStructure), newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return false
        },
        getLoci(pickingId: PickingId) {
            return getLinkLoci(pickingId, currentStructure, renderObject.id)
        },
        mark(loci: Loci, action: MarkerAction) {
            markLink(loci, action, currentStructure, renderObject.values)
        },
        destroy() {
            // TODO
        }
    }
}

function createCarbohydrateLinkIterator(structure: Structure): LocationIterator {
    const { elements, links } = structure.carbohydrates
    const elementCount = links.length
    const instanceCount = 1
    const location = Link.Location()
    const getLocation = (elementIndex: number, instanceIndex: number) => {
        const link = links[elementIndex]
        const carbA = elements[link.carbohydrateIndexA]
        const carbB = elements[link.carbohydrateIndexB]
        const indexA = OrderedSet.findPredecessorIndex(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.findPredecessorIndex(carbB.unit.elements, carbB.anomericCarbon)
        location.aUnit = carbA.unit
        location.aIndex = indexA as StructureElement.UnitIndex
        location.bUnit = carbB.unit
        location.bIndex = indexB as StructureElement.UnitIndex
        return location
    }
    return LocationIterator(elementCount, instanceCount, getLocation)
}

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const { links, elements } = structure.carbohydrates
        const l = links[elementId]
        const carbA = elements[l.carbohydrateIndexA]
        const carbB = elements[l.carbohydrateIndexB]
        const indexA = OrderedSet.findPredecessorIndex(carbA.unit.elements, carbA.anomericCarbon)
        const indexB = OrderedSet.findPredecessorIndex(carbB.unit.elements, carbB.anomericCarbon)
        return Link.Loci([
            Link.Location(
                carbA.unit, indexA as StructureElement.UnitIndex,
                carbB.unit, indexB as StructureElement.UnitIndex
            )
        ])
    }
    return EmptyLoci
}

function markLink(loci: Loci, action: MarkerAction, structure: Structure, values: MarkerData) {
    const tMarker = values.tMarker

    const { getLinkIndex } = structure.carbohydrates
    const elementCount = structure.carbohydrates.elements.length

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        if (applyMarkerAction(array, 0, elementCount, action)) {
            changed = true
        }
    } else if (Link.isLoci(loci)) {
        for (const l of loci.links) {
            const idx = getLinkIndex(l.aUnit, l.aUnit.elements[l.aIndex], l.bUnit, l.bUnit.elements[l.bIndex])
            if (idx !== undefined) {
                if (applyMarkerAction(array, idx, idx + 1, action) && !changed) {
                    changed = true
                }
            }
        }
    } else {
        return
    }
    if (changed) {
        ValueCell.update(tMarker, tMarker.ref.value)
    }
}