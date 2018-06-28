/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Link, Structure } from 'mol-model/structure';
import { DefaultStructureProps, StructureVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { LinkCylinderProps, DefaultLinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../../util/color-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { createIdentityTransform } from './util/common';
import { updateMeshValues, updateRenderableState, createMeshValues, createRenderableState } from '../../util';
// import { chainIdLinkColorData } from '../../../theme/structure/color/chain-id';

async function createInterUnitLinkCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links

    if (!bondCount) return Mesh.createEmpty(mesh)

    const builderProps = {
        linkCount: bondCount,
        referencePosition: (edgeIndex: number) => null, // TODO
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = bonds[edgeIndex]
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
        },
        order: (edgeIndex: number) => bonds[edgeIndex].order,
        flags: (edgeIndex: number) => bonds[edgeIndex].flag
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultInterUnitLinkProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
}
export type InterUnitLinkProps = Partial<typeof DefaultInterUnitLinkProps>

export function InterUnitLinkVisual(): StructureVisual<InterUnitLinkProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultInterUnitLinkProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, structure: Structure, props: InterUnitLinkProps = {}) {
            currentProps = Object.assign({}, DefaultInterUnitLinkProps, props)
            currentStructure = structure

            const elementCount = structure.links.bondCount
            const instanceCount = 1

            mesh = await createInterUnitLinkCylinderMesh(ctx, structure, currentProps)

            const transforms = createIdentityTransform()
            const color = createUniformColor({ value: 0x999911 }) // TODO
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts),
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: InterUnitLinkProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            // TODO create in-place
            if (currentProps.radialSegments !== newProps.radialSegments) return false

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

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

function getLinkLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, elementId } = pickingId
    if (id === objectId) {
        const bond = structure.links.bonds[elementId]
        return Link.Loci([{
            aUnit: bond.unitA,
            aIndex: bond.indexA,
            bUnit: bond.unitB,
            bIndex: bond.indexB
        }])
    }
    return EmptyLoci
}

function markLink(loci: Loci, action: MarkerAction, structure: Structure, values: MarkerData) {
    const tMarker = values.tMarker

    const links = structure.links
    const elementCount = links.bondCount
    const instanceCount = 1

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        applyMarkerAction(array, 0, elementCount * instanceCount, action)
        changed = true
    } else if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const _idx = structure.links.getBondIndex(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (_idx !== -1) {
                const idx = _idx
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