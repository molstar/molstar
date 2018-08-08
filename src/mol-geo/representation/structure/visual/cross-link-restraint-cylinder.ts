/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Link, Structure, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, StructureVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { LinkCylinderProps, DefaultLinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { createValueColor } from '../../../util/color-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { createIdentityTransform } from './util/common';
import { updateMeshValues, updateRenderableState, createMeshValues, createRenderableState } from '../../util';
// import { chainIdLinkColorData } from '../../../theme/structure/color/chain-id';

async function createCrossLinkRestraintCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {

    const crossLinks = structure.crossLinkRestraints
    if (!crossLinks.count) return Mesh.createEmpty(mesh)

    const builderProps = {
        linkCount: crossLinks.count,
        referencePosition: (edgeIndex: number) => null,
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            const b = crossLinks.pairs[edgeIndex]
            // console.log(b)
            const uA = b.unitA, uB = b.unitB
            uA.conformation.position(uA.elements[b.indexA], posA)
            uB.conformation.position(uB.elements[b.indexB], posB)
            // console.log(posA, posB)
        },
        order: (edgeIndex: number) => 1,
        flags: (edgeIndex: number) => 0
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultCrossLinkRestraintProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
    flipSided: false,
    flatShaded: false,
}
export type CrossLinkRestraintProps = Partial<typeof DefaultCrossLinkRestraintProps>

export function CrossLinkRestraintVisual(): StructureVisual<CrossLinkRestraintProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultCrossLinkRestraintProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, structure: Structure, props: CrossLinkRestraintProps = {}) {
            currentProps = Object.assign({}, DefaultCrossLinkRestraintProps, props)
            currentStructure = structure

            const elementCount = structure.crossLinkRestraints.count
            const instanceCount = 1

            mesh = await createCrossLinkRestraintCylinderMesh(ctx, structure, currentProps)

            const transforms = createIdentityTransform()
            const color = createValueColor(0x119911) // TODO
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
        async update(ctx: RuntimeContext, props: CrossLinkRestraintProps) {
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
        const pair = structure.crossLinkRestraints.pairs[elementId]
        if (pair) {
            return Link.Loci([
                Link.Location(
                    pair.unitA, pair.indexA as StructureElement.UnitIndex,
                    pair.unitB, pair.indexB as StructureElement.UnitIndex
                )
            ])
        }
    }
    return EmptyLoci
}

function markLink(loci: Loci, action: MarkerAction, structure: Structure, values: MarkerData) {
    const tMarker = values.tMarker

    const crossLinks = structure.crossLinkRestraints
    const elementCount = crossLinks.count
    const instanceCount = 1

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        applyMarkerAction(array, 0, elementCount * instanceCount, action)
        changed = true
    } else if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const indices = crossLinks.getPairIndices(b.aIndex, b.aUnit, b.bIndex, b.bUnit)
            if (indices) {
                for (let i = 0, il = indices.length; i < il; ++i) {
                    const idx = indices[i]
                    if (applyMarkerAction(array, idx, idx + 1, action) && !changed) {
                        changed = true
                    }
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