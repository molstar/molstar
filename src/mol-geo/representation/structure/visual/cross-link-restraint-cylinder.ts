/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
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
    const renderObjects: RenderObject[] = []
    let cylinders: MeshRenderObject
    let currentProps: typeof DefaultCrossLinkRestraintProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        renderObjects,
        async create(ctx: RuntimeContext, structure: Structure, props: CrossLinkRestraintProps = {}) {
            currentProps = Object.assign({}, DefaultCrossLinkRestraintProps, props)

            renderObjects.length = 0 // clear
            currentStructure = structure

            const elementCount = structure.crossLinkRestraints.count
            const instanceCount = 1

            mesh = await createCrossLinkRestraintCylinderMesh(ctx, structure, currentProps)

            if (ctx.shouldUpdate) await ctx.update('Computing link transforms');
            const transforms = createIdentityTransform()

            if (ctx.shouldUpdate) await ctx.update('Computing link colors');
            const color = createUniformColor({ value: 0x119911 }) // TODO

            if (ctx.shouldUpdate) await ctx.update('Computing link marks');
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

            cylinders = createMeshRenderObject(values, state)
            console.log(values, instanceCount, elementCount)
            renderObjects.push(cylinders)
        },
        async update(ctx: RuntimeContext, props: CrossLinkRestraintProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!cylinders) return false

            // TODO create in-place
            if (currentProps.radialSegments !== newProps.radialSegments) return false

            updateMeshValues(cylinders.values, newProps)
            updateRenderableState(cylinders.state, newProps)

            return false
        },
        getLoci(pickingId: PickingId) {
            return getLinkLoci(pickingId, currentStructure, cylinders.id)
        },
        mark(loci: Loci, action: MarkerAction) {
            markLink(loci, action, currentStructure, cylinders.values)
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
            return Link.Loci([{
                aUnit: pair.unitA,
                aIndex: pair.indexA,
                bUnit: pair.unitB,
                bIndex: pair.indexB
            }])
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