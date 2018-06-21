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
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../../util/color-data';
import { defaults } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { createIdentityTransform } from './util/common';
// import { chainIdLinkColorData } from '../../../theme/structure/color/chain-id';

async function createInterUnitLinkCylinderMesh(ctx: RuntimeContext, structure: Structure, props: LinkCylinderProps, mesh?: Mesh) {
    const links = structure.links
    const { bondCount, bonds } = links

    if (!bondCount) return Mesh.createEmpty(mesh)

    function referencePosition(edgeIndex: number): Vec3 | null {
        // TODO
        return null
    }

    function position(posA: Vec3, posB: Vec3, edgeIndex: number): void {
        const b = bonds[edgeIndex]
        const uA = b.unitA, uB = b.unitB
        uA.conformation.position(uA.elements[b.indexA], posA)
        uB.conformation.position(uB.elements[b.indexB], posB)
    }

    function order(edgeIndex: number): number {
        return bonds[edgeIndex].order
    }

    const linkBuilder = { linkCount: bondCount, referencePosition, position, order }

    return createLinkCylinderMesh(ctx, linkBuilder, props, mesh)
}

export const DefaultInterUnitLinkProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
    flipSided: false,
    flatShaded: false,
}
export type InterUnitLinkProps = Partial<typeof DefaultInterUnitLinkProps>

export function InterUnitLinkVisual(): StructureVisual<InterUnitLinkProps> {
    const renderObjects: RenderObject[] = []
    let cylinders: MeshRenderObject
    let currentProps: typeof DefaultInterUnitLinkProps
    let mesh: Mesh
    let currentStructure: Structure

    return {
        renderObjects,
        async create(ctx: RuntimeContext, structure: Structure, props: InterUnitLinkProps = {}) {
            currentProps = Object.assign({}, DefaultInterUnitLinkProps, props)

            renderObjects.length = 0 // clear
            currentStructure = structure

            const elementCount = structure.links.bondCount
            const instanceCount = 1

            mesh = await createInterUnitLinkCylinderMesh(ctx, structure, currentProps)

            if (ctx.shouldUpdate) await ctx.update('Computing link transforms');
            const transforms = createIdentityTransform()

            if (ctx.shouldUpdate) await ctx.update('Computing link colors');
            const color = createUniformColor({ value: 0x999911 })
            // const color = chainIdLinkColorData({ group, elementCount })

            if (ctx.shouldUpdate) await ctx.update('Computing link marks');
            const marker = createMarkers(instanceCount * elementCount)

            const values: MeshValues = {
                ...getMeshData(mesh),
                aTransform: transforms,
                aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                ...color,
                ...marker,

                uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                uInstanceCount: ValueCell.create(instanceCount),
                uElementCount: ValueCell.create(elementCount),

                elements: mesh.indexBuffer,

                drawCount: ValueCell.create(mesh.triangleCount * 3),
                instanceCount: ValueCell.create(instanceCount),

                dDoubleSided: ValueCell.create(defaults(props.doubleSided, true)),
                dFlatShaded: ValueCell.create(defaults(props.flatShaded, false)),
                dFlipSided: ValueCell.create(defaults(props.flipSided, false)),
                dUseFog: ValueCell.create(defaults(props.useFog, true)),
            }
            const state: RenderableState = {
                depthMask: defaults(props.depthMask, true),
                visible: defaults(props.visible, true)
            }

            cylinders = createMeshRenderObject(values, state)
            renderObjects.push(cylinders)
        },
        async update(ctx: RuntimeContext, props: InterUnitLinkProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!cylinders) return false
            // TODO

            ValueCell.updateIfChanged(cylinders.values.uAlpha, newProps.alpha)
            ValueCell.updateIfChanged(cylinders.values.dDoubleSided, newProps.doubleSided)
            ValueCell.updateIfChanged(cylinders.values.dFlipSided, newProps.flipSided)
            ValueCell.updateIfChanged(cylinders.values.dFlatShaded, newProps.flatShaded)

            cylinders.state.visible = newProps.visible
            cylinders.state.depthMask = newProps.depthMask

            return true
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