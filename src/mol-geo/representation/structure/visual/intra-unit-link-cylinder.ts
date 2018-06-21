/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Link } from 'mol-model/structure';
import { UnitsVisual, DefaultStructureProps } from '../index';
import { RuntimeContext } from 'mol-task'
import { DefaultLinkCylinderProps, LinkCylinderProps, createLinkCylinderMesh } from './util/link';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Vec3 } from 'mol-math/linear-algebra';
// import { createUniformColor } from '../../../util/color-data';
import { defaults } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { chainIdLinkColorData } from '../../../theme/structure/color/chain-id';
import { createTransforms } from './util/common';

async function createIntraUnitLinkCylinderMesh(ctx: RuntimeContext, unit: Unit, props: LinkCylinderProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const elements = unit.elements;
    const links = unit.links
    const { edgeCount, a, b, edgeProps, offset } = links
    const { order: _order, flags: _flags } = edgeProps

    if (!edgeCount) return Mesh.createEmpty(mesh)

    const vRef = Vec3.zero()
    const pos = unit.conformation.invariantPosition

    const builderProps = {
        linkCount: edgeCount * 2,
        referencePosition: (edgeIndex: number) => {
            let aI = a[edgeIndex], bI = b[edgeIndex];
            if (aI > bI) [aI, bI] = [bI, aI]
            for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
                if (b[i] !== bI) return pos(elements[b[i]], vRef)
            }
            for (let i = offset[bI], il = offset[bI + 1]; i < il; ++i) {
                if (a[i] !== aI) return pos(elements[a[i]], vRef)
            }
            return null
        },
        position: (posA: Vec3, posB: Vec3, edgeIndex: number) => {
            pos(elements[a[edgeIndex]], posA)
            pos(elements[b[edgeIndex]], posB)
        },
        order: (edgeIndex: number) => _order[edgeIndex],
        flags: (edgeIndex: number) => _flags[edgeIndex]
    }

    return createLinkCylinderMesh(ctx, builderProps, props, mesh)
}

export const DefaultIntraUnitLinkProps = {
    ...DefaultStructureProps,
    ...DefaultLinkCylinderProps,
    sizeTheme: { name: 'physical', factor: 0.3 } as SizeTheme,
    flipSided: false,
    flatShaded: false,
}
export type IntraUnitLinkProps = Partial<typeof DefaultIntraUnitLinkProps>

export function IntraUnitLinkVisual(): UnitsVisual<IntraUnitLinkProps> {
    const renderObjects: RenderObject[] = []
    let cylinders: MeshRenderObject
    let currentProps: typeof DefaultIntraUnitLinkProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        renderObjects,
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: IntraUnitLinkProps = {}) {
            currentProps = Object.assign({}, DefaultIntraUnitLinkProps, props)

            renderObjects.length = 0 // clear
            currentGroup = group

            const unit = group.units[0]
            const elementCount = Unit.isAtomic(unit) ? unit.links.edgeCount * 2 : 0
            const instanceCount = group.units.length

            mesh = await createIntraUnitLinkCylinderMesh(ctx, unit, currentProps)

            if (ctx.shouldUpdate) await ctx.update('Computing link transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing link colors');
            // const color = createUniformColor({ value: 0xFF0000 })
            const color = chainIdLinkColorData({ group, elementCount })

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
        async update(ctx: RuntimeContext, props: IntraUnitLinkProps) {
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
            return getLinkLoci(pickingId, currentGroup, cylinders.id)
        },
        mark(loci: Loci, action: MarkerAction) {
            markLink(loci, action, currentGroup, cylinders.values)
        },
        destroy() {
            // TODO
        }
    }
}

function getLinkLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number) {
    const { objectId, instanceId, elementId } = pickingId
    const unit = group.units[instanceId]
    if (id === objectId && Unit.isAtomic(unit)) {
        return Link.Loci([{
            aUnit: unit,
            aIndex: unit.links.a[elementId],
            bUnit: unit,
            bIndex: unit.links.b[elementId]
        }])
    }
    return EmptyLoci
}

function markLink(loci: Loci, action: MarkerAction, group: Unit.SymmetryGroup, values: MarkerData) {
    const tMarker = values.tMarker
    const unit = group.units[0]
    if (!Unit.isAtomic(unit)) return

    const elementCount = unit.links.edgeCount * 2
    const instanceCount = group.units.length

    let changed = false
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        applyMarkerAction(array, 0, elementCount * instanceCount, action)
        changed = true
    } else if (Link.isLoci(loci)) {
        for (const b of loci.links) {
            const unitIdx = Unit.findUnitById(b.aUnit.id, group.units)
            if (unitIdx !== -1) {
                const _idx = unit.links.getDirectedEdgeIndex(b.aIndex, b.bIndex)
                if (_idx !== -1) {
                    const idx = _idx
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