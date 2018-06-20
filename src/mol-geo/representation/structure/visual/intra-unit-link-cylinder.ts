/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO multiple cylinders for higher bond orders

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Link } from 'mol-model/structure';
import { UnitsVisual, DefaultStructureProps } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms } from '../utils';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
// import { createUniformColor } from '../../../util/color-data';
import { defaults } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers, MarkerData } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { chainIdLinkColorData } from '../../../theme/structure/color/chain-id';

const DefaultLinkCylinderProps = {
    linkScale: 0.4,
    linkSpacing: 1,
    linkRadius: 0.25,
    radialSegments: 16
}
type LinkCylinderProps = typeof DefaultLinkCylinderProps

async function createLinkCylinderMesh(ctx: RuntimeContext, unit: Unit, props: LinkCylinderProps, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const elements = unit.elements;
    const links = unit.links
    const { edgeCount, a, b, edgeProps, offset } = links
    const orders = edgeProps.order

    if (!edgeCount) return Mesh.createEmpty(mesh)

    // approximate vertextCount, exact calculation would need to take link orders into account
    const vertexCount = 32 * edgeCount
    const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vd = Vec3.zero()
    const vc = Vec3.zero()
    const m = Mat4.identity()
    const mt = Mat4.identity()

    const vShift = Vec3.zero()
    const vCenter = Vec3.zero()
    const vRef = Vec3.zero()

    const { linkScale, linkSpacing, linkRadius, radialSegments } = props

    const cylinderParams = {
        height: 1,
        radiusTop: linkRadius,
        radiusBottom: linkRadius,
        radialSegments,
        openEnded: true
    }

    const pos = unit.conformation.invariantPosition
    // const l = Element.Location()
    // l.unit = unit

    // assumes aI < bI
    function getRefPos(aI: number, bI: number) {
        if (aI > bI) console.log('aI > bI')
        for (let i = offset[aI], il = offset[aI + 1]; i < il; ++i) {
            if (b[i] !== bI) return pos(elements[b[i]], vRef)
        }
        for (let i = offset[bI], il = offset[bI + 1]; i < il; ++i) {
            if (a[i] !== aI) return pos(elements[a[i]], vRef)
        }
        // console.log('no ref', aI, bI, unit.model.atomicHierarchy.atoms.auth_atom_id.value(aI), unit.model.atomicHierarchy.atoms.auth_atom_id.value(bI), offset[aI], offset[aI + 1], offset[bI], offset[bI + 1], offset)
        return null
    }

    for (let edgeIndex = 0, _eI = edgeCount * 2; edgeIndex < _eI; ++edgeIndex) {
        const aI = elements[a[edgeIndex]], bI = elements[b[edgeIndex]];
        // Each edge is included twice to allow for coloring/picking
        // the half closer to the first vertex, i.e. vertex a.
        pos(aI, va)
        pos(bI, vb)
        const d = Vec3.distance(va, vb)

        Vec3.sub(vd, vb, va)
        Vec3.scale(vd, Vec3.normalize(vd, vd), d / 4)
        Vec3.add(vc, va, vd)
        // ensure both edge halfs are pointing in the the same direction so the triangles align
        if (aI > bI) Vec3.scale(vd, vd, -1)
        Vec3.makeRotation(m, Vec3.create(0, 1, 0), vd)

        const order = orders[edgeIndex]
        meshBuilder.setId(edgeIndex)
        cylinderParams.height = d / 2

        if (order === 2 || order === 3) {
            const multiRadius = linkRadius * (linkScale / (0.5 * order))
            const absOffset = (linkRadius - multiRadius) * linkSpacing

            if (aI < bI) {
                calculateShiftDir(vShift, va, vb, getRefPos(a[edgeIndex], b[edgeIndex]))
            } else {
                calculateShiftDir(vShift, vb, va, getRefPos(b[edgeIndex], a[edgeIndex]))
            }
            Vec3.setMagnitude(vShift, vShift, absOffset)

            cylinderParams.radiusTop = multiRadius
            cylinderParams.radiusBottom = multiRadius

            if (order === 3) {
                Mat4.fromTranslation(mt, vc)
                Mat4.mul(mt, mt, m)
                meshBuilder.addCylinder(mt, cylinderParams)
            }

            Vec3.add(vCenter, vc, vShift)
            Mat4.fromTranslation(mt, vCenter)
            Mat4.mul(mt, mt, m)
            meshBuilder.addCylinder(mt, cylinderParams)

            Vec3.sub(vCenter, vc, vShift)
            Mat4.fromTranslation(mt, vCenter)
            Mat4.mul(mt, mt, m)
            meshBuilder.addCylinder(mt, cylinderParams)
        } else {
            cylinderParams.radiusTop = linkRadius
            cylinderParams.radiusBottom = linkRadius

            Mat4.setTranslation(m, vc)
            meshBuilder.addCylinder(m, cylinderParams)
        }

        if (edgeIndex % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Cylinder mesh', current: edgeIndex, max: edgeCount });
        }
    }

    return meshBuilder.getMesh()
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

            mesh = await createLinkCylinderMesh(ctx, unit, currentProps)

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

const tmpShiftV12 = Vec3.zero()
const tmpShiftV13 = Vec3.zero()

/** Calculate 'shift' direction that is perpendiculat to v1 - v2 and goes through v3 */
function calculateShiftDir (out: Vec3, v1: Vec3, v2: Vec3, v3: Vec3 | null) {
    Vec3.sub(tmpShiftV12, v1, v2)

    if (v3 !== null) {
        Vec3.sub(tmpShiftV13, v1, v3)
    } else {
        Vec3.copy(tmpShiftV13, v1)  // no reference point, use v1
    }
    Vec3.normalize(tmpShiftV13, tmpShiftV13)

    // ensure v13 and v12 are not colinear
    let dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
    if (1 - Math.abs(dp) < 1e-5) {
        Vec3.set(tmpShiftV13, 1, 0, 0)
        dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
        if (1 - Math.abs(dp) < 1e-5) {
            Vec3.set(tmpShiftV13, 0, 1, 0)
            dp = Vec3.dot(tmpShiftV12, tmpShiftV13)
        }
    }

    Vec3.setMagnitude(tmpShiftV12, tmpShiftV12, dp)
    Vec3.sub(tmpShiftV13, tmpShiftV13, tmpShiftV12)
    return Vec3.normalize(out, tmpShiftV13)
}