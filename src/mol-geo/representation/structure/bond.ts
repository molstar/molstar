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
import { UnitsVisual, DefaultStructureProps } from './index';
import { RuntimeContext } from 'mol-task'
import { createTransforms } from './utils';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../util/mesh-data';
import { Mesh } from '../../shape/mesh';
import { PickingId } from '../../util/picking';
import { MeshBuilder } from '../../shape/mesh-builder';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { createUniformColor } from '../../util/color-data';
import { defaults } from 'mol-util';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction, createMarkers } from '../../util/marker-data';

async function createBondMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh)

    const elements = unit.elements;
    const bonds = unit.links
    const { edgeCount, a, b } = bonds

    if (!edgeCount) return Mesh.createEmpty(mesh)

    // TODO calculate vertextCount properly
    const vertexCount = 32 * edgeCount
    const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

    const va = Vec3.zero()
    const vb = Vec3.zero()
    const vt = Vec3.zero()
    const m = Mat4.identity()

    const pos = unit.conformation.invariantPosition
    // const l = Element.Location()
    // l.unit = unit

    for (let edgeIndex = 0, _eI = edgeCount * 2; edgeIndex < _eI; ++edgeIndex) {
        const aI = elements[a[edgeIndex]], bI = elements[b[edgeIndex]];
        // each edge is included twice because of the "adjacency list" structure
        // keep only the 1st occurence.
        if (aI >= bI) continue;
        pos(aI, va)
        pos(bI, vb)

        Vec3.scale(vt, Vec3.add(vt, va, vb), 0.5)
        Vec3.makeRotation(m, Vec3.create(0, 1, 0), Vec3.sub(vb, vb, va))
        Mat4.setTranslation(m, vt)
        
        meshBuilder.setId(edgeIndex)
        meshBuilder.addCylinder(m, { radiusTop: 0.2, radiusBottom: 0.2 })

        if (edgeIndex % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Cylinder mesh', current: edgeIndex, max: edgeCount });
        }
    }

    return meshBuilder.getMesh()
}

export const DefaultBondProps = {
    ...DefaultStructureProps,
    flipSided: false,
    flatShaded: false,
}
export type BondProps = Partial<typeof DefaultBondProps>

export default function IntraUnitBondVisual(): UnitsVisual<BondProps> {
    const renderObjects: RenderObject[] = []
    let cylinders: MeshRenderObject
    let currentProps: typeof DefaultBondProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    // let vertexMap: VertexMap

    return {
        renderObjects,
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: BondProps = {}) {
            currentProps = Object.assign({}, DefaultBondProps, props)

            renderObjects.length = 0 // clear
            currentGroup = group

            const unit = group.units[0]
            const elementCount = Unit.isAtomic(unit) ? unit.links.edgeCount * 2 : 0
            const instanceCount = group.units.length

            mesh = await createBondMesh(ctx, unit)

            // console.log(mesh)
            // vertexMap = VertexMap.fromMesh(mesh)

            if (ctx.shouldUpdate) await ctx.update('Computing bond transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing bond colors');
            const color = createUniformColor({ value: 0xFF0000 })

            if (ctx.shouldUpdate) await ctx.update('Computing bond marks');
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
        async update(ctx: RuntimeContext, props: BondProps) {
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
            const { objectId, instanceId, elementId } = pickingId
            const unit = currentGroup.units[instanceId]
            if (cylinders.id === objectId && Unit.isAtomic(unit)) {
                return Link.Loci([{
                    aUnit: unit,
                    aIndex: unit.links.a[elementId],
                    bUnit: unit,
                    bIndex: unit.links.b[elementId]
                }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            const group = currentGroup
            const tMarker = cylinders.values.tMarker
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
                        const _idx = unit.links.getEdgeIndex(b.aIndex, b.bIndex)
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
        },
        destroy() {
            // TODO
        }
    }
}
