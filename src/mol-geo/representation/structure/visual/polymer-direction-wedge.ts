/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { MeshRenderObject } from 'mol-gl/render-object'
import { Unit, StructureElement } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '..';
import { RuntimeContext } from 'mol-task'
import { createColors, createUnitsMeshRenderObject } from './util/common';
import { markElement } from './util/element';
import { deepEqual } from 'mol-util';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { OrderedSet } from 'mol-data/int';
import { MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { updateMeshValues, updateRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerTraceIterator, createCurveSegmentState, interpolateCurveSegment } from './util/polymer';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';
import { SecondaryStructureType, MoleculeType } from 'mol-model/structure/model/types';
import { StructureElementIterator } from './util/location-iterator';

const t = Mat4.identity()
const sVec = Vec3.zero()
const n0 = Vec3.zero()
const n1 = Vec3.zero()
const upVec = Vec3.zero()

async function createPolymerDirectionWedgeMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount direction', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegments = 1

    const state = createCurveSegmentState(linearSegments)
    const { normalVectors, binormalVectors } = state

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()
        builder.setId(v.center.element)

        const isNucleic = v.moleculeType === MoleculeType.DNA || v.moleculeType === MoleculeType.RNA
        const isSheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta)
        const tension = (isNucleic || isSheet) ? 0.5 : 0.9

        interpolateCurveSegment(state, v, tension)

        if ((isSheet && !v.secStrucChange) || !isSheet) {

            let width = 0.5, height = 1.2, depth = 0.6
            if (isNucleic) {
                Vec3.fromArray(n0, binormalVectors, 0)
                Vec3.fromArray(n1, binormalVectors, 3)
                Vec3.normalize(upVec, Vec3.add(upVec, n0, n1))
                depth = 0.9
            } else {
                Vec3.fromArray(n0, normalVectors, 0)
                Vec3.fromArray(n1, normalVectors, 3)
                Vec3.normalize(upVec, Vec3.add(upVec, n0, n1))
            }

            Mat4.targetTo(t, v.p3, v.p1, upVec)
            Mat4.mul(t, t, Mat4.rotY90Z180)
            Mat4.scale(t, t, Vec3.set(sVec, height, width, depth))
            Mat4.setTranslation(t, v.p2)
            builder.addWedge(t)
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Polymer direction mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerDirectionProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type PolymerDirectionProps = Partial<typeof DefaultPolymerDirectionProps>

export function PolymerDirectionVisual(): UnitsVisual<PolymerDirectionProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultPolymerDirectionProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PolymerDirectionProps = {}) {
            currentProps = Object.assign({}, DefaultPolymerDirectionProps, props)
            currentGroup = group

            const { unitKinds } = { ...DefaultPolymerDirectionProps, ...props }
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createPolymerDirectionWedgeMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)

            const locationIt = StructureElementIterator.fromGroup(group)
            renderObject = createUnitsMeshRenderObject(group, mesh, locationIt, currentProps)
        },
        async update(ctx: RuntimeContext, props: PolymerDirectionProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createPolymerDirectionWedgeMesh(ctx, unit, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                if (ctx.shouldUpdate) await ctx.update('Computing direction colors');
                createColors(StructureElementIterator.fromGroup(currentGroup), newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (renderObject.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as StructureElement.UnitIndex);
                return StructureElement.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(renderObject.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
