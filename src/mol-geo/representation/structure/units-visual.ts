/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { RepresentationProps, Visual } from '..';
import { DefaultStructureMeshProps, MeshUpdateState } from '.';
import { RuntimeContext } from 'mol-task';
import { PickingId } from '../../util/picking';
import { LocationIterator } from '../../util/location-iterator';
import { Mesh } from '../../mesh/mesh';
import { MarkerAction, applyMarkerAction } from '../../util/marker-data';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MeshRenderObject } from 'mol-gl/render-object';
import { createUnitsMeshRenderObject, createColors } from './visual/util/common';
import { deepEqual, ValueCell } from 'mol-util';
import { updateMeshValues, updateRenderableState } from '../util';
import { Interval } from 'mol-data/int';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<Unit.SymmetryGroup, P> { }

export const DefaultUnitsMeshProps = {
    ...DefaultStructureMeshProps,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type UnitsMeshProps = typeof DefaultUnitsMeshProps

export interface UnitsMeshVisualBuilder<P extends UnitsMeshProps> {
    defaultProps: P
    createMesh(ctx: RuntimeContext, unit: Unit, props: P, mesh?: Mesh): Promise<Mesh>
    createLocationIterator(group: Unit.SymmetryGroup): LocationIterator
    getLoci(pickingId: PickingId, group: Unit.SymmetryGroup, id: number): Loci
    mark(loci: Loci, group: Unit.SymmetryGroup, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: MeshUpdateState, newProps: P, currentProps: P): void
}

export function UnitsMeshVisual<P extends UnitsMeshProps>(builder: UnitsMeshVisualBuilder<P>): UnitsVisual<P> {
    const { defaultProps, createMesh, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = MeshUpdateState.create()

    let renderObject: MeshRenderObject | undefined
    let currentProps: P
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    let locationIt: LocationIterator

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: Partial<P> = {}) {
            currentProps = Object.assign({}, defaultProps, props)
            currentGroup = group

            const unit = group.units[0]
            mesh = currentProps.unitKinds.includes(unit.kind)
                ? await createMesh(ctx, unit, currentProps, mesh)
                : Mesh.createEmpty(mesh)

            locationIt = createLocationIterator(group)
            renderObject = createUnitsMeshRenderObject(group, mesh, locationIt, currentProps)
        },
        async update(ctx: RuntimeContext, props: Partial<P>) {
            const newProps = Object.assign({}, currentProps, props)
            const unit = currentGroup.units[0]

            if (!renderObject) return false

            locationIt.reset()
            MeshUpdateState.reset(updateState)
            setUpdateState(updateState, newProps, currentProps)

            if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) {
                updateState.createMesh = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateState.updateColor = true
            }

            //

            if (updateState.createMesh) {
                mesh = await createMesh(ctx, unit, newProps, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateState.updateColor = true
            }

            if (updateState.updateColor) {
                createColors(locationIt, newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentGroup, renderObject.id) : EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!renderObject) return
            const { tMarker } = renderObject.values
            const { groupCount, instanceCount } = locationIt

            function apply(interval: Interval) {
                const start = Interval.start(interval)
                const end = Interval.end(interval)
                return applyMarkerAction(tMarker.ref.value.array, start, end, action)
            }

            let changed = false
            if (isEveryLoci(loci)) {
                apply(Interval.ofBounds(0, groupCount * instanceCount))
                changed = true
            } else {
                changed = mark(loci, currentGroup, apply)
            }
            if (changed) {
                ValueCell.update(tMarker, tMarker.ref.value)
            }
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}