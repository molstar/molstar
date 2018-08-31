/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { Visual } from '..';
import { MeshRenderObject } from 'mol-gl/render-object';
import { Mesh } from '../../mesh/mesh';
import { RuntimeContext } from 'mol-task';
import { LocationIterator } from '../../util/location-iterator';
import { createComplexMeshRenderObject, createColors } from './visual/util/common';
import { StructureProps, DefaultStructureMeshProps, MeshUpdateState } from '.';
import { deepEqual, ValueCell } from 'mol-util';
import { updateMeshValues, updateRenderableState } from '../util';
import { PickingId } from '../../util/picking';
import { Loci, isEveryLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../util/marker-data';
import { Interval } from 'mol-data/int';

export interface  ComplexVisual<P extends StructureProps> extends Visual<Structure, P> { }

export const DefaultComplexMeshProps = {
    ...DefaultStructureMeshProps
}
export type ComplexMeshProps = typeof DefaultComplexMeshProps

export interface ComplexMeshVisualBuilder<P extends ComplexMeshProps> {
    defaultProps: P
    createMesh(ctx: RuntimeContext, structure: Structure, props: P, mesh?: Mesh): Promise<Mesh>
    createLocationIterator(structure: Structure): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    mark(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean): boolean,
    setUpdateState(state: MeshUpdateState, newProps: P, currentProps: P): void
}

export function ComplexMeshVisual<P extends ComplexMeshProps>(builder: ComplexMeshVisualBuilder<P>): ComplexVisual<P> {
    const { defaultProps, createMesh, createLocationIterator, getLoci, mark, setUpdateState } = builder
    const updateState = MeshUpdateState.create()

    let renderObject: MeshRenderObject | undefined
    let currentProps: P
    let mesh: Mesh
    let currentStructure: Structure
    let locationIt: LocationIterator
    let conformationHashCode: number

    async function create(ctx: RuntimeContext, structure: Structure, props: Partial<P> = {}) {
        currentProps = Object.assign({}, defaultProps, props)
        currentStructure = structure

        conformationHashCode = Structure.conformationHash(currentStructure)
        mesh = await createMesh(ctx, currentStructure, currentProps, mesh)

        locationIt = createLocationIterator(structure)
        renderObject = await createComplexMeshRenderObject(ctx, structure, mesh, locationIt, currentProps)
    }

    async function update(ctx: RuntimeContext, props: Partial<P>) {
        const newProps = Object.assign({}, currentProps, props)

        if (!renderObject) return false

        locationIt.reset()
        MeshUpdateState.reset(updateState)
        setUpdateState(updateState, newProps, currentProps)

        const newConformationHashCode = Structure.conformationHash(currentStructure)
        if (newConformationHashCode !== conformationHashCode) {
            conformationHashCode = newConformationHashCode
            updateState.createMesh = true
        }

        if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateState.createMesh = true
        if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateState.updateColor = true
        // if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) updateState.createMesh = true // TODO

        //

        if (updateState.createMesh) {
            mesh = await createMesh(ctx, currentStructure, newProps, mesh)
            ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
            updateState.updateColor = true
        }

        if (updateState.updateColor) {
            await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
        }

        updateMeshValues(renderObject.values, newProps)
        updateRenderableState(renderObject.state, newProps)

        currentProps = newProps
        return true
    }

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: Partial<P> = {}, structure?: Structure) {
            if (!structure && !currentStructure) {
                throw new Error('missing structure')
            } else if (structure && (!currentStructure || !renderObject)) {
                await create(ctx, structure, props)
            } else if (structure && structure.hashCode !== currentStructure.hashCode) {
                await create(ctx, structure, props)
            } else {
                if (structure && Structure.conformationHash(structure) !== Structure.conformationHash(currentStructure)) {
                    currentStructure = structure
                }
                await update(ctx, props)
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentStructure, renderObject.id) : EmptyLoci
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
                changed = mark(loci, currentStructure, apply)
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