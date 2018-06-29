/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element, StructureProperties } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from './util/common';
import { deepEqual } from 'mol-util';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { Segmentation } from 'mol-data/int';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { Vec3 } from 'mol-math/linear-algebra';
import { getPolymerElementCount } from './util/polymer';
import { MoleculeType } from 'mol-model/structure/model/types';

async function createPolymerBackboneCylinderMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)

    const chemCompMap = unit.model.properties.chemicalComponentMap
    const { elements } = unit
    const curV = Vec3.zero()
    const prevV = Vec3.zero()
    const pos = unit.conformation.invariantPosition
    const l = Element.Location(unit)

    if (Unit.isAtomic(unit)) {
        const { polymerSegments, residueSegments } = unit.model.atomicHierarchy
        const polymerIt = Segmentation.transientSegments(polymerSegments, elements);
        const residuesIt = Segmentation.transientSegments(residueSegments, elements);

        let i = 0
        let first = true

        while (polymerIt.hasNext) {
            const polymerSegment = polymerIt.move();
            residuesIt.setSegment(polymerSegment);
            first = true
            while (residuesIt.hasNext) {
                const residueSegment = residuesIt.move();
                l.element = elements[residueSegment.start];

                const compId = StructureProperties.residue.label_comp_id(l)
                const cc = chemCompMap.get(compId)
                const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown
                let traceName = ''
                if (moleculeType === MoleculeType.protein) {
                    traceName = 'CA'
                } else if (moleculeType === MoleculeType.DNA || moleculeType === MoleculeType.RNA) {
                    traceName = 'P'
                }

                for (let j = residueSegment.start, _j = residueSegment.end; j < _j; j++) {
                    l.element = elements[j];
                    if (StructureProperties.atom.label_atom_id(l) === traceName) break
                }
                pos(l.element, curV)

                if (!first) {
                    builder.setId(residueSegment.start)
                    builder.addCylinder(prevV, curV, 1, { radiusTop: 0.2, radiusBottom: 0.2 })
                } else {
                    first = false
                }

                Vec3.copy(prevV, curV)

                if (i % 10000 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
                }
                ++i
            }
        }
    } else if (Unit.isSpheres(unit)) {
        let prevSeqIdEnd = -1
        for (let i = 0, il = elements.length; i < il; ++i) {
            l.element = elements[i]
            if (StructureProperties.entity.type(l) !== 'polymer') continue;

            pos(elements[i], curV)
            const seqIdBegin = StructureProperties.coarse.seq_id_begin(l)
            const seqIdEnd = StructureProperties.coarse.seq_id_end(l)

            pos(elements[i], curV)

            if (seqIdBegin - 1 === prevSeqIdEnd) {
                builder.setId(i)
                builder.addCylinder(prevV, curV, 1, { radiusTop: 0.2, radiusBottom: 0.2 })
            }

            Vec3.copy(prevV, curV)
            prevSeqIdEnd = seqIdEnd

            if (i % 10000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
            }
        }
    }

    return builder.getMesh()
}

export const DefaultPolymerBackboneProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type PolymerBackboneProps = Partial<typeof DefaultPolymerBackboneProps>

export function PolymerBackboneVisual(): UnitsVisual<PolymerBackboneProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultPolymerBackboneProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PolymerBackboneProps = {}) {
            currentProps = Object.assign({}, DefaultPolymerBackboneProps, props)
            currentGroup = group

            const { colorTheme, unitKinds } = { ...DefaultPolymerBackboneProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createPolymerBackboneCylinderMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)

            if (ctx.shouldUpdate) await ctx.update('Computing trace transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing trace colors');
            const color = createColors(group, elementCount, colorTheme)

            if (ctx.shouldUpdate) await ctx.update('Computing trace marks');
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
        async update(ctx: RuntimeContext, props: PolymerBackboneProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createPolymerBackboneCylinderMesh(ctx, unit, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                const elementCount = currentGroup.elements.length
                if (ctx.shouldUpdate) await ctx.update('Computing trace colors');
                createColors(currentGroup, elementCount, newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
        },
        destroy() {
            // TODO
        }
    }
}
