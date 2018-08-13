/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Unit, Structure } from 'mol-model/structure';
import { Mat4 } from 'mol-math/linear-algebra'

import { createUniformColor, ColorData, createElementColor, createElementInstanceColor, createInstanceColor } from '../../../../util/color-data';
import { createUniformSize, SizeData, createElementSize, createElementInstanceSize, createInstanceSize } from '../../../../util/size-data';
import { ColorThemeProps, SizeThemeProps } from '../../../../theme';
import { ColorTheme } from '../../../../theme/structure/color';
import { ValueCell } from 'mol-util';
import { LocationIterator } from './location-iterator';
import { Mesh } from '../../../../shape/mesh';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../../util/mesh-data';
import { MeshProps, createMeshValues, createRenderableState } from '../../../util';
import { StructureProps } from '../..';
import { createMarkers } from '../../../../util/marker-data';
import { createMeshRenderObject } from 'mol-gl/render-object';
import { SizeTheme } from '../../../../theme/structure/size';

export function createTransforms({ units }: Unit.SymmetryGroup, transforms?: ValueCell<Float32Array>) {
    const unitCount = units.length
    const n = unitCount * 16
    const array = transforms && transforms.ref.value.length >= n ? transforms.ref.value : new Float32Array(n)
    for (let i = 0; i < unitCount; i++) {
        Mat4.toArray(units[i].conformation.operator.matrix, array, i * 16)
    }
    return transforms ? ValueCell.update(transforms, array) : ValueCell.create(array)
}

const identityTransform = new Float32Array(16)
Mat4.toArray(Mat4.identity(), identityTransform, 0)
export function createIdentityTransform(transforms?: ValueCell<Float32Array>) {
    return transforms ? ValueCell.update(transforms, identityTransform) : ValueCell.create(identityTransform)
}

export function createColors(locationIt: LocationIterator, props: ColorThemeProps, colorData?: ColorData) {
    const colorTheme = ColorTheme(props)
    switch (colorTheme.kind) {
        case 'uniform': return createUniformColor(locationIt, colorTheme.color, colorData)
        case 'element': return createElementColor(locationIt, colorTheme.color, colorData)
        case 'elementInstance': return createElementInstanceColor(locationIt, colorTheme.color, colorData)
        case 'instance': return createInstanceColor(locationIt, colorTheme.color, colorData)
    }
}

export function createSizes(locationIt: LocationIterator, props: SizeThemeProps, sizeData?: SizeData): SizeData {
    const sizeTheme = SizeTheme(props)
    switch (sizeTheme.kind) {
        case 'uniform': return createUniformSize(locationIt, sizeTheme.size, sizeData)
        case 'element': return createElementSize(locationIt, sizeTheme.size, sizeData)
        case 'elementInstance': return createElementInstanceSize(locationIt, sizeTheme.size, sizeData)
        case 'instance': return createInstanceSize(locationIt, sizeTheme.size, sizeData)
    }
}

type StructureMeshProps = Required<MeshProps & StructureProps>

function _createMeshValues(transforms: ValueCell<Float32Array>, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): MeshValues {
    const { instanceCount, elementCount } = locationIt
    const color = createColors(locationIt, props.colorTheme)
    const marker = createMarkers(instanceCount * elementCount)

    const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

    return {
        ...getMeshData(mesh),
        ...color,
        ...marker,
        aTransform: transforms,
        elements: mesh.indexBuffer,
        ...createMeshValues(props, counts)
    }
}

export function createComplexMeshValues(structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): MeshValues {
    const transforms = createIdentityTransform()
    return _createMeshValues(transforms, mesh, locationIt, props)
}

export function createUnitsMeshValues(group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): MeshValues {
    const transforms = createTransforms(group)
    return _createMeshValues(transforms, mesh, locationIt, props)
}

export function createComplexMeshRenderObject(structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const values = createComplexMeshValues(structure, mesh, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export function createUnitsMeshRenderObject(group: Unit.SymmetryGroup, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps) {
    const values = createUnitsMeshValues(group, mesh, locationIt, props)
    const state = createRenderableState(props)
    return createMeshRenderObject(values, state)
}

export function updateComplexMeshRenderObject(structure: Structure, mesh: Mesh, locationIt: LocationIterator, props: StructureMeshProps): MeshValues {
    const transforms = createIdentityTransform()
    return _createMeshValues(transforms, mesh, locationIt, props)
}