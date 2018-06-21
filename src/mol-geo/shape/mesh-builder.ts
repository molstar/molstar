/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4, Mat3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

import Box, { BoxProps } from '../primitive/box';
import Cylinder, { CylinderProps } from '../primitive/cylinder';
import Icosahedron, { IcosahedronProps } from '../primitive/icosahedron';
import { Mesh } from './mesh';
import { getNormalMatrix } from '../util';

type Primitive = {
    vertices: Float32Array
    normals: Float32Array
    indices: Uint32Array
}

export interface MeshBuilder {
    add(t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices?: Uint32Array): void
    addBox(t: Mat4, props?: BoxProps): void
    addCylinder(start: Vec3, end: Vec3, lengthScale: number, props: CylinderProps): void
    addDoubleCylinder(start: Vec3, end: Vec3, lengthScale: number, shift: Vec3, props: CylinderProps): void
    addFixedCountDashedCylinder(start: Vec3, end: Vec3, lengthScale: number, segmentCount: number, props: CylinderProps): void
    addIcosahedron(center: Vec3, radius: number, detail: number): void
    setId(id: number): void
    getMesh(): Mesh
}

const cylinderMap = new Map<string, Primitive>()
const icosahedronMap = new Map<string, Primitive>()

const up = Vec3.create(0, 1, 0)
const tmpV = Vec3.zero()
const tmpMat3 = Mat3.zero()

const tmpCylinderDir = Vec3.zero()
const tmpCylinderMatDir = Vec3.zero()
const tmpCylinderCenter = Vec3.zero()
const tmpCylinderMat = Mat4.zero()
// const tmpCylinderMatTrans = Mat4.zero()
const tmpCylinderStart = Vec3.zero()

function setCylinderMat(m: Mat4, start: Vec3, dir: Vec3, length: number) {
    Vec3.setMagnitude(tmpCylinderMatDir, dir, length / 2)
    Vec3.add(tmpCylinderCenter, start, tmpCylinderMatDir)
    // ensure the direction use to create the rotation is always pointing in the same
    // direction so the triangles of adjacent cylinder will line up
    if (Vec3.dot(tmpCylinderMatDir, up) < 0) Vec3.scale(tmpCylinderMatDir, tmpCylinderMatDir, -1)
    Vec3.makeRotation(m, up, tmpCylinderMatDir)
    // Mat4.fromTranslation(tmpCylinderMatTrans, tmpCylinderCenter)
    // Mat4.mul(m, tmpCylinderMatTrans, m)
    Mat4.setTranslation(m, tmpCylinderCenter)
    return m
}

function getCylinder(props: CylinderProps) {
    const key = JSON.stringify(props)
    let cylinder = cylinderMap.get(key)
    if (cylinder === undefined) {
        cylinder = Cylinder(props)
        cylinderMap.set(key, cylinder)
    }
    return cylinder
}

const tmpIcosahedronMat = Mat4.identity()

function setIcosahedronMat(m: Mat4, center: Vec3) {
    Mat4.setTranslation(m, center)
    return m
}

function getIcosahedron(props: IcosahedronProps) {
    const key = JSON.stringify(props)
    let icosahedron = icosahedronMap.get(key)
    if (icosahedron === undefined) {
        icosahedron = Icosahedron(props)
        icosahedronMap.set(key, icosahedron)
    }
    return icosahedron
}

// TODO cache primitives based on props

export namespace MeshBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, mesh?: Mesh): MeshBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.vertexBuffer.ref.value : initialCount);
        const normals = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.normalBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, mesh ? mesh.indexBuffer.ref.value : initialCount * 3);

        const ids = ChunkedArray.create(Float32Array, 1, chunkSize, mesh ? mesh.idBuffer.ref.value : initialCount);
        const offsets = ChunkedArray.create(Uint32Array, 1, chunkSize, mesh ? mesh.offsetBuffer.ref.value : initialCount);

        let currentId = -1

        function add(t: Mat4, _vertices: Float32Array, _normals: Float32Array, _indices: Uint32Array) {
            const { elementCount } = vertices
            const n = getNormalMatrix(tmpMat3, t)
            for (let i = 0, il = _vertices.length; i < il; i += 3) {
                // position
                Vec3.fromArray(tmpV, _vertices, i)
                Vec3.transformMat4(tmpV, tmpV, t)
                ChunkedArray.add3(vertices, tmpV[0], tmpV[1], tmpV[2]);
                // normal
                Vec3.fromArray(tmpV, _normals, i)
                Vec3.transformMat3(tmpV, tmpV, n)
                ChunkedArray.add3(normals, tmpV[0], tmpV[1], tmpV[2]);
                // id
                ChunkedArray.add(ids, currentId);
            }
            for (let i = 0, il = _indices.length; i < il; i += 3) {
                ChunkedArray.add3(indices, _indices[i] + elementCount, _indices[i + 1] + elementCount, _indices[i + 2] + elementCount);
            }
        }

        return {
            add,
            addBox: (t: Mat4, props?: BoxProps) => {
                const box = Box(props)
                add(t, box.vertices, box.normals, box.indices)
            },
            addCylinder: (start: Vec3, end: Vec3, lengthScale: number, props: CylinderProps) => {
                const d = Vec3.distance(start, end) * lengthScale
                props.height = d
                const { vertices, normals, indices } = getCylinder(props)
                Vec3.sub(tmpCylinderDir, end, start)
                setCylinderMat(tmpCylinderMat, start, tmpCylinderDir, d)
                add(tmpCylinderMat, vertices, normals, indices)
            },
            addDoubleCylinder: (start: Vec3, end: Vec3, lengthScale: number, shift: Vec3, props: CylinderProps) => {
                const d = Vec3.distance(start, end) * lengthScale
                props.height = d
                const { vertices, normals, indices } = getCylinder(props)
                Vec3.sub(tmpCylinderDir, end, start)
                // positivly shifted cylinder
                Vec3.add(tmpCylinderStart, start, shift)
                setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d)
                add(tmpCylinderMat, vertices, normals, indices)
                // negativly shifted cylinder
                Vec3.sub(tmpCylinderStart, start, shift)
                setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d)
                add(tmpCylinderMat, vertices, normals, indices)
            },
            addFixedCountDashedCylinder: (start: Vec3, end: Vec3, lengthScale: number, segmentCount: number, props: CylinderProps) => {
                const s = Math.floor(segmentCount / 2)
                const step = 1 / segmentCount

                // automatically adjust length so links/bonds that are rendered as two half cylinders
                // have evenly spaced dashed cylinders
                if (lengthScale < 1) {
                    const bias = lengthScale / 2 / segmentCount
                    lengthScale += segmentCount % 2 === 1 ? bias : -bias
                }

                const d = Vec3.distance(start, end) * lengthScale
                props.height = d * step
                const { vertices, normals, indices } = getCylinder(props)
                Vec3.sub(tmpCylinderDir, end, start)

                for (let j = 0; j < s; ++j) {
                    const f = step * (j * 2 + 1)
                    Vec3.setMagnitude(tmpCylinderDir, tmpCylinderDir, d * f)
                    Vec3.add(tmpCylinderStart, start, tmpCylinderDir)
                    setCylinderMat(tmpCylinderMat, tmpCylinderStart, tmpCylinderDir, d * step)
                    add(tmpCylinderMat, vertices, normals, indices)
                }
            },
            addIcosahedron: (center: Vec3, radius: number, detail: number) => {
                const { vertices, normals, indices } = getIcosahedron({ radius, detail })
                setIcosahedronMat(tmpIcosahedronMat, center)
                add(tmpIcosahedronMat, vertices, normals, indices)
            },
            setId: (id: number) => {
                if (currentId !== id) {
                    currentId = id
                    ChunkedArray.add(offsets, vertices.elementCount)
                }
            },
            getMesh: () => {
                ChunkedArray.add(offsets, vertices.elementCount)
                const vb = ChunkedArray.compact(vertices, true) as Float32Array
                const ib = ChunkedArray.compact(indices, true) as Uint32Array
                const nb = ChunkedArray.compact(normals, true) as Float32Array
                const idb = ChunkedArray.compact(ids, true) as Float32Array
                const ob = ChunkedArray.compact(offsets, true) as Uint32Array
                return {
                    vertexCount: vertices.elementCount,
                    triangleCount: indices.elementCount,
                    offsetCount: offsets.elementCount,
                    vertexBuffer: mesh ? ValueCell.update(mesh.vertexBuffer, vb) : ValueCell.create(vb),
                    indexBuffer: mesh ? ValueCell.update(mesh.indexBuffer, ib) : ValueCell.create(ib),
                    normalBuffer: mesh ? ValueCell.update(mesh.normalBuffer, nb) : ValueCell.create(nb),
                    idBuffer: mesh ? ValueCell.update(mesh.idBuffer, idb) : ValueCell.create(idb),
                    offsetBuffer: mesh ? ValueCell.update(mesh.offsetBuffer, ob) : ValueCell.create(ob),
                    normalsComputed: true,
                    offsetsComputed: true,
                }
            }
        }
    }
}