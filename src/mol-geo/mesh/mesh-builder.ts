/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Vec3, Mat4, Mat3 } from 'mol-math/linear-algebra';
import { ChunkedArray } from 'mol-data/util';

import { Plane, PlaneProps } from '../primitive/plane';
import { Cylinder, CylinderProps } from '../primitive/cylinder';
import { Sphere, SphereProps } from '../primitive/sphere';
import { Mesh } from './mesh';
import { getNormalMatrix } from '../util';
import { addSheet } from './sheet';
import { addTube } from './tube';
import { StarProps, Star } from '../primitive/star';
import { Octahedron, PerforatedOctahedron } from '../primitive/octahedron';
import { Primitive } from '../primitive/primitive';
import { DiamondPrism, PentagonalPrism, HexagonalPrism } from '../primitive/prism';
import { OctagonalPyramide, PerforatedOctagonalPyramid } from '../primitive/pyramid';
import { PerforatedBox, Box } from '../primitive/box';
import { Wedge } from '../primitive/wedge';

export interface MeshBuilderState {
    vertices: ChunkedArray<number, 3>
    normals: ChunkedArray<number, 3>
    indices: ChunkedArray<number, 3>
}

export interface MeshBuilder {
    add(t: Mat4, _vertices: ArrayLike<number>, _normals: ArrayLike<number>, _indices?: ArrayLike<number>): void
    addBox(t: Mat4): void
    addPerforatedBox(t: Mat4): void
    addPlane(t: Mat4, props?: PlaneProps): void
    addWedge(t: Mat4): void
    addDiamondPrism(t: Mat4): void
    addPentagonalPrism(t: Mat4): void
    addHexagonalPrism(t: Mat4): void
    addOctagonalPyramid(t: Mat4): void
    addPerforatedOctagonalPyramid(t: Mat4): void
    addStar(t: Mat4, props?: StarProps): void
    addOctahedron(t: Mat4): void
    addPerforatedOctahedron(t: Mat4): void
    addCylinder(start: Vec3, end: Vec3, lengthScale: number, props: CylinderProps): void
    addDoubleCylinder(start: Vec3, end: Vec3, lengthScale: number, shift: Vec3, props: CylinderProps): void
    addFixedCountDashedCylinder(start: Vec3, end: Vec3, lengthScale: number, segmentCount: number, props: CylinderProps): void
    addSphere(center: Vec3, radius: number, detail: number): void
    addTube(centers: ArrayLike<number>, normals: ArrayLike<number>, binormals: ArrayLike<number>, linearSegments: number, radialSegments: number, width: number, height: number, waveFactor: number, startCap: boolean, endCap: boolean): void
    addSheet(centers: ArrayLike<number>, normals: ArrayLike<number>, binormals: ArrayLike<number>, linearSegments: number, width: number, height: number, arrowHeight: number, startCap: boolean, endCap: boolean): void
    setGroup(id: number): void
    getMesh(): Mesh
}

const cylinderMap = new Map<string, Primitive>()
const sphereMap = new Map<string, Primitive>()

const up = Vec3.create(0, 1, 0)
const tmpV = Vec3.zero()
const tmpMat3 = Mat3.zero()

const tmpCylinderDir = Vec3.zero()
const tmpCylinderMatDir = Vec3.zero()
const tmpCylinderCenter = Vec3.zero()
const tmpCylinderMat = Mat4.zero()
// const tmpCylinderMatTrans = Mat4.zero()
const tmpCylinderStart = Vec3.zero()
const tmpUp = Vec3.zero()

function setCylinderMat(m: Mat4, start: Vec3, dir: Vec3, length: number) {
    Vec3.setMagnitude(tmpCylinderMatDir, dir, length / 2)
    Vec3.add(tmpCylinderCenter, start, tmpCylinderMatDir)
    // ensure the direction used to create the rotation is always pointing in the same
    // direction so the triangles of adjacent cylinder will line up
    Vec3.copy(tmpUp, up)
    if (Vec3.dot(tmpCylinderMatDir, tmpUp) < 0) Vec3.scale(tmpUp, tmpUp, -1)
    Vec3.makeRotation(m, tmpUp, tmpCylinderMatDir)
    return Mat4.setTranslation(m, tmpCylinderCenter)
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

const tmpSphereMat = Mat4.identity()

function setSphereMat(m: Mat4, center: Vec3) {
    return Mat4.setTranslation(m, center)
}

function getSphere(props: SphereProps) {
    const key = JSON.stringify(props)
    let sphere = sphereMap.get(key)
    if (sphere === undefined) {
        sphere = Sphere(props)
        sphereMap.set(key, sphere)
    }
    return sphere
}

// TODO cache primitives based on props

export namespace MeshBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, mesh?: Mesh): MeshBuilder {
        const vertices = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.vertexBuffer.ref.value : initialCount);
        const normals = ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.normalBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize * 3, mesh ? mesh.indexBuffer.ref.value : initialCount * 3);
        const state: MeshBuilderState = { vertices, normals, indices };

        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, mesh ? mesh.groupBuffer.ref.value : initialCount);

        let currentGroup = -1

        function add(t: Mat4, _vertices: ArrayLike<number>, _normals: ArrayLike<number>, _indices: ArrayLike<number>) {
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
                // group
                ChunkedArray.add(groups, currentGroup);
            }
            for (let i = 0, il = _indices.length; i < il; i += 3) {
                ChunkedArray.add3(indices, _indices[i] + elementCount, _indices[i + 1] + elementCount, _indices[i + 2] + elementCount);
            }
        }

        return {
            add,
            addBox: (t: Mat4) => {
                const { vertices, normals, indices } = Box()
                add(t, vertices, normals, indices)
            },
            addPerforatedBox: (t: Mat4) => {
                const { vertices, normals, indices } = PerforatedBox()
                add(t, vertices, normals, indices)
            },
            addPlane: (t: Mat4, props?: PlaneProps) => {
                const { vertices, normals, indices } = Plane(props)
                add(t, vertices, normals, indices)
            },
            addWedge: (t: Mat4) => {
                const { vertices, normals, indices } = Wedge()
                add(t, vertices, normals, indices)
            },
            addDiamondPrism: (t: Mat4) => {
                const { vertices, normals, indices } = DiamondPrism()
                add(t, vertices, normals, indices)
            },
            addPentagonalPrism: (t: Mat4) => {
                const { vertices, normals, indices } = PentagonalPrism()
                add(t, vertices, normals, indices)
            },
            addHexagonalPrism: (t: Mat4) => {
                const { vertices, normals, indices } = HexagonalPrism()
                add(t, vertices, normals, indices)
            },
            addOctagonalPyramid: (t: Mat4) => {
                const { vertices, normals, indices } = OctagonalPyramide()
                add(t, vertices, normals, indices)
            },
            addPerforatedOctagonalPyramid: (t: Mat4) => {
                const { vertices, normals, indices } = PerforatedOctagonalPyramid()
                add(t, vertices, normals, indices)
            },
            addStar: (t: Mat4, props?: StarProps) => {
                const { vertices, normals, indices } = Star(props)
                add(t, vertices, normals, indices)
            },
            addOctahedron: (t: Mat4) => {
                const { vertices, normals, indices } = Octahedron()
                add(t, vertices, normals, indices)
            },
            addPerforatedOctahedron: (t: Mat4) => {
                const { vertices, normals, indices } = PerforatedOctahedron()
                add(t, vertices, normals, indices)
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
            addSphere: (center: Vec3, radius: number, detail: number) => {
                const { vertices, normals, indices } = getSphere({ radius, detail })
                setSphereMat(tmpSphereMat, center)
                add(tmpSphereMat, vertices, normals, indices)
            },
            addTube: (centers: ArrayLike<number>, normals: ArrayLike<number>, binormals: ArrayLike<number>, linearSegments: number, radialSegments: number, width: number, height: number, waveFactor: number, startCap: boolean, endCap: boolean) => {
                const addedVertexCount = addTube(centers, normals, binormals, linearSegments, radialSegments, width, height, waveFactor, startCap, endCap, state)
                for (let i = 0, il = addedVertexCount; i < il; ++i) ChunkedArray.add(groups, currentGroup);
            },
            addSheet: (controls: ArrayLike<number>, normals: ArrayLike<number>, binormals: ArrayLike<number>, linearSegments: number, width: number, height: number, arrowHeight: number, startCap: boolean, endCap: boolean) => {
                const addedVertexCount = addSheet(controls, normals, binormals, linearSegments, width, height, arrowHeight, startCap, endCap, state)
                for (let i = 0, il = addedVertexCount; i < il; ++i) ChunkedArray.add(groups, currentGroup);
            },
            setGroup: (group: number) => {
                currentGroup = group
            },
            getMesh: () => {
                const vb = ChunkedArray.compact(vertices, true) as Float32Array
                const ib = ChunkedArray.compact(indices, true) as Uint32Array
                const nb = ChunkedArray.compact(normals, true) as Float32Array
                const idb = ChunkedArray.compact(groups, true) as Float32Array
                return {
                    vertexCount: vertices.elementCount,
                    triangleCount: indices.elementCount,
                    vertexBuffer: mesh ? ValueCell.update(mesh.vertexBuffer, vb) : ValueCell.create(vb),
                    indexBuffer: mesh ? ValueCell.update(mesh.indexBuffer, ib) : ValueCell.create(ib),
                    normalBuffer: mesh ? ValueCell.update(mesh.normalBuffer, nb) : ValueCell.create(nb),
                    groupBuffer: mesh ? ValueCell.update(mesh.groupBuffer, idb) : ValueCell.create(idb),
                    normalsComputed: true,
                }
            }
        }
    }
}