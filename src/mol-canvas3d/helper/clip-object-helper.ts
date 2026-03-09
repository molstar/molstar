/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject, GraphicsRenderObject, getNextMaterialId } from '../../mol-gl/render-object';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { Clip } from '../../mol-util/clip';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Quat } from '../../mol-math/linear-algebra/3d/quat';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Box } from '../../mol-geo/primitive/box';
import { Plane } from '../../mol-geo/primitive/plane';
import { Cylinder } from '../../mol-geo/primitive/cylinder';
import { Sphere } from '../../mol-geo/primitive/sphere';

export const ClipObjectHelperParams = {
    clipObjects: PD.Boolean(false, { description: 'Show clip-objects of visible render objects.' }),
};
export type ClipObjectHelperParams = typeof ClipObjectHelperParams;
export type ClipObjectHelperProps = PD.Values<ClipObjectHelperParams>;

//

/** Serializes clip object params to a string key for deduplication */
function clipObjectKey(type: number, invert: boolean, position: ArrayLike<number>, posOffset: number, rotation: ArrayLike<number>, rotOffset: number, scale: ArrayLike<number>, scaleOffset: number, transform: ArrayLike<number>, transformOffset: number): string {
    // Round floats to 5 decimal places to avoid floating point noise
    const r = (v: number) => Math.round(v * 100000) / 100000;
    const parts = [
        type, invert ? 1 : 0,
        r(position[posOffset]), r(position[posOffset + 1]), r(position[posOffset + 2]),
        r(rotation[rotOffset]), r(rotation[rotOffset + 1]), r(rotation[rotOffset + 2]), r(rotation[rotOffset + 3]),
        r(scale[scaleOffset]), r(scale[scaleOffset + 1]), r(scale[scaleOffset + 2]),
    ];
    for (let j = 0; j < 16; ++j) {
        parts.push(r(transform[transformOffset + j]));
    }
    return parts.join(',');
}

type ClipObjectData = {
    key: string,
    renderObject: GraphicsRenderObject,
    indicatorRenderObject: GraphicsRenderObject,
    mesh: Mesh,
}

const clipObjectColors: Record<number, Color> = {
    [Clip.Type.plane]: ColorNames.orange,
    [Clip.Type.sphere]: ColorNames.green,
    [Clip.Type.cube]: ColorNames.dodgerblue,
    [Clip.Type.cylinder]: ColorNames.gold,
    [Clip.Type.infiniteCone]: ColorNames.crimson,
};

const clipMaterialId = getNextMaterialId();
const indicatorMaterialId = getNextMaterialId();

// Pre-rotation matrices for aligning primitives to GLSL SDF local frames
// Plane: Rx(-90°) maps primitive Z-normal to GLSL Y-normal
const preRotPlaneQuat = Quat.setAxisAngle(Quat(), Vec3.create(1, 0, 0), -Math.PI / 2);
const preRotPlaneMat = Mat4.fromQuat(Mat4(), preRotPlaneQuat);
// Cone: Rx(+90°) maps primitive Y-axis to GLSL Z-axis
const preRotConeQuat = Quat.setAxisAngle(Quat(), Vec3.create(1, 0, 0), Math.PI / 2);
const preRotConeMat = Mat4.fromQuat(Mat4(), preRotConeQuat);

// Temp variables for constructing transforms
const _position = Vec3();
const _rotation = Quat();
const _scale = Vec3();
const _clipTransform = Mat4();
const _invClipTransform = Mat4();
const _rotMat = Mat4();
const _translateMat = Mat4();
const _baseMat = Mat4();
const _tmpMat = Mat4();
const _axisEnd = Vec3();
const _yAxis = Vec3.create(0, 1, 0);
const _zAxis = Vec3.create(0, 0, 1);
const _indicatorPos = Vec3();

export class ClipObjectHelper {
    readonly scene: Scene;

    private readonly parent: Scene;
    private _props: ClipObjectHelperProps;
    private objectsData = new Map<string, ClipObjectData>();

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<ClipObjectHelperProps>) {
        this.scene = Scene.create(ctx, 'blended');
        this.parent = parent;
        this._props = { ...PD.getDefaultValues(ClipObjectHelperParams), ...props };
    }

    update() {
        const currentKeys = new Set<string>();
        const sceneRadius = this.parent.boundingSphereVisible.radius || 50;

        this.parent.forEach((r, ro) => {
            if (!ro.state.visible) return;

            const count = ro.values.dClipObjectCount.ref.value;
            if (count === 0) return;

            const types = ro.values.uClipObjectType.ref.value;
            const inverts = ro.values.uClipObjectInvert.ref.value;
            const positions = ro.values.uClipObjectPosition.ref.value;
            const rotations = ro.values.uClipObjectRotation.ref.value;
            const scales = ro.values.uClipObjectScale.ref.value;
            const transforms = ro.values.uClipObjectTransform.ref.value;

            for (let i = 0; i < count; ++i) {
                const type = types[i];
                if (type === Clip.Type.none) continue;

                const key = clipObjectKey(
                    type, inverts[i],
                    positions, i * 3,
                    rotations, i * 4,
                    scales, i * 3,
                    transforms, i * 16
                );

                currentKeys.add(key);

                if (this.objectsData.has(key)) continue;

                // Extract per-object params
                Vec3.fromArray(_position, positions, i * 3);
                Quat.fromArray(_rotation, rotations, i * 4);
                Quat.normalize(_rotation, _rotation); // ensure unit quaternion for proper rotation
                Vec3.fromArray(_scale, scales, i * 3);
                Mat4.fromArray(_clipTransform, transforms, i * 16);

                // Build base transform (translate * rotate) without scale,
                // so each shape can insert pre-rotations before scale.
                Mat4.fromQuat(_rotMat, _rotation);
                Mat4.fromTranslation(_translateMat, _position);
                Mat4.mul(_baseMat, _translateMat, _rotMat);

                // apply inverse of clip transform
                if (!Mat4.isIdentity(_clipTransform)) {
                    Mat4.invert(_invClipTransform, _clipTransform);
                    Mat4.mul(_baseMat, _invClipTransform, _baseMat);
                }

                const mesh = createClipObjectMesh(type, _baseMat, _scale, sceneRadius);
                const color = clipObjectColors[type] || ColorNames.white;
                const renderObject = createClipObjectRenderObject(mesh, color, clipMaterialId, type);

                // Create position/rotation indicator mesh
                const invert = inverts[i];
                const indicatorMesh = createIndicatorMesh(_position, _rotation, _clipTransform, _scale, type, invert);
                const indicatorRenderObject = createIndicatorRenderObject(indicatorMesh, indicatorMaterialId);

                this.scene.add(renderObject);
                this.scene.add(indicatorRenderObject);
                this.objectsData.set(key, { key, renderObject, indicatorRenderObject, mesh });
            }
        });

        // Remove clip objects no longer present
        this.objectsData.forEach((data, key) => {
            if (!currentKeys.has(key)) {
                this.scene.remove(data.renderObject);
                this.scene.remove(data.indicatorRenderObject);
                this.objectsData.delete(key);
            }
        });

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    syncVisibility() {
        const visible = this._props.clipObjects;
        this.objectsData.forEach(data => {
            data.renderObject.state.visible = visible;
            data.indicatorRenderObject.state.visible = visible;
        });
    }

    clear() {
        this.objectsData.clear();
        this.scene.clear();
    }

    get isEnabled() {
        return this._props.clipObjects;
    }

    get props() { return this._props as Readonly<ClipObjectHelperProps>; }

    setProps(props: Partial<ClipObjectHelperProps>) {
        Object.assign(this._props, props);
        if (this.isEnabled) this.update();
    }
}

//

function createClipObjectMesh(type: number, baseMat: Mat4, scale: Vec3, sceneRadius: number): Mesh {
    switch (type) {
        case Clip.Type.plane: return createPlaneMesh(baseMat, sceneRadius);
        case Clip.Type.sphere: return createSphereMesh(baseMat, scale);
        case Clip.Type.cube: return createCubeMesh(baseMat, scale);
        case Clip.Type.cylinder: return createCylinderMesh(baseMat, scale);
        case Clip.Type.infiniteCone: return createConeMesh(baseMat, scale, sceneRadius);
        default: return createSphereMesh(baseMat, scale); // fallback
    }
}

/**
 * Plane: GLSL normal is quaternionTransform(rotation, vec3(0,1,0)) — Y-up default.
 * Plane() primitive lies in XY with normal (0,0,1) along Z.
 * Pre-rotate Rx(-90°) to align primitive Z-normal to GLSL Y-normal.
 * Sized to cover the scene bounding sphere. Clip scale is ignored (plane is infinite in GLSL).
 */
function createPlaneMesh(baseMat: Mat4, sceneRadius: number): Mesh {
    const size = Math.max(sceneRadius * 2, 10);
    // baseMat * preRotPlane * uniformScale(size)
    Mat4.mul(_tmpMat, baseMat, preRotPlaneMat);
    Mat4.scale(_tmpMat, _tmpMat, Vec3.create(size, size, 1));

    const plane = Plane();
    const builderState = MeshBuilder.createState(256, 128);
    MeshBuilder.addPrimitive(builderState, _tmpMat, plane);
    // Add flipped backface for double-sided visibility
    MeshBuilder.addPrimitiveFlipped(builderState, _tmpMat, plane);
    return MeshBuilder.getMesh(builderState);
}

/**
 * Sphere: SDF uses scale * 0.5 as the radii (ellipsoid).
 * Sphere primitive has radius 1.
 * Transform: baseMat * scale * 0.5
 */
function createSphereMesh(baseMat: Mat4, scale: Vec3): Mesh {
    const detail = 2;
    const sphere = getSphereForHelper(detail);
    // baseMat * scale(scale * 0.5)
    Mat4.scale(_tmpMat, baseMat, Vec3.create(scale[0] * 0.5, scale[1] * 0.5, scale[2] * 0.5));

    const vertexCount = 10 * Math.pow(2, 2 * detail) + 2;
    const builderState = MeshBuilder.createState(vertexCount * 3, vertexCount);
    MeshBuilder.addPrimitive(builderState, _tmpMat, sphere);
    return MeshBuilder.getMesh(builderState);
}

let _helperSphere: ReturnType<typeof Sphere> | undefined;
function getSphereForHelper(detail: number) {
    if (!_helperSphere) _helperSphere = Sphere(detail);
    return _helperSphere;
}

/**
 * Cube: SDF uses scale * 0.5 as half-extents.
 * Box() primitive is ±0.5 (unit cube), so scaling by `scale` gives half-extents of scale*0.5.
 */
function createCubeMesh(baseMat: Mat4, scale: Vec3): Mesh {
    // baseMat * scale(scale)
    Mat4.scale(_tmpMat, baseMat, scale);

    const box = Box();
    const builderState = MeshBuilder.createState(256, 128);
    MeshBuilder.addPrimitive(builderState, _tmpMat, box);
    return MeshBuilder.getMesh(builderState);
}

/**
 * Cylinder: SDF axis along Y, radius = scale.x * 0.5, half-height = scale.y * 0.5.
 * Cylinder primitive: axis along Y, radius=1 in XZ, half-height=0.5 in Y.
 * Need: X/Z *= scale.x * 0.5 (radius 1 → scale.x*0.5), Y *= scale.y (half-height 0.5 → scale.y*0.5).
 */
function createCylinderMesh(baseMat: Mat4, scale: Vec3): Mesh {
    const cyl = Cylinder({ radiusTop: 1, radiusBottom: 1, height: 1, radialSegments: 16, heightSegments: 1, topCap: true, bottomCap: true });
    // baseMat * scale(scale.x * 0.5, scale.y, scale.x * 0.5) — use scale.x for both radial axes
    Mat4.scale(_tmpMat, baseMat, Vec3.create(scale[0] * 0.5, scale[1], scale[0] * 0.5));

    const vertexCount = cyl.vertices.length / 3;
    const builderState = MeshBuilder.createState(vertexCount * 3, vertexCount);
    MeshBuilder.addPrimitive(builderState, _tmpMat, cyl);
    return MeshBuilder.getMesh(builderState);
}

/**
 * InfiniteCone: GLSL SDF axis along Z, radial in XY.
 *   surface: size.x * length(t.xy) + size.y * t.z = 0  (size = scale * 0.5)
 *   half-angle: tan(θ) = scale.y / scale.x
 *   apex at clip position (origin), opens in -Z direction.
 *
 * Cylinder primitive (radiusTop=0, radiusBottom=1, height=1):
 *   axis along Y, tip at Y=+0.5, base at Y=-0.5, base radius=1.
 *
 * Transform chain (right-to-left):
 *   1. Scale(baseRadius, coneLength, baseRadius): stretch primitive to correct proportions
 *   2. Translate(0, -0.5*coneLength, 0): move tip from Y=+0.5·cL to Y=0 (apex at origin)
 *      (after scale, tip is at Y=+0.5·cL; shifting by -0.5·cL puts it at Y=0)
 *   3. preRotCone Rx(+90°): map prim-Y→Z, so cone axis becomes Z, opening in -Z
 *   4. baseMat: position + rotation of clip object
 */
function createConeMesh(baseMat: Mat4, scale: Vec3, sceneRadius: number): Mesh {
    const cone = Cylinder({ radiusTop: 0, radiusBottom: 1, height: 1, radialSegments: 16, heightSegments: 1, topCap: false, bottomCap: true });

    // Visible length of the (infinite) cone, and base radius matching the GLSL half-angle
    const coneLength = Math.max(sceneRadius * 2, 10);
    const tanHalfAngle = (scale[1] || 1) / (scale[0] || 1); // tan(θ) = scaleY / scaleX
    const baseRadius = coneLength * tanHalfAngle;

    // baseMat * preRotCone * Translate(0, -coneLength/2, 0) * Scale(baseRadius, coneLength, baseRadius)
    const scaleMat = Mat4.fromScaling(Mat4(), Vec3.create(baseRadius, coneLength, baseRadius));
    const translateMat = Mat4.fromTranslation(Mat4(), Vec3.create(0, -coneLength * 0.5, 0));
    Mat4.mul(_tmpMat, translateMat, scaleMat);
    Mat4.mul(_tmpMat, preRotConeMat, _tmpMat);
    Mat4.mul(_tmpMat, baseMat, _tmpMat);

    const vertexCount = cone.vertices.length / 3;
    const builderState = MeshBuilder.createState(vertexCount * 3, vertexCount);
    MeshBuilder.addPrimitive(builderState, _tmpMat, cone);
    return MeshBuilder.getMesh(builderState);
}

function createClipObjectRenderObject(mesh: Mesh, color: Color, materialId: number, type: number) {
    const alpha = type === Clip.Type.plane ? 0.25 : 0.15;
    const values = Mesh.Utils.createValuesSimple(mesh, { alpha, doubleSided: false, cellSize: 0, batchSize: 0 }, color, 1);
    return createRenderObject('mesh', values, { disposed: false, visible: true, alphaFactor: 1, pickable: false, colorOnly: false, opaque: false, writeDepth: false }, materialId);
}

/**
 * Create a mesh with a sphere at the clip object position and a cylinder
 * along the characteristic axis to indicate orientation.
 *
 * - Plane/sphere/cube/cylinder: axis = rotated Y (matches GLSL normal/axis direction)
 * - InfiniteCone: axis = rotated Z (cone axis is Z in local frame)
 * - Plane with invert: direction is flipped
 */
function createIndicatorMesh(position: Vec3, rotation: Quat, clipTransform: Mat4, scale: Vec3, type: number, invert: boolean): Mesh {
    const objectSize = Math.max(scale[0], scale[1], scale[2]);
    const sphereRadius = Math.max(objectSize * 0.004, 0.01);
    const cylinderRadius = sphereRadius * 0.4;
    const axisLength = Math.max(objectSize * 0.1, 2);

    // Transform position by inverse of clipTransform if non-identity
    Vec3.copy(_indicatorPos, position);
    if (!Mat4.isIdentity(clipTransform)) {
        Mat4.invert(_invClipTransform, clipTransform);
        Vec3.transformMat4(_indicatorPos, _indicatorPos, _invClipTransform);
    }

    // Choose the local-frame axis based on clip type
    const localAxis = type === Clip.Type.infiniteCone ? _zAxis : _yAxis;
    Vec3.transformQuat(_axisEnd, localAxis, rotation);

    // Cone opens in -Z locally, so negate to point along the cone opening
    if (type === Clip.Type.infiniteCone) {
        Vec3.negate(_axisEnd, _axisEnd);
    }

    // For planes, the normal points toward the clipped (removed) side.
    // Flip so the indicator points toward the non-clipped (kept) geometry.
    // When inverted, the kept side is the normal side, so don't flip.
    if (type === Clip.Type.plane && !invert) {
        Vec3.negate(_axisEnd, _axisEnd);
    }

    // If clipTransform is non-identity, also transform the axis direction
    if (!Mat4.isIdentity(clipTransform)) {
        // Transform direction (not position) by inverse clipTransform
        const endWorld = Vec3();
        Vec3.add(endWorld, position, Vec3.scale(Vec3(), _axisEnd, axisLength));
        Vec3.transformMat4(endWorld, endWorld, _invClipTransform);
        Vec3.sub(_axisEnd, endWorld, _indicatorPos);
        Vec3.normalize(_axisEnd, _axisEnd);
    }

    // Axis cylinder endpoint
    const axisEndPoint = Vec3();
    Vec3.scaleAndAdd(axisEndPoint, _indicatorPos, _axisEnd, axisLength);

    const builderState = MeshBuilder.createState(512, 256);
    // Position sphere
    addSphere(builderState, _indicatorPos, sphereRadius, 1);
    // Rotation axis cylinder
    addCylinder(builderState, _indicatorPos, axisEndPoint, 1, { radiusTop: cylinderRadius, radiusBottom: cylinderRadius, radialSegments: 8 });
    // Small sphere at tip of axis
    addSphere(builderState, axisEndPoint, cylinderRadius * 1.5, 1);
    return MeshBuilder.getMesh(builderState);
}

function createIndicatorRenderObject(mesh: Mesh, materialId: number) {
    const values = Mesh.Utils.createValuesSimple(mesh, { alpha: 0.7, doubleSided: false, cellSize: 0, batchSize: 0 }, ColorNames.white, 1);
    return createRenderObject('mesh', values, { disposed: false, visible: true, alphaFactor: 1, pickable: false, colorOnly: false, opaque: false, writeDepth: false }, materialId);
}
