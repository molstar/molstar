/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject, GraphicsRenderObject, getNextMaterialId } from '../../mol-gl/render-object';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ColorNames } from '../../mol-util/color/names';
import { Color } from '../../mol-util/color';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { Quat } from '../../mol-math/linear-algebra/3d/quat';
import { ImageValues } from '../../mol-gl/renderable/image';
import { Clip } from '../../mol-util/clip';
import { addSphere as addLinesSphere } from '../../mol-geo/geometry/lines/builder/sphere';
import { addBox } from '../../mol-geo/geometry/lines/builder/box';
import { addPlane } from '../../mol-geo/geometry/lines/builder/plane';

export const ImageHelperParams = {
    imageEdges: PD.Boolean(false, { description: 'Show edges of visible image render objects.' }),
};
export type ImageHelperParams = typeof ImageHelperParams;
export type ImageHelperProps = PD.Values<ImageHelperParams>;

const imageEdgeMaterialId = getNextMaterialId();
const imageTrimMaterialId = getNextMaterialId();

// Temp vectors
const _trimPos = Vec3();
const _trimScale = Vec3();
const _trimRot = Quat();
const _trimTransform = Mat4();
const _tmpMat = Mat4();

export class ImageHelper {
    readonly scene: Scene;

    private readonly parent: Scene;
    private _props: ImageHelperProps;
    private renderObjects = new Map<number, { roList: GraphicsRenderObject[], version: number }>();

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<ImageHelperProps>) {
        this.scene = Scene.create(ctx, 'blended');
        this.parent = parent;
        this._props = { ...PD.getDefaultValues(ImageHelperParams), ...props };
    }

    update() {
        const previousIds = new Set(this.renderObjects.keys());

        this.parent.forEach((r, ro) => {
            if (!ro.state.visible) return;
            if (ro.type !== 'image') return;

            const values = ro.values as ImageValues;
            const version = values.aPosition.ref.version
                + values.uTrimType.ref.version + values.uTrimCenter.ref.version
                + values.uTrimRotation.ref.version + values.uTrimScale.ref.version
                + values.uTrimTransform.ref.version + values.aTransform.ref.version;

            const existing = this.renderObjects.get(ro.id);
            if (existing && existing.version === version) {
                previousIds.delete(ro.id);
                return;
            }

            // Remove old entries if version changed
            if (existing) {
                for (const oldRO of existing.roList) this.scene.remove(oldRO);
                this.renderObjects.delete(ro.id);
            }

            const roList: GraphicsRenderObject[] = [];

            const edgeLines = createImageEdgeLines(values);
            if (edgeLines) {
                const edgeRO = createLinesRenderObject(edgeLines, imageEdgeMaterialId, ColorNames.cyan, 0.8);
                this.scene.add(edgeRO);
                roList.push(edgeRO);
            }

            const trimLines = createTrimEdgeLines(values);
            if (trimLines) {
                const trimRO = createLinesRenderObject(trimLines, imageTrimMaterialId, ColorNames.yellow, 0.7);
                this.scene.add(trimRO);
                roList.push(trimRO);
            }

            if (roList.length > 0) {
                this.renderObjects.set(ro.id, { roList, version });
            }
            previousIds.delete(ro.id);
        });

        for (const id of previousIds) {
            const entry = this.renderObjects.get(id);
            if (entry) {
                for (const ro of entry.roList) this.scene.remove(ro);
                this.renderObjects.delete(id);
            }
        }

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    syncVisibility() {
        const visible = this._props.imageEdges;
        this.renderObjects.forEach(entry => {
            for (const ro of entry.roList) ro.state.visible = visible;
        });
    }

    clear() {
        this.renderObjects.clear();
        this.scene.clear();
    }

    get isEnabled() {
        return this._props.imageEdges;
    }

    get props() { return this._props as Readonly<ImageHelperProps>; }

    setProps(props: Partial<ImageHelperProps>) {
        Object.assign(this._props, props);
        if (this.isEnabled) this.update();
    }
}

//

/**
 * Image quad vertex layout (from image.ts):
 *   Vertex 0: UV (0,1) — top-left
 *   Vertex 1: UV (0,0) — bottom-left
 *   Vertex 2: UV (1,1) — top-right
 *   Vertex 3: UV (1,0) — bottom-right
 *
 * addPlane expects corners in winding order (0→1→2→3→0),
 * so we reorder to: top-left, bottom-left, bottom-right, top-right.
 */
const _planeCorners = new Float32Array(12);

function createImageEdgeLines(values: ImageValues): Lines | undefined {
    const positions = values.aPosition.ref.value;
    const transforms = values.aTransform.ref.value;
    const instanceCount = values.uInstanceCount.ref.value;

    if (positions.length < 12) return undefined; // need 4 vertices × 3 components

    // Reorder from [TL, BL, TR, BR] to winding order [TL, BL, BR, TR]
    // V0 (TL) → slot 0
    _planeCorners[0] = positions[0]; _planeCorners[1] = positions[1]; _planeCorners[2] = positions[2];
    // V1 (BL) → slot 1
    _planeCorners[3] = positions[3]; _planeCorners[4] = positions[4]; _planeCorners[5] = positions[5];
    // V3 (BR) → slot 2
    _planeCorners[6] = positions[9]; _planeCorners[7] = positions[10]; _planeCorners[8] = positions[11];
    // V2 (TR) → slot 3
    _planeCorners[9] = positions[6]; _planeCorners[10] = positions[7]; _planeCorners[11] = positions[8];

    const builder = LinesBuilder.create(4 * instanceCount);

    for (let inst = 0; inst < instanceCount; ++inst) {
        const transform = Mat4();
        Mat4.fromArray(transform, transforms, inst * 16);
        addPlane(builder, _planeCorners, transform, 0);
    }

    return builder.getLines();
}

function createTrimEdgeLines(values: ImageValues): Lines | undefined {
    const trimType = values.uTrimType.ref.value as number;
    if (trimType === 0) return undefined; // no trim

    const transforms = values.aTransform.ref.value;
    const instanceCount = values.uInstanceCount.ref.value;

    Vec3.copy(_trimPos, values.uTrimCenter.ref.value);
    Quat.copy(_trimRot, values.uTrimRotation.ref.value);
    Vec3.copy(_trimScale, values.uTrimScale.ref.value);
    Mat4.copy(_trimTransform, values.uTrimTransform.ref.value);

    if (trimType === Clip.Type.cube) {
        return createCubeTrimLines(transforms, instanceCount);
    } else if (trimType === Clip.Type.sphere) {
        return createSphereTrimLines(transforms, instanceCount);
    }

    // For other trim types (plane, cylinder, cone), draw a cube outline as a fallback
    // using the trim center/scale/rotation
    return createCubeTrimLines(transforms, instanceCount);
}

function createCubeTrimLines(transforms: Float32Array, instanceCount: number): Lines | undefined {
    // Build cube transform: translate * rotate * scale
    const rotMat = Mat4.fromQuat(Mat4(), _trimRot);
    const translateMat = Mat4.fromTranslation(Mat4(), _trimPos);
    const scaleMat = Mat4.fromScaling(Mat4(), _trimScale);
    Mat4.mul(_tmpMat, translateMat, rotMat);
    Mat4.mul(_tmpMat, _tmpMat, scaleMat);

    // Apply inverse of trim transform
    if (!Mat4.isIdentity(_trimTransform)) {
        const invTrimTransform = Mat4.invert(Mat4(), _trimTransform);
        Mat4.mul(_tmpMat, invTrimTransform, _tmpMat);
    }

    // addBox uses [0,1]^3, trim cube uses [-0.5,0.5]^3 — prepend offset
    const offset = Mat4.fromTranslation(Mat4(), Vec3.create(-0.5, -0.5, -0.5));
    Mat4.mul(_tmpMat, _tmpMat, offset);

    const builder = LinesBuilder.create(12 * instanceCount);

    for (let inst = 0; inst < instanceCount; ++inst) {
        const instTransform = Mat4();
        Mat4.fromArray(instTransform, transforms, inst * 16);

        const combined = Mat4.mul(Mat4(), instTransform, _tmpMat);
        addBox(builder, combined, 0);
    }

    return builder.getLines();
}

function createSphereTrimLines(transforms: Float32Array, instanceCount: number): Lines | undefined {
    const radius = Math.max(_trimScale[0] * 0.5, _trimScale[1] * 0.5, _trimScale[2] * 0.5);

    const rotMat = Mat4.fromQuat(Mat4(), _trimRot);
    const translateMat = Mat4.fromTranslation(Mat4(), _trimPos);
    Mat4.mul(_tmpMat, translateMat, rotMat);

    if (!Mat4.isIdentity(_trimTransform)) {
        const invTrimTransform = Mat4.invert(Mat4(), _trimTransform);
        Mat4.mul(_tmpMat, invTrimTransform, _tmpMat);
    }

    const segments = 32;
    const circlesPerDimension = 3;
    const builder = LinesBuilder.create(segments * 3 * circlesPerDimension * instanceCount);

    for (let inst = 0; inst < instanceCount; ++inst) {
        const instTransform = Mat4();
        Mat4.fromArray(instTransform, transforms, inst * 16);
        const combined = Mat4.mul(Mat4(), instTransform, _tmpMat);

        addLinesSphere(builder, radius, combined, 0, { segments, circlesPerDimension });
    }

    return builder.getLines();
}

function createLinesRenderObject(lines: Lines, materialId: number, color: Color, alpha: number): GraphicsRenderObject {
    const props = { ...PD.getDefaultValues(Lines.Params), sizeFactor: 1, alpha };
    const values = Lines.Utils.createValuesSimple(lines, props, color, 1);
    const state = Lines.Utils.createRenderableState(props);
    state.pickable = false;
    return createRenderObject('lines', values, state, materialId);
}
