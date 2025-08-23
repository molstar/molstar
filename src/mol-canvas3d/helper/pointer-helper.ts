/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { Scene } from '../../mol-gl/scene';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ColorNames } from '../../mol-util/color/names';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Geometry } from '../../mol-geo/geometry/geometry';
import { addCylinderFromRay3D } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { Camera, ICamera } from '../camera';
import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Viewport } from '../camera/util';
import { Shape } from '../../mol-model/shape/shape';

export const PointerHelperParams = {
    ...Mesh.Params,
    enabled: PD.Select('off', PD.arrayToOptions(['on', 'off']), { isEssential: true }),
    ignoreLight: { ...Mesh.Params.ignoreLight, defaultValue: true },
    color: PD.Color(ColorNames.grey, { isEssential: true }),
    hitColor: PD.Color(ColorNames.pink, { isEssential: true }),
};
export type PointerHelperParams = typeof PointerHelperParams
export type PointerHelperProps = PD.Values<PointerHelperParams>

export class PointerHelper {
    readonly scene: Scene;
    readonly camera: Camera;
    readonly props: PointerHelperProps;

    private renderObject: GraphicsRenderObject<'mesh'>;
    private shape: Shape<Mesh>;

    private pointers: Ray3D[] = [];
    private points: Vec3[] = [];
    private hit: Vec3 | undefined = undefined;

    setProps(props: Partial<PointerHelperProps>) {
        Object.assign(this.props, props);
        if (this.isEnabled) this.update(this.pointers, this.points, this.hit);
    }

    ensureEnabled() {
        if (this.props.enabled !== 'on') this.props.enabled = 'on';
    }

    get isEnabled() {
        return this.props.enabled === 'on';
    }

    setCamera(camera: ICamera) {
        Viewport.copy(this.camera.viewport, camera.viewport);
        Mat4.copy(this.camera.view, camera.view);
        Mat4.copy(this.camera.projection, camera.projection);
        Mat4.copy(this.camera.projectionView, camera.projectionView);
        Mat4.copy(this.camera.headRotation, camera.headRotation);
        Camera.copyViewOffset(this.camera.viewOffset, camera.viewOffset);
        this.camera.far = camera.far;
        this.camera.near = camera.near;
        this.camera.fogFar = camera.fogFar;
        this.camera.fogNear = camera.fogNear;

        this.camera.forceFull = camera.forceFull;
        this.camera.scale = 1;
    }

    update(pointers: Ray3D[], points: Vec3[], hit: Vec3 | undefined) {
        this.pointers = pointers;
        this.points = points;
        this.hit = hit;

        const p = this.props;
        if (p.enabled !== 'on') {
            if (this.renderObject) this.renderObject.state.visible = false;
            return;
        }

        this.shape = getPointerMeshShape(this.getData(), this.props, this.shape);

        ValueCell.updateIfChanged(this.renderObject.values.drawCount, Geometry.getDrawCount(this.shape.geometry));
        ValueCell.updateIfChanged(this.renderObject.values.uVertexCount, Geometry.getVertexCount(this.shape.geometry));
        ValueCell.updateIfChanged(this.renderObject.values.uGroupCount, 2);
        Mesh.Utils.updateBoundingSphere(this.renderObject.values, this.shape.geometry);
        Mesh.Utils.updateValues(this.renderObject.values, this.props);
        Mesh.Utils.updateRenderableState(this.renderObject.state, this.props);

        this.renderObject.state.visible = true;

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    private getData() {
        return {
            pointers: this.pointers,
            points: this.points,
            hit: this.hit,
        };
    }

    constructor(webgl: WebGLContext, props: Partial<PointerHelperProps> = {}) {
        this.scene = Scene.create(webgl, 'blended');
        this.props = { ...PD.getDefaultValues(PointerHelperParams), ...props };

        this.camera = new Camera();

        this.shape = getPointerMeshShape(this.getData(), this.props, this.shape);
        this.renderObject = createMeshRenderObject(this.shape, this.props);
        this.scene.add(this.renderObject);
    }
}

type PointerData = {
    pointers: Ray3D[]
    points: Vec3[]
    hit?: Vec3
}

export enum PointerHelperGroup {
    None = 0,
    Hit,
}

function createPointerMesh(data: PointerData, mesh?: Mesh) {
    const scale = 1; // 1 / 0.01;
    const state = MeshBuilder.createState(512, 256, mesh);
    const radius = 0.0005 * scale;
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments: 32 };

    state.currentGroup = PointerHelperGroup.None;
    for (const pointer of data.pointers) {
        addCylinderFromRay3D(state, pointer, 0.2 * scale, cylinderProps);
        addSphere(state, pointer.origin, 0.001 * scale, 1);
    }
    for (const point of data.points) {
        addSphere(state, point, 0.0025 * scale, 1);
    }

    if (data.hit) {
        state.currentGroup = PointerHelperGroup.Hit;
        addSphere(state, data.hit, 0.0025 * scale, 1);
    }

    return MeshBuilder.getMesh(state);
}

function getPointerMeshShape(data: PointerData, props: PointerHelperProps, shape?: Shape<Mesh>) {
    const mesh = createPointerMesh(data, shape?.geometry);
    const getColor = (groupId: number) => {
        switch (groupId) {
            case PointerHelperGroup.Hit: return props.hitColor;
            default: return props.color;
        }
    };
    return Shape.create('pointer-mesh', data, mesh, getColor, () => 1, () => '', undefined, 2);
}

function createMeshRenderObject(shape: Shape<Mesh>, props: PointerHelperProps) {
    return Shape.createRenderObject(shape, {
        ...PD.getDefaultValues(Mesh.Params),
        ...props,
        ignoreLight: props.ignoreLight,
        cellSize: 0,
    }) as GraphicsRenderObject<'mesh'>;
}
