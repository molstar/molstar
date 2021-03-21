/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import produce from 'immer';
import { Interval } from '../../mol-data/int/interval';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { PickingId } from '../../mol-geo/geometry/picking';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Sphere3D } from '../../mol-math/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { DataLoci, EmptyLoci, Loci } from '../../mol-model/loci';
import { Shape } from '../../mol-model/shape';
import { Visual } from '../../mol-repr/visual';
import { ColorNames } from '../../mol-util/color/names';
import { MarkerAction, MarkerActions } from '../../mol-util/marker-action';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera, ICamera } from '../camera';
import { Viewport } from '../camera/util';

// TODO add scale line/grid

const AxesParams = {
    ...Mesh.Params,
    alpha: { ...Mesh.Params.alpha, defaultValue: 0.51 },
    ignoreLight: { ...Mesh.Params.ignoreLight, defaultValue: true },
    colorX: PD.Color(ColorNames.red, { isEssential: true }),
    colorY: PD.Color(ColorNames.green, { isEssential: true }),
    colorZ: PD.Color(ColorNames.blue, { isEssential: true }),
    scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.1 }, { isEssential: true }),
};
type AxesParams = typeof AxesParams
type AxesProps = PD.Values<AxesParams>

export const CameraHelperParams = {
    axes: PD.MappedStatic('on', {
        on: PD.Group(AxesParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Show camera orientation axes' }),
};
export type CameraHelperParams = typeof CameraHelperParams
export type CameraHelperProps = PD.Values<CameraHelperParams>

export class CameraHelper {
    scene: Scene
    camera: Camera
    props: CameraHelperProps = {
        axes: { name: 'off', params: {} }
    }

    private renderObject: GraphicsRenderObject | undefined

    constructor(private webgl: WebGLContext, props: Partial<CameraHelperProps> = {}) {
        this.scene = Scene.create(webgl);

        this.camera = new Camera();
        Vec3.set(this.camera.up, 0, 1, 0);
        Vec3.set(this.camera.target, 0, 0, 0);

        this.setProps(props);
    }

    setProps(props: Partial<CameraHelperProps>) {
        this.props = produce(this.props, p => {
            if (props.axes !== undefined) {
                p.axes.name = props.axes.name;
                if (props.axes.name === 'on') {
                    this.scene.clear();
                    const params = { ...props.axes.params, scale: props.axes.params.scale * this.webgl.pixelRatio };
                    this.renderObject = createAxesRenderObject(params);
                    this.renderObject.state.noClip = true;
                    this.scene.add(this.renderObject);
                    this.scene.commit();

                    Vec3.set(this.camera.position, 0, 0, params.scale * 200);
                    Mat4.lookAt(this.camera.view, this.camera.position, this.camera.target, this.camera.up);

                    p.axes.params = { ...props.axes.params };
                }
            }
        });
    }

    get isEnabled() {
        return this.props.axes.name === 'on';
    }

    getLoci(pickingId: PickingId) {
        const { objectId, groupId, instanceId } = pickingId;
        if (!this.renderObject || objectId !== this.renderObject.id || groupId === CameraHelperAxis.None) return EmptyLoci;
        return CameraAxesLoci(this, groupId, instanceId);
    }

    private eachGroup = (loci: Loci, apply: (interval: Interval) => boolean): boolean => {
        if (!this.renderObject) return false;
        if (!isCameraAxesLoci(loci)) return false;
        let changed = false;
        const groupCount = this.renderObject.values.uGroupCount.ref.value;
        const { elements } = loci;
        for (const { groupId, instanceId } of elements) {
            const idx = instanceId * groupCount + groupId;
            if (apply(Interval.ofSingleton(idx))) changed = true;
        }
        return changed;
    }

    mark(loci: Loci, action: MarkerAction) {
        if (!MarkerActions.is(MarkerActions.Highlighting, action)) return false;
        if (!isCameraAxesLoci(loci)) return false;
        if (loci.data !== this) return false;
        return Visual.mark(this.renderObject, loci, action, this.eachGroup);
    }

    update(camera: ICamera) {
        if (!this.renderObject) return;

        updateCamera(this.camera, camera.viewport, camera.viewOffset);
        Mat4.extractRotation(this.scene.view, camera.view);

        const r = this.renderObject.values.boundingSphere.ref.value.radius;
        Mat4.setTranslation(this.scene.view, Vec3.create(
            -camera.viewport.width / 2 + r,
            -camera.viewport.height / 2 + r,
            0
        ));
    }
}

export const enum CameraHelperAxis {
    None = 0,
    X,
    Y,
    Z,
    XY,
    XZ,
    YZ
}

function getAxisLabel(axis: number) {
    switch (axis) {
        case CameraHelperAxis.X: return 'X Axis';
        case CameraHelperAxis.Y: return 'Y Axis';
        case CameraHelperAxis.Z: return 'Z Axis';
        case CameraHelperAxis.XY: return 'XY Plane';
        case CameraHelperAxis.XZ: return 'XZ Plane';
        case CameraHelperAxis.YZ: return 'YZ Plane';
        default: return 'Axes';
    }
}

function CameraAxesLoci(cameraHelper: CameraHelper, groupId: number, instanceId: number) {
    return DataLoci('camera-axes', cameraHelper, [{ groupId, instanceId }],
        void 0 /** bounding sphere */,
        () => getAxisLabel(groupId));
}
export type CameraAxesLoci = ReturnType<typeof CameraAxesLoci>
export function isCameraAxesLoci(x: Loci): x is CameraAxesLoci {
    return x.kind === 'data-loci' && x.tag === 'camera-axes';
}

function updateCamera(camera: Camera, viewport: Viewport, viewOffset: Camera.ViewOffset) {
    const { near, far } = camera;

    const fullLeft = -viewport.width / 2;
    const fullRight = viewport.width / 2;
    const fullTop = viewport.height / 2;
    const fullBottom = -viewport.height / 2;

    const dx = (fullRight - fullLeft) / 2;
    const dy = (fullTop - fullBottom) / 2;
    const cx = (fullRight + fullLeft) / 2;
    const cy = (fullTop + fullBottom) / 2;

    let left = cx - dx;
    let right = cx + dx;
    let top = cy + dy;
    let bottom = cy - dy;

    if (viewOffset.enabled) {
        const scaleW = (fullRight - fullLeft) / viewOffset.width;
        const scaleH = (fullTop - fullBottom) / viewOffset.height;
        left += scaleW * viewOffset.offsetX;
        right = left + scaleW * viewOffset.width;
        top -= scaleH * viewOffset.offsetY;
        bottom = top - scaleH * viewOffset.height;
    }

    Mat4.ortho(camera.projection, left, right, top, bottom, near, far);
}

function createAxesMesh(scale: number, mesh?: Mesh) {
    const state = MeshBuilder.createState(512, 256, mesh);
    const radius = 0.075 * scale;
    const x = Vec3.scale(Vec3(), Vec3.unitX, scale);
    const y = Vec3.scale(Vec3(), Vec3.unitY, scale);
    const z = Vec3.scale(Vec3(), Vec3.unitZ, scale);
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments: 32 };

    state.currentGroup = CameraHelperAxis.None;
    addSphere(state, Vec3.origin, radius, 2);

    state.currentGroup = CameraHelperAxis.X;
    addSphere(state, x, radius, 2);
    addCylinder(state, Vec3.origin, x, 1, cylinderProps);

    state.currentGroup = CameraHelperAxis.Y;
    addSphere(state, y, radius, 2);
    addCylinder(state, Vec3.origin, y, 1, cylinderProps);

    state.currentGroup = CameraHelperAxis.Z;
    addSphere(state, z, radius, 2);
    addCylinder(state, Vec3.origin, z, 1, cylinderProps);

    Vec3.scale(x, x, 0.5);
    Vec3.scale(y, y, 0.5);
    Vec3.scale(z, z, 0.5);

    state.currentGroup = CameraHelperAxis.XY;
    MeshBuilder.addTriangle(state, Vec3.origin, x, y);
    MeshBuilder.addTriangle(state, Vec3.origin, y, x);
    const xy = Vec3.add(Vec3(), x, y);
    MeshBuilder.addTriangle(state, xy, x, y);
    MeshBuilder.addTriangle(state, xy, y, x);

    state.currentGroup = CameraHelperAxis.XZ;
    MeshBuilder.addTriangle(state, Vec3.origin, x, z);
    MeshBuilder.addTriangle(state, Vec3.origin, z, x);
    const xz = Vec3.add(Vec3(), x, z);
    MeshBuilder.addTriangle(state, xz, x, z);
    MeshBuilder.addTriangle(state, xz, z, x);

    state.currentGroup = CameraHelperAxis.YZ;
    MeshBuilder.addTriangle(state, Vec3.origin, y, z);
    MeshBuilder.addTriangle(state, Vec3.origin, z, y);
    const yz = Vec3.add(Vec3(), y, z);
    MeshBuilder.addTriangle(state, yz, y, z);
    MeshBuilder.addTriangle(state, yz, z, y);

    return MeshBuilder.getMesh(state);
}

function getAxesShape(props: AxesProps, shape?: Shape<Mesh>) {
    const scale = 100 * props.scale;
    const mesh = createAxesMesh(scale, shape?.geometry);
    mesh.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale + scale / 4));
    const getColor = (groupId: number) => {
        switch (groupId) {
            case 1: return props.colorX;
            case 2: return props.colorY;
            case 3: return props.colorZ;
            default: return ColorNames.grey;
        }
    };
    return Shape.create('axes', {}, mesh, getColor, () => 1, () => '');
}

function createAxesRenderObject(props: AxesProps) {
    const shape = getAxesShape(props);
    return Shape.createRenderObject(shape, props);
}