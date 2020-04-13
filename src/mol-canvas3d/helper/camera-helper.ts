/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import Scene from '../../mol-gl/scene';
import { Camera } from '../camera';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat4 } from '../../mol-math/linear-algebra';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ColorNames } from '../../mol-util/color/names';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Viewport } from '../camera/util';
import { ValueCell } from '../../mol-util';
import { Sphere3D } from '../../mol-math/geometry';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import produce from 'immer';
import { Shape } from '../../mol-model/shape';

// TODO add scale line/grid

const AxesParams = {
    ...Mesh.Params,
    alpha: { ...Mesh.Params.alpha, defaultValue: 0.33 },
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

    update(camera: Camera) {
        if (!this.renderObject) return;

        updateCamera(this.camera, camera.viewport);

        const m = this.renderObject.values.aTransform.ref.value as unknown as Mat4;
        Mat4.extractRotation(m, camera.view);

        const r = this.renderObject.values.boundingSphere.ref.value.radius;
        Mat4.setTranslation(m, Vec3.create(
            -camera.viewport.width / 2 + r,
            -camera.viewport.height / 2 + r,
            0
        ));

        ValueCell.update(this.renderObject.values.aTransform, this.renderObject.values.aTransform.ref.value);
        this.scene.update([this.renderObject], true);
    }
}

function updateCamera(camera: Camera, viewport: Viewport) {
    const { near, far } = camera;

    const fullLeft = -(viewport.width - viewport.x) / 2;
    const fullRight = (viewport.width - viewport.x) / 2;
    const fullTop = (viewport.height - viewport.y) / 2;
    const fullBottom = -(viewport.height - viewport.y) / 2;

    const dx = (fullRight - fullLeft) / 2;
    const dy = (fullTop - fullBottom) / 2;
    const cx = (fullRight + fullLeft) / 2;
    const cy = (fullTop + fullBottom) / 2;

    const left = cx - dx;
    const right = cx + dx;
    const top = cy + dy;
    const bottom = cy - dy;

    Mat4.ortho(camera.projection, left, right, top, bottom, near, far);
}

function createAxesMesh(scale: number, mesh?: Mesh) {
    const state = MeshBuilder.createState(512, 256, mesh);
    const radius = 0.05 * scale;
    const x = Vec3.scale(Vec3(), Vec3.unitX, scale);
    const y = Vec3.scale(Vec3(), Vec3.unitY, scale);
    const z = Vec3.scale(Vec3(), Vec3.unitZ, scale);
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments: 32 };

    state.currentGroup = 0;
    addSphere(state, Vec3.origin, radius, 2);

    state.currentGroup = 1;
    addSphere(state, x, radius, 2);
    addCylinder(state, Vec3.origin, x, 1, cylinderProps);

    state.currentGroup = 2;
    addSphere(state, y, radius, 2);
    addCylinder(state, Vec3.origin, y, 1, cylinderProps);

    state.currentGroup = 3;
    addSphere(state, z, radius, 2);
    addCylinder(state, Vec3.origin, z, 1, cylinderProps);

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