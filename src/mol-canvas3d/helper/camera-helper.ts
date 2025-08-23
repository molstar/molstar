/**
 * Copyright (c) 2020-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { produce } from '../../mol-util/produce';
import { Interval } from '../../mol-data/int/interval';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Text } from '../../mol-geo/geometry/text/text';
import { TextBuilder } from '../../mol-geo/geometry/text/text-builder';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Sphere3D } from '../../mol-math/geometry';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { DataLoci, EmptyLoci, isEveryLoci, Loci } from '../../mol-model/loci';
import { Shape } from '../../mol-model/shape';
import { Visual } from '../../mol-repr/visual';
import { ColorNames } from '../../mol-util/color/names';
import { MarkerAction, MarkerActions } from '../../mol-util/marker-action';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { assertUnreachable } from '../../mol-util/type-helpers';
import { Camera, ICamera } from '../camera';
import { Viewport } from '../camera/util';

// TODO add scale line/grid

const AxesParams = {
    alpha: PD.Numeric(0.51, { min: 0, max: 1, step: 0.01 }, { isEssential: true, label: 'Opacity' }),
    colorX: PD.Color(ColorNames.red, { isEssential: true }),
    colorY: PD.Color(ColorNames.green, { isEssential: true }),
    colorZ: PD.Color(ColorNames.blue, { isEssential: true }),
    scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.1 }, { isEssential: true }),
    location: PD.Select('bottom-left', PD.arrayToOptions(['bottom-left', 'bottom-right', 'top-left', 'top-right'] as const)),
    locationOffsetX: PD.Numeric(0),
    locationOffsetY: PD.Numeric(0),
    originColor: PD.Color(ColorNames.grey),
    radiusScale: PD.Numeric(0.075, { min: 0.01, max: 0.3, step: 0.001 }),
    showPlanes: PD.Boolean(true),
    planeColorXY: PD.Color(ColorNames.grey, { label: 'Plane Color XY' }),
    planeColorXZ: PD.Color(ColorNames.grey, { label: 'Plane Color XZ' }),
    planeColorYZ: PD.Color(ColorNames.grey, { label: 'Plane Color YZ' }),
    showLabels: PD.Boolean(false),
    labelX: PD.Text('X'),
    labelY: PD.Text('Y'),
    labelZ: PD.Text('Z'),
    labelColorX: PD.Color(ColorNames.grey),
    labelColorY: PD.Color(ColorNames.grey),
    labelColorZ: PD.Color(ColorNames.grey),
    labelOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    labelScale: PD.Numeric(0.25, { min: 0.1, max: 1.0, step: 0.01 }),
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
    scene: Scene;
    camera: Camera;
    props: CameraHelperProps = {
        axes: { name: 'off', params: {} }
    };

    private meshRenderObject: GraphicsRenderObject | undefined;
    private textRenderObject: GraphicsRenderObject | undefined;
    private pixelRatio = 1;

    constructor(private webgl: WebGLContext, props: Partial<CameraHelperProps> = {}) {
        this.scene = Scene.create(webgl, 'blended');

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
                    this.pixelRatio = this.webgl.pixelRatio;
                    const params = {
                        ...props.axes.params,
                        scale: props.axes.params.scale * this.pixelRatio,
                        labelScale: props.axes.params.labelScale * this.pixelRatio,
                    };
                    this.meshRenderObject = createMeshRenderObject(params);
                    this.scene.add(this.meshRenderObject);
                    if (props.axes.params.showLabels) {
                        this.textRenderObject = createTextRenderObject(params);
                        this.scene.add(this.textRenderObject);
                    } else {
                        this.textRenderObject = undefined;
                    }
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
        if ((
            (!this.meshRenderObject || objectId !== this.meshRenderObject.id) &&
            (!this.textRenderObject || objectId !== this.textRenderObject.id)
        ) || groupId === CameraHelperAxis.None) return EmptyLoci;
        return CameraAxesLoci(this, groupId, instanceId);
    }

    private eachGroup = (loci: Loci, apply: (interval: Interval) => boolean): boolean => {
        if (!isCameraAxesLoci(loci)) return false;
        let changed = false;
        if (this.meshRenderObject) {
            const groupCount = this.meshRenderObject.values.uGroupCount.ref.value;
            for (const { groupId, instanceId } of loci.elements) {
                const idx = instanceId * groupCount + groupId;
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
        if (this.textRenderObject) {
            const groupCount = this.textRenderObject.values.uGroupCount.ref.value;
            for (const { groupId, instanceId } of loci.elements) {
                const idx = instanceId * groupCount + groupId;
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
        return changed;
    };

    mark(loci: Loci, action: MarkerAction) {
        if (!MarkerActions.is(MarkerActions.Highlighting, action)) return false;
        if (!isEveryLoci(loci)) {
            if (!isCameraAxesLoci(loci)) return false;
            if (loci.data !== this) return false;
        }
        return (
            Visual.mark(this.meshRenderObject, loci, action, this.eachGroup) ||
            Visual.mark(this.textRenderObject, loci, action, this.eachGroup)
        );
    }

    update(camera: ICamera) {
        if (!this.meshRenderObject || this.props.axes.name === 'off') return;

        if (this.pixelRatio !== this.webgl.pixelRatio) {
            this.setProps(this.props);
        }

        updateCamera(this.camera, camera.viewport, camera.viewOffset);
        Mat4.extractRotation(this.scene.view, camera.view);

        const r = this.textRenderObject
            ? this.textRenderObject.values.boundingSphere.ref.value.radius
            : this.meshRenderObject.values.boundingSphere.ref.value.radius;
        const l = this.props.axes.params.location;
        const ox = this.props.axes.params.locationOffsetX * this.pixelRatio;
        const oy = this.props.axes.params.locationOffsetY * this.pixelRatio;
        if (l === 'bottom-left') {
            Mat4.setTranslation(this.scene.view, Vec3.create(
                -camera.viewport.width / 2 + r + ox,
                -camera.viewport.height / 2 + r + oy,
                0
            ));
        } else if (l === 'bottom-right') {
            Mat4.setTranslation(this.scene.view, Vec3.create(
                camera.viewport.width / 2 - r - ox,
                -camera.viewport.height / 2 + r + oy,
                0
            ));
        } else if (l === 'top-left') {
            Mat4.setTranslation(this.scene.view, Vec3.create(
                -camera.viewport.width / 2 + r + ox,
                camera.viewport.height / 2 - r - oy,
                0
            ));
        } else if (l === 'top-right') {
            Mat4.setTranslation(this.scene.view, Vec3.create(
                camera.viewport.width / 2 - r - ox,
                camera.viewport.height / 2 - r - oy,
                0
            ));
        } else {
            assertUnreachable(l);
        }
    }
}

export enum CameraHelperAxis {
    None = 0,
    X,
    Y,
    Z,
    XY,
    XZ,
    YZ,
    Origin
}

function getAxisLabel(axis: number, cameraHelper: CameraHelper) {
    const a = cameraHelper.props.axes;
    const x = a.name === 'on' ? a.params.labelX : 'X';
    const y = a.name === 'on' ? a.params.labelY : 'Y';
    const z = a.name === 'on' ? a.params.labelZ : 'Z';
    switch (axis) {
        case CameraHelperAxis.X: return `${x} Axis`;
        case CameraHelperAxis.Y: return `${y} Axis`;
        case CameraHelperAxis.Z: return `${z} Axis`;
        case CameraHelperAxis.XY: return `${x}${y} Plane`;
        case CameraHelperAxis.XZ: return `${x}${z} Plane`;
        case CameraHelperAxis.YZ: return `${y}${z} Plane`;
        case CameraHelperAxis.Origin: return 'Origin';
        default: return 'Axes';
    }
}

function CameraAxesLoci(cameraHelper: CameraHelper, groupId: number, instanceId: number) {
    return DataLoci('camera-axes', cameraHelper, [{ groupId, instanceId }],
        void 0 /** bounding sphere */,
        () => getAxisLabel(groupId, cameraHelper));
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

function createAxesMesh(props: AxesProps, mesh?: Mesh) {
    const state = MeshBuilder.createState(512, 256, mesh);
    const scale = 100 * props.scale;
    const radius = props.radiusScale * scale;
    const textScale = props.showLabels ? 100 * props.labelScale / 3 : 0;
    const x = Vec3.scale(Vec3(), Vec3.unitX, scale - textScale);
    const y = Vec3.scale(Vec3(), Vec3.unitY, scale - textScale);
    const z = Vec3.scale(Vec3(), Vec3.unitZ, scale - textScale);
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments: 32 };

    state.currentGroup = CameraHelperAxis.Origin;
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

    if (props.showPlanes) {
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
    }

    return MeshBuilder.getMesh(state);
}

function getAxesMeshShape(props: AxesProps, shape?: Shape<Mesh>) {
    const scale = 100 * props.scale;
    const mesh = createAxesMesh(props, shape?.geometry);
    mesh.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale + scale / 4));
    const getColor = (groupId: number) => {
        switch (groupId) {
            case CameraHelperAxis.X: return props.colorX;
            case CameraHelperAxis.Y: return props.colorY;
            case CameraHelperAxis.Z: return props.colorZ;
            case CameraHelperAxis.XY: return props.planeColorXY;
            case CameraHelperAxis.XZ: return props.planeColorXZ;
            case CameraHelperAxis.YZ: return props.planeColorYZ;
            case CameraHelperAxis.Origin: return props.originColor;
            default: return ColorNames.grey;
        }
    };
    return Shape.create('axes-mesh', {}, mesh, getColor, () => 1, () => '');
}

function createMeshRenderObject(props: AxesProps) {
    const shape = getAxesMeshShape(props);
    return Shape.createRenderObject(shape, {
        ...PD.getDefaultValues(Mesh.Params),
        ...props,
        ignoreLight: true,
        cellSize: 0,
    });
}

//

function createAxesText(props: AxesProps, text?: Text) {
    const builder = TextBuilder.create(props, 8, 8, text);
    const scale = 100 * props.scale;

    const x = Vec3.scale(Vec3(), Vec3.unitX, scale);
    const y = Vec3.scale(Vec3(), Vec3.unitY, scale);
    const z = Vec3.scale(Vec3(), Vec3.unitZ, scale);

    const textScale = 100 * props.labelScale;
    builder.add(props.labelX, x[0], x[1], x[2], 0.0, textScale, CameraHelperAxis.X);
    builder.add(props.labelY, y[0], y[1], y[2], 0.0, textScale, CameraHelperAxis.Y);
    builder.add(props.labelZ, z[0], z[1], z[2], 0.0, textScale, CameraHelperAxis.Z);

    return builder.getText();
}

function getAxesTextShape(props: AxesProps, shape?: Shape<Text>) {
    const scale = 100 * props.scale;
    const text = createAxesText(props, shape?.geometry);
    text.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale));
    const getColor = (groupId: number) => {
        switch (groupId) {
            case CameraHelperAxis.X: return props.labelColorX;
            case CameraHelperAxis.Y: return props.labelColorY;
            case CameraHelperAxis.Z: return props.labelColorZ;
            default: return ColorNames.grey;
        }
    };
    return Shape.create('axes-text', {}, text, getColor, () => 1, () => '');
}

function createTextRenderObject(props: AxesProps) {
    const shape = getAxesTextShape(props);
    return Shape.createRenderObject(shape, {
        ...PD.getDefaultValues(Text.Params),
        ...props,
        alpha: props.labelOpacity,
        cellSize: 0,
    });
}
