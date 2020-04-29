/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import Scene from '../../mol-gl/scene';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3, Mat4, Mat3 } from '../../mol-math/linear-algebra';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { GraphicsRenderObject } from '../../mol-gl/render-object';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ColorNames } from '../../mol-util/color/names';
import { addCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { ValueCell } from '../../mol-util';
import { Sphere3D } from '../../mol-math/geometry';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import produce from 'immer';
import { Shape } from '../../mol-model/shape';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Camera } from '../camera';
import { DataLoci, EmptyLoci, Loci } from '../../mol-model/loci';
import { MarkerAction, MarkerActions } from '../../mol-util/marker-action';
import { Visual } from '../../mol-repr/visual';
import { Interval } from '../../mol-data/int';

const HandleParams = {
    ...Mesh.Params,
    alpha: { ...Mesh.Params.alpha, defaultValue: 1 },
    ignoreLight: { ...Mesh.Params.ignoreLight, defaultValue: true },
    colorX: PD.Color(ColorNames.red, { isEssential: true }),
    colorY: PD.Color(ColorNames.green, { isEssential: true }),
    colorZ: PD.Color(ColorNames.blue, { isEssential: true }),
    scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.1 }, { isEssential: true }),
};
type HandleParams = typeof HandleParams
type HandleProps = PD.Values<HandleParams>

export const HandleHelperParams = {
    handle: PD.MappedStatic('off', {
        on: PD.Group(HandleParams),
        off: PD.Group({})
    }, { cycle: true, description: 'Show handle tool' }),
};
export type HandleHelperParams = typeof HandleHelperParams
export type HandleHelperProps = PD.Values<HandleHelperParams>

export class HandleHelper {
    scene: Scene
    props: HandleHelperProps = {
        handle: { name: 'off', params: {} }
    }

    private renderObject: GraphicsRenderObject | undefined

    private _transform = Mat4();
    getBoundingSphere(out: Sphere3D, instanceId: number) {
        if (this.renderObject) {
            Sphere3D.copy(out, this.renderObject.values.invariantBoundingSphere.ref.value);
            Mat4.fromArray(this._transform, this.renderObject.values.aTransform.ref.value, instanceId * 16);
            Sphere3D.transform(out, out, this._transform);
        }
        return out;
    }

    setProps(props: Partial<HandleHelperProps>) {
        this.props = produce(this.props, p => {
            if (props.handle !== undefined) {
                p.handle.name = props.handle.name;
                if (props.handle.name === 'on') {
                    this.scene.clear();
                    const params = { ...props.handle.params, scale: props.handle.params.scale * this.webgl.pixelRatio };
                    this.renderObject = createHandleRenderObject(params);
                    this.scene.add(this.renderObject);
                    this.scene.commit();

                    p.handle.params = { ...props.handle.params };
                }
            }
        });
    }

    get isEnabled() {
        return this.props.handle.name === 'on';
    }

    // TODO could be a lists of position/rotation if we want to show more than one handle tool,
    //      they would be distingishable by their instanceId
    update(camera: Camera, position: Vec3, rotation: Mat3) {
        if (!this.renderObject) return;

        Mat4.setTranslation(this.renderObject.values.aTransform.ref.value as unknown as Mat4, position);
        Mat4.fromMat3(this.renderObject.values.aTransform.ref.value as unknown as Mat4, rotation);

        // TODO make invariant to camera scaling by adjusting renderObject transform

        ValueCell.update(this.renderObject.values.aTransform, this.renderObject.values.aTransform.ref.value);
        this.scene.update([this.renderObject], true);
    }

    getLoci(pickingId: PickingId) {
        const { objectId, groupId, instanceId } = pickingId;
        if (!this.renderObject || objectId !== this.renderObject.id) return EmptyLoci;
        return HandleLoci(this, groupId, instanceId);
    }

    private eachGroup = (loci: Loci, apply: (interval: Interval) => boolean): boolean => {
        if (!this.renderObject) return false;
        if (!isHandleLoci(loci)) return false;
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
        if (!isHandleLoci(loci)) return false;
        if (loci.data !== this) return false;
        return Visual.mark(this.renderObject, loci, action, this.eachGroup);
    }

    constructor(private webgl: WebGLContext, props: Partial<HandleHelperProps> = {}) {
        this.scene = Scene.create(webgl);
        this.setProps(props);
    }
}

function createHandleMesh(scale: number, mesh?: Mesh) {
    const state = MeshBuilder.createState(512, 256, mesh);
    const radius = 0.05 * scale;
    const x = Vec3.scale(Vec3(), Vec3.unitX, scale);
    const y = Vec3.scale(Vec3(), Vec3.unitY, scale);
    const z = Vec3.scale(Vec3(), Vec3.unitZ, scale);
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments: 32 };

    state.currentGroup = HandleGroup.TranslateScreenXY;
    addSphere(state, Vec3.origin, radius * 3, 2);

    state.currentGroup = HandleGroup.TranslateObjectX;
    addSphere(state, x, radius, 2);
    addCylinder(state, Vec3.origin, x, 1, cylinderProps);

    state.currentGroup = HandleGroup.TranslateObjectY;
    addSphere(state, y, radius, 2);
    addCylinder(state, Vec3.origin, y, 1, cylinderProps);

    state.currentGroup = HandleGroup.TranslateObjectZ;
    addSphere(state, z, radius, 2);
    addCylinder(state, Vec3.origin, z, 1, cylinderProps);

    // TODO add more helper geometries for the other HandleGroup options
    // TODO add props to create subset of geometries

    return MeshBuilder.getMesh(state);
}

export const HandleGroup = {
    None: 0,
    TranslateScreenXY: 1,
    // TranslateScreenZ: 2,
    TranslateObjectX: 3,
    TranslateObjectY: 4,
    TranslateObjectZ: 5,
    // TranslateObjectXY: 6,
    // TranslateObjectXZ: 7,
    // TranslateObjectYZ: 8,

    // RotateScreenZ: 9,
    // RotateObjectX: 10,
    // RotateObjectY: 11,
    // RotateObjectZ: 12,
} as const;

function HandleLoci(handleHelper: HandleHelper, groupId: number, instanceId: number) {
    return DataLoci('handle', handleHelper, [{ groupId, instanceId }],
        (boundingSphere: Sphere3D) => handleHelper.getBoundingSphere(boundingSphere, instanceId),
        () => `Handle Helper | Group Id ${groupId} | Instance Id ${instanceId}`);
}
export type HandleLoci = ReturnType<typeof HandleLoci>
export function isHandleLoci(x: Loci): x is HandleLoci {
    return x.kind === 'data-loci' && x.tag === 'handle';
}

function getHandleShape(props: HandleProps, shape?: Shape<Mesh>) {
    const scale = 10 * props.scale;
    const mesh = createHandleMesh(scale, shape?.geometry);
    mesh.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale + scale / 4));
    const getColor = (groupId: number) => {
        switch (groupId) {
            case HandleGroup.TranslateObjectX: return props.colorX;
            case HandleGroup.TranslateObjectY: return props.colorY;
            case HandleGroup.TranslateObjectZ: return props.colorZ;
            default: return ColorNames.grey;
        }
    };
    return Shape.create('handle', {}, mesh, getColor, () => 1, () => '');
}

function createHandleRenderObject(props: HandleProps) {
    const shape = getHandleShape(props);
    return Shape.createRenderObject(shape, props);
}