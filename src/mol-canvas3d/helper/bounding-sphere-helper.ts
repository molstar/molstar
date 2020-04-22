/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject, GraphicsRenderObject, getNextMaterialId } from '../../mol-gl/render-object';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from '../../mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import Scene from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Sphere3D } from '../../mol-math/geometry';
import { Color } from '../../mol-util/color';
import { ColorNames } from '../../mol-util/color/names';
import { TransformData } from '../../mol-geo/geometry/transform-data';
import { sphereVertexCount } from '../../mol-geo/primitive/sphere';
import { ValueCell } from '../../mol-util';
import { Geometry } from '../../mol-geo/geometry/geometry';

export const DebugHelperParams = {
    sceneBoundingSpheres: PD.Boolean(false, { description: 'Show full scene bounding spheres.' }),
    visibleSceneBoundingSpheres: PD.Boolean(false, { description: 'Show visible scene bounding spheres.' }),
    objectBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of visible render objects.' }),
    instanceBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of visible instances.' }),
};
export type DebugHelperParams = typeof DebugHelperParams
export type DebugHelperProps = PD.Values<DebugHelperParams>

type BoundingSphereData = { boundingSphere: Sphere3D, renderObject: GraphicsRenderObject, mesh: Mesh }

export class BoundingSphereHelper {
    readonly scene: Scene

    private readonly parent: Scene
    private _props: DebugHelperProps
    private objectsData = new Map<GraphicsRenderObject, BoundingSphereData>()
    private instancesData = new Map<GraphicsRenderObject, BoundingSphereData>()
    private sceneData: BoundingSphereData | undefined
    private visibleSceneData: BoundingSphereData | undefined

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<DebugHelperProps>) {
        this.scene = Scene.create(ctx);
        this.parent = parent;
        this._props = { ...PD.getDefaultValues(DebugHelperParams), ...props };
    }

    update() {
        const newSceneData = updateBoundingSphereData(this.scene, this.parent.boundingSphere, this.sceneData, ColorNames.lightgrey, sceneMaterialId);
        if (newSceneData) this.sceneData = newSceneData;

        const newVisibleSceneData = updateBoundingSphereData(this.scene, this.parent.boundingSphereVisible, this.visibleSceneData, ColorNames.black, visibleSceneMaterialId);
        if (newVisibleSceneData) this.visibleSceneData = newVisibleSceneData;

        this.parent.forEach((r, ro) => {
            const objectData = this.objectsData.get(ro);
            const newObjectData = updateBoundingSphereData(this.scene, r.values.boundingSphere.ref.value, objectData, ColorNames.tomato, objectMaterialId);
            if (newObjectData) this.objectsData.set(ro, newObjectData);

            const instanceData = this.instancesData.get(ro);
            const newInstanceData = updateBoundingSphereData(this.scene, r.values.invariantBoundingSphere.ref.value, instanceData, ColorNames.skyblue, instanceMaterialId, {
                aTransform: ro.values.aTransform,
                matrix: ro.values.matrix,
                transform: ro.values.transform,
                extraTransform: ro.values.extraTransform,
                uInstanceCount: ro.values.uInstanceCount,
                instanceCount: ro.values.instanceCount,
                aInstance: ro.values.aInstance,
            });
            if (newInstanceData) this.instancesData.set(ro, newInstanceData);
        });

        this.objectsData.forEach((objectData, ro) => {
            if (!this.parent.has(ro)) {
                this.scene.remove(objectData.renderObject);
                this.objectsData.delete(ro);
            }
        });
        this.instancesData.forEach((instanceData, ro) => {
            if (!this.parent.has(ro)) {
                this.scene.remove(instanceData.renderObject);
                this.instancesData.delete(ro);
            }
        });

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    syncVisibility() {
        if (this.sceneData) {
            this.sceneData.renderObject.state.visible = this._props.sceneBoundingSpheres;
        }

        if (this.visibleSceneData) {
            this.visibleSceneData.renderObject.state.visible = this._props.visibleSceneBoundingSpheres;
        }

        this.parent.forEach((_, ro) => {
            const objectData = this.objectsData.get(ro);
            if (objectData) objectData.renderObject.state.visible = ro.state.visible && this._props.objectBoundingSpheres;

            const instanceData = this.instancesData.get(ro);
            if (instanceData) instanceData.renderObject.state.visible = ro.state.visible && this._props.instanceBoundingSpheres;
        });
    }

    clear() {
        this.sceneData = undefined;
        this.objectsData.clear();
        this.scene.clear();
    }

    get isEnabled() {
        return (
            this._props.sceneBoundingSpheres || this._props.visibleSceneBoundingSpheres ||
            this._props.objectBoundingSpheres || this._props.instanceBoundingSpheres
        );
    }
    get props() { return this._props as Readonly<DebugHelperProps>; }

    setProps (props: Partial<DebugHelperProps>) {
        Object.assign(this._props, props);
        if (this.isEnabled) this.update();
    }
}

function updateBoundingSphereData(scene: Scene, boundingSphere: Sphere3D, data: BoundingSphereData | undefined, color: Color, materialId: number, transform?: TransformData) {
    if (!data || !Sphere3D.equals(data.boundingSphere, boundingSphere)) {
        const mesh = createBoundingSphereMesh(boundingSphere, data && data.mesh);
        const renderObject = data ? data.renderObject : createBoundingSphereRenderObject(mesh, color, materialId, transform);
        if (data) {
            ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(mesh));
        } else {
            scene.add(renderObject);
        }
        return { boundingSphere: Sphere3D.clone(boundingSphere), renderObject, mesh };
    }
}

function createBoundingSphereMesh(boundingSphere: Sphere3D, mesh?: Mesh) {
    const detail = 2;
    const vertexCount = sphereVertexCount(detail);
    const builderState = MeshBuilder.createState(vertexCount, vertexCount / 2, mesh);
    if (boundingSphere.radius) {
        addSphere(builderState, boundingSphere.center, boundingSphere.radius, detail);
        if (Sphere3D.hasExtrema(boundingSphere)) {
            for (const e of boundingSphere.extrema) addSphere(builderState, e, 1.0, 0);
        }
    }
    return MeshBuilder.getMesh(builderState);
}

const sceneMaterialId = getNextMaterialId();
const visibleSceneMaterialId = getNextMaterialId();
const objectMaterialId = getNextMaterialId();
const instanceMaterialId = getNextMaterialId();

function createBoundingSphereRenderObject(mesh: Mesh, color: Color, materialId: number, transform?: TransformData) {
    const values = Mesh.Utils.createValuesSimple(mesh, { alpha: 0.1, doubleSided: false }, color, 1, transform);
    return createRenderObject('mesh', values, { visible: true, alphaFactor: 1, pickable: false, opaque: false, writeDepth: false }, materialId);
}