/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createMeshRenderObject, RenderObject } from 'mol-gl/render-object'
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import Scene from 'mol-gl/scene';
import { WebGLContext } from 'mol-gl/webgl/context';
import { Sphere3D } from 'mol-math/geometry';
import { Vec3, Mat4 } from 'mol-math/linear-algebra';

export const DebugHelperParams = {
    sceneBoundingSpheres: PD.Boolean(false, { description: 'Show scene bounding spheres.' }),
    objectBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of render objects.' }),
    instanceBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of instances.' }),
}
export type DebugHelperParams = typeof DebugHelperParams
export type DebugHelperProps = PD.Values<DebugHelperParams>

type BoundingSphereData = { boundingSphere: Sphere3D, renderObject: RenderObject }

export class BoundingSphereHelper {
    readonly scene: Scene
    private readonly parent: Scene
    private _props: DebugHelperProps
    private objectsData = new Map<RenderObject, BoundingSphereData>()
    private instancesData = new Map<RenderObject, BoundingSphereData>()
    private sceneData: BoundingSphereData | undefined

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<DebugHelperProps>) {
        this.scene = Scene.create(ctx)
        this.parent = parent
        this._props = { ...PD.getDefaultValues(DebugHelperParams), ...props }
    }

    update() {
        const newSceneData = updateBoundingSphereData(this.scene, this.parent.boundingSphere, this.sceneData)
        if (newSceneData) this.sceneData = newSceneData

        const oldRO = new Set<RenderObject>()
        this.parent.forEach((r, ro) => {
            const objectData = this.objectsData.get(ro)
            const newObjectData = updateBoundingSphereData(this.scene, r.boundingSphere, objectData)
            if (newObjectData) this.objectsData.set(ro, newObjectData)

            if (ro.type === 'mesh' || ro.type === 'lines' || ro.type === 'points') {
                const instanceData = this.instancesData.get(ro)
                const newInstanceData = updateBoundingSphereData(this.scene, ro.values.invariantBoundingSphere.ref.value, instanceData, ro.values.aTransform.ref.value, ro.values.instanceCount.ref.value)
                if (newInstanceData) this.instancesData.set(ro, newInstanceData)
            }

            oldRO.delete(ro)
        })
        oldRO.forEach(ro => {
            const objectData = this.objectsData.get(ro)
            if (objectData) {
                this.scene.remove(objectData.renderObject)
                this.objectsData.delete(ro)
            }
        })
    }

    syncVisibility() {
        if(this.sceneData) {
            this.sceneData.renderObject.state.visible = this._props.sceneBoundingSpheres
        }

        this.parent.forEach((_, ro) => {
            const objectData = this.objectsData.get(ro)
            if (objectData) objectData.renderObject.state.visible = ro.state.visible && this._props.objectBoundingSpheres

            const instanceData = this.instancesData.get(ro)
            if (instanceData) instanceData.renderObject.state.visible = ro.state.visible && this._props.instanceBoundingSpheres
        })
    }

    clear() {
        this.sceneData = undefined
        this.objectsData.clear()
        this.scene.clear()
    }

    get isEnabled() {
        return this._props.sceneBoundingSpheres || this._props.objectBoundingSpheres || this._props.instanceBoundingSpheres
    }
    get props() { return this._props as Readonly<DebugHelperProps> }

    setProps (props: Partial<DebugHelperProps>) {
        Object.assign(this._props, props)
        if (this.isEnabled) this.update()
    }
}

function updateBoundingSphereData(scene: Scene, boundingSphere: Sphere3D, data: BoundingSphereData | undefined, transform?: Float32Array, transformCount?: number) {
    if (!data || !Sphere3D.exactEquals(data.boundingSphere, boundingSphere)) {
        if (data) scene.remove(data.renderObject)
        const renderObject = createBoundingSphereRenderObject(boundingSphere, transform, transformCount)
        scene.add(renderObject)
        return { boundingSphere, renderObject }
    }
}

const tmpCenter = Vec3.zero()
const tmpM = Mat4.identity()
function createBoundingSphereRenderObject(boundingSphere: Sphere3D, transform?: Float32Array, transformCount?: number) {
    const builderState = MeshBuilder.createState(1024, 512)
    if (transform && transformCount) {
        // TODO create instanced mesh?
        for (let i = 0, _i = transformCount; i < _i; ++i) {
            Mat4.fromArray(tmpM, transform, i * 16)
            Vec3.transformMat4(tmpCenter, boundingSphere.center, tmpM)
            if (boundingSphere.radius) addSphere(builderState, tmpCenter, boundingSphere.radius, 1)
        }
    } else {
        if (boundingSphere.radius) addSphere(builderState, boundingSphere.center, boundingSphere.radius, 2)
    }
    const mesh = MeshBuilder.getMesh(builderState)
    const values = Mesh.createValuesSimple(mesh, { alpha: 0.1, doubleSided: false })
    return createMeshRenderObject(values, { visible: true, pickable: false, opaque: false })
}