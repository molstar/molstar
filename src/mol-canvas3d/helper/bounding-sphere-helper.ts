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

export const DebugHelperParams = {
    sceneBoundingSpheres: PD.Boolean(false, { description: 'Show scene bounding spheres.' }),
    objectBoundingSpheres: PD.Boolean(false, { description: 'Show bounding spheres of render objects.' }),
}
export type DebugHelperParams = typeof DebugHelperParams
export type DebugHelperProps = PD.Values<DebugHelperParams>

type BoundingSphereData = { boundingSphere: Sphere3D, renderObject: RenderObject }

// TODO per-object-transform bounding spheres

export class BoundingSphereHelper {
    readonly scene: Scene
    private readonly parent: Scene
    private _props: DebugHelperProps
    private objectsData = new Map<RenderObject, BoundingSphereData>()
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
            let objectData = this.objectsData.get(ro)
            const newObjectData = updateBoundingSphereData(this.scene, r.boundingSphere, objectData)
            if (newObjectData) this.objectsData.set(ro, newObjectData)
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
            else console.error('expected to have debug render object')
        })
    }

    clear() {
        this.sceneData = undefined
        this.objectsData.clear()
        this.scene.clear()
    }

    get isEnabled() {
        return this._props.sceneBoundingSpheres || this._props.objectBoundingSpheres
    }
    get props() { return this._props as Readonly<DebugHelperProps> }

    setProps (props: Partial<DebugHelperProps>) {
        Object.assign(this._props, props)
        if (this.isEnabled) this.update()
    }
}

function updateBoundingSphereData(scene: Scene, boundingSphere: Sphere3D, data?: BoundingSphereData) {
    if (!data || !Sphere3D.exactEquals(data.boundingSphere, boundingSphere)) {
        if (data) scene.remove(data.renderObject)
        const renderObject = createBoundingSphereRenderObject(boundingSphere)
        scene.add(renderObject)
        return { boundingSphere, renderObject }
    }
}

function createBoundingSphereRenderObject(boundingSphere: Sphere3D) {
    const builderState = MeshBuilder.createState(1024, 512)
    if (boundingSphere.radius) {
        addSphere(builderState, boundingSphere.center, boundingSphere.radius, 2)
    } else if (isNaN(boundingSphere.radius)) {
        console.warn('boundingSphere.radius is NaN')
    }
    const mesh = MeshBuilder.getMesh(builderState)
    const values = Mesh.createValuesSimple(mesh, { alpha: 0.1, doubleSided: false })
    return createMeshRenderObject(values, { visible: true, pickable: false, opaque: false })
}