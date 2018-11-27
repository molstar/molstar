/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Geometry } from 'mol-geo/geometry/geometry';
import { ValueCell } from 'mol-util';
import Scene from 'mol-gl/scene';

export class BoundingSphereHelper {
    private mesh: Mesh
    private renderObject: MeshRenderObject

    constructor(private scene: Scene, visible: boolean) {
        this.mesh = MeshBuilder.getMesh(MeshBuilder.createState(1024, 512))
        const values = Mesh.createValuesSimple(this.mesh, { alpha: 0.1, doubleSided: false })
        this.renderObject = createMeshRenderObject(values, { visible, pickable: false, opaque: false })
        scene.add(this.renderObject)
    }

    update() {
        const builderState = MeshBuilder.createState(1024, 512, this.mesh)
        if (this.scene.boundingSphere.radius) {
            addSphere(builderState, this.scene.boundingSphere.center, this.scene.boundingSphere.radius, 2)
        }
        this.scene.forEach(r => {
            if (r.boundingSphere.radius) {
                addSphere(builderState, r.boundingSphere.center, r.boundingSphere.radius, 2)
            }
        })
        this.mesh = MeshBuilder.getMesh(builderState)
        ValueCell.update(this.renderObject.values.drawCount, Geometry.getDrawCount(this.mesh))
    }

    get visible() { return this.renderObject.state.visible }
    set visible(value: boolean) { this.renderObject.state.visible = value }
}