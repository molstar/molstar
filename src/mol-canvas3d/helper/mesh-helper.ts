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
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../mol-geo/geometry/lines/lines-builder';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { MeshValues } from '../../mol-gl/renderable/mesh';

export const MeshHelperParams = {
    meshNormals: PD.Boolean(false, { description: 'Show normals of visible mesh render objects.' }),
};
export type MeshHelperParams = typeof MeshHelperParams;
export type MeshHelperProps = PD.Values<MeshHelperParams>;

const meshHelperMaterialId = getNextMaterialId();

const _v = Vec3();
const _n = Vec3();
const _start = Vec3();
const _end = Vec3();

export class MeshHelper {
    readonly scene: Scene;

    private readonly parent: Scene;
    private _props: MeshHelperProps;
    private renderObjects = new Map<number, GraphicsRenderObject>();

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<MeshHelperProps>) {
        this.scene = Scene.create(ctx, 'blended');
        this.parent = parent;
        this._props = { ...PD.getDefaultValues(MeshHelperParams), ...props };
    }

    update() {
        const previousIds = new Set(this.renderObjects.keys());
        const currentIds = new Set<number>();

        this.parent.forEach((r, ro) => {
            if (!ro.state.visible) return;
            if (ro.type !== 'mesh') return;

            currentIds.add(ro.id);

            // Skip if we already have normals for this render object
            if (this.renderObjects.has(ro.id)) {
                previousIds.delete(ro.id);
                return;
            }

            const values = ro.values as MeshValues;
            const lines = createNormalLines(values);
            if (!lines) return;

            const linesRO = createNormalLinesRenderObject(lines, meshHelperMaterialId);
            this.scene.add(linesRO);
            this.renderObjects.set(ro.id, linesRO);
        });

        // Remove normals for render objects no longer present
        for (const id of previousIds) {
            const linesRO = this.renderObjects.get(id);
            if (linesRO) {
                this.scene.remove(linesRO);
                this.renderObjects.delete(id);
            }
        }

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    syncVisibility() {
        const visible = this._props.meshNormals;
        this.renderObjects.forEach(ro => {
            ro.state.visible = visible;
        });
    }

    clear() {
        this.renderObjects.clear();
        this.scene.clear();
    }

    get isEnabled() {
        return this._props.meshNormals;
    }

    get props() { return this._props as Readonly<MeshHelperProps>; }

    setProps(props: Partial<MeshHelperProps>) {
        Object.assign(this._props, props);
        if (this.isEnabled) this.update();
    }
}

//

function createNormalLines(values: MeshValues): Lines | undefined {
    const positions = values.aPosition.ref.value;
    const normals = values.aNormal.ref.value;
    const indices = values.elements.ref.value;
    const transforms = values.aTransform.ref.value;
    const instanceCount = values.uInstanceCount.ref.value;

    const vertexCount = positions.length / 3;
    if (vertexCount === 0) return undefined;

    // Determine normal line length: proportional to bounding sphere radius
    const bs = values.boundingSphere.ref.value;
    const normalLength = Math.max(bs.radius * 0.01, 0.1);

    // Count unique vertices referenced by indices
    const indexCount = values.drawCount.ref.value;

    const builder = LinesBuilder.create(indexCount * instanceCount);

    for (let inst = 0; inst < instanceCount; ++inst) {
        const tOffset = inst * 16;
        const transform = Mat4();
        Mat4.fromArray(transform, transforms, tOffset);

        // Use a set to avoid drawing duplicate normals for shared vertices
        const visited = new Set<number>();

        for (let i = 0; i < indexCount; ++i) {
            const vi = indices[i];
            if (visited.has(vi)) continue;
            visited.add(vi);

            const vo = vi * 3;
            Vec3.set(_v, positions[vo], positions[vo + 1], positions[vo + 2]);
            Vec3.set(_n, normals[vo], normals[vo + 1], normals[vo + 2]);

            // Transform vertex position and normal direction by instance transform
            Vec3.transformMat4(_start, _v, transform);
            Vec3.transformDirection(_end, _n, transform);
            Vec3.normalize(_end, _end);
            Vec3.scaleAndAdd(_end, _start, _end, normalLength);

            builder.addVec(_start, _end, 0);
        }
    }

    return builder.getLines();
}

function createNormalLinesRenderObject(lines: Lines, materialId: number): GraphicsRenderObject {
    const props = { ...PD.getDefaultValues(Lines.Params), sizeFactor: 1, alpha: 0.7 };
    const values = Lines.Utils.createValuesSimple(lines, props, ColorNames.magenta, 1);
    const state = Lines.Utils.createRenderableState(props);
    state.pickable = false;
    return createRenderObject('lines', values, state, materialId);
}
