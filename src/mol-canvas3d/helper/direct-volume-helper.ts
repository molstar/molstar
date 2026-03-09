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
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { DirectVolumeValues } from '../../mol-gl/renderable/direct-volume';
import { addBox } from '../../mol-geo/geometry/lines/builder/box';

export const DirectVolumeHelperParams = {
    directVolumeEdges: PD.Boolean(false, { description: 'Show edges of visible direct-volume render objects.' }),
};
export type DirectVolumeHelperParams = typeof DirectVolumeHelperParams;
export type DirectVolumeHelperProps = PD.Values<DirectVolumeHelperParams>;

const directVolumeMaterialId = getNextMaterialId();

type TrackedEntry = { ro: GraphicsRenderObject, version: number };

export class DirectVolumeHelper {
    readonly scene: Scene;

    private readonly parent: Scene;
    private _props: DirectVolumeHelperProps;
    private renderObjects = new Map<number, TrackedEntry>();

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<DirectVolumeHelperProps>) {
        this.scene = Scene.create(ctx, 'blended');
        this.parent = parent;
        this._props = { ...PD.getDefaultValues(DirectVolumeHelperParams), ...props };
    }

    update() {
        const previousIds = new Set(this.renderObjects.keys());

        this.parent.forEach((r, ro) => {
            if (!ro.state.visible) return;
            if (ro.type !== 'direct-volume') return;

            const values = ro.values as DirectVolumeValues;
            const version = values.uUnitToCartn.ref.version + values.uGridDim.ref.version + values.aTransform.ref.version;

            const existing = this.renderObjects.get(ro.id);
            if (existing && existing.version === version) {
                previousIds.delete(ro.id);
                return;
            }

            // Remove old entry if version changed
            if (existing) {
                this.scene.remove(existing.ro);
                this.renderObjects.delete(ro.id);
            }

            const lines = createVolumeEdgeLines(values);
            if (!lines) return;

            const linesRO = createLinesRenderObject(lines, directVolumeMaterialId);
            this.scene.add(linesRO);
            this.renderObjects.set(ro.id, { ro: linesRO, version });
            previousIds.delete(ro.id);
        });

        for (const id of previousIds) {
            const entry = this.renderObjects.get(id);
            if (entry) {
                this.scene.remove(entry.ro);
                this.renderObjects.delete(id);
            }
        }

        this.scene.update(void 0, false);
        this.scene.commit();
    }

    syncVisibility() {
        const visible = this._props.directVolumeEdges;
        this.renderObjects.forEach(entry => {
            entry.ro.state.visible = visible;
        });
    }

    clear() {
        this.renderObjects.clear();
        this.scene.clear();
    }

    get isEnabled() {
        return this._props.directVolumeEdges;
    }

    get props() { return this._props as Readonly<DirectVolumeHelperProps>; }

    setProps(props: Partial<DirectVolumeHelperProps>) {
        Object.assign(this._props, props);
        if (this.isEnabled) this.update();
    }
}

//

/**
 * The volume proxy box in the shader uses aPosition in [-0.5, 0.5]^3,
 * shifted to [0,1]^3 (unitCoord = aPosition + 0.5), then transformed by:
 *   uUnitToCartn → Cartesian space
 *   aTransform → instance space
 *
 * We replicate this pipeline to get the correct world-space edges.
 * Grid ticks are placed at 1/gridDim intervals along each edge.
 */
function createVolumeEdgeLines(values: DirectVolumeValues): Lines | undefined {
    const unitToCartn = values.uUnitToCartn.ref.value;
    const transforms = values.aTransform.ref.value;
    const instanceCount = values.uInstanceCount.ref.value;

    const bs = values.boundingSphere.ref.value;
    if (bs.radius < 1e-6) return undefined;

    const builder = LinesBuilder.create(128 * instanceCount);
    for (let inst = 0; inst < instanceCount; ++inst) {
        const instTransform = Mat4();
        Mat4.fromArray(instTransform, transforms, inst * 16);
        // Combined transform: aTransform * uUnitToCartn
        const combined = Mat4.mul(Mat4(), instTransform, unitToCartn);
        addBox(builder, combined, 0);
    }
    return builder.getLines();
}

function createLinesRenderObject(lines: Lines, materialId: number): GraphicsRenderObject {
    const props = { ...PD.getDefaultValues(Lines.Params), sizeFactor: 1, alpha: 0.8 };
    const values = Lines.Utils.createValuesSimple(lines, props, ColorNames.orange, 1);
    const state = Lines.Utils.createRenderableState(props);
    state.pickable = false;
    return createRenderObject('lines', values, state, materialId);
}
