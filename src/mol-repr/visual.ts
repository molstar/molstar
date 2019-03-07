/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from 'mol-task'
import { GraphicsRenderObject } from 'mol-gl/render-object'
import { PickingId } from '../mol-geo/geometry/picking';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../mol-geo/geometry/marker-data';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { WebGLContext } from 'mol-gl/webgl/context';
import { Theme } from 'mol-theme/theme';
import { Mat4 } from 'mol-math/linear-algebra';
import { updateTransformData, fillIdentityTransform } from 'mol-geo/geometry/transform-data';
import { calculateTransformBoundingSphere } from 'mol-gl/renderable/util';
import { ValueCell } from 'mol-util';
import { Overpaint } from 'mol-theme/overpaint';

export interface VisualContext {
    readonly runtime: RuntimeContext
    readonly webgl?: WebGLContext
}
// export type VisualFactory<D, P extends PD.Params> = (ctx: VisualContext) => Visual<D, P>

export { Visual }
interface Visual<D, P extends PD.Params> {
    /** Number of addressable groups in all instances of the visual */
    readonly groupCount: number
    readonly renderObject: GraphicsRenderObject | undefined
    createOrUpdate: (ctx: VisualContext, theme: Theme, props?: Partial<PD.Values<P>>, data?: D) => Promise<void> | void
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    setVisibility: (visible: boolean) => void
    setAlphaFactor: (alphaFactor: number) => void
    setPickable: (pickable: boolean) => void
    setTransform: (matrix?: Mat4, instanceMatrices?: Float32Array | null) => void
    setOverpaint: (layers: Overpaint.Layers) => void
    destroy: () => void
}
namespace Visual {
    export function setVisibility(renderObject: GraphicsRenderObject | undefined, visible: boolean) {
        if (renderObject) renderObject.state.visible = visible
    }

    export function setAlphaFactor(renderObject: GraphicsRenderObject | undefined, alphaFactor: number) {
        if (renderObject) renderObject.state.alphaFactor = alphaFactor
    }

    export function setPickable(renderObject: GraphicsRenderObject | undefined, pickable: boolean) {
        if (renderObject) renderObject.state.pickable = pickable
    }

    export function setTransform(renderObject: GraphicsRenderObject | undefined, transform?: Mat4, instanceTransforms?: Float32Array | null) {
        if (renderObject && (transform || instanceTransforms)) {
            const { values } = renderObject
            if (transform) {
                Mat4.copy(values.matrix.ref.value, transform)
                ValueCell.update(values.matrix, values.matrix.ref.value)
            }
            if (instanceTransforms) {
                values.extraTransform.ref.value.set(instanceTransforms)
                ValueCell.update(values.extraTransform, values.extraTransform.ref.value)
            } else if (instanceTransforms === null) {
                fillIdentityTransform(values.extraTransform.ref.value, values.instanceCount.ref.value)
                ValueCell.update(values.extraTransform, values.extraTransform.ref.value)
            }
            updateTransformData(values)
            const boundingSphere = calculateTransformBoundingSphere(values.invariantBoundingSphere.ref.value, values.aTransform.ref.value, values.instanceCount.ref.value)
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
    }
}