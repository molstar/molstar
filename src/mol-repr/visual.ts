/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext } from '../mol-task'
import { GraphicsRenderObject } from '../mol-gl/render-object'
import { PickingId } from '../mol-geo/geometry/picking';
import { Loci, isEmptyLoci } from '../mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../mol-geo/geometry/marker-data';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { WebGLContext } from '../mol-gl/webgl/context';
import { Theme } from '../mol-theme/theme';
import { Mat4 } from '../mol-math/linear-algebra';
import { updateTransformData, fillIdentityTransform } from '../mol-geo/geometry/transform-data';
import { calculateTransformBoundingSphere } from '../mol-gl/renderable/util';
import { ValueCell } from '../mol-util';
import { Overpaint } from '../mol-theme/overpaint';
import { createOverpaint, clearOverpaint, applyOverpaintColor } from '../mol-geo/geometry/overpaint-data';
import { Interval } from '../mol-data/int';
import { Transparency } from '../mol-theme/transparency';
import { createTransparency, clearTransparency, applyTransparencyValue } from '../mol-geo/geometry/transparency-data';

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
    setOverpaint: (overpaint: Overpaint) => void
    setTransparency: (transparency: Transparency) => void
    destroy: () => void
}
namespace Visual {
    export type LociApply = (loci: Loci, apply: (interval: Interval) => boolean) => boolean

    export function setVisibility(renderObject: GraphicsRenderObject | undefined, visible: boolean) {
        if (renderObject) renderObject.state.visible = visible
    }

    export function setAlphaFactor(renderObject: GraphicsRenderObject | undefined, alphaFactor: number) {
        if (renderObject) renderObject.state.alphaFactor = alphaFactor
    }

    export function setPickable(renderObject: GraphicsRenderObject | undefined, pickable: boolean) {
        if (renderObject) renderObject.state.pickable = pickable
    }

    export function mark(renderObject: GraphicsRenderObject | undefined, loci: Loci, action: MarkerAction, lociApply: LociApply) {
        if (!renderObject) return false

        const { tMarker } = renderObject.values

        function apply(interval: Interval) {
            const start = Interval.start(interval)
            const end = Interval.end(interval)
            return applyMarkerAction(tMarker.ref.value.array, start, end, action)
        }

        const changed = lociApply(loci, apply)
        if (changed) ValueCell.update(tMarker, tMarker.ref.value)
        return changed
    }

    export function setOverpaint(renderObject: GraphicsRenderObject | undefined, overpaint: Overpaint, lociApply: LociApply, clear: boolean) {
        if (!renderObject) return

        const { tOverpaint, uGroupCount, instanceCount } = renderObject.values
        const count = uGroupCount.ref.value * instanceCount.ref.value

        // ensure texture has right size
        createOverpaint(overpaint.layers.length ? count : 0, renderObject.values)

        // clear if requested
        if (clear) clearOverpaint(tOverpaint.ref.value.array, 0, count)

        for (let i = 0, il = overpaint.layers.length; i < il; ++i) {
            const { loci, color } = overpaint.layers[i]
            const apply = (interval: Interval) => {
                const start = Interval.start(interval)
                const end = Interval.end(interval)
                return applyOverpaintColor(tOverpaint.ref.value.array, start, end, color, overpaint.alpha)
            }
            lociApply(loci, apply)
        }
        ValueCell.update(tOverpaint, tOverpaint.ref.value)
    }

    export function setTransparency(renderObject: GraphicsRenderObject | undefined, transparency: Transparency, lociApply: LociApply, clear: boolean) {
        if (!renderObject) return

        const { tTransparency, uGroupCount, instanceCount } = renderObject.values
        const count = uGroupCount.ref.value * instanceCount.ref.value

        const { loci, value, variant } = transparency

        // ensure texture has right size and variant
        createTransparency(value && !isEmptyLoci(loci) ? count : 0, variant, renderObject.values)

        // clear if requested
        if (clear) clearTransparency(tTransparency.ref.value.array, 0, count)

        const apply = (interval: Interval) => {
            const start = Interval.start(interval)
            const end = Interval.end(interval)
            return applyTransparencyValue(tTransparency.ref.value.array, start, end, value)
        }
        lociApply(loci, apply)

        ValueCell.update(tTransparency, tTransparency.ref.value)
    }

    export function setTransform(renderObject: GraphicsRenderObject | undefined, transform?: Mat4, instanceTransforms?: Float32Array | null) {
        if (!renderObject || (!transform && !instanceTransforms)) return

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