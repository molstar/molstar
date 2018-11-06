/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object'
import { PickingId } from '../mol-geo/geometry/picking';
import { Loci, isEmptyLoci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../mol-geo/geometry/marker-data';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { WebGLContext } from 'mol-gl/webgl/context';
// import { ColorTheme } from 'mol-theme/color';

// export interface RepresentationProps {
//     visuals?: string[]
// }
export type RepresentationProps = { [k: string]: any }

export interface RepresentationContext {
    webgl?: WebGLContext
}

export interface Representation<D, P extends RepresentationProps = {}> {
    readonly label: string
    readonly params: PD.Params
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<P>
    createOrUpdate: (ctx: RepresentationContext, props?: Partial<P>, data?: D) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}

export namespace Representation {
    export function createMulti<D, P extends RepresentationProps = {}>(label: string, params: PD.Params, defaultProps: P, reprList: Representation<D, P>[]): Representation<D, P> {
        let currentProps: P
        let currentData: D

        const visualsOptions: [string, string][] = []
        for (let i = 0, il = reprList.length; i < il; ++i) {
            visualsOptions.push([ i.toString(), reprList[i].label ])
        }
        params['visuals'] = PD.MultiSelectParam<string>('Visuals', '', ['surface'], visualsOptions)

        if (!defaultProps.visuals) {
            defaultProps.visuals = reprList.map((r, i) => i.toString())
        }

        return {
            label,
            params,
            get renderObjects() {
                const { visuals } = currentProps
                const renderObjects: RenderObject[] = []
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    if (!visuals || visuals.includes(i.toString())) {
                        renderObjects.push(...reprList[i].renderObjects)
                    }
                }
                return renderObjects
            },
            get props() {
                const props = {}
                reprList.forEach(r => Object.assign(props, r.props))
                return props as P
            },
            createOrUpdate: (ctx: RepresentationContext, props: Partial<P> = {}, data?: D) => {
                if (data) currentData = data
                // const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
                // currentProps = Object.assign({}, DefaultCartoonProps, currentProps, props, qualityProps)
                currentProps = Object.assign({}, defaultProps, currentProps, props)

                const { visuals } = currentProps
                return Task.create(`Creating '${label}' representation`, async runtime => {
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(i.toString())) {
                            await reprList[i].createOrUpdate(ctx, currentProps, currentData).runInContext(runtime)
                        }
                    }
                })
            },
            getLoci: (pickingId: PickingId) => {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    const loci = reprList[i].getLoci(pickingId)
                    if (!isEmptyLoci(loci)) return loci
                }
                return EmptyLoci
            },
            mark: (loci: Loci, action: MarkerAction) => {
                let marked = false
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    marked = reprList[i].mark(loci, action) || marked
                }
                return marked
            },
            destroy() {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].destroy()
                }
            }
        }
    }
}

//

export interface VisualContext extends RepresentationContext {
    runtime: RuntimeContext,
    // TODO
    // colorTheme: ColorTheme,
}

export interface Visual<D, P extends RepresentationProps> {
    readonly renderObject: RenderObject | undefined
    createOrUpdate: (ctx: VisualContext, props?: Partial<P>, data?: D) => Promise<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}