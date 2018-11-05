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
import { Params, MultiSelectParam } from 'mol-util/parameter';

// export interface RepresentationProps {
//     visuals?: string[]
// }
export type RepresentationProps = { [k: string]: any }

export interface Representation<D, P extends RepresentationProps = {}> {
    readonly label: string
    readonly params: Params
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<P>
    createOrUpdate: (props?: Partial<P>, data?: D) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}

export namespace Representation {
    export function createMulti<D, P extends RepresentationProps = {}>(label: string, params: Params, defaultProps: P, reprList: Representation<D, P>[]): Representation<D, P> {
        let currentProps: P
        let currentData: D

        const visualsOptions: [string, string][] = []
        for (let i = 0, il = reprList.length; i < il; ++i) {
            visualsOptions.push([ i.toString(), reprList[i].label ])
        }
        params['visuals'] = MultiSelectParam<string>('Visuals', '', ['surface'], visualsOptions)

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
            createOrUpdate: (props: Partial<P> = {}, data?: D) => {
                if (data) currentData = data
                // const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
                // currentProps = Object.assign({}, DefaultCartoonProps, currentProps, props, qualityProps)
                currentProps = Object.assign({}, defaultProps, currentProps, props)

                const { visuals } = currentProps
                return Task.create(`Creating '${label}' representation`, async ctx => {
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(i.toString())) {
                            await reprList[i].createOrUpdate(currentProps, currentData).runInContext(ctx)
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

export interface Visual<D, P extends RepresentationProps> {
    readonly renderObject: RenderObject | undefined
    createOrUpdate: (ctx: RuntimeContext, props?: Partial<P>, data?: D) => Promise<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}