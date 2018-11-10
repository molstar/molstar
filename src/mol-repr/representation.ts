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
import { getQualityProps } from './util';
import { ColorTheme } from 'mol-theme/color';
import { SizeTheme } from 'mol-theme/size';
import { ThemeProps, Theme, ThemeRegistryContext } from 'mol-theme/theme';

// export interface RepresentationProps {
//     visuals?: string[]
// }
export type RepresentationProps = { [k: string]: any }

//

export interface RepresentationProvider<D, P extends PD.Params> {
    readonly factory: (defaultProps: PD.DefaultValues<P>) => Representation<D, PD.DefaultValues<P>>
    readonly params: (ctx: ThemeRegistryContext, data: D) => P
}

export class RepresentationRegistry<D> {
    private _list: { name: string, provider: RepresentationProvider<D, any> }[] = []
    private _map = new Map<string, RepresentationProvider<D, any>>()

    constructor() {};

    add<P extends PD.Params>(name: string, factory: RepresentationProvider<D, P>['factory'], params: RepresentationProvider<D, P>['params']) {
        const provider = { factory, params } as RepresentationProvider<D, P>
        this._list.push({ name, provider })
        this._map.set(name, provider)
    }

    get(id: string) {
        return this._map.get(id)
    }

    create(id: string, ctx: ThemeRegistryContext, data: D, props = {}): Representation<D, any> {
        const provider = this.get(id)
        return provider ? provider.factory({ ...PD.getDefaultValues(provider.params(ctx, data)), ...props }) : Representation.Empty
    }

    get list() {
        return this._list
    }
}

//

export interface RepresentationContext {
    webgl?: WebGLContext
    colorThemeRegistry: ColorTheme.Registry
    sizeThemeRegistry: SizeTheme.Registry
}

export { Representation }
interface Representation<D, P extends RepresentationProps = {}> {
    readonly label: string
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<P>
    createOrUpdate: (ctx: RepresentationContext, props?: Partial<P>, themeProps?: ThemeProps, data?: D) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}
namespace Representation {
    export type Any = Representation<any>
    export const Empty: Representation<any> = {
        label: '', renderObjects: [], props: {},
        createOrUpdate: () => Task.constant('', undefined),
        getLoci: () => EmptyLoci,
        mark: () => false,
        destroy: () => {}
    }

    export type Def<P extends RepresentationProps = {}> = { [k: string]: (defaultProps: P) => Representation<any, P> }

    export function createMulti<D, P extends RepresentationProps = {}>(label: string, defaultProps: P, reprDefs: Def<P>): Representation<D, P> {
        let currentProps: P = Object.assign({}, defaultProps)
        let currentData: D

        const reprMap: { [k: number]: string } = {}
        const reprList: Representation<D, P>[] = Object.keys(reprDefs).map((name, i) => {
            reprMap[i] = name
            return reprDefs[name](defaultProps)
        })

        return {
            label,
            get renderObjects() {
                const { visuals } = currentProps
                const renderObjects: RenderObject[] = []
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    if (!visuals || visuals.includes(reprMap[i])) {
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
            createOrUpdate: (ctx: RepresentationContext, props: Partial<P> = {}, themeProps: ThemeProps = {}, data?: D) => {
                if (data) currentData = data
                const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentData)
                Object.assign(currentProps, props, qualityProps)

                const { visuals } = currentProps
                return Task.create(`Creating '${label}' representation`, async runtime => {
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            await reprList[i].createOrUpdate(ctx, currentProps, themeProps, currentData).runInContext(runtime)
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

export interface VisualContext {
    webgl?: WebGLContext
    runtime: RuntimeContext,
}

export interface Visual<D, P extends RepresentationProps> {
    readonly renderObject: RenderObject | undefined
    createOrUpdate: (ctx: VisualContext, theme: Theme, props?: Partial<P>, data?: D) => Promise<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    destroy: () => void
}