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
import { Theme, ThemeRegistryContext } from 'mol-theme/theme';
import { Subject } from 'rxjs';

// export interface RepresentationProps {
//     visuals?: string[]
// }
export type RepresentationProps = { [k: string]: any }

export type RepresentationParamsGetter<D, P extends PD.Params> = (ctx: ThemeRegistryContext, data: D) => P

//

export interface RepresentationProvider<D, P extends PD.Params> {
    readonly label: string
    readonly description: string
    readonly factory: (getParams: RepresentationParamsGetter<D, P>) => Representation<D, P>
    readonly getParams: (ctx: ThemeRegistryContext, data: D) => P
    readonly defaultValues: PD.Values<P>
}

export type AnyRepresentationProvider = RepresentationProvider<any, {}>

export const EmptyRepresentationProvider = {
    label: '',
    description: '',
    factory: () => Representation.Empty,
    getParams: () => ({}),
    defaultValues: {}
}

export class RepresentationRegistry<D> {
    private _list: { name: string, provider: RepresentationProvider<D, any> }[] = []
    private _map = new Map<string, RepresentationProvider<D, any>>()

    get default() { return this._list[0]; }
    get types(): [string, string][] {
        return this._list.map(e => [e.name, e.provider.label] as [string, string]);
    }

    constructor() {};

    add<P extends PD.Params>(name: string, provider: RepresentationProvider<D, P>) {
        this._list.push({ name, provider })
        this._map.set(name, provider)
    }

    get<P extends PD.Params>(name: string): RepresentationProvider<D, P> {
        return this._map.get(name) || EmptyRepresentationProvider as unknown as RepresentationProvider<D, P>
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
interface Representation<D, P extends PD.Params = {}> {
    readonly label: string
    readonly updated: Subject<number>
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<PD.Values<P>>
    readonly params: Readonly<P>
    createOrUpdate: (ctx: RepresentationContext, props?: Partial<PD.Values<P>>, data?: D) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    setVisibility: (value: boolean) => void
    setPickable: (value: boolean) => void
    destroy: () => void
}
namespace Representation {
    export type Any = Representation<any>
    export const Empty: Any = {
        label: '', renderObjects: [], props: {}, params: {}, updated: new Subject(),
        createOrUpdate: () => Task.constant('', undefined),
        getLoci: () => EmptyLoci,
        mark: () => false,
        setVisibility: () => {},
        setPickable: () => {},
        destroy: () => {}
    }

    export type Def<D, P extends PD.Params = {}> = { [k: string]: (getParams: RepresentationParamsGetter<D, P>) => Representation<any, P> }

    export function createMulti<D, P extends PD.Params = {}>(label: string, getParams: RepresentationParamsGetter<D, P>, reprDefs: Def<D, P>): Representation<D, P> {
        let version = 0
        const updated = new Subject<number>()

        let currentParams: P
        let currentProps: PD.Values<P>
        let currentData: D

        const reprMap: { [k: number]: string } = {}
        const reprList: Representation<D, P>[] = Object.keys(reprDefs).map((name, i) => {
            reprMap[i] = name
            return reprDefs[name](getParams)
        })

        return {
            label,
            updated,
            get renderObjects() {
                const renderObjects: RenderObject[] = []
                if (currentProps) {
                    const { visuals } = currentProps
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            renderObjects.push(...reprList[i].renderObjects)
                        }
                    }
                }
                return renderObjects
            },
            get props() {
                const props = {}
                reprList.forEach(r => Object.assign(props, r.props))
                return props as P
            },
            get params() { return currentParams },
            createOrUpdate: (ctx: RepresentationContext, props: Partial<P> = {}, data?: D) => {
                if (data && data !== currentData) {
                    currentParams = getParams(ctx, data)
                    currentData = data
                    if (!currentProps) currentProps = PD.getDefaultValues(currentParams) as P
                }
                const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentData)
                Object.assign(currentProps, props, qualityProps)

                const { visuals } = currentProps
                return Task.create(`Creating '${label}' representation`, async runtime => {
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            await reprList[i].createOrUpdate(ctx, currentProps, currentData).runInContext(runtime)
                        }
                    }
                    updated.next(version++)
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
            setVisibility: (value: boolean) => {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].setVisibility(value)
                }
            },
            setPickable: (value: boolean) => {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].setPickable(value)
                }
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

export interface Visual<D, P extends PD.Params> {
    readonly renderObject: RenderObject | undefined
    createOrUpdate: (ctx: VisualContext, theme: Theme, props?: Partial<PD.Values<P>>, data?: D) => Promise<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => boolean
    setVisibility: (value: boolean) => void
    setPickable: (value: boolean) => void
    destroy: () => void
}