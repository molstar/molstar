/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../mol-util/param-definition';
import { WebGLContext } from '../mol-gl/webgl/context';
import { ColorTheme } from '../mol-theme/color';
import { SizeTheme } from '../mol-theme/size';
import { ThemeRegistryContext, Theme } from '../mol-theme/theme';
import { Subject } from 'rxjs';
import { GraphicsRenderObject } from '../mol-gl/render-object';
import { Task } from '../mol-task';
import { PickingId } from '../mol-geo/geometry/picking';
import { MarkerAction, MarkerActions } from '../mol-util/marker-action';
import { Loci as ModelLoci, EmptyLoci, isEmptyLoci } from '../mol-model/loci';
import { Overpaint } from '../mol-theme/overpaint';
import { Transparency } from '../mol-theme/transparency';
import { Mat4 } from '../mol-math/linear-algebra';
import { getQualityProps } from './util';
import { BaseGeometry } from '../mol-geo/geometry/base';
import { Visual } from './visual';
import { CustomProperty } from '../mol-model-props/common/custom-property';
import { Clipping } from '../mol-theme/clipping';

// export interface RepresentationProps {
//     visuals?: string[]
// }
export type RepresentationProps = { [k: string]: any }

export interface RepresentationContext {
    readonly webgl?: WebGLContext
    readonly colorThemeRegistry: ColorTheme.Registry
    readonly sizeThemeRegistry: SizeTheme.Registry
}

export type RepresentationParamsGetter<D, P extends PD.Params> = (ctx: ThemeRegistryContext, data: D) => P
export type RepresentationFactory<D, P extends PD.Params, S extends Representation.State> = (ctx: RepresentationContext, getParams: RepresentationParamsGetter<D, P>) => Representation<D, P, S>

//

export interface RepresentationProvider<D = any, P extends PD.Params = any, S extends Representation.State = any, Id extends string = string> {
    readonly name: Id,
    readonly label: string
    readonly description: string
    readonly factory: RepresentationFactory<D, P, S>
    readonly getParams: RepresentationParamsGetter<D, P>
    readonly defaultValues: PD.Values<P>
    readonly defaultColorTheme: { name: string, props?: {} }
    readonly defaultSizeTheme: { name: string, props?: {} }
    readonly isApplicable: (data: D) => boolean
    readonly ensureCustomProperties?: {
        attach: (ctx: CustomProperty.Context, data: D) => Promise<void>,
        detach: (data: D) => void
    }
}

export namespace RepresentationProvider {
    export type ParamValues<R extends RepresentationProvider<any, any, any>> = R extends RepresentationProvider<any, infer P, any> ? PD.Values<P> : never;

    export function getDetaultParams<R extends RepresentationProvider<D, any, any>, D>(r: R, ctx: ThemeRegistryContext, data: D) {
        return PD.getDefaultValues(r.getParams(ctx, data));
    }
}

export type AnyRepresentationProvider = RepresentationProvider<any, {}, Representation.State>

export const EmptyRepresentationProvider = {
    label: '',
    description: '',
    factory: () => Representation.Empty,
    getParams: () => ({}),
    defaultValues: {}
};

function getTypes(list: { name: string, provider: RepresentationProvider<any, any, any> }[]) {
    return list.map(e => [e.name, e.provider.label] as [string, string]);
}

export class RepresentationRegistry<D, S extends Representation.State> {
    private _list: { name: string, provider: RepresentationProvider<D, any, any> }[] = []
    private _map = new Map<string, RepresentationProvider<D, any, any>>()
    private _name = new Map<RepresentationProvider<D, any, any>, string>()

    get default() { return this._list[0]; }
    get types(): [string, string][] { return getTypes(this._list); }

    constructor() {};

    add<P extends PD.Params>(provider: RepresentationProvider<D, P, S>) {
        if (this._map.has(provider.name)) {
            throw new Error(`${provider.name} already registered.`);
        }

        this._list.push({ name: provider.name, provider });
        this._map.set(provider.name, provider);
        this._name.set(provider, provider.name);
    }

    getName(provider: RepresentationProvider<D, any, any>): string {
        if (!this._name.has(provider)) throw new Error(`'${provider.label}' is not a registered represenatation provider.`);
        return this._name.get(provider)!;
    }

    remove(provider: RepresentationProvider<D, any, any>) {
        const name = provider.name;

        this._list.splice(this._list.findIndex(e => e.name === name), 1);
        const p = this._map.get(name);
        if (p) {
            this._map.delete(name);
            this._name.delete(p);
        }
    }

    get<P extends PD.Params>(name: string): RepresentationProvider<D, P, S> {
        return this._map.get(name) || EmptyRepresentationProvider as unknown as RepresentationProvider<D, P, S>;
    }

    get list() {
        return this._list;
    }

    getApplicableList(data: D) {
        return this._list.filter(e => e.provider.isApplicable(data));
    }

    getApplicableTypes(data: D) {
        return getTypes(this.getApplicableList(data));
    }
}

//

export { Representation };
interface Representation<D, P extends PD.Params = {}, S extends Representation.State = Representation.State> {
    readonly label: string
    readonly updated: Subject<number>
    /** Number of addressable groups in all visuals of the representation */
    readonly groupCount: number
    readonly renderObjects: ReadonlyArray<GraphicsRenderObject>
    readonly props: Readonly<PD.Values<P>>
    readonly params: Readonly<P>
    readonly state: Readonly<S>
    readonly theme: Readonly<Theme>
    createOrUpdate: (props?: Partial<PD.Values<P>>, data?: D) => Task<void>
    setState: (state: Partial<S>) => void
    setTheme: (theme: Theme) => void
    /** If no pickingId is given, returns a Loci for the whole representation */
    getLoci: (pickingId?: PickingId) => ModelLoci
    mark: (loci: ModelLoci, action: MarkerAction) => boolean
    destroy: () => void
}
namespace Representation {
    export interface Loci<T extends ModelLoci = ModelLoci> { loci: T, repr?: Representation.Any }

    export namespace Loci {
        export function areEqual(a: Loci, b: Loci) {
            return a.repr === b.repr && ModelLoci.areEqual(a.loci, b.loci);
        }

        export const Empty: Loci = { loci: EmptyLoci };
    }

    export interface State {
        /** Controls if the representation's renderobjects are rendered or not */
        visible: boolean
        /** A factor applied to alpha value of the representation's renderobjects */
        alphaFactor: number
        /** Controls if the representation's renderobjects are pickable or not */
        pickable: boolean
        /** Overpaint applied to the representation's renderobjects */
        overpaint: Overpaint
        /** Per group transparency applied to the representation's renderobjects */
        transparency: Transparency
        /** Bit mask of per group clipping applied to the representation's renderobjects */
        clipping: Clipping
        /** Controls if the representation's renderobjects are synced automatically with GPU or not */
        syncManually: boolean
        /** A transformation applied to the representation's renderobjects */
        transform: Mat4
        /** Bit mask of allowed marker actions */
        markerActions: MarkerActions
    }
    export function createState(): State {
        return { visible: true, alphaFactor: 1, pickable: true, syncManually: false, transform: Mat4.identity(), overpaint: Overpaint.Empty, transparency: Transparency.Empty, clipping: Clipping.Empty, markerActions: MarkerActions.All };
    }
    export function updateState(state: State, update: Partial<State>) {
        if (update.visible !== undefined) state.visible = update.visible;
        if (update.alphaFactor !== undefined) state.alphaFactor = update.alphaFactor;
        if (update.pickable !== undefined) state.pickable = update.pickable;
        if (update.overpaint !== undefined) state.overpaint = update.overpaint;
        if (update.transparency !== undefined) state.transparency = update.transparency;
        if (update.clipping !== undefined) state.clipping = update.clipping;
        if (update.syncManually !== undefined) state.syncManually = update.syncManually;
        if (update.transform !== undefined) Mat4.copy(state.transform, update.transform);
        if (update.markerActions !== undefined) state.markerActions = update.markerActions;
    }
    export interface StateBuilder<S extends State> {
        create(): S
        update(state: S, update: Partial<S>): void
    }
    export const StateBuilder: StateBuilder<State> = { create: createState, update: updateState };

    export type Any = Representation<any, any, any>
    export const Empty: Any = {
        label: '', groupCount: 0, renderObjects: [], props: {}, params: {}, updated: new Subject(), state: createState(), theme: Theme.createEmpty(),
        createOrUpdate: () => Task.constant('', undefined),
        setState: () => {},
        setTheme: () => {},
        getLoci: () => EmptyLoci,
        mark: () => false,
        destroy: () => {}
    };

    export type Def<D, P extends PD.Params = {}, S extends State = State> = { [k: string]: RepresentationFactory<D, P, S> }

    export function createMulti<D, P extends PD.Params = {}, S extends State = State>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<D, P>, stateBuilder: StateBuilder<S>, reprDefs: Def<D, P>): Representation<D, P, S> {
        let version = 0;
        const updated = new Subject<number>();
        const currentState = stateBuilder.create();
        let currentTheme = Theme.createEmpty();

        let currentParams: P;
        let currentProps: PD.Values<P>;
        let currentData: D;

        const reprMap: { [k: number]: string } = {};
        const reprList: Representation<D, P>[] = Object.keys(reprDefs).map((name, i) => {
            reprMap[i] = name;
            const repr = reprDefs[name](ctx, getParams);
            repr.setState(currentState);
            return repr;
        });

        return {
            label,
            updated,
            get groupCount() {
                let groupCount = 0;
                if (currentProps) {
                    const { visuals } = currentProps;
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            groupCount += reprList[i].groupCount;
                        }
                    }
                }
                return groupCount;
            },
            get renderObjects() {
                const renderObjects: GraphicsRenderObject[] = [];
                if (currentProps) {
                    const { visuals } = currentProps;
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            renderObjects.push(...reprList[i].renderObjects);
                        }
                    }
                }
                return renderObjects;
            },
            get props() { return currentProps; },
            get params() { return currentParams; },
            createOrUpdate: (props: Partial<P> = {}, data?: D) => {
                if (data && data !== currentData) {
                    currentParams = getParams(ctx, data);
                    currentData = data;
                    if (!currentProps) currentProps = PD.getDefaultValues(currentParams) as P;
                }
                const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentData);
                Object.assign(currentProps, props, qualityProps);

                const { visuals } = currentProps;
                return Task.create(`Creating or updating '${label}' representation`, async runtime => {
                    for (let i = 0, il = reprList.length; i < il; ++i) {
                        if (!visuals || visuals.includes(reprMap[i])) {
                            await reprList[i].createOrUpdate(currentProps, currentData).runInContext(runtime);
                        }
                    }
                    updated.next(version++);
                });
            },
            get state() { return currentState; },
            get theme() { return currentTheme; },
            getLoci: (pickingId?: PickingId) => {
                const { visuals } = currentProps;
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    if (!visuals || visuals.includes(reprMap[i])) {
                        const loci = reprList[i].getLoci(pickingId);
                        if (!isEmptyLoci(loci)) return loci;
                    }
                }
                return EmptyLoci;
            },
            mark: (loci: ModelLoci, action: MarkerAction) => {
                let marked = false;
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    marked = reprList[i].mark(loci, action) || marked;
                }
                return marked;
            },
            setState: (state: Partial<S>) => {
                stateBuilder.update(currentState, state);
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].setState(currentState);
                }
            },
            setTheme: (theme: Theme) => {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].setTheme(theme);
                }
            },
            destroy() {
                for (let i = 0, il = reprList.length; i < il; ++i) {
                    reprList[i].destroy();
                }
            }
        };
    }

    export function fromRenderObject(label: string, renderObject: GraphicsRenderObject): Representation<GraphicsRenderObject, BaseGeometry.Params> {
        let version = 0;
        const updated = new Subject<number>();
        const currentState = Representation.createState();
        const currentTheme = Theme.createEmpty();

        const currentParams = PD.clone(BaseGeometry.Params);
        const currentProps = PD.getDefaultValues(BaseGeometry.Params);

        return {
            label,
            updated,
            get groupCount() { return renderObject.values.uGroupCount.ref.value; },
            get renderObjects() { return [renderObject]; },
            get props() { return currentProps; },
            get params() { return currentParams; },
            createOrUpdate: (props: Partial<PD.Values<BaseGeometry.Params>> = {}) => {
                const qualityProps = getQualityProps(Object.assign({}, currentProps, props));
                Object.assign(currentProps, props, qualityProps);

                return Task.create(`Updating '${label}' representation`, async runtime => {
                    // TODO
                    updated.next(version++);
                });
            },
            get state() { return currentState; },
            get theme() { return currentTheme; },
            getLoci: () => {
                // TODO
                return EmptyLoci;
            },
            mark: (loci: ModelLoci, action: MarkerAction) => {
                // TODO
                return false;
            },
            setState: (state: Partial<State>) => {
                if (state.visible !== undefined) Visual.setVisibility(renderObject, state.visible);
                if (state.alphaFactor !== undefined) Visual.setAlphaFactor(renderObject, state.alphaFactor);
                if (state.pickable !== undefined) Visual.setPickable(renderObject, state.pickable);
                if (state.overpaint !== undefined) {
                    // TODO
                }
                if (state.transparency !== undefined) {
                    // TODO
                }
                if (state.transform !== undefined) Visual.setTransform(renderObject, state.transform);

                Representation.updateState(currentState, state);
            },
            setTheme: () => { },
            destroy() { }
        };
    }
}