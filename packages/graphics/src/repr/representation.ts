/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import type { ColorTheme } from '@molstar/graphics/theme/color';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { type ThemeRegistryContext, Theme } from '@molstar/graphics/theme/theme';
import { Subject } from 'rxjs';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Task } from '@molstar/core/task';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { Loci as ModelLoci, EmptyLoci, isEmptyLoci } from '@molstar/model/model/loci';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { Transparency } from '@molstar/graphics/theme/transparency';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { type LocationCallback, getQualityProps } from './util.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { Visual } from './visual.js';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { SetUtils } from '@molstar/core/util/set';
import { cantorPairing } from '@molstar/core/data/util';
import { Substance } from '@molstar/graphics/theme/substance';
import { Emissive } from '@molstar/graphics/theme/emissive';
import { Wiggle } from '@molstar/graphics/theme/wiggle';
import type { Location } from '@molstar/model/model/location';

export type RepresentationProps = { [k: string]: any };

export interface RepresentationContext {
  readonly webgl?: WebGLContext;
  readonly colorThemeRegistry: ColorTheme.Registry;
  readonly sizeThemeRegistry: SizeTheme.Registry;
}

export type RepresentationParamsGetter<D, P extends PD.Params> = (ctx: ThemeRegistryContext, data: D) => P;
export type RepresentationFactory<D, P extends PD.Params, S extends Representation.State> = (
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<D, P>,
) => Representation<D, P, S>;

//

export interface RepresentationProvider<
  D = any,
  P extends PD.Params = any,
  S extends Representation.State = any,
  Id extends string = string,
> {
  readonly name: Id;
  readonly label: string;
  readonly description: string;
  readonly factory: RepresentationFactory<D, P, S>;
  readonly getParams: RepresentationParamsGetter<D, P>;
  readonly defaultValues: PD.Values<P>;
  readonly defaultColorTheme: { name: string; props?: {} };
  readonly defaultSizeTheme: { name: string; props?: {} };
  readonly isApplicable: (data: D) => boolean;
  readonly ensureCustomProperties?: {
    attach: (ctx: CustomProperty.Context, data: D) => Promise<void>;
    detach: (data: D) => void;
  };
  readonly getData?: (data: D, props: PD.Values<P>) => D;
  readonly mustRecreate?: (oldProps: PD.Values<P>, newProps: PD.Values<P>) => boolean;
  readonly locationKinds?: ReadonlyArray<Location['kind']>;
}

export namespace RepresentationProvider {
  export type ParamValues<R extends RepresentationProvider<any, any, any>> =
    R extends RepresentationProvider<any, infer P, any> ? PD.Values<P> : never;

  export function getDefaultParams<R extends RepresentationProvider<D, any, any>, D>(
    r: R,
    ctx: ThemeRegistryContext,
    data: D,
  ) {
    return PD.getDefaultValues(r.getParams(ctx, data));
  }
}

export type AnyRepresentationProvider = RepresentationProvider<any, {}, Representation.State>;

export const EmptyRepresentationProvider: RepresentationProvider = {
  name: '',
  label: '',
  description: '',
  factory: () => Representation.Empty,
  getParams: () => ({}),
  defaultValues: {},
  defaultColorTheme: { name: '' },
  defaultSizeTheme: { name: '' },
  isApplicable: () => true,
};

function getTypes(list: { name: string; provider: RepresentationProvider<any, any, any> }[]) {
  return list.map((e) => [e.name, e.provider.label] as [string, string]);
}

export class RepresentationRegistry<D, S extends Representation.State> {
  private _list: { name: string; provider: RepresentationProvider<D, any, any> }[] = [];
  private _map = new Map<string, RepresentationProvider<D, any, any>>();
  private _name = new Map<RepresentationProvider<D, any, any>, string>();
  private _count = new Map<RepresentationProvider<D, any, any>, number>();

  /** The first registered entry, or `undefined` when the registry is empty. */
  get default(): { name: string; provider: RepresentationProvider<D, any, any> } | undefined {
    return this._list[0];
  }
  get types(): [string, string][] {
    return getTypes(this._list);
  }

  constructor() {}

  /** Returns the message `add(provider)` would throw, or `undefined` when it would succeed. Changes nothing. */
  findConflict(provider: RepresentationProvider<D, any, any>): string | undefined {
    const existing = this._map.get(provider.name);
    if (existing && existing !== provider) {
      return `Representation '${provider.name}' is already registered with a different provider.`;
    }
    return undefined;
  }

  /** Registering the same object again increments a count; a different object under the same name throws. */
  add<P extends PD.Params>(provider: RepresentationProvider<D, P, S>) {
    const conflict = this.findConflict(provider);
    if (conflict) throw new Error(conflict);

    const count = this._count.get(provider);
    if (count !== undefined) {
      this._count.set(provider, count + 1);
      return;
    }

    this._list.push({ name: provider.name, provider });
    this._map.set(provider.name, provider);
    this._name.set(provider, provider.name);
    this._count.set(provider, 1);
  }

  getName(provider: RepresentationProvider<D, any, any>): string {
    if (!this._name.has(provider)) throw new Error(`'${provider.label}' is not a registered represenatation provider.`);
    return this._name.get(provider)!;
  }

  /** Decrements the count and removes the provider at zero. Removing an unknown provider is a no-op. */
  remove(provider: RepresentationProvider<D, any, any>) {
    const count = this._count.get(provider);
    if (count === undefined) return;
    if (count > 1) {
      this._count.set(provider, count - 1);
      return;
    }

    this._count.delete(provider);
    this._map.delete(provider.name);
    this._name.delete(provider);
    const i = this._list.findIndex((e) => e.provider === provider);
    if (i >= 0) this._list.splice(i, 1);
  }

  has(nameOrProvider: string | RepresentationProvider<D, any, any>): boolean {
    return typeof nameOrProvider === 'string' ? this._map.has(nameOrProvider) : this._count.has(nameOrProvider);
  }

  get<P extends PD.Params>(name: string): RepresentationProvider<D, P, S> {
    return this._map.get(name) || EmptyRepresentationProvider;
  }

  get list() {
    return this._list;
  }

  getApplicableList(data: D) {
    return this._list.filter((e) => e.provider.isApplicable(data));
  }

  getApplicableTypes(data: D) {
    return getTypes(this.getApplicableList(data));
  }

  clear() {
    this._list.length = 0;
    this._map.clear();
    this._name.clear();
    this._count.clear();
  }
}

//

export { Representation };
interface Representation<D, P extends PD.Params = PD.Params, S extends Representation.State = Representation.State> {
  readonly label: string;
  readonly updated: Subject<number>;
  /** Number of addressable groups in all visuals of the representation */
  readonly groupCount: number;
  readonly renderObjects: ReadonlyArray<GraphicsRenderObject>;
  readonly geometryVersion: number;
  readonly props: Readonly<PD.Values<P>>;
  readonly params: Readonly<P>;
  readonly state: Readonly<S>;
  readonly theme: Readonly<Theme>;
  createOrUpdate: (props?: Partial<PD.Values<P>>, data?: D) => Task<void>;
  setState: (state: Partial<S>) => void;
  setTheme: (theme: Theme) => void;
  getLoci: (pickingId: PickingId) => ModelLoci;
  getAllLoci: () => ModelLoci[];
  eachLocation: (cb: LocationCallback) => void;
  mark: (loci: ModelLoci, action: MarkerAction) => boolean;
  destroy: () => void;
}
namespace Representation {
  export interface Loci<T extends ModelLoci = ModelLoci> {
    loci: T;
    repr?: Representation.Any;
  }

  export namespace Loci {
    export function areEqual(a: Loci, b: Loci) {
      return a.repr === b.repr && ModelLoci.areEqual(a.loci, b.loci);
    }

    export function isEmpty(a: Loci) {
      return ModelLoci.isEmpty(a.loci);
    }

    export const Empty: Loci = { loci: EmptyLoci };
  }

  export interface State {
    /** Controls if the representation's renderobjects are rendered or not */
    visible: boolean;
    /** A factor applied to alpha value of the representation's renderobjects */
    alphaFactor: number;
    /** Controls if the representation's renderobjects are pickable or not */
    pickable: boolean;
    /** Controls if the representation's renderobjects is rendered in color pass (i.e., not pick and depth) or not */
    colorOnly: boolean;
    /** Overpaint applied to the representation's renderobjects */
    overpaint: Overpaint;
    /** Per group transparency applied to the representation's renderobjects */
    transparency: Transparency;
    /** Per group emissive applied to the representation's renderobjects */
    emissive: Emissive;
    /** Per group material applied to the representation's renderobjects */
    substance: Substance;
    /** Per group wiggle applied to the representation's renderobjects */
    wiggle: Wiggle;
    /** Bit mask of per group clipping applied to the representation's renderobjects */
    clipping: Clipping;
    /** Strength of the representations overpaint, transparency, emissive, substance, wiggle */
    themeStrength: { overpaint: number; transparency: number; emissive: number; substance: number; wiggle: number };
    /** Controls if the representation's renderobjects are synced automatically with GPU or not */
    syncManually: boolean;
    /** A transformation applied to the representation's renderobjects */
    transform: Mat4;
    /** Bit mask of allowed marker actions */
    markerActions: MarkerActions;
  }
  export function createState(): State {
    return {
      visible: true,
      alphaFactor: 1,
      pickable: true,
      colorOnly: false,
      syncManually: false,
      transform: Mat4.identity(),
      overpaint: Overpaint.Empty,
      transparency: Transparency.Empty,
      emissive: Emissive.Empty,
      substance: Substance.Empty,
      wiggle: Wiggle.Empty,
      clipping: Clipping.Empty,
      themeStrength: { overpaint: 1, transparency: 1, emissive: 1, substance: 1, wiggle: 1 },
      markerActions: MarkerActions.All,
    };
  }
  export function updateState(state: State, update: Partial<State>) {
    if (update.visible !== undefined) state.visible = update.visible;
    if (update.alphaFactor !== undefined) state.alphaFactor = update.alphaFactor;
    if (update.pickable !== undefined) state.pickable = update.pickable;
    if (update.colorOnly !== undefined) state.colorOnly = update.colorOnly;
    if (update.overpaint !== undefined) state.overpaint = update.overpaint;
    if (update.transparency !== undefined) state.transparency = update.transparency;
    if (update.emissive !== undefined) state.emissive = update.emissive;
    if (update.substance !== undefined) state.substance = update.substance;
    if (update.wiggle !== undefined) state.wiggle = update.wiggle;
    if (update.clipping !== undefined) state.clipping = update.clipping;
    if (update.themeStrength !== undefined) state.themeStrength = update.themeStrength;
    if (update.syncManually !== undefined) state.syncManually = update.syncManually;
    if (update.transform !== undefined) Mat4.copy(state.transform, update.transform);
    if (update.markerActions !== undefined) state.markerActions = update.markerActions;
  }
  export interface StateBuilder<S extends State> {
    create(): S;
    update(state: S, update: Partial<S>): void;
  }
  export const StateBuilder: StateBuilder<State> = { create: createState, update: updateState };

  export type Any<P extends PD.Params = PD.Params, S extends State = State> = Representation<any, P, S>;

  export declare const Empty: Any;

  export function createEmpty(): Any {
    return {
      label: '',
      groupCount: 0,
      renderObjects: [],
      geometryVersion: -1,
      props: {},
      params: {},
      updated: new Subject(),
      state: createState(),
      theme: Theme.createEmpty(),
      createOrUpdate: () => Task.constant('', undefined),
      setState: () => {},
      setTheme: () => {},
      getLoci: () => EmptyLoci,
      getAllLoci: () => [],
      eachLocation: () => {},
      mark: () => false,
      destroy: () => {},
    };
  }

  export type Def<D, P extends PD.Params = PD.Params, S extends State = State> = {
    [k: string]: RepresentationFactory<D, P, S>;
  };

  export class GeometryState {
    private curr = new Set<number>();
    private next = new Set<number>();

    private _version = -1;
    get version() {
      return this._version;
    }

    add(id: number, version: number) {
      this.next.add(cantorPairing(id, version));
    }

    snapshot() {
      if (!SetUtils.areEqual(this.curr, this.next)) {
        this._version += 1;
      }
      [this.curr, this.next] = [this.next, this.curr];
      this.next.clear();
    }
  }

  export function createMulti<D, P extends PD.Params = PD.Params, S extends State = State>(
    label: string,
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<D, P>,
    stateBuilder: StateBuilder<S>,
    reprDefs: Def<D, P>,
  ): Representation<D, P, S> {
    let version = 0;
    const updated = new Subject<number>();
    const geometryState = new GeometryState();
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
      get geometryVersion() {
        return geometryState.version;
      },
      get props() {
        return currentProps;
      },
      get params() {
        return currentParams;
      },
      createOrUpdate: (props: Partial<P> = {}, data?: D) => {
        if (data && data !== currentData) {
          currentParams = getParams(ctx, data);
          currentData = data;
          if (!currentProps) currentProps = PD.getDefaultValues(currentParams) as P;
        }
        const qualityProps = getQualityProps(Object.assign({}, currentProps, props), currentData);
        Object.assign(currentProps, props, qualityProps);

        const { visuals } = currentProps;
        return Task.create(`Creating or updating '${label}' representation`, async (runtime) => {
          for (let i = 0, il = reprList.length; i < il; ++i) {
            if (!visuals || visuals.includes(reprMap[i])) {
              await reprList[i].createOrUpdate(currentProps, currentData).runInContext(runtime);
            }
            geometryState.add(i, reprList[i].geometryVersion);
          }
          geometryState.snapshot();
          updated.next(version++);
        });
      },
      get state() {
        return currentState;
      },
      get theme() {
        return currentTheme;
      },
      getLoci: (pickingId: PickingId) => {
        const { visuals } = currentProps;
        for (let i = 0, il = reprList.length; i < il; ++i) {
          if (!visuals || visuals.includes(reprMap[i])) {
            const loci = reprList[i].getLoci(pickingId);
            if (!isEmptyLoci(loci)) return loci;
          }
        }
        return EmptyLoci;
      },
      getAllLoci: () => {
        const loci: ModelLoci[] = [];
        const { visuals } = currentProps;
        for (let i = 0, il = reprList.length; i < il; ++i) {
          if (!visuals || visuals.includes(reprMap[i])) {
            loci.push(...reprList[i].getAllLoci());
          }
        }
        return loci;
      },
      eachLocation: (cb: LocationCallback) => {
        const { visuals } = currentProps;
        for (let i = 0, il = reprList.length; i < il; ++i) {
          if (!visuals || visuals.includes(reprMap[i])) {
            reprList[i].eachLocation(cb);
          }
        }
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
          reprList[i].setState(state); // only set the new (partial) state
        }
      },
      setTheme: (theme: Theme) => {
        currentTheme = theme;
        for (let i = 0, il = reprList.length; i < il; ++i) {
          reprList[i].setTheme(theme);
        }
      },
      destroy() {
        for (let i = 0, il = reprList.length; i < il; ++i) {
          reprList[i].destroy();
        }
      },
    };
  }

  export function fromRenderObject(
    label: string,
    renderObject: GraphicsRenderObject,
  ): Representation<GraphicsRenderObject, BaseGeometry.Params> {
    let version = 0;
    const updated = new Subject<number>();
    const geometryState = new GeometryState();
    const currentState = Representation.createState();
    const currentTheme = Theme.createEmpty();

    const currentParams = PD.clone(BaseGeometry.Params);
    const currentProps = PD.getDefaultValues(BaseGeometry.Params);

    return {
      label,
      updated,
      get groupCount() {
        return renderObject.values.uGroupCount.ref.value;
      },
      get renderObjects() {
        return [renderObject];
      },
      get geometryVersion() {
        return geometryState.version;
      },
      get props() {
        return currentProps;
      },
      get params() {
        return currentParams;
      },
      createOrUpdate: (props: Partial<PD.Values<BaseGeometry.Params>> = {}) => {
        const qualityProps = getQualityProps(Object.assign({}, currentProps, props));
        Object.assign(currentProps, props, qualityProps);

        return Task.create(`Updating '${label}' representation`, async (runtime) => {
          // TODO
          geometryState.add(0, renderObject.id);
          geometryState.snapshot();
          updated.next(version++);
        });
      },
      get state() {
        return currentState;
      },
      get theme() {
        return currentTheme;
      },
      getLoci: () => {
        // TODO
        return EmptyLoci;
      },
      getAllLoci: () => {
        // TODO
        return [];
      },
      eachLocation: () => {
        // TODO
      },
      mark: (loci: ModelLoci, action: MarkerAction) => {
        // TODO
        return false;
      },
      setState: (state: Partial<State>) => {
        if (state.visible !== undefined) Visual.setVisibility(renderObject, state.visible);
        if (state.alphaFactor !== undefined) Visual.setAlphaFactor(renderObject, state.alphaFactor);
        if (state.pickable !== undefined) Visual.setPickable(renderObject, state.pickable);
        if (state.colorOnly !== undefined) Visual.setColorOnly(renderObject, state.colorOnly);
        if (state.overpaint !== undefined) {
          // TODO
        }
        if (state.transparency !== undefined) {
          // TODO
        }
        if (state.emissive !== undefined) {
          // TODO
        }
        if (state.substance !== undefined) {
          // TODO
        }
        if (state.wiggle !== undefined) {
          // TODO
        }
        if (state.clipping !== undefined) {
          // TODO
        }
        if (state.themeStrength !== undefined) Visual.setThemeStrength(renderObject, state.themeStrength);
        if (state.transform !== undefined) Visual.setTransform(renderObject, state.transform);

        Representation.updateState(currentState, state);
      },
      setTheme: () => {},
      destroy() {},
    };
  }
}

let _EmptyRepresentation: Representation.Any | undefined = undefined;
Object.defineProperty(Representation, 'Empty', {
  get: () => {
    return (_EmptyRepresentation ??= Representation.createEmpty());
  },
});
