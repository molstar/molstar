/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * After the proof of concept by Paul Pillot (molstar/molstar#1939).
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Model } from '../../../model.js';
import type { Structure } from '../../structure.js';
import type { Unit } from '../../unit.js';
import type { BondComputationProps } from './common.js';
import type { InterUnitBonds, IntraUnitBonds } from './data.js';
import type { InterBondComputationProps } from './inter-compute.js';

/**
 * A kind of bonds: named, with params, making instances (bond-providers plan §2). The state tree stores the name and
 * the params; the instance is what units and structures hold.
 */
export interface BondProvider<P extends PD.Params = PD.Params> {
  /** The registry key, and what the state tree stores. */
  readonly name: string;
  readonly label: string;
  readonly description?: string;
  /**
   * The params for a model; called without one when a param definition is built without data, and must not throw.
   * `ctx` is there when a registry builds them, for providers that compose others.
   */
  getParams(model: Model | undefined, ctx?: BondProvider.Context): P;
  /** A provider not applicable to a model is not offered for it, and a stored choice of it falls back to the default. */
  isApplicable?(model: Model): boolean;
  create(props: PD.Values<P>, ctx?: BondProvider.Context): BondProvider.Instance;
}

export namespace BondProvider {
  /**
   * What a `Unit.Atomic` and a `Structure` hold, fixed at their construction: two reads of `unit.bonds` return the same
   * object. Not bound to a model or a structure, so it survives `remapModel` (bond-providers plan §2).
   */
  export interface Instance {
    readonly provider: BondProvider;
    /**
     * Equal keys give equal bonds for an equal model and element set: the key of the model's intra-unit bond cache
     * (`ElementSetIntraBondCache`). `undefined` caches nothing.
     */
    readonly cacheKey: string | undefined;
    /**
     * How a child unit gets its bonds when its parent has none yet (unit.md §4): `'self'` (the default) computes the
     * child's own, right for decisions per atom pair; `'root'` computes the root unit's and restricts them, right for
     * decisions over rings or conjugated systems, at the cost of the whole root unit once.
     */
    readonly childBonds?: 'self' | 'root';
    /** The bonds of an atomic unit; without it, the built-in computation. */
    intra?(unit: Unit.Atomic, ctx: IntraContext): IntraUnitBonds;
    /** The bonds between the units of a structure; without it, the built-in computation. */
    inter?(structure: Structure, ctx: InterContext): InterUnitBonds;
  }

  export interface IntraContext {
    /** The built-in computation for the unit (`computeIntraUnitBonds`). */
    builtIn(props?: Partial<BondComputationProps>): IntraUnitBonds;
  }

  export interface InterContext {
    /** The built-in computation for the structure, with the props the structure would use (`computeInterUnitBonds`). */
    builtIn(props?: Partial<InterBondComputationProps>): InterUnitBonds;
  }

  /** A choice of provider as the state tree stores it. */
  export interface Choice {
    name: string;
    params: any;
  }

  /** What a registry offers the providers it builds, so a provider can compose the others (`'per-unit'`). */
  export interface Context {
    /** The select over the registry's providers for a model, as the structure's `bonds` param, without `exclude`. */
    choiceParam(model: Model | undefined, options?: { exclude?: readonly string[] }): PD.Mapped<PD.NamedParams<any>>;
    /** The instance of a choice; the default instance for none, an unknown or a non-applicable one. */
    resolve(choice: Choice | undefined, model?: Model): Instance;
  }

  export const DefaultName = 'default';

  /** The built-in bonds: the default and the fallback. Never registered: it is always there. */
  export const Default: BondProvider<{}> = {
    name: DefaultName,
    label: 'Default',
    description: 'Bonds from the file and computed from distances and the component dictionary',
    getParams: () => ({}),
    create: () => DefaultInstance,
  };

  /** What every unit and structure gets unless told otherwise; its intra-unit cache is the model's existing one. */
  export const DefaultInstance: Instance = { provider: Default, cacheKey: DefaultName };

  export function isDefault(instance: Instance) {
    return instance === DefaultInstance;
  }
}
