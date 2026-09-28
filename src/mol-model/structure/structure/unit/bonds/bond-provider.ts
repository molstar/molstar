/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import type { Model } from '../../../model/model';
import type { Unit } from '../../unit';
import type { IntraUnitBonds } from './data';

/**
 * A model-scoped source of atomic intra-unit bonds.
 *
 * The selected provider owns the full computation. It can delegate to the
 * exported default computation, replace it, or build a graph from scratch.
 */
export interface BondProvider {
    readonly name: string
    getBonds(unit: Unit.Atomic): IntraUnitBonds | undefined
}

export class BondProviderRegistry {
    private static readonly PropertyName = '__BondProviderRegistry__';

    static get(model: Model): BondProviderRegistry {
        let registry = model._dynamicPropertyData[BondProviderRegistry.PropertyName] as BondProviderRegistry | undefined;
        if (!registry) {
            registry = new BondProviderRegistry();
            model._dynamicPropertyData[BondProviderRegistry.PropertyName] = registry;
        }
        return registry;
    }

    private readonly providers = new Map<string, BondProvider>();
    private activeName: string | undefined;

    get list(): ReadonlyArray<BondProvider> {
        return Array.from(this.providers.values());
    }

    get active(): BondProvider | undefined {
        return this.activeName ? this.providers.get(this.activeName) : undefined;
    }

    add(provider: BondProvider) {
        if (this.providers.has(provider.name)) return;
        this.providers.set(provider.name, provider);
    }

    remove(provider: BondProvider) {
        if (this.providers.get(provider.name) !== provider) return false;
        if (this.activeName === provider.name) this.select(undefined);
        this.providers.delete(provider.name);
        return true;
    }

    get(name: string): BondProvider | undefined {
        return this.providers.get(name);
    }

    has(name: string): boolean {
        return this.providers.has(name);
    }

    /**
     * Select a registered provider, or restore default computation with
     * `undefined`.
     */
    select(name: string | undefined) {
        if (name !== undefined && !this.providers.has(name)) {
            throw new Error(`Bond provider '${name}' is not registered.`);
        }
        this.activeName = name;
    }
}
