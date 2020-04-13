/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MoleculeType, isPolymer } from '../../../mol-model/structure/model/types';
import { Column, Table } from '../../../mol-data/db';
import { BasicSchema } from '../basic/schema';

export type EntityCompound = { chains: string[], description: string }

// TODO add support for `branched`
type EntityType = 'water' | 'polymer' | 'non-polymer'

export class EntityBuilder {
    private count = 0
    private ids: string[] = []
    private types: EntityType[] = []
    private descriptions: string[][] = []

    private compoundsMap = new Map<string, string>()
    private namesMap = new Map<string, string>()
    private heteroMap = new Map<string, string>()
    private chainMap = new Map<string, string>()
    private waterId?: string

    private set(type: EntityType, description: string) {
        this.count += 1;
        this.ids.push(`${this.count}`);
        this.types.push(type);
        this.descriptions.push([description]);
    }

    getEntityId(compId: string, moleculeType: MoleculeType, chainId: string): string {
        if (moleculeType === MoleculeType.Water) {
            if (this.waterId === undefined) {
                this.set('water', 'Water');
                this.waterId = `${this.count}`;
            }
            return this.waterId;
        } else if (isPolymer(moleculeType)) {
            if (this.compoundsMap.has(chainId)) {
                return this.compoundsMap.get(chainId)!;
            } else {
                if (!this.chainMap.has(chainId)) {
                    this.set('polymer', `Polymer ${this.chainMap.size + 1}`);
                    this.chainMap.set(chainId, `${this.count}`);
                }
                return this.chainMap.get(chainId)!;
            }
        } else {
            if (!this.heteroMap.has(compId)) {
                this.set('non-polymer', this.namesMap.get(compId) || compId);
                this.heteroMap.set(compId, `${this.count}`);
            }
            return this.heteroMap.get(compId)!;
        }
    }

    getEntityTable() {
        return Table.ofPartialColumns(BasicSchema.entity, {
            id: Column.ofStringArray(this.ids),
            type: Column.ofStringAliasArray(this.types),
            pdbx_description: Column.ofStringListArray(this.descriptions),
        }, this.count);
    }

    setCompounds(compounds: EntityCompound[]) {
        for (let i = 0, il = compounds.length; i < il; ++i) {
            const { chains, description } = compounds[i];
            this.set('polymer', description);
            for (let j = 0, jl = chains.length; j < jl; ++j) {
                this.compoundsMap.set(chains[j], `${this.count}`);
            }
        }
    }

    setNames(names: [string, string][]) {
        names.forEach(n => this.namesMap.set(n[0], n[1]));
    }
}