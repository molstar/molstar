/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MoleculeType, isPolymer } from '../../../mol-model/structure/model/types';
import { Column, Table } from '../../../mol-data/db';
import { BasicSchema } from '../basic/schema';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';

export type EntityCompound = { chains: string[], description: string }

type EntityType = mmCIF_Schema['entity']['type']['T']

export class EntityBuilder {
    private count = 0;
    private ids: string[] = [];
    private types: EntityType[] = [];
    private descriptions: string[][] = [];

    private compoundsMap = new Map<string, string>();
    private seqresMap = new Map<string, { key: string, residues: Set<string> }>();
    private namesMap = new Map<string, string>();
    private heteroMap = new Map<string, string>();
    private chainMap = new Map<string, string>();
    private sequenceMap = new Map<string, string>();
    private waterId?: string;

    private polymerCount = 0;

    private set(type: EntityType, description: string) {
        this.count += 1;
        this.ids.push(`${this.count}`);
        this.types.push(type);
        this.descriptions.push([description]);
    }

    private addPolymer(map: Map<string, string>, key: string, options?: { customName?: string }) {
        if (!map.has(key)) {
            this.polymerCount += 1;
            this.set('polymer', options?.customName || `Polymer ${this.polymerCount}`);
            map.set(key, `${this.count}`);
        }
        return map.get(key)!;
    }

    private addNonPolymer(map: Map<string, string>, key: string, moleculeType: MoleculeType, options?: { customName?: string }) {
        if (!map.has(key)) {
            const type = moleculeType === MoleculeType.Saccharide ? 'branched' : 'non-polymer';
            this.set(type, options?.customName || this.namesMap.get(key) || key);
            map.set(key, `${this.count}`);
        }
        return map.get(key)!;
    }

    getEntityId(compId: string, moleculeType: MoleculeType, chainId: string, options?: { customName?: string }): string {
        if (moleculeType === MoleculeType.Water) {
            if (this.waterId === undefined) {
                this.set('water', options?.customName || 'Water');
                this.waterId = `${this.count}`;
            }
            return this.waterId;
        } else if (isPolymer(moleculeType)) {
            if (this.compoundsMap.has(chainId)) {
                return this.compoundsMap.get(chainId)!;
            } else {
                if (this.seqresMap.has(chainId)) {
                    const { key, residues } = this.seqresMap.get(chainId)!;
                    if (residues.has(compId)) {
                        return this.addPolymer(this.sequenceMap, key, options);
                    }
                }
                return this.addPolymer(this.chainMap, chainId, options);
            }
        } else {
            if (this.seqresMap.has(chainId)) {
                const { key, residues } = this.seqresMap.get(chainId)!;
                if (residues.has(compId)) {
                    return this.addNonPolymer(this.sequenceMap, key, moleculeType, options);
                }
            }
            return this.addNonPolymer(this.heteroMap, compId, moleculeType, options);
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

    setSeqres(seqresMap: Map<string, string[]>) {
        for (const [chainId, residues] of seqresMap) {
            this.seqresMap.set(chainId, {
                key: residues.join('-'),
                residues: new Set(residues)
            });
        }
    }
}
