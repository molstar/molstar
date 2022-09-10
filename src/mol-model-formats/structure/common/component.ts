/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table, Column } from '../../../mol-data/db';
import { WaterNames, PolymerNames } from '../../../mol-model/structure/model/types';
import { SetUtils } from '../../../mol-util/set';
import { BasicSchema } from '../basic/schema';
import { mmCIF_chemComp_schema } from '../../../mol-io/reader/cif/schema/mmcif-extras';
import { SaccharideCompIdMap } from '../../../mol-model/structure/structure/carbohydrates/constants';

type Component = Table.Row<Pick<mmCIF_chemComp_schema, 'id' | 'name' | 'type'>>

const ProteinAtomIdsList = [
    new Set(['CA']),
    new Set(['C']),
    new Set(['N'])
];
const RnaAtomIdsList = [
    new Set(['P', 'O3\'', 'O3*']),
    new Set(['C4\'', 'C4*']),
    new Set(['O2\'', 'O2*', 'F2\'', 'F2*'])
];
const DnaAtomIdsList = [
    new Set(['P', 'O3\'', 'O3*']),
    new Set(['C3\'', 'C3*']),
    new Set(['O2\'', 'O2*', 'F2\'', 'F2*'])
];

/** Used to reduce false positives for atom name-based type guessing */
const NonPolymerNames = new Set([
    'FMN', 'NCN', 'FNS', 'FMA', 'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', // Mononucleotides
    'LIG'
]);

const StandardComponents = (function () {
    const map = new Map<string, Component>();
    const components: Component[] = [
        { id: 'HIS', name: 'HISTIDINE', type: 'l-peptide linking' },
        { id: 'ARG', name: 'ARGININE', type: 'l-peptide linking' },
        { id: 'LYS', name: 'LYSINE', type: 'l-peptide linking' },
        { id: 'ILE', name: 'ISOLEUCINE', type: 'l-peptide linking' },
        { id: 'PHE', name: 'PHENYLALANINE', type: 'l-peptide linking' },
        { id: 'LEU', name: 'LEUCINE', type: 'l-peptide linking' },
        { id: 'TRP', name: 'TRYPTOPHAN', type: 'l-peptide linking' },
        { id: 'ALA', name: 'ALANINE', type: 'l-peptide linking' },
        { id: 'MET', name: 'METHIONINE', type: 'l-peptide linking' },
        { id: 'CYS', name: 'CYSTEINE', type: 'l-peptide linking' },
        { id: 'ASN', name: 'ASPARAGINE', type: 'l-peptide linking' },
        { id: 'VAL', name: 'VALINE', type: 'l-peptide linking' },
        { id: 'GLY', name: 'GLYCINE', type: 'peptide linking' },
        { id: 'SER', name: 'SERINE', type: 'l-peptide linking' },
        { id: 'GLN', name: 'GLUTAMINE', type: 'l-peptide linking' },
        { id: 'TYR', name: 'TYROSINE', type: 'l-peptide linking' },
        { id: 'ASP', name: 'ASPARTIC ACID', type: 'l-peptide linking' },
        { id: 'GLU', name: 'GLUTAMIC ACID', type: 'l-peptide linking' },
        { id: 'THR', name: 'THREONINE', type: 'l-peptide linking' },
        { id: 'PRO', name: 'PROLINE', type: 'l-peptide linking' },
        { id: 'SEC', name: 'SELENOCYSTEINE', type: 'l-peptide linking' },
        { id: 'PYL', name: 'PYRROLYSINE', type: 'l-peptide linking' },

        { id: 'MSE', name: 'SELENOMETHIONINE', type: 'l-peptide linking' },
        { id: 'SEP', name: 'PHOSPHOSERINE', type: 'l-peptide linking' },
        { id: 'TPO', name: 'PHOSPHOTHREONINE', type: 'l-peptide linking' },
        { id: 'PTR', name: 'O-PHOSPHOTYROSINE', type: 'l-peptide linking' },
        { id: 'PCA', name: 'PYROGLUTAMIC ACID', type: 'l-peptide linking' },

        { id: 'A', name: 'ADENOSINE-5\'-MONOPHOSPHATE', type: 'rna linking' },
        { id: 'C', name: 'CYTIDINE-5\'-MONOPHOSPHATE', type: 'rna linking' },
        { id: 'T', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'rna linking' },
        { id: 'G', name: 'GUANOSINE-5\'-MONOPHOSPHATE', type: 'rna linking' },
        { id: 'I', name: 'INOSINIC ACID', type: 'rna linking' },
        { id: 'U', name: 'URIDINE-5\'-MONOPHOSPHATE', type: 'rna linking' },

        { id: 'DA', name: '2\'-DEOXYADENOSINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
        { id: 'DC', name: '2\'-DEOXYCYTIDINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
        { id: 'DT', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
        { id: 'DG', name: '2\'-DEOXYGUANOSINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
        { id: 'DI', name: '2\'-DEOXYINOSINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
        { id: 'DU', name: '2\'-DEOXYURIDINE-5\'-MONOPHOSPHATE', type: 'dna linking' },
    ];
    components.forEach(c => map.set(c.id, c));
    return map;
})();

const CharmmIonComponents = (function () {
    const map = new Map<string, Component>();
    const components: Component[] = [
        { id: 'ZN2', name: 'ZINC ION', type: 'ion' },
        { id: 'SOD', name: 'SODIUM ION', type: 'ion' },
        { id: 'CES', name: 'CESIUM ION', type: 'ion' },
        { id: 'CLA', name: 'CHLORIDE ION', type: 'ion' },
        { id: 'CAL', name: 'CALCIUM ION', type: 'ion' },
        { id: 'POT', name: 'POTASSIUM ION', type: 'ion' },
    ];
    components.forEach(c => map.set(c.id, c));
    return map;
})();

export class ComponentBuilder {
    private namesMap = new Map<string, string>();
    private comps = new Map<string, Component>();
    private ids: string[] = [];
    private names: string[] = [];
    private types: mmCIF_chemComp_schema['type']['T'][] = [];
    private mon_nstd_flags: mmCIF_chemComp_schema['mon_nstd_flag']['T'][] = [];

    private set(c: Component) {
        this.comps.set(c.id, c);
        this.ids.push(c.id);
        this.names.push(c.name);
        this.types.push(c.type);
        this.mon_nstd_flags.push(PolymerNames.has(c.id) ? 'y' : 'n');
    }

    private getAtomIds(index: number) {
        const atomIds = new Set<string>();
        const prevSeqId = this.seqId.value(index);
        while (index < this.seqId.rowCount) {
            const seqId = this.seqId.value(index);
            if (seqId !== prevSeqId) break;
            atomIds.add(this.atomId.value(index));
            prevSeqId - seqId;
            index += 1;
        }
        return atomIds;
    }

    private hasAtomIds(atomIds: Set<string>, atomIdsList: Set<string>[]) {
        for (let i = 0, il = atomIdsList.length; i < il; ++i) {
            if (!SetUtils.areIntersecting(atomIds, atomIdsList[i])) {
                return false;
            }
        }
        return true;
    }

    private getType(atomIds: Set<string>): Component['type'] {
        if (this.hasAtomIds(atomIds, ProteinAtomIdsList)) {
            return 'peptide linking';
        } else if (this.hasAtomIds(atomIds, RnaAtomIdsList)) {
            return 'rna linking';
        } else if (this.hasAtomIds(atomIds, DnaAtomIdsList)) {
            return 'dna linking';
        } else {
            return 'other';
        }
    }

    has(compId: string) { return this.comps.has(compId); }
    get(compId: string) { return this.comps.get(compId); }

    add(compId: string, index: number) {
        if (!this.has(compId)) {
            if (StandardComponents.has(compId)) {
                this.set(StandardComponents.get(compId)!);
            } else if (WaterNames.has(compId)) {
                this.set({ id: compId, name: 'WATER', type: 'non-polymer' });
            } else if (NonPolymerNames.has(compId.toUpperCase())) {
                this.set({ id: compId, name: this.namesMap.get(compId) || compId, type: 'non-polymer' });
            } else if (SaccharideCompIdMap.has(compId.toUpperCase())) {
                this.set({ id: compId, name: this.namesMap.get(compId) || compId, type: 'saccharide' });
            } else {
                const atomIds = this.getAtomIds(index);
                if (atomIds.size === 1 && CharmmIonComponents.has(compId)) {
                    this.set(CharmmIonComponents.get(compId)!);
                } else {
                    const type = this.getType(atomIds);
                    this.set({ id: compId, name: this.namesMap.get(compId) || compId, type });
                }
            }
        }
        return this.get(compId)!;
    }

    getChemCompTable() {
        return Table.ofPartialColumns(BasicSchema.chem_comp, {
            id: Column.ofStringArray(this.ids),
            name: Column.ofStringArray(this.names),
            type: Column.ofStringAliasArray(this.types),
            mon_nstd_flag: Column.ofStringAliasArray(this.mon_nstd_flags),
        }, this.ids.length);
    }

    setNames(names: [string, string][]) {
        names.forEach(n => this.namesMap.set(n[0], n[1]));
    }

    constructor(private seqId: Column<number>, private atomId: Column<string>) {

    }
}