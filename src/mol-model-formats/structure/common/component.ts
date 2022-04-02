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
]);

const StandardComponents = (function () {
    const map = new Map<string, Component>();
    const components: Component[] = [
        { id: 'HIS', name: 'HISTIDINE', type: 'L-PEPTIDE LINKING' },
        { id: 'ARG', name: 'ARGININE', type: 'L-PEPTIDE LINKING' },
        { id: 'LYS', name: 'LYSINE', type: 'L-PEPTIDE LINKING' },
        { id: 'ILE', name: 'ISOLEUCINE', type: 'L-PEPTIDE LINKING' },
        { id: 'PHE', name: 'PHENYLALANINE', type: 'L-PEPTIDE LINKING' },
        { id: 'LEU', name: 'LEUCINE', type: 'L-PEPTIDE LINKING' },
        { id: 'TRP', name: 'TRYPTOPHAN', type: 'L-PEPTIDE LINKING' },
        { id: 'ALA', name: 'ALANINE', type: 'L-PEPTIDE LINKING' },
        { id: 'MET', name: 'METHIONINE', type: 'L-PEPTIDE LINKING' },
        { id: 'CYS', name: 'CYSTEINE', type: 'L-PEPTIDE LINKING' },
        { id: 'ASN', name: 'ASPARAGINE', type: 'L-PEPTIDE LINKING' },
        { id: 'VAL', name: 'VALINE', type: 'L-PEPTIDE LINKING' },
        { id: 'GLY', name: 'GLYCINE', type: 'PEPTIDE LINKING' },
        { id: 'SER', name: 'SERINE', type: 'L-PEPTIDE LINKING' },
        { id: 'GLN', name: 'GLUTAMINE', type: 'L-PEPTIDE LINKING' },
        { id: 'TYR', name: 'TYROSINE', type: 'L-PEPTIDE LINKING' },
        { id: 'ASP', name: 'ASPARTIC ACID', type: 'L-PEPTIDE LINKING' },
        { id: 'GLU', name: 'GLUTAMIC ACID', type: 'L-PEPTIDE LINKING' },
        { id: 'THR', name: 'THREONINE', type: 'L-PEPTIDE LINKING' },
        { id: 'PRO', name: 'PROLINE', type: 'L-PEPTIDE LINKING' },
        { id: 'SEC', name: 'SELENOCYSTEINE', type: 'L-PEPTIDE LINKING' },
        { id: 'PYL', name: 'PYRROLYSINE', type: 'L-PEPTIDE LINKING' },

        { id: 'MSE', name: 'SELENOMETHIONINE', type: 'L-PEPTIDE LINKING' },
        { id: 'SEP', name: 'PHOSPHOSERINE', type: 'L-PEPTIDE LINKING' },
        { id: 'TPO', name: 'PHOSPHOTHREONINE', type: 'L-PEPTIDE LINKING' },
        { id: 'PTR', name: 'O-PHOSPHOTYROSINE', type: 'L-PEPTIDE LINKING' },
        { id: 'PCA', name: 'PYROGLUTAMIC ACID', type: 'L-PEPTIDE LINKING' },

        { id: 'A', name: 'ADENOSINE-5\'-MONOPHOSPHATE', type: 'RNA LINKING' },
        { id: 'C', name: 'CYTIDINE-5\'-MONOPHOSPHATE', type: 'RNA LINKING' },
        { id: 'T', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'RNA LINKING' },
        { id: 'G', name: 'GUANOSINE-5\'-MONOPHOSPHATE', type: 'RNA LINKING' },
        { id: 'I', name: 'INOSINIC ACID', type: 'RNA LINKING' },
        { id: 'U', name: 'URIDINE-5\'-MONOPHOSPHATE', type: 'RNA LINKING' },

        { id: 'DA', name: '2\'-DEOXYADENOSINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
        { id: 'DC', name: '2\'-DEOXYCYTIDINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
        { id: 'DT', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
        { id: 'DG', name: '2\'-DEOXYGUANOSINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
        { id: 'DI', name: '2\'-DEOXYINOSINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
        { id: 'DU', name: '2\'-DEOXYURIDINE-5\'-MONOPHOSPHATE', type: 'DNA LINKING' },
    ];
    components.forEach(c => map.set(c.id, c));
    return map;
})();

const CharmmIonComponents = (function () {
    const map = new Map<string, Component>();
    const components: Component[] = [
        { id: 'ZN2', name: 'ZINC ION', type: 'ION' },
        { id: 'SOD', name: 'SODIUM ION', type: 'ION' },
        { id: 'CES', name: 'CESIUM ION', type: 'ION' },
        { id: 'CLA', name: 'CHLORIDE ION', type: 'ION' },
        { id: 'CAL', name: 'CALCIUM ION', type: 'ION' },
        { id: 'POT', name: 'POTASSIUM ION', type: 'ION' },
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
        this.mon_nstd_flags.push(PolymerNames.has(c.id) ? 'Y' : 'N');
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
            return 'PEPTIDE LINKING';
        } else if (this.hasAtomIds(atomIds, RnaAtomIdsList)) {
            return 'RNA LINKING';
        } else if (this.hasAtomIds(atomIds, DnaAtomIdsList)) {
            return 'DNA LINKING';
        } else {
            return 'OTHER';
        }
    }

    has(compId: string) { return this.comps.has(compId); }
    get(compId: string) { return this.comps.get(compId); }

    add(compId: string, index: number) {
        if (!this.has(compId)) {
            if (StandardComponents.has(compId)) {
                this.set(StandardComponents.get(compId)!);
            } else if (WaterNames.has(compId)) {
                this.set({ id: compId, name: 'WATER', type: 'NON-POLYMER' });
            } else if (NonPolymerNames.has(compId.toUpperCase())) {
                this.set({ id: compId, name: this.namesMap.get(compId) || compId, type: 'NON-POLYMER' });
            } else if (SaccharideCompIdMap.has(compId.toUpperCase())) {
                this.set({ id: compId, name: this.namesMap.get(compId) || compId, type: 'SACCHARIDE' });
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