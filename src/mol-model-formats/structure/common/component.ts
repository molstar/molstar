/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Table, Column } from '../../../mol-data/db';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { WaterNames, PolymerNames } from '../../../mol-model/structure/model/types';
import { SetUtils } from '../../../mol-util/set';
import { BasicSchema } from '../basic/schema';

type Component = Table.Row<Pick<mmCIF_Schema['chem_comp'], 'id' | 'name' | 'type'>>

const ProteinAtomIdsList = [
    new Set([ 'CA' ]),
    new Set([ 'C' ]),
    new Set([ 'N' ])
];
const RnaAtomIdsList = [
    new Set([ 'P', 'O3\'', 'O3*' ]),
    new Set([ 'C4\'', 'C4*' ]),
    new Set([ 'O2\'', 'O2*', 'F2\'', 'F2*' ])
];
const DnaAtomIdsList = [
    new Set([ 'P', 'O3\'', 'O3*' ]),
    new Set([ 'C3\'', 'C3*' ]),
    new Set([ 'O2\'', 'O2*', 'F2\'', 'F2*' ])
];

const StandardComponents = (function() {
    const map = new Map<string, Component>();
    const components: Component[] = [
        { id: 'HIS', name: 'HISTIDINE', type: 'L-peptide linking' },
        { id: 'ARG', name: 'ARGININE', type: 'L-peptide linking' },
        { id: 'LYS', name: 'LYSINE', type: 'L-peptide linking' },
        { id: 'ILE', name: 'ISOLEUCINE', type: 'L-peptide linking' },
        { id: 'PHE', name: 'PHENYLALANINE', type: 'L-peptide linking' },
        { id: 'LEU', name: 'LEUCINE', type: 'L-peptide linking' },
        { id: 'TRP', name: 'TRYPTOPHAN', type: 'L-peptide linking' },
        { id: 'ALA', name: 'ALANINE', type: 'L-peptide linking' },
        { id: 'MET', name: 'METHIONINE', type: 'L-peptide linking' },
        { id: 'CYS', name: 'CYSTEINE', type: 'L-peptide linking' },
        { id: 'ASN', name: 'ASPARAGINE', type: 'L-peptide linking' },
        { id: 'VAL', name: 'VALINE', type: 'L-peptide linking' },
        { id: 'GLY', name: 'GLYCINE', type: 'peptide linking' },
        { id: 'SER', name: 'SERINE', type: 'L-peptide linking' },
        { id: 'GLN', name: 'GLUTAMINE', type: 'L-peptide linking' },
        { id: 'TYR', name: 'TYROSINE', type: 'L-peptide linking' },
        { id: 'ASP', name: 'ASPARTIC ACID', type: 'L-peptide linking' },
        { id: 'GLU', name: 'GLUTAMIC ACID', type: 'L-peptide linking' },
        { id: 'THR', name: 'THREONINE', type: 'L-peptide linking' },
        { id: 'SEC', name: 'SELENOCYSTEINE', type: 'L-peptide linking' },
        { id: 'PYL', name: 'PYRROLYSINE', type: 'L-peptide linking' },

        { id: 'A', name: 'ADENOSINE-5\'-MONOPHOSPHATE', type: 'RNA linking' },
        { id: 'C', name: 'CYTIDINE-5\'-MONOPHOSPHATE', type: 'RNA linking' },
        { id: 'T', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'RNA linking' },
        { id: 'G', name: 'GUANOSINE-5\'-MONOPHOSPHATE', type: 'RNA linking' },
        { id: 'I', name: 'INOSINIC ACID', type: 'RNA linking' },
        { id: 'U', name: 'URIDINE-5\'-MONOPHOSPHATE', type: 'RNA linking' },

        { id: 'DA', name: '2\'-DEOXYADENOSINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
        { id: 'DC', name: '2\'-DEOXYCYTIDINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
        { id: 'DT', name: 'THYMIDINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
        { id: 'DG', name: '2\'-DEOXYGUANOSINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
        { id: 'DI', name: '2\'-DEOXYINOSINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
        { id: 'DU', name: '2\'-DEOXYURIDINE-5\'-MONOPHOSPHATE', type: 'DNA linking' },
    ];
    components.forEach(c => map.set(c.id, c));
    return map;
})();

export class ComponentBuilder {
    private namesMap = new Map<string, string>()
    private comps = new Map<string, Component>()
    private ids: string[] = []
    private names: string[] = []
    private types: mmCIF_Schema['chem_comp']['type']['T'][] = []
    private mon_nstd_flags: mmCIF_Schema['chem_comp']['mon_nstd_flag']['T'][] = []

    private set(c: Component) {
        this.comps.set(c.id, c);
        this.ids.push(c.id);
        this.names.push(c.name);
        this.types.push(c.type);
        this.mon_nstd_flags.push(PolymerNames.has(c.id) ? 'y' : 'n');
    }

    private getAtomIds(index: number) {
        const atomIds = new Set<string>();
        let prevSeqId = this.seqId.value(index);
        while (index < this.seqId.rowCount) {
            const seqId = this.seqId.value(index);
            if (seqId !== prevSeqId) break;
            atomIds.add(this.atomId.value(index));
            prevSeqId - seqId;
            index += 1;
        }
        return atomIds;
    }

    private hasAtomIds (atomIds: Set<string>, atomIdsList: Set<string>[]) {
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
            return 'RNA linking';
        } else if (this.hasAtomIds(atomIds, DnaAtomIdsList)) {
            return 'DNA linking';
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
            } else {
                const type = this.getType(this.getAtomIds(index));
                this.set({ id: compId, name: this.namesMap.get(compId) || compId, type });
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