/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from 'mol-data/db';
import { Interval, Segmentation } from 'mol-data/int';
import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif';
import UUID from 'mol-util/uuid';
import { ElementIndex } from '../../../../structure';
import Format from '../../format';
import { Model } from '../../model';
import { AtomicConformation, AtomicData, AtomicHierarchy, AtomicSegments, AtomsSchema, ChainsSchema, ResiduesSchema } from '../../properties/atomic';
import { getAtomicIndex } from '../../properties/utils/atomic-index';
import { ElementSymbol } from '../../types';
import { Entities } from '../../properties/common';

import mmCIF_Format = Format.mmCIF
import { getAtomicRanges } from '../../properties/utils/atomic-ranges';
import { FormatData } from '../mmcif';
import { getAtomicDerivedData } from '../../properties/utils/atomic-derived';

type AtomSite = mmCIF_Database['atom_site']

function findHierarchyOffsets(atom_site: AtomSite) {
    if (atom_site._rowCount === 0) return { residues: [], chains: [] };

    const start = 0, end = atom_site._rowCount;
    const residues = [start as ElementIndex], chains = [start as ElementIndex];

    const { label_entity_id, label_asym_id, label_seq_id, auth_seq_id, pdbx_PDB_ins_code, label_comp_id } = atom_site;

    for (let i = start + 1 as ElementIndex; i < end; i++) {
        const newChain = !label_entity_id.areValuesEqual(i - 1, i) || !label_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !label_seq_id.areValuesEqual(i - 1, i)
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i)
            || !label_comp_id.areValuesEqual(i - 1, i);

        if (newResidue) residues[residues.length] = i as ElementIndex;
        if (newChain) chains[chains.length] = i as ElementIndex;
    }
    return { residues, chains };
}

function createHierarchyData(atom_site: AtomSite, offsets: { residues: ArrayLike<number>, chains: ArrayLike<number> }): AtomicData {
    const atoms = Table.ofColumns(AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(atom_site.type_symbol, ElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: atom_site.label_atom_id,
        auth_atom_id: atom_site.auth_atom_id,
        label_alt_id: atom_site.label_alt_id,
        pdbx_formal_charge: atom_site.pdbx_formal_charge
    });
    const residues = Table.view(atom_site, ResiduesSchema, offsets.residues);
    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);
    const chains = Table.view(atom_site, ChainsSchema, offsets.chains);
    return { atoms, residues, chains };
}

function getConformation(atom_site: AtomSite): AtomicConformation {
    return {
        id: UUID.create22(),
        atomId: atom_site.id,
        occupancy: atom_site.occupancy,
        B_iso_or_equiv: atom_site.B_iso_or_equiv,
        x: atom_site.Cartn_x.toArray({ array: Float32Array }),
        y: atom_site.Cartn_y.toArray({ array: Float32Array }),
        z: atom_site.Cartn_z.toArray({ array: Float32Array }),
    }
}

function isHierarchyDataEqual(a: AtomicData, b: AtomicData) {
    // TODO need to cast because of how TS handles type resolution for interfaces https://github.com/Microsoft/TypeScript/issues/15300
    return Table.areEqual(a.chains as Table<ChainsSchema>, b.chains as Table<ChainsSchema>)
        && Table.areEqual(a.residues as Table<ResiduesSchema>, b.residues as Table<ResiduesSchema>)
        && Table.areEqual(a.atoms as Table<AtomsSchema>, b.atoms as Table<AtomsSchema>)
}

export function getAtomicHierarchyAndConformation(format: mmCIF_Format, atom_site: AtomSite, entities: Entities, formatData: FormatData, previous?: Model) {
    const hierarchyOffsets = findHierarchyOffsets(atom_site);
    const hierarchyData = createHierarchyData(atom_site, hierarchyOffsets);

    if (previous && isHierarchyDataEqual(previous.atomicHierarchy, hierarchyData)) {
        return {
            sameAsPrevious: true,
            hierarchy: previous.atomicHierarchy,
            conformation: getConformation(atom_site)
        };
    }

    const conformation = getConformation(atom_site)

    const hierarchySegments: AtomicSegments = {
        residueAtomSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, Interval.ofBounds(0, atom_site._rowCount)),
        chainAtomSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, Interval.ofBounds(0, atom_site._rowCount)),
    }

    const index = getAtomicIndex(hierarchyData, entities, hierarchySegments);
    console.time('derived')
    const derived = getAtomicDerivedData(hierarchyData, index, formatData.chemicalComponentMap);
    console.timeEnd('derived')
    console.log(derived)
    const hierarchyRanges = getAtomicRanges(hierarchyData, hierarchySegments, conformation, formatData.chemicalComponentMap);
    const hierarchy: AtomicHierarchy = { ...hierarchyData, ...hierarchySegments, ...hierarchyRanges, index, derived };
    return { sameAsPrevious: false, hierarchy, conformation };
}