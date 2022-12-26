/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../../mol-data/db';
import { Interval, Segmentation } from '../../../mol-data/int';
import { toDatabase } from '../../../mol-io/reader/cif/schema';
import { SymmetryOperator } from '../../../mol-math/geometry';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { ChainIndex, ElementIndex } from '../../../mol-model/structure';
import { AtomSiteOperatorMappingSchema } from '../../../mol-model/structure/export/categories/atom_site_operator_mapping';
import { Model } from '../../../mol-model/structure/model/model';
import { AtomicConformation, AtomicData, AtomicHierarchy, AtomicSegments, AtomsSchema, ChainsSchema, ResiduesSchema } from '../../../mol-model/structure/model/properties/atomic';
import { Entities } from '../../../mol-model/structure/model/properties/common';
import { getAtomicDerivedData } from '../../../mol-model/structure/model/properties/utils/atomic-derived';
import { getAtomicIndex } from '../../../mol-model/structure/model/properties/utils/atomic-index';
import { ElementSymbol } from '../../../mol-model/structure/model/types';
import { UUID } from '../../../mol-util/uuid';
import { ModelFormat } from '../../format';
import { MmcifFormat } from '../mmcif';
import { AtomSite } from './schema';

function findHierarchyOffsets(atom_site: AtomSite) {
    if (atom_site._rowCount === 0) return { residues: [], chains: [] };

    const start = 0, end = atom_site._rowCount;
    const residues = [start as ElementIndex], chains = [start as ElementIndex];

    const { label_entity_id, label_asym_id, label_seq_id, auth_seq_id, pdbx_PDB_ins_code } = atom_site;

    for (let i = start + 1 as ElementIndex; i < end; i++) {
        const newChain = !label_entity_id.areValuesEqual(i - 1, i) || !label_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !label_seq_id.areValuesEqual(i - 1, i)
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i);
        // not checking label_comp_id to allow for MICROHETEROGENEITY

        if (newResidue) residues[residues.length] = i as ElementIndex;
        if (newChain) chains[chains.length] = i as ElementIndex;
    }
    return { residues, chains };
}

function createHierarchyData(atom_site: AtomSite, sourceIndex: Column<number>, offsets: { residues: ArrayLike<number>, chains: ArrayLike<number> }): AtomicData {
    const atoms = Table.ofColumns(AtomsSchema, {
        type_symbol: Column.ofArray({ array: Column.mapToArray(atom_site.type_symbol, ElementSymbol), schema: Column.Schema.Aliased<ElementSymbol>(Column.Schema.str) }),
        label_atom_id: atom_site.label_atom_id,
        auth_atom_id: atom_site.auth_atom_id,
        label_alt_id: atom_site.label_alt_id,
        label_comp_id: atom_site.label_comp_id,
        auth_comp_id: atom_site.auth_comp_id,
        pdbx_formal_charge: atom_site.pdbx_formal_charge
    });

    const residues = Table.view(atom_site, ResiduesSchema, offsets.residues);
    const chains = Table.view(atom_site, ChainsSchema, offsets.chains);

    if (!residues.label_seq_id.isDefined) {
        const seqIds = new Int32Array(residues.label_seq_id.rowCount);
        const { residues: residueOffsets, chains: chainOffsets } = offsets;
        let cI = 0;
        let seqId = 0;
        for (let i = 0, il = seqIds.length; i < il; ++i) {
            if (residueOffsets[i] >= chainOffsets[cI + 1]) {
                cI += 1;
                seqId = 0;
            }
            seqIds[i] = ++seqId; // start id on one
        }
        residues.label_seq_id = Column.ofIntArray(seqIds);
    }

    // Optimize the numeric columns
    Table.columnToArray(residues, 'label_seq_id', Int32Array);
    Table.columnToArray(residues, 'auth_seq_id', Int32Array);

    return { atoms, residues, chains, atomSourceIndex: sourceIndex };
}

function getConformation(atom_site: AtomSite): AtomicConformation {
    return {
        id: UUID.create22(),
        atomId: atom_site.id,
        occupancy: atom_site.occupancy.isDefined ? atom_site.occupancy : Column.ofConst(1, atom_site._rowCount, Column.Schema.float),
        B_iso_or_equiv: atom_site.B_iso_or_equiv,
        xyzDefined: atom_site.Cartn_x.isDefined && atom_site.Cartn_y.isDefined && atom_site.Cartn_z.isDefined,
        x: atom_site.Cartn_x.toArray({ array: Float32Array }),
        y: atom_site.Cartn_y.toArray({ array: Float32Array }),
        z: atom_site.Cartn_z.toArray({ array: Float32Array }),
    };
}

function isHierarchyDataEqual(a: AtomicData, b: AtomicData) {
    return Table.areEqual(a.chains, b.chains)
        && Table.areEqual(a.residues, b.residues)
        && Table.areEqual(a.atoms, b.atoms);
}

function createChainOperatorMappingAndSubstituteNames(hierarchy: AtomicData, format: ModelFormat) {
    const mapping = new Map<ChainIndex, SymmetryOperator>();
    if (!MmcifFormat.is(format)) return mapping;

    const { molstar_atom_site_operator_mapping: entries } = toDatabase(AtomSiteOperatorMappingSchema, format.data.frame);
    if (entries._rowCount === 0) return mapping;

    const labelMap = new Map<string, { name: string, operator: SymmetryOperator }>();
    const authMap = new Map<string, string>();

    for (let i = 0; i < entries._rowCount; i++) {
        const assembly: SymmetryOperator['assembly'] = entries.assembly_operator_id.valueKind(i) === Column.ValueKinds.Present
            ? { id: entries.assembly_id.value(i), operList: [], operId: entries.assembly_operator_id.value(i) }
            : void 0;

        const operator = SymmetryOperator.create(entries.operator_name.value(i), Mat4.identity(), {
            assembly,
            spgrOp: entries.symmetry_operator_index.valueKind(i) === Column.ValueKinds.Present ? entries.symmetry_operator_index.value(i) : void 0,
            hkl: Vec3.ofArray(entries.symmetry_hkl.value(i)),
            ncsId: entries.ncs_id.value(i)
        });

        const suffix = entries.suffix.value(i);
        const label = entries.label_asym_id.value(i);
        labelMap.set(`${label}${suffix}`, { name: label, operator });
        const auth = entries.auth_asym_id.value(i);
        authMap.set(`${auth}${suffix}`, auth);
    }

    const { label_asym_id, auth_asym_id } = hierarchy.chains;
    const mappedLabel: string[] = new Array(label_asym_id.rowCount);
    const mappedAuth: string[] = new Array(label_asym_id.rowCount);

    for (let i = 0 as ChainIndex; i < label_asym_id.rowCount; i++) {
        const label = label_asym_id.value(i), auth = auth_asym_id.value(i);
        if (!labelMap.has(label)) {
            mappedLabel[i] = label;
            mappedAuth[i] = auth;
            continue;
        }

        const { name, operator } = labelMap.get(label)!;
        mapping.set(i, operator);

        mappedLabel[i] = name;
        mappedAuth[i] = authMap.get(auth) || auth;
    }

    hierarchy.chains.label_asym_id = Column.ofArray({ array: mappedLabel, valueKind: hierarchy.chains.label_asym_id.valueKind, schema: hierarchy.chains.label_asym_id.schema });
    hierarchy.chains.auth_asym_id = Column.ofArray({ array: mappedAuth, valueKind: hierarchy.chains.auth_asym_id.valueKind, schema: hierarchy.chains.auth_asym_id.schema });

    return mapping;
}

function getAtomicHierarchy(atom_site: AtomSite, sourceIndex: Column<number>, entities: Entities, chemicalComponentMap: Model['properties']['chemicalComponentMap'], format: ModelFormat, previous?: Model) {
    const hierarchyOffsets = findHierarchyOffsets(atom_site);
    const hierarchyData = createHierarchyData(atom_site, sourceIndex, hierarchyOffsets);
    const chainOperatorMapping = createChainOperatorMappingAndSubstituteNames(hierarchyData, format);

    if (previous && isHierarchyDataEqual(previous.atomicHierarchy, hierarchyData)) {
        return {
            sameAsPrevious: true,
            hierarchy: previous.atomicHierarchy,
            chainOperatorMapping
        };
    }

    const hierarchySegments: AtomicSegments = {
        residueAtomSegments: Segmentation.ofOffsets(hierarchyOffsets.residues, Interval.ofBounds(0, atom_site._rowCount)),
        chainAtomSegments: Segmentation.ofOffsets(hierarchyOffsets.chains, Interval.ofBounds(0, atom_site._rowCount)),
    };

    const index = getAtomicIndex(hierarchyData, entities, hierarchySegments);
    const derived = getAtomicDerivedData(hierarchyData, hierarchySegments, index, chemicalComponentMap);
    const hierarchy: AtomicHierarchy = { ...hierarchyData, ...hierarchySegments, index, derived };
    return { sameAsPrevious: false, hierarchy, chainOperatorMapping };
}

export function getAtomicHierarchyAndConformation(atom_site: AtomSite, sourceIndex: Column<number>, entities: Entities, chemicalComponentMap: Model['properties']['chemicalComponentMap'], format: ModelFormat, previous?: Model) {
    const { sameAsPrevious, hierarchy, chainOperatorMapping } = getAtomicHierarchy(atom_site, sourceIndex, entities, chemicalComponentMap, format, previous);
    const conformation = getConformation(atom_site);
    return { sameAsPrevious, hierarchy, conformation, chainOperatorMapping };
}