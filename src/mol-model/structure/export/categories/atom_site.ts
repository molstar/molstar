/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../../mol-data/db';
import { mmCIF_Database } from '../../../../mol-io/reader/cif/schema/mmcif';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { SIFTSMapping } from '../../../../mol-model-props/sequence/sifts-mapping';
import { StructureElement, Structure, StructureProperties as P } from '../../structure';
import { CifExportContext } from '../mmcif';
import CifField = CifWriter.Field
import CifCategory = CifWriter.Category
import E = CifWriter.Encodings

function atom_site_label_asym_id(e: StructureElement.Location) {
    const l = P.chain.label_asym_id(e);
    const suffix = e.unit.conformation.operator.suffix;
    if (!suffix) return l;
    return l + suffix;
}

function atom_site_auth_asym_id(e: StructureElement.Location) {
    const l = P.chain.auth_asym_id(e);
    const suffix = e.unit.conformation.operator.suffix;
    if (!suffix) return l;
    return l + suffix;
}


const atom_site_pdbx_label_index = {
    shouldInclude(s: AtomSiteData) {
        return !!s.atom_site?.pdbx_label_index.isDefined;
    },
    value(e: StructureElement.Location, d: AtomSiteData) {
        const srcIndex = d.sourceIndex.value(e.element);
        return d.atom_site!.pdbx_label_index.value(srcIndex);
    },
};

const SIFTS = {
    shouldInclude(s: AtomSiteData) {
        return SIFTSMapping.isAvailable(s.structure.models[0]);
    },
    pdbx_sifts_xref_db_name: {
        value(e: StructureElement.Location, d: AtomSiteData) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_name.value(srcIndex);
        },
        valueKind(e: StructureElement.Location, d: any) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_name.valueKind(srcIndex);
        },
    },
    pdbx_sifts_xref_db_acc: {
        value(e: StructureElement.Location, d: AtomSiteData) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_acc.value(srcIndex);
        },
        valueKind(e: StructureElement.Location, d: any) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_acc.valueKind(srcIndex);
        },
    },
    pdbx_sifts_xref_db_num: {
        value(e: StructureElement.Location, d: AtomSiteData) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_num.value(srcIndex);
        },
        valueKind(e: StructureElement.Location, d: any) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_num.valueKind(srcIndex);
        },
    },
    pdbx_sifts_xref_db_res: {
        value(e: StructureElement.Location, d: AtomSiteData) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_res.value(srcIndex);
        },
        valueKind(e: StructureElement.Location, d: any) {
            const srcIndex = d.sourceIndex.value(e.element);
            return d.atom_site!.pdbx_sifts_xref_db_res.valueKind(srcIndex);
        },
    }
};

const atom_site_fields = () => CifWriter.fields<StructureElement.Location, AtomSiteData>()
    .str('group_PDB', P.residue.group_PDB)
    .index('id')
    .str('type_symbol', P.atom.type_symbol as any)
    .str('label_atom_id', P.atom.label_atom_id)

    .str('label_comp_id', P.atom.label_comp_id)
    .int('label_seq_id', P.residue.label_seq_id, {
        encoder: E.deltaRLE,
        valueKind: (k, d) => {
            const m = k.unit.model;
            return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[k.element]);
        }
    })
    .str('label_alt_id', P.atom.label_alt_id)
    .str('pdbx_PDB_ins_code', P.residue.pdbx_PDB_ins_code)

    .str('label_asym_id', atom_site_label_asym_id)
    .str('label_entity_id', P.chain.label_entity_id)

    .float('Cartn_x', P.atom.x, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_y', P.atom.y, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('Cartn_z', P.atom.z, { digitCount: 3, encoder: E.fixedPoint3 })
    .float('occupancy', P.atom.occupancy, { digitCount: 2, encoder: E.fixedPoint2 })
    .float('B_iso_or_equiv', P.atom.B_iso_or_equiv, { digitCount: 2, encoder: E.fixedPoint2 })
    .int('pdbx_formal_charge', P.atom.pdbx_formal_charge, {
        encoder: E.deltaRLE,
        valueKind: (k, d) => k.unit.model.atomicHierarchy.atoms.pdbx_formal_charge.valueKind(k.element)
    })

    .str('auth_atom_id', P.atom.auth_atom_id)
    .str('auth_comp_id', P.atom.auth_comp_id)
    .int('auth_seq_id', P.residue.auth_seq_id, { encoder: E.deltaRLE })
    .str('auth_asym_id', atom_site_auth_asym_id)

    .int('pdbx_PDB_model_num', P.unit.model_num, { encoder: E.deltaRLE })

    .int('pdbx_label_index', atom_site_pdbx_label_index.value, { shouldInclude: atom_site_pdbx_label_index.shouldInclude })

    // SIFTS
    .str('pdbx_sifts_xref_db_name', SIFTS.pdbx_sifts_xref_db_name.value, { shouldInclude: SIFTS.shouldInclude, valueKind: SIFTS.pdbx_sifts_xref_db_name.valueKind })
    .str('pdbx_sifts_xref_db_acc', SIFTS.pdbx_sifts_xref_db_acc.value, { shouldInclude: SIFTS.shouldInclude, valueKind: SIFTS.pdbx_sifts_xref_db_acc.valueKind })
    .str('pdbx_sifts_xref_db_num', SIFTS.pdbx_sifts_xref_db_num.value, { shouldInclude: SIFTS.shouldInclude, valueKind: SIFTS.pdbx_sifts_xref_db_num.valueKind })
    .str('pdbx_sifts_xref_db_res', SIFTS.pdbx_sifts_xref_db_res.value, { shouldInclude: SIFTS.shouldInclude, valueKind: SIFTS.pdbx_sifts_xref_db_res.valueKind })

    // .str('operator_name', P.unit.operator_name, {
    //     shouldInclude: structure => structure.units.some(u => !u.conformation.operator.isIdentity)
    // })
    .getFields();

interface AtomSiteData {
    structure: Structure,
    sourceIndex: Column<number>,
    atom_site?: mmCIF_Database['atom_site']
}

export const _atom_site: CifCategory<CifExportContext> = {
    name: 'atom_site',
    instance({ structures }: CifExportContext) {
        return {
            fields: atom_site_fields(),
            source: structures.map(s => ({
                data: {
                    structure: s,
                    sourceIndex: s.model.atomicHierarchy.atomSourceIndex,
                    atom_site: MmcifFormat.is(s.model.sourceData) ? s.model.sourceData.data.db.atom_site : void 0
                } as AtomSiteData,
                rowCount: s.elementCount,
                keys: () => s.elementLocations()
            }))
        };
    }
};

function prepostfixed(prefix: string | undefined, name: string) {
    if (prefix) return `${prefix}_${name}`;
    return name;
}

function prefixedInsCode(prefix: string | undefined) {
    if (!prefix) return 'pdbx_PDB_ins_code';
    return `pdbx_${prefix}_PDB_ins_code`;
}

function mappedProp<K, D>(loc: (key: K, data: D) => StructureElement.Location, prop: (e: StructureElement.Location) => any) {
    return (k: K, d: D) => prop(loc(k, d));
}

function addModelNum<K, D>(fields: CifWriter.Field.Builder<K, D>, getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions) {
    if (options && options.includeModelNum) {
        fields.int('pdbx_PDB_model_num', mappedProp(getLocation, P.unit.model_num));
    }
}

export interface IdFieldsOptions {
    prefix?: string,
    includeModelNum?: boolean
}

export function residueIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix;
    const ret = CifWriter.fields<K, D>()
        .str(prepostfixed(prefix, `label_comp_id`), mappedProp(getLocation, P.atom.label_comp_id))
        .int(prepostfixed(prefix, `label_seq_id`), mappedProp(getLocation, P.residue.label_seq_id), {
            encoder: E.deltaRLE,
            valueKind: (k, d) => {
                const e = getLocation(k, d);
                const m = e.unit.model;
                return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[e.element]);
            }
        })
        .str(prefixedInsCode(prefix), mappedProp(getLocation, P.residue.pdbx_PDB_ins_code))

        .str(prepostfixed(prefix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))

        .str(prepostfixed(prefix, `auth_comp_id`), mappedProp(getLocation, P.atom.auth_comp_id))
        .int(prepostfixed(prefix, `auth_seq_id`), mappedProp(getLocation, P.residue.auth_seq_id), { encoder: E.deltaRLE })
        .str(prepostfixed(prefix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function chainIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix;
    const ret = CifField.build<K, D>()
        .str(prepostfixed(prefix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))
        .str(prepostfixed(prefix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function entityIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix;
    const ret = CifField.build<K, D>()
        .str(prepostfixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}

export function atomIdFields<K, D>(getLocation: (key: K, data: D) => StructureElement.Location, options?: IdFieldsOptions): CifField<K, D>[] {
    const prefix = options && options.prefix;
    const ret = CifWriter.fields<K, D>()
        .str(prepostfixed(prefix, `label_atom_id`), mappedProp(getLocation, P.atom.label_atom_id))
        .str(prepostfixed(prefix, `label_comp_id`), mappedProp(getLocation, P.atom.label_comp_id))
        .int(prepostfixed(prefix, `label_seq_id`), mappedProp(getLocation, P.residue.label_seq_id), {
            encoder: E.deltaRLE,
            valueKind: (k, d) => {
                const e = getLocation(k, d);
                const m = e.unit.model;
                return m.atomicHierarchy.residues.label_seq_id.valueKind(m.atomicHierarchy.residueAtomSegments.index[e.element]);
            }
        })
        .str(prepostfixed(prefix, `label_alt_id`), mappedProp(getLocation, P.atom.label_alt_id))
        .str(prefixedInsCode(prefix), mappedProp(getLocation, P.residue.pdbx_PDB_ins_code))

        .str(prepostfixed(prefix, `label_asym_id`), mappedProp(getLocation, P.chain.label_asym_id))
        .str(prepostfixed(prefix, `label_entity_id`), mappedProp(getLocation, P.chain.label_entity_id))

        .str(prepostfixed(prefix, `auth_atom_id`), mappedProp(getLocation, P.atom.auth_atom_id))
        .str(prepostfixed(prefix, `auth_comp_id`), mappedProp(getLocation, P.atom.auth_comp_id))
        .int(prepostfixed(prefix, `auth_seq_id`), mappedProp(getLocation, P.residue.auth_seq_id), { encoder: E.deltaRLE })
        .str(prepostfixed(prefix, `auth_asym_id`), mappedProp(getLocation, P.chain.auth_asym_id));

    addModelNum(ret, getLocation, options);
    return ret.getFields();
}