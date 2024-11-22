/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yana Rose <yana.v.rose@gmail.com>
 */

import { substringStartsWith } from '../../../mol-util/string';
import { CifCategory, CifField, CifFrame } from '../../../mol-io/reader/cif';
import { Tokenizer } from '../../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../../mol-io/reader/pdb/schema';
import { parseCryst1, parseRemark350, parseMtrix } from './assembly';
import { parseHelix, parseSheet } from './secondary-structure';
import { parseCmpnd, parseHetnam } from './entity';
import { ComponentBuilder } from '../common/component';
import { EntityBuilder } from '../common/entity';
import { Column } from '../../../mol-data/db';
import { getMoleculeType } from '../../../mol-model/structure/model/types';
import { getAtomSiteTemplate, addAtom, getAtomSite, LabelAsymIdHelper } from './atom-site';
import { addAnisotropic, getAnisotropicTemplate, getAnisotropic } from './anisotropic';
import { parseConect } from './conect';
import { isDebugMode } from '../../../mol-util/debug';
import { PdbHeaderData, addHeader } from './header';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';

export async function pdbToMmCif(pdb: PdbFile): Promise<CifFrame> {
    const { lines } = pdb;
    const { data, indices } = lines;
    const tokenizer = Tokenizer(data);
    const isPdbqt = !!pdb.isPdbqt;

    // Count the atoms
    let atomCount = 0;
    let anisotropicCount = 0;
    for (let i = 0, _i = lines.count; i < _i; i++) {
        const s = indices[2 * i], e = indices[2 * i + 1];
        switch (data[s]) {
            case 'A':
                if (substringStartsWith(data, s, e, 'ATOM  ')) atomCount++;
                else if (substringStartsWith(data, s, e, 'ANISOU')) anisotropicCount++;
                break;
            case 'H':
                if (substringStartsWith(data, s, e, 'HETATM')) atomCount++;
                break;
        }
    }
    const header: PdbHeaderData = {};
    const atomSite = getAtomSiteTemplate(data, atomCount);
    const anisotropic = getAnisotropicTemplate(data, anisotropicCount);
    const entityBuilder = new EntityBuilder();
    const helperCategories: CifCategory[] = [];
    const heteroNames: [string, string][] = [];

    let modelNum = 0, modelStr = '';
    let conectRange: [number, number] | undefined = undefined;
    let hasAssemblies = false;
    const terIndices = new Set<number>();

    for (let i = 0, _i = lines.count; i < _i; i++) {
        let s = indices[2 * i], e = indices[2 * i + 1];
        switch (data[s]) {
            case 'A':
                if (substringStartsWith(data, s, e, 'ATOM  ')) {
                    if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                    addAtom(atomSite, modelStr, tokenizer, s, e, isPdbqt);
                } else if (substringStartsWith(data, s, e, 'ANISOU')) {
                    addAnisotropic(anisotropic, modelStr, tokenizer, s, e);
                }
                break;
            case 'C':
                if (substringStartsWith(data, s, e, 'CRYST1')) {
                    helperCategories.push(...parseCryst1(pdb.id || '?', data.substring(s, e)));
                } else if (substringStartsWith(data, s, e, 'CONECT')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'CONECT')) break;
                        j++;
                    }
                    if (conectRange) {
                        if (isDebugMode) {
                            console.log('only single CONECT block allowed, ignoring others');
                        }
                    } else {
                        conectRange = [i, j];
                    }
                    i = j - 1;
                } else if (substringStartsWith(data, s, e, 'COMPND')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'COMPND')) break;
                        j++;
                    }
                    entityBuilder.setCompounds(parseCmpnd(lines, i, j));
                    i = j - 1;
                }
                break;
            case 'H':
                if (substringStartsWith(data, s, e, 'HEADER')) {
                    addHeader(data, s, e, header);
                } else if (substringStartsWith(data, s, e, 'HETATM')) {
                    if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                    addAtom(atomSite, modelStr, tokenizer, s, e, isPdbqt);
                } else if (substringStartsWith(data, s, e, 'HELIX')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'HELIX')) break;
                        j++;
                    }
                    helperCategories.push(parseHelix(lines, i, j));
                    i = j - 1;
                } else if (substringStartsWith(data, s, e, 'HETNAM')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'HETNAM')) break;
                        j++;
                    }
                    heteroNames.push(...Array.from(parseHetnam(lines, i, j).entries()));
                    i = j - 1;
                }
                break;
            case 'M':
                if (substringStartsWith(data, s, e, 'MODEL ')) {
                    modelNum++;
                    modelStr = '' + modelNum;
                }
                if (substringStartsWith(data, s, e, 'MTRIX')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'MTRIX')) break;
                        j++;
                    }
                    helperCategories.push(...parseMtrix(lines, i, j));
                    i = j - 1;
                }
                // TODO: MODRES records => pdbx_struct_mod_residue
                break;
            case 'O':
                // TODO: ORIGX record => cif.database_PDB_matrix.origx, cif.database_PDB_matrix.origx_vector
                break;
            case 'R':
                if (substringStartsWith(data, s, e, 'REMARK 350')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'REMARK 350')) break;
                        j++;
                    }
                    helperCategories.push(...parseRemark350(lines, i, j));
                    i = j - 1;
                    hasAssemblies = true;
                }
                break;
            case 'S':
                if (substringStartsWith(data, s, e, 'SHEET')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'SHEET')) break;
                        j++;
                    }
                    helperCategories.push(parseSheet(lines, i, j));
                    i = j - 1;
                }
                // TODO: SCALE record => cif.atom_sites.fract_transf_matrix, cif.atom_sites.fract_transf_vector
                break;
            case 'T':
                if (substringStartsWith(data, s, e, 'TER')) {
                    terIndices.add(atomSite.index);
                }
        }
    }

    // build entry, struct_keywords and pdbx_database_status
    if (header.id_code) {
        const entry: CifCategory.SomeFields<mmCIF_Schema['entry']> = {
            id: CifField.ofString(header.id_code)
        };
        helperCategories.push(CifCategory.ofFields('entry', entry));
    }
    if (header.classification) {
        const struct_keywords: CifCategory.SomeFields<mmCIF_Schema['struct_keywords']> = {
            pdbx_keywords: CifField.ofString(header.classification)
        };
        helperCategories.push(CifCategory.ofFields('struct_keywords', struct_keywords));
    }
    if (header.dep_date) {
        const pdbx_database_status: CifCategory.SomeFields<mmCIF_Schema['pdbx_database_status']> = {
            recvd_initial_deposition_date: CifField.ofString(header.dep_date)
        };
        helperCategories.push(CifCategory.ofFields('pdbx_database_status', pdbx_database_status));
    }

    // build entity and chem_comp categories
    const seqIds = Column.ofIntTokens(atomSite.auth_seq_id);
    const atomIds = Column.ofStringTokens(atomSite.auth_atom_id);
    const compIds = Column.ofStringTokens(atomSite.auth_comp_id);
    const asymIds = Column.ofStringTokens(atomSite.auth_asym_id);
    const labelAsymIdHelper = new LabelAsymIdHelper(asymIds, atomSite.pdbx_PDB_model_num, terIndices, hasAssemblies);
    const componentBuilder = new ComponentBuilder(seqIds, atomIds);
    componentBuilder.setNames(heteroNames);
    entityBuilder.setNames(heteroNames);
    for (let i = 0, il = compIds.rowCount; i < il; ++i) {
        const compId = compIds.value(i);
        const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId);
        const asymId = labelAsymIdHelper.get(i);
        atomSite.label_entity_id[i] = entityBuilder.getEntityId(compId, moleculeType, asymId);
    }
    const atom_site = getAtomSite(atomSite, labelAsymIdHelper, { hasAssemblies });
    if (!isPdbqt) delete atom_site.partial_charge;

    if (conectRange) {
        helperCategories.push(parseConect(lines, conectRange[0], conectRange[1], atom_site));
    }

    const categories = {
        entity: CifCategory.ofTable('entity', entityBuilder.getEntityTable()),
        chem_comp: CifCategory.ofTable('chem_comp', componentBuilder.getChemCompTable()),
        atom_site: CifCategory.ofFields('atom_site', atom_site),
        atom_site_anisotrop: CifCategory.ofFields('atom_site_anisotrop', getAnisotropic(anisotropic))
    } as any;

    for (const c of helperCategories) {
        categories[c.name] = c;
    }

    return {
        header: pdb.id || 'PDB',
        categoryNames: Object.keys(categories),
        categories
    };
}