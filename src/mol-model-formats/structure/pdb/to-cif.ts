/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { substringStartsWith } from '../../../mol-util/string';
import { CifCategory, CifFrame } from '../../../mol-io/reader/cif';
import { Tokenizer } from '../../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../../mol-io/reader/pdb/schema';
import { parseCryst1, parseRemark350, parseMtrix } from './assembly';
import { parseHelix, parseSheet } from './secondary-structure';
import { parseCmpnd, parseHetnam } from './entity';
import { ComponentBuilder } from '../common/component';
import { EntityBuilder } from '../common/entity';
import { Column } from '../../../mol-data/db';
import { getMoleculeType } from '../../../mol-model/structure/model/types';
import { getAtomSiteTemplate, addAtom, getAtomSite } from './atom-site';
import { addAnisotropic, getAnisotropicTemplate, getAnisotropic } from './anisotropic';

export async function pdbToMmCif(pdb: PdbFile): Promise<CifFrame> {
    const { lines } = pdb;
    const { data, indices } = lines;
    const tokenizer = Tokenizer(data);

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

    const atomSite = getAtomSiteTemplate(data, atomCount);
    const anisotropic = getAnisotropicTemplate(data, anisotropicCount);
    const entityBuilder = new EntityBuilder();
    const helperCategories: CifCategory[] = [];
    const heteroNames: [string, string][] = [];

    let modelNum = 0, modelStr = '';

    for (let i = 0, _i = lines.count; i < _i; i++) {
        let s = indices[2 * i], e = indices[2 * i + 1];
        switch (data[s]) {
            case 'A':
                if (substringStartsWith(data, s, e, 'ATOM  ')) {
                    if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                    addAtom(atomSite, modelStr, tokenizer, s, e);
                } else if (substringStartsWith(data, s, e, 'ANISOU')) {
                    addAnisotropic(anisotropic, modelStr, tokenizer, s, e);
                }
                break;
            case 'C':
                if (substringStartsWith(data, s, e, 'CRYST1')) {
                    helperCategories.push(...parseCryst1(pdb.id || '?', data.substring(s, e)));
                } else if (substringStartsWith(data, s, e, 'CONNECT')) {
                    // TODO: CONNECT records => struct_conn
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
                if (substringStartsWith(data, s, e, 'HETATM')) {
                    if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                    addAtom(atomSite, modelStr, tokenizer, s, e);
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
        }
    }

    // build entity and chem_comp categories
    const seqIds = Column.ofIntTokens(atomSite.auth_seq_id);
    const atomIds = Column.ofStringTokens(atomSite.auth_atom_id);
    const compIds = Column.ofStringTokens(atomSite.auth_comp_id);
    const asymIds = Column.ofStringTokens(atomSite.auth_asym_id);
    const componentBuilder = new ComponentBuilder(seqIds, atomIds);
    componentBuilder.setNames(heteroNames);
    entityBuilder.setNames(heteroNames);
    for (let i = 0, il = compIds.rowCount; i < il; ++i) {
        const compId = compIds.value(i);
        const moleculeType = getMoleculeType(componentBuilder.add(compId, i).type, compId);
        atomSite.label_entity_id[i] = entityBuilder.getEntityId(compId, moleculeType, asymIds.value(i));
    }

    const categories = {
        entity: CifCategory.ofTable('entity', entityBuilder.getEntityTable()),
        chem_comp: CifCategory.ofTable('chem_comp', componentBuilder.getChemCompTable()),
        atom_site: CifCategory.ofFields('atom_site', getAtomSite(atomSite)),
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