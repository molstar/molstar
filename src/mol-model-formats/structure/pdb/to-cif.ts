/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Yana Rose <yana.v.rose@gmail.com>
 */

import { CifCategory, CifField, CifFrame } from '../../../mol-io/reader/cif';
import { Tokenizer } from '../../../mol-io/reader/common/text/tokenizer';
import { PdbFile } from '../../../mol-io/reader/pdb/schema';
import { parseCryst1, parseRemark350, parseMtrix } from './assembly';
import { parseHelix, parseSheet } from './secondary-structure';
import { parseCmpnd, parseHetnam, parseSeqres } from './entity';
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
import { StringLike } from '../../../mol-io/common/string-like';

function substringStartsWith(str: StringLike, start: number, end: number, target: string) {
    const len = target.length;
    if (len > end - start) return false;
    for (let i = 0; i < len; i++) {
        if (str.charCodeAt(start + i) !== target.charCodeAt(i)) return false;
    }
    return true;
}

export async function pdbToMmCif(pdb: PdbFile): Promise<CifFrame> {
    const { lines, variant } = pdb;
    const { data, indices } = lines;
    const tokenizer = Tokenizer(data);

    // Count the atoms
    let atomCount = 0;
    let anisotropicCount = 0;
    for (let i = 0, _i = lines.count; i < _i; i++) {
        const s = indices[2 * i], e = indices[2 * i + 1];
        switch (data.charAt(s)) {
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
    let hasSeqRes = false;
    let seqresMap: Map<string, string[]> | undefined = undefined;
    const terIndices = new Set<number>();

    for (let i = 0, _i = lines.count; i < _i; i++) {
        let s = indices[2 * i], e = indices[2 * i + 1];
        switch (data.charAt(s)) {
            case 'A':
                if (substringStartsWith(data, s, e, 'ATOM  ')) {
                    if (!modelNum) { modelNum++; modelStr = '' + modelNum; }
                    addAtom(atomSite, modelStr, tokenizer, s, e, variant);
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
                    addAtom(atomSite, modelStr, tokenizer, s, e, variant);
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
                } else if (substringStartsWith(data, s, e, 'SEQRES')) {
                    let j = i + 1;
                    while (true) {
                        s = indices[2 * j]; e = indices[2 * j + 1];
                        if (!substringStartsWith(data, s, e, 'SEQRES')) break;
                        j++;
                    }
                    seqresMap = parseSeqres(lines, i, j);
                    entityBuilder.setSeqres(seqresMap);
                    i = j - 1;
                    hasSeqRes = true;
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
    const atom_site = getAtomSite(atomSite, labelAsymIdHelper, { hasAssemblies, hasSeqRes, seqresMap });
    if (variant === 'pdb') delete atom_site.partial_charge;

    if (conectRange) {
        helperCategories.push(parseConect(lines, conectRange[0], conectRange[1], atom_site));
    }

    const categories = {
        entity: CifCategory.ofTable('entity', entityBuilder.getEntityTable()),
        chem_comp: CifCategory.ofTable('chem_comp', componentBuilder.getChemCompTable()),
        atom_site: CifCategory.ofFields('atom_site', atom_site),
        atom_site_anisotrop: CifCategory.ofFields('atom_site_anisotrop', getAnisotropic(anisotropic))
    } as any;

    // Build entity_poly_seq from SEQRES
    if (seqresMap) {
        const epsEntityIds: string[] = [];
        const epsNums: number[] = [];
        const epsMonIds: string[] = [];
        const epsHeteros: string[] = [];
        const processedEntities = new Set<string>();

        for (const [chainId, residues] of seqresMap) {
            const entityId = entityBuilder.getEntityIdForChain(chainId);
            if (!entityId || processedEntities.has(entityId)) continue;
            processedEntities.add(entityId);

            for (let j = 0; j < residues.length; j++) {
                epsEntityIds.push(entityId);
                epsNums.push(j + 1);
                epsMonIds.push(residues[j]);
                epsHeteros.push('no');
            }
        }

        if (epsEntityIds.length > 0) {
            const entity_poly_seq: CifCategory.SomeFields<mmCIF_Schema['entity_poly_seq']> = {
                entity_id: CifField.ofStrings(epsEntityIds),
                num: CifField.ofNumbers(epsNums),
                mon_id: CifField.ofStrings(epsMonIds),
                hetero: CifField.ofStrings(epsHeteros),
            };
            categories.entity_poly_seq = CifCategory.ofFields('entity_poly_seq', entity_poly_seq);
        }

        // Build pdbx_unobs_or_zero_occ_residues by comparing SEQRES with observed ATOM records
        // Collect observed (label_asym_id, label_seq_id) pairs per model, and auth_asym_id -> label_asym_id mapping.
        // Only include atoms belonging to the polymer entity for each SEQRES chain.
        // Non-polymer HETATMs (ions, ligands, water) share auth_asym_id in PDB format
        // and their sequential label_seq_id values can collide with unobserved SEQRES positions.
        const polymerEntityIds = new Map<string, string>();
        for (const chainId of seqresMap.keys()) {
            const entityId = entityBuilder.getEntityIdForChain(chainId);
            if (entityId) polymerEntityIds.set(chainId, entityId);
        }

        const observedResidues = new Set<string>();
        const authToLabelAsym = new Map<string, string>();
        const rowCount = atom_site.label_asym_id!.rowCount;
        for (let i = 0; i < rowCount; ++i) {
            const authAsym = atom_site.auth_asym_id!.str(i);
            const entityId = atom_site.label_entity_id!.str(i);

            // Skip non-polymer atoms: their label_seq_id can collide with SEQRES positions
            const polymerEntityId = polymerEntityIds.get(authAsym);
            if (polymerEntityId && entityId !== polymerEntityId) continue;

            const labelAsym = atom_site.label_asym_id!.str(i);
            const labelSeq = atom_site.label_seq_id!.int(i);
            const modelN = atom_site.pdbx_PDB_model_num!.int(i);
            observedResidues.add(`${modelN}|${labelAsym}|${labelSeq}`);
            if (!authToLabelAsym.has(authAsym)) {
                authToLabelAsym.set(authAsym, labelAsym);
            }
        }

        const unobsIds: number[] = [];
        const unobsModelNums: number[] = [];
        const unobsLabelAsymIds: string[] = [];
        const unobsLabelCompIds: string[] = [];
        const unobsLabelSeqIds: number[] = [];
        const unobsAuthAsymIds: string[] = [];
        const unobsAuthCompIds: string[] = [];
        const unobsAuthSeqIds: number[] = [];
        const unobsPolymerFlags: string[] = [];
        const unobsOccupancyFlags: number[] = [];
        let unobsId = 0;
        const modelCount = Math.max(modelNum, 1);

        for (let m = 1; m <= modelCount; ++m) {
            for (const [chainId, residues] of seqresMap) {
                const labelAsymId = authToLabelAsym.get(chainId);
                if (!labelAsymId) continue;

                for (let j = 0; j < residues.length; j++) {
                    const seqId = j + 1; // 1-based label_seq_id
                    if (!observedResidues.has(`${m}|${labelAsymId}|${seqId}`)) {
                        unobsId++;
                        unobsIds.push(unobsId);
                        unobsModelNums.push(m);
                        unobsLabelAsymIds.push(labelAsymId);
                        unobsLabelCompIds.push(residues[j]);
                        unobsLabelSeqIds.push(seqId);
                        unobsAuthAsymIds.push(chainId);
                        unobsAuthCompIds.push(residues[j]);
                        unobsAuthSeqIds.push(seqId);
                        unobsPolymerFlags.push('y');
                        unobsOccupancyFlags.push(1);
                    }
                }
            }
        }

        if (unobsIds.length > 0) {
            const pdbx_unobs: CifCategory.SomeFields<mmCIF_Schema['pdbx_unobs_or_zero_occ_residues']> = {
                id: CifField.ofNumbers(unobsIds),
                PDB_model_num: CifField.ofNumbers(unobsModelNums),
                polymer_flag: CifField.ofStrings(unobsPolymerFlags),
                occupancy_flag: CifField.ofNumbers(unobsOccupancyFlags),
                label_asym_id: CifField.ofStrings(unobsLabelAsymIds),
                label_comp_id: CifField.ofStrings(unobsLabelCompIds),
                label_seq_id: CifField.ofNumbers(unobsLabelSeqIds),
                auth_asym_id: CifField.ofStrings(unobsAuthAsymIds),
                auth_comp_id: CifField.ofStrings(unobsAuthCompIds),
                auth_seq_id: CifField.ofNumbers(unobsAuthSeqIds),
            };
            categories.pdbx_unobs_or_zero_occ_residues = CifCategory.ofFields('pdbx_unobs_or_zero_occ_residues', pdbx_unobs);
        }
    }

    for (const c of helperCategories) {
        categories[c.name] = c;
    }

    return {
        header: pdb.id || variant.toUpperCase(),
        categoryNames: Object.keys(categories),
        categories
    };
}