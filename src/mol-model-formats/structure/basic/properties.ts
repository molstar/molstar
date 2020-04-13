/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Model } from '../../../mol-model/structure/model/model';
import { ChemicalComponent, MissingResidue, StructAsym } from '../../../mol-model/structure/model/properties/common';
import { getMoleculeType, MoleculeType, getDefaultChemicalComponent } from '../../../mol-model/structure/model/types';
import { SaccharideComponentMap, SaccharideComponent, SaccharidesSnfgMap, SaccharideCompIdMap, UnknownSaccharideComponent } from '../../../mol-model/structure/structure/carbohydrates/constants';
import { memoize1 } from '../../../mol-util/memoize';
import { BasicData } from './schema';
import { Table } from '../../../mol-data/db';

function getMissingResidues(data: BasicData): Model['properties']['missingResidues'] {
    const map = new Map<string, MissingResidue>();
    const getKey = (model_num: number, asym_id: string, seq_id: number) => {
        return `${model_num}|${asym_id}|${seq_id}`;
    };

    const c = data.pdbx_unobs_or_zero_occ_residues;
    for (let i = 0, il = c._rowCount; i < il; ++i) {
        const key = getKey(c.PDB_model_num.value(i), c.label_asym_id.value(i), c.label_seq_id.value(i));
        map.set(key, { polymer_flag: c.polymer_flag.value(i), occupancy_flag: c.occupancy_flag.value(i) });
    }

    return {
        has: (model_num: number, asym_id: string, seq_id: number) => {
            return map.has(getKey(model_num, asym_id, seq_id));
        },
        get: (model_num: number, asym_id: string, seq_id: number) => {
            return map.get(getKey(model_num, asym_id, seq_id));
        },
        size: map.size
    };
}

function getChemicalComponentMap(data: BasicData): Model['properties']['chemicalComponentMap'] {
    const map = new Map<string, ChemicalComponent>();

    if (data.chem_comp._rowCount > 0) {
        const { id } = data.chem_comp;
        for (let i = 0, il = id.rowCount; i < il; ++i) {
            map.set(id.value(i), Table.getRow(data.chem_comp, i));
        }
    } else {
        const uniqueNames = getUniqueComponentNames(data);
        uniqueNames.forEach(n => {
            map.set(n, getDefaultChemicalComponent(n));
        });
    }
    return map;
}

function getSaccharideComponentMap(data: BasicData): SaccharideComponentMap {
    const map = new Map<string, SaccharideComponent>();

    if (data.pdbx_chem_comp_identifier._rowCount > 0) {
        // note that `pdbx_chem_comp_identifier` does not contain
        // a 'SNFG CARBOHYDRATE SYMBOL' entry for 'Unknown' saccharide components
        // so we always need to check `chem_comp` for those
        const { comp_id, type, identifier } = data.pdbx_chem_comp_identifier;
        for (let i = 0, il = comp_id.rowCount; i < il; ++i) {
            if (type.value(i) === 'SNFG CARBOHYDRATE SYMBOL' ||
                type.value(i) === 'SNFG CARB SYMBOL' // legacy, to be removed from mmCIF dictionary
            ) {
                const snfgName = identifier.value(i);
                const saccharideComp = SaccharidesSnfgMap.get(snfgName);
                if (saccharideComp) {
                    map.set(comp_id.value(i), saccharideComp);
                } else {
                    console.warn(`Unknown SNFG name '${snfgName}'`);
                }
            }
        }
    }

    if (data.chem_comp._rowCount > 0) {
        const { id, type  } = data.chem_comp;
        for (let i = 0, il = id.rowCount; i < il; ++i) {
            const _id = id.value(i);
            if (map.has(_id)) continue;
            const _type = type.value(i);
            if (SaccharideCompIdMap.has(_id)) {
                map.set(_id, SaccharideCompIdMap.get(_id)!);
            } else if (getMoleculeType(_type, _id) === MoleculeType.Saccharide) {
                map.set(_id, UnknownSaccharideComponent);
            }
        }
    } else {
        const uniqueNames = getUniqueComponentNames(data);
        SaccharideCompIdMap.forEach((v, k) => {
            if (!map.has(k) && uniqueNames.has(k)) map.set(k, v);
        });
    }
    return map;
}

const getUniqueComponentNames = memoize1((data: BasicData) => {
    const uniqueNames = new Set<string>();
    const { label_comp_id, auth_comp_id } = data.atom_site;
    const comp_id = label_comp_id.isDefined ? label_comp_id : auth_comp_id;
    for (let i = 0, il = comp_id.rowCount; i < il; ++i) {
        uniqueNames.add(comp_id.value(i));
    }
    return uniqueNames;
});


function getStructAsymMap(data: BasicData): Model['properties']['structAsymMap'] {
    const map = new Map<string, StructAsym>();

    const { label_asym_id, auth_asym_id, label_entity_id } = data.atom_site;
    for (let i = 0, il = label_asym_id.rowCount; i < il; ++i) {
        const id = label_asym_id.value(i);
        if (!map.has(id)) {
            map.set(id, {
                id,
                auth_id: auth_asym_id.value(i),
                entity_id: label_entity_id.value(i)
            });
        }
    }

    if (data.struct_asym._rowCount > 0) {
        const { id, entity_id } = data.struct_asym;
        for (let i = 0, il = id.rowCount; i < il; ++i) {
            const _id = id.value(i);
            if (!map.has(_id)) {
                map.set(_id, {
                    id: _id,
                    auth_id: '',
                    entity_id: entity_id.value(i)
                });
            }
        }
    }
    return map;
}

export function getProperties(data: BasicData): Model['properties'] {
    return {
        missingResidues: getMissingResidues(data),
        chemicalComponentMap: getChemicalComponentMap(data),
        saccharideComponentMap: getSaccharideComponentMap(data),
        structAsymMap: getStructAsymMap(data)
    };
}