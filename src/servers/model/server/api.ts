/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Queries, Structure, StructureQuery, StructureSymmetry } from 'mol-model/structure';
import { getAtomsTests } from '../query/atoms';

export enum QueryParamType {
    JSON,
    String,
    Integer,
    Float
}

export interface QueryParamInfo {
    name: string,
    type: QueryParamType,
    description?: string,
    required?: boolean,
    defaultValue?: any,
    exampleValues?: string[],
    validation?: (v: any) => void
}

export interface QueryDefinition {
    name: string,
    niceName: string,
    exampleId: string, // default is 1cbs
    query: (params: any, structure: Structure) => StructureQuery,
    description: string,
    params: QueryParamInfo[],
    structureTransform?: (params: any, s: Structure) => Promise<Structure>
}

// const AtomSiteParams = {
//     entity_id: <QueryParamInfo>{ name: 'entity_id', type: QueryParamType.String, description: 'Corresponds to the \'_entity.id\' or \'*.label_entity_id\' field, depending on the context.' },

//     label_asym_id: <QueryParamInfo>{ name: 'label_asym_id', type: QueryParamType.String, description: 'Corresponds to the \'_atom_site.label_asym_id\' field.' },
//     auth_asym_id: <QueryParamInfo>{ name: 'auth_asym_id', type: QueryParamType.String, exampleValue: 'A', description: 'Corresponds to the \'_atom_site.auth_asym_id\' field.' },

//     label_seq_id: <QueryParamInfo>{ name: 'label_seq_id', type: QueryParamType.Integer, description: 'Residue seq. number. Corresponds to the \'_atom_site.label_seq_id\' field.' },
//     auth_seq_id: <QueryParamInfo>{ name: 'auth_seq_id', type: QueryParamType.Integer, exampleValue: '200', description: 'Author residue seq. number. Corresponds to the \'_atom_site.auth_seq_id\' field.' },
//     label_comp_id: <QueryParamInfo>{ name: 'label_comp_id', type: QueryParamType.String, description: 'Residue name. Corresponds to the \'_atom_site.label_comp_id\' field.' },
//     auth_comp_id: <QueryParamInfo>{ name: 'auth_comp_id', type: QueryParamType.String, exampleValue: 'REA', description: 'Author residue name. Corresponds to the \'_atom_site.auth_comp_id\' field.' },
//     pdbx_PDB_ins_code: <QueryParamInfo>{ name: 'pdbx_PDB_ins_code', type: QueryParamType.String, description: 'Corresponds to the \'_atom_site.pdbx_PDB_ins_code\' field.' },
// };

const AtomSiteTestParams: QueryParamInfo = {
    name: 'atom_site',
    type: QueryParamType.JSON,
    description: 'Object or array of objects describing atom properties. Name are same as in wwPDB mmCIF dictionary of the atom_site category.',
    exampleValues: [`{ label_comp_id: 'ALA' }`, `{ label_seq_id: 123, label_asym_id: 'A' }`]
};

const RadiusParam: QueryParamInfo = {
    name: 'radius',
    type: QueryParamType.Float,
    defaultValue: 5,
    exampleValues: ['5'],
    description: 'Value in Angstroms.',
    validation(v: any) {
        if (v < 1 || v > 10) {
            throw `Invalid radius for residue interaction query (must be a value between 1 and 10).`;
        }
    }
};

const QueryMap: { [id: string]: Partial<QueryDefinition> } = {
    'full': { niceName: 'Full Structure', query: () => Queries.generators.all, description: 'The full structure.' },
    'atoms': {
        niceName: 'Atoms',
        description: 'Atoms satisfying the given criteria.',
        query: p => Queries.combinators.merge(getAtomsTests(p.atom_site).map(test => Queries.generators.atoms(test))),
        params: [ AtomSiteTestParams ]
    },
    'symmetryMates': {
        niceName: 'Symmetry Mates',
        description: 'Computes crystal symmetry mates within the specified radius.',
        query: () => Queries.generators.all,
        structureTransform(p, s) {
            return StructureSymmetry.builderSymmetryMates(s, p.radius).run();
        },
    },
    'assembly': {
        niceName: 'Assembly',
        description: 'Computes structural assembly.',
        query: () => Queries.generators.all,
        structureTransform(p, s) {
            return StructureSymmetry.buildAssembly(s, '' + (p.name || '1')).run();
        },
        params: [{
            name: 'name',
            type: QueryParamType.String,
            defaultValue: '1',
            exampleValues: ['1'],
            description: 'Assembly name.'
        }]
    },
    'residueInteraction': {
        niceName: 'Residue Interaction',
        description: 'Identifies all residues within the given radius from the source residue. Takes crystal symmetry into account.',
        query(p) {
            const tests = getAtomsTests(p.atom_site);
            const center = Queries.combinators.merge(tests.map(test => Queries.generators.atoms({
                ...test,
                entityTest: test.entityTest
                    ? ctx => test.entityTest!(ctx) && ctx.element.unit.conformation.operator.isIdentity
                    : ctx => ctx.element.unit.conformation.operator.isIdentity
            })));
            return Queries.modifiers.includeSurroundings(center, { radius: p.radius, wholeResidues: true });
        },
        structureTransform(p, s) {
            return StructureSymmetry.builderSymmetryMates(s, p.radius).run();
        },
        params: [ AtomSiteTestParams, RadiusParam ]
    },
};

export function getQueryByName(name: string): QueryDefinition {
    return QueryMap[name] as QueryDefinition;
}

export const QueryList = (function () {
    const list: { name: string, definition: QueryDefinition }[] = [];
    for (const k of Object.keys(QueryMap)) list.push({ name: k, definition: <QueryDefinition>QueryMap[k] });
    list.sort(function (a, b) { return a.name < b.name ? -1 : a.name > b.name ? 1 : 0 });
    return list;
})();

// normalize the queries
(function () {
    for (let q of QueryList) {
        const m = q.definition;
        m.name = q.name;
        m.params = m.params || [];
    }
})();

// function _normalizeQueryParams(params: { [p: string]: string }, paramList: QueryParamInfo[]): { [p: string]: string | number | boolean } {
//     const ret: any = {};
//     for (const p of paramList) {
//         const key = p.name;
//         const value = params[key];
//         if (typeof value === 'undefined' || (typeof value !== 'undefined' && value !== null && value['length'] === 0)) {
//             if (p.required) {
//                 throw `The parameter '${key}' is required.`;
//             }
//             if (typeof p.defaultValue !== 'undefined') ret[key] = p.defaultValue;
//         } else {
//             switch (p.type) {
//                 case QueryParamType.JSON: ret[key] = JSON.parse(value); break;
//                 case QueryParamType.String: ret[key] = value; break;
//                 case QueryParamType.Integer: ret[key] = parseInt(value); break;
//                 case QueryParamType.Float: ret[key] = parseFloat(value); break;
//             }

//             if (p.validation) p.validation(ret[key]);
//         }
//     }

//     return ret;
// }

export function normalizeQueryParams(query: QueryDefinition, params: any) {
    return params;
    //return _normalizeQueryParams(params, query.params);
}