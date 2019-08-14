/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Queries, Structure, StructureQuery, StructureSymmetry } from '../../../mol-model/structure';
import { getAtomsTests } from '../query/atoms';
import { CifWriter } from '../../../mol-io/writer/cif';
import { QuerySchemas } from '../query/schemas';

export enum QueryParamType {
    JSON,
    String,
    Integer,
    Float,
    Boolean
}

export interface QueryParamInfo<T extends string | number = string | number> {
    name: string,
    type: QueryParamType,
    description?: string,
    required?: boolean,
    defaultValue?: any,
    exampleValues?: any[],
    validation?: (v: T) => void,
    supportedValues?: string[]
}

export interface QueryDefinition<Params = any> {
    name: string,
    niceName: string,
    exampleId: string, // default is 1cbs
    query: (params: any, structure: Structure) => StructureQuery,
    description: string,
    jsonParams: QueryParamInfo[],
    restParams: QueryParamInfo[],
    structureTransform?: (params: any, s: Structure) => Promise<Structure>,
    filter?: CifWriter.Category.Filter,
    '@params': Params
}

export const CommonQueryParamsInfo: QueryParamInfo[] = [
    { name: 'model_nums', type: QueryParamType.String, description: `A comma-separated list of model ids (i.e. 1,2). If set, only include atoms with the corresponding '_atom_site.pdbx_PDB_model_num' field.` },
    { name: 'encoding', type: QueryParamType.String, defaultValue: 'cif', description: `Determines the output encoding (text based 'CIF' or binary 'BCIF').`, supportedValues: ['cif', 'bcif'] },
    { name: 'data_Source', type: QueryParamType.String, defaultValue: '', description: 'Allows to control how the provided data source ID maps to input file (as specified by the server instance config).' }
];

export interface CommonQueryParamsInfo {
    model_nums?: number[],
    encoding?: 'cif' | 'bcif',
    data_source?: string
}

export const AtomSiteSchemaElement = {
    label_entity_id: { type: QueryParamType.String },

    label_asym_id: { type: QueryParamType.String },
    auth_asym_id: { type: QueryParamType.String },

    label_comp_id: { type: QueryParamType.String },
    auth_comp_id: { type: QueryParamType.String },
    label_seq_id: { type: QueryParamType.Integer },
    auth_seq_id: { type: QueryParamType.Integer },
    pdbx_PDB_ins_code: { type: QueryParamType.String },

    label_atom_id: { type: QueryParamType.String },
    auth_atom_id: { type: QueryParamType.String },
    type_symbol: { type: QueryParamType.String }
}

export interface AtomSiteSchemaElement {
    label_entity_id?: string,

    label_asym_id?: string,
    auth_asym_id?: string,

    label_comp_id?: string,
    auth_comp_id?: string,
    label_seq_id?: number,
    auth_seq_id?: number,
    pdbx_PDB_ins_code?: string,

    label_atom_id?: string,
    auth_atom_id?: string,
    type_symbol?: string
}

export type AtomSiteSchema = AtomSiteSchemaElement | AtomSiteSchemaElement[]

const AtomSiteTestJsonParam: QueryParamInfo = {
    name: 'atom_site',
    type: QueryParamType.JSON,
    description: 'Object or array of objects describing atom properties. Names are same as in wwPDB mmCIF dictionary of the atom_site category.',
    exampleValues: [{ label_comp_id: 'ALA' }, { label_seq_id: 123, label_asym_id: 'A' }]
};

export const AtomSiteTestRestParams = (function() {
    const params: QueryParamInfo[] = [];
    for (const k of Object.keys(AtomSiteSchemaElement)) {
        const p = (AtomSiteSchemaElement as any)[k] as QueryParamInfo;
        p.name = k;
        params.push(p);
    }
    return params;
})();

const RadiusParam: QueryParamInfo = {
    name: 'radius',
    type: QueryParamType.Float,
    defaultValue: 5,
    exampleValues: [5],
    description: 'Value in Angstroms.',
    validation(v: any) {
        if (v < 1 || v > 10) {
            throw `Invalid radius for residue interaction query (must be a value between 1 and 10).`;
        }
    }
};

function Q<Params = any>(definition: Partial<QueryDefinition<Params>>) {
    return definition;
}

const QueryMap = {
    'full': Q<{} | undefined>({ niceName: 'Full Structure', query: () => Queries.generators.all, description: 'The full structure.' }),
    'atoms': Q<{ atom_site: AtomSiteSchema }>({
        niceName: 'Atoms',
        description: 'Atoms satisfying the given criteria.',
        query: p => Queries.combinators.merge(getAtomsTests(p).map(test => Queries.generators.atoms(test))),
        jsonParams: [ AtomSiteTestJsonParam ],
        restParams: AtomSiteTestRestParams
    }),
    'symmetryMates': Q<{ radius: number }>({
        niceName: 'Symmetry Mates',
        description: 'Computes crystal symmetry mates within the specified radius.',
        query: () => Queries.generators.all,
        structureTransform(p, s) {
            return StructureSymmetry.builderSymmetryMates(s, p.radius).run();
        },
        jsonParams: [ RadiusParam ]
    }),
    'assembly': Q<{ name: string }>({
        niceName: 'Assembly',
        description: 'Computes structural assembly.',
        query: () => Queries.generators.all,
        structureTransform(p, s) {
            return StructureSymmetry.buildAssembly(s, '' + (p.name || '1')).run();
        },
        jsonParams: [{
            name: 'name',
            type: QueryParamType.String,
            defaultValue: '1',
            exampleValues: ['1'],
            description: 'Assembly name.'
        }]
    }),
    'residueInteraction': Q<AtomSiteSchema & { radius: number }>({
        niceName: 'Residue Interaction',
        description: 'Identifies all residues within the given radius from the source residue. Takes crystal symmetry into account.',
        query(p) {
            const tests = getAtomsTests(p);
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
        jsonParams: [ AtomSiteTestJsonParam, RadiusParam ],
        restParams: [ ...AtomSiteTestRestParams, RadiusParam ],
        filter: QuerySchemas.interaction
    }),
    'residueSurroundings': Q<AtomSiteSchema & { radius: number }>({
        niceName: 'Residue Surroundings',
        description: 'Identifies all residues within the given radius from the source residue.',
        query(p) {
            const tests = getAtomsTests(p);
            const center = Queries.combinators.merge(tests.map(test => Queries.generators.atoms(test)));
            return Queries.modifiers.includeSurroundings(center, { radius: p.radius, wholeResidues: true });
        },
        jsonParams: [ AtomSiteTestJsonParam, RadiusParam ],
        restParams: [ ...AtomSiteTestRestParams, RadiusParam ],
        filter: QuerySchemas.interaction
    })
};

export type QueryName = keyof typeof QueryMap
export type QueryParams<Q extends QueryName> = Partial<(typeof QueryMap)[Q]['@params']>

export function getQueryByName(name: QueryName): QueryDefinition {
    return QueryMap[name] as QueryDefinition;
}

export const QueryList = (function () {
    const list: { name: string, definition: QueryDefinition }[] = [];
    for (const k of Object.keys(QueryMap)) list.push({ name: k, definition: <QueryDefinition>QueryMap[k as QueryName] });
    list.sort(function (a, b) { return a.name < b.name ? -1 : a.name > b.name ? 1 : 0 });
    return list;
})();

// normalize the queries
(function () {
    for (let q of QueryList) {
        const m = q.definition;
        m.name = q.name;
        m.jsonParams = m.jsonParams || [];
        m.restParams = m.restParams || m.jsonParams;
    }
})();

function _normalizeQueryParams(params: { [p: string]: string }, paramList: QueryParamInfo[]): { [p: string]: string | number | boolean } {
    const ret: any = {};
    for (const p of paramList) {
        const key = p.name;
        const value = params[key];
        if (typeof value === 'undefined' || (typeof value !== 'undefined' && value !== null && value['length'] === 0)) {
            if (p.required) {
                throw `The parameter '${key}' is required.`;
            }
            if (typeof p.defaultValue !== 'undefined') ret[key] = p.defaultValue;
        } else {
            switch (p.type) {
                case QueryParamType.JSON: ret[key] = JSON.parse(value); break;
                case QueryParamType.String: ret[key] = value; break;
                case QueryParamType.Integer: ret[key] = parseInt(value); break;
                case QueryParamType.Float: ret[key] = parseFloat(value); break;
            }

            if (p.validation) p.validation(ret[key]);
        }
    }

    return ret;
}

export function normalizeRestQueryParams(query: QueryDefinition, params: any) {
    // return params;
    return _normalizeQueryParams(params, query.restParams);
}

export function normalizeRestCommonParams(params: any): CommonQueryParamsInfo {
    return {
        model_nums: params.model_nums ? ('' + params.model_nums).split(',').map(n => n.trim()).filter(n => !!n).map(n => +n) : void 0,
        data_source: params.data_source,
        encoding: ('' + params.encoding).toLocaleLowerCase() === 'bcif' ? 'bcif' : 'cif'
    };
}