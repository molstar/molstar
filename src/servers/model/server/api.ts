/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Queries, Structure, StructureQuery, StructureSymmetry, StructureProperties } from '../../../mol-model/structure';
import { getAtomsTests } from '../query/atoms';
import { CifWriter } from '../../../mol-io/writer/cif';
import { QuerySchemas } from '../query/schemas';
import { Mat4 } from '../../../mol-math/linear-algebra';

export enum QueryParamType {
    JSON,
    String,
    Integer,
    Boolean,
    Float
}

export interface QueryParamInfo<T extends string | number = string | number> {
    name: string,
    type: QueryParamType,
    description?: string,
    required?: boolean,
    defaultValue?: any,
    exampleValues?: any[],
    validation?: (v: T) => void,
    supportedValues?: string[],
    groupName?: string
}

export interface QueryDefinition<Params = any> {
    name: string,
    niceName: string,
    exampleId: string, // default is 1cbs
    query: (params: Params, structure: Structure, numModels: number[]) => StructureQuery,
    description: string,
    jsonParams: QueryParamInfo[],
    restParams: QueryParamInfo[],
    structureTransform?: (params: any, s: Structure) => Promise<Structure>,
    filter?: CifWriter.Category.Filter,
    '@params': Params
}

export const CommonQueryParamsInfo: QueryParamInfo[] = [
    { name: 'model_nums', type: QueryParamType.String, description: `A comma-separated list of model ids (i.e. 1,2). If set, only include atoms with the corresponding '_atom_site.pdbx_PDB_model_num' field.` },
    { name: 'encoding', type: QueryParamType.String, defaultValue: 'cif', description: `Determines the output encoding (text based 'CIF' or binary 'BCIF'). Ligands can also be exported as 'SDF', 'MOL', or 'MOL2'.`, supportedValues: ['cif', 'bcif', 'sdf', 'mol', 'mol2'] },
    { name: 'copy_all_categories', type: QueryParamType.Boolean, defaultValue: false, description: 'If true, copy all categories from the input file.' },
    { name: 'data_source', type: QueryParamType.String, defaultValue: '', description: 'Allows to control how the provided data source ID maps to input file (as specified by the server instance config).' },
    { name: 'transform', type: QueryParamType.String, description: `Transformation to apply to coordinates in '_atom_site'. Accepts a 4x4 transformation matrix, provided as array of 16 float values.` },
    { name: 'download', type: QueryParamType.Boolean, defaultValue: false, description: 'If true, browser will download text files.' },
    { name: 'filename', type: QueryParamType.String, defaultValue: '', description: `Controls the filename for downloaded files. Will force download if specified.` }
];

export type Encoding = 'cif' | 'bcif' | 'sdf' | 'mol' | 'mol2';
export interface CommonQueryParamsInfo {
    model_nums?: number[],
    encoding?: Encoding,
    copy_all_categories?: boolean
    data_source?: string,
    transform?: Mat4,
    download?: boolean,
    filename?: string
}

export const AtomSiteSchemaElement = {
    label_entity_id: { type: QueryParamType.String, groupName: 'atom_site' },

    label_asym_id: { type: QueryParamType.String, groupName: 'atom_site' },
    auth_asym_id: { type: QueryParamType.String, groupName: 'atom_site' },

    label_comp_id: { type: QueryParamType.String, groupName: 'atom_site' },
    auth_comp_id: { type: QueryParamType.String, groupName: 'atom_site' },
    label_seq_id: { type: QueryParamType.Integer, groupName: 'atom_site' },
    auth_seq_id: { type: QueryParamType.Integer, groupName: 'atom_site' },
    pdbx_PDB_ins_code: { type: QueryParamType.String, groupName: 'atom_site' },

    label_atom_id: { type: QueryParamType.String, groupName: 'atom_site' },
    auth_atom_id: { type: QueryParamType.String, groupName: 'atom_site' },
    type_symbol: { type: QueryParamType.String, groupName: 'atom_site' }
};

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
    exampleValues: [[{ label_seq_id: 30, label_asym_id: 'A' }, { label_seq_id: 31, label_asym_id: 'A' }], { label_comp_id: 'ALA' }]
};

export const AtomSiteTestRestParams = (function () {
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
            throw new Error('Invalid radius for residue interaction query (must be a value between 1 and 10).');
        }
    }
};

const AssemblyNameParam: QueryParamInfo = {
    name: 'assembly_name',
    type: QueryParamType.String,
    description: 'Assembly name. If none is provided, crystal symmetry (where available) or deposited model is used.'
};

const OmitWaterParam: QueryParamInfo = {
    name: 'omit_water',
    type: QueryParamType.Boolean,
    required: false,
    defaultValue: false
};

function Q<Params = any>(definition: Partial<QueryDefinition<Params>>) {
    return definition;
}

const QueryMap = {
    'full': Q<{} | undefined>({ niceName: 'Full Structure', query: () => Queries.generators.all, description: 'The full structure.' }),
    'ligand': Q<{ atom_site: AtomSiteSchema }>({
        niceName: 'Ligand',
        description: 'Coordinates of the first group satisfying the given criteria.',
        query: (p, _s, numModels) => {
            const tests = getAtomsTests(p.atom_site);
            const ligands = Queries.combinators.merge(tests.map(test => Queries.generators.atoms({
                ...test,
                unitTest: ctx => StructureProperties.unit.model_num(ctx.element) === numModels[0],
                groupBy: ctx => StructureProperties.residue.key(ctx.element)
            })));
            return Queries.filters.first(ligands);
        },
        jsonParams: [AtomSiteTestJsonParam],
        restParams: AtomSiteTestRestParams
    }),
    'atoms': Q<{ atom_site: AtomSiteSchema }>({
        niceName: 'Atoms',
        description: 'Atoms satisfying the given criteria.',
        query: p => {
            return Queries.combinators.merge(getAtomsTests(p.atom_site).map(test => Queries.generators.atoms(test)));
        },
        jsonParams: [AtomSiteTestJsonParam],
        restParams: AtomSiteTestRestParams
    }),
    'symmetryMates': Q<{ radius: number }>({
        niceName: 'Symmetry Mates',
        description: 'Computes crystal symmetry mates within the specified radius.',
        query: () => Queries.generators.all,
        structureTransform(p, s) {
            return StructureSymmetry.builderSymmetryMates(s, p.radius).run();
        },
        jsonParams: [RadiusParam],
        filter: QuerySchemas.assembly
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
        }],
        filter: QuerySchemas.assembly
    }),
    'residueInteraction': Q<{ atom_site: AtomSiteSchema, radius: number, assembly_name: string }>({
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
            return Queries.modifiers.includeSurroundings(center, { radius: p.radius !== void 0 ? p.radius : 5, wholeResidues: true });
        },
        structureTransform(p, s) {
            if (p.assembly_name) return StructureSymmetry.buildAssembly(s, '' + p.assembly_name).run();
            return StructureSymmetry.builderSymmetryMates(s, p.radius !== void 0 ? p.radius : 5).run();
        },
        jsonParams: [AtomSiteTestJsonParam, RadiusParam, AssemblyNameParam],
        restParams: [...AtomSiteTestRestParams, RadiusParam, AssemblyNameParam],
        filter: QuerySchemas.interaction
    }),
    'residueSurroundings': Q<{ atom_site: AtomSiteSchema, radius: number }>({
        niceName: 'Residue Surroundings',
        description: 'Identifies all residues within the given radius from the source residue.',
        query(p) {
            const center = Queries.combinators.merge(getAtomsTests(p.atom_site).map(test => Queries.generators.atoms(test)));
            return Queries.modifiers.includeSurroundings(center, { radius: p.radius, wholeResidues: true });
        },
        jsonParams: [AtomSiteTestJsonParam, RadiusParam],
        restParams: [...AtomSiteTestRestParams, RadiusParam],
        filter: QuerySchemas.interaction
    }),
    'surroundingLigands': Q<{ atom_site: AtomSiteSchema, radius: number, assembly_name: string, omit_water: boolean }>({
        niceName: 'Surrounding Ligands',
        description: 'Identifies (complete) ligands within the given radius from the source atom set. Takes crystal symmetry into account.',
        query(p) {
            const tests = getAtomsTests(p.atom_site);
            const center = Queries.combinators.merge(tests.map(test => Queries.generators.atoms({
                ...test,
                entityTest: test.entityTest
                    ? ctx => test.entityTest!(ctx) && ctx.element.unit.conformation.operator.isIdentity
                    : ctx => ctx.element.unit.conformation.operator.isIdentity
            })));
            return Queries.modifiers.surroundingLigands({ query: center, radius: p.radius !== void 0 ? p.radius : 5, includeWater: !p.omit_water });
        },
        structureTransform(p, s) {
            if (p.assembly_name) return StructureSymmetry.buildAssembly(s, '' + p.assembly_name).run();
            return StructureSymmetry.builderSymmetryMates(s, p.radius !== void 0 ? p.radius : 5).run();
        },
        jsonParams: [AtomSiteTestJsonParam, RadiusParam, OmitWaterParam, AssemblyNameParam],
        restParams: [...AtomSiteTestRestParams, RadiusParam, OmitWaterParam, AssemblyNameParam],
        filter: QuerySchemas.interaction
    }),
};

export type QueryName = keyof typeof QueryMap
export type QueryParams<Q extends QueryName> = Partial<(typeof QueryMap)[Q]['@params']>

export function getQueryByName(name: QueryName): QueryDefinition {
    return QueryMap[name] as QueryDefinition;
}

export const QueryList = (function () {
    const list: { name: string, definition: QueryDefinition }[] = [];
    for (const k of Object.keys(QueryMap)) list.push({ name: k, definition: <QueryDefinition>QueryMap[k as QueryName] });
    list.sort(function (a, b) { return a.name < b.name ? -1 : a.name > b.name ? 1 : 0; });
    return list;
})();

// normalize the queries
(function () {
    for (const q of QueryList) {
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
        let el: any;
        if (typeof value === 'undefined' || (typeof value !== 'undefined' && value !== null && value['length'] === 0)) {
            if (p.required) {
                throw new Error(`The parameter '${key}' is required.`);
            }
            if (typeof p.defaultValue !== 'undefined') el = p.defaultValue;
        } else {
            switch (p.type) {
                case QueryParamType.JSON: el = JSON.parse(value); break;
                case QueryParamType.String: el = value; break;
                case QueryParamType.Integer: el = parseInt(value); break;
                case QueryParamType.Float: el = parseFloat(value); break;
                case QueryParamType.Boolean: el = Boolean(+value); break;
            }

            if (p.validation) p.validation(el);
        }
        if (typeof el === 'undefined') continue;

        if (p.groupName) {
            if (typeof ret[p.groupName] === 'undefined') ret[p.groupName] = {};
            ret[p.groupName][key] = el;
        } else {
            ret[key] = el;
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
        copy_all_categories: isTrue(params.copy_all_categories),
        encoding: mapEncoding(('' + params.encoding).toLocaleLowerCase()),
        transform: params.transform ? ('' + params.transform).split(',').map(n => n.trim()).map(n => +n) as Mat4 : Mat4.identity(),
        download: isTrue(params.download) || !!params.filename,
        filename: params.filename
    };
}

function isTrue(val: any): boolean {
    const b = Boolean(val);
    if (!b) return false;
    if (typeof val === 'string') return val !== '0' && val.toLowerCase() !== 'false';
    return b;
}

function mapEncoding(value: string) {
    switch (value) {
        case 'bcif':
            return 'bcif';
        case 'mol':
            return 'mol';
        case 'mol2':
            return 'mol2';
        case 'sdf':
            return 'sdf';
        default:
            return 'cif';
    }
}
