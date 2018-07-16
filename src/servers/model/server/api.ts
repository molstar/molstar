/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery, Queries, Structure, StructureElement, StructureSymmetry, StructureProperties as Props, QueryPredicate } from 'mol-model/structure';

export enum QueryParamType {
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
    exampleValue?: string,
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

const AtomSiteParameters = {
    entity_id: <QueryParamInfo>{ name: 'entity_id', type: QueryParamType.String, description: 'Corresponds to the \'_entity.id\' or \'*.label_entity_id\' field, depending on the context.' },

    label_asym_id: <QueryParamInfo>{ name: 'label_asym_id', type: QueryParamType.String, description: 'Corresponds to the \'_atom_site.label_asym_id\' field.' },
    auth_asym_id: <QueryParamInfo>{ name: 'auth_asym_id', type: QueryParamType.String, exampleValue: 'A', description: 'Corresponds to the \'_atom_site.auth_asym_id\' field.' },

    label_seq_id: <QueryParamInfo>{ name: 'label_seq_id', type: QueryParamType.Integer, description: 'Residue seq. number. Corresponds to the \'_atom_site.label_seq_id\' field.' },
    auth_seq_id: <QueryParamInfo>{ name: 'auth_seq_id', type: QueryParamType.Integer, exampleValue: '200', description: 'Author residue seq. number. Corresponds to the \'_atom_site.auth_seq_id\' field.' },
    label_comp_id: <QueryParamInfo>{ name: 'label_comp_id', type: QueryParamType.String, description: 'Residue name. Corresponds to the \'_atom_site.label_comp_id\' field.' },
    auth_comp_id: <QueryParamInfo>{ name: 'auth_comp_id', type: QueryParamType.String, exampleValue: 'REA', description: 'Author residue name. Corresponds to the \'_atom_site.auth_comp_id\' field.' },
    pdbx_PDB_ins_code: <QueryParamInfo>{ name: 'pdbx_PDB_ins_code', type: QueryParamType.String, description: 'Corresponds to the \'_atom_site.pdbx_PDB_ins_code\' field.' },
};

// function entityTest(params: any): Element.Predicate | undefined {
//     if (typeof params.entity_id === 'undefined') return void 0;
//     const p = Props.entity.id, id = '' + params.entityId;
//     return Element.property(l => p(l) === id);
// }

function entityTest1_555(params: any): QueryPredicate | undefined {
    if (typeof params.entity_id === 'undefined') return ctx => ctx.element.unit.conformation.operator.isIdentity;
    const p = Props.entity.id, id = '' + params.entityId;
    return ctx => ctx.element.unit.conformation.operator.isIdentity && p(ctx.element) === id;
}

function chainTest(params: any): QueryPredicate | undefined {
    if (typeof params.label_asym_id !== 'undefined') {
        const p = Props.chain.label_asym_id, id = '' + params.label_asym_id;
        return ctx => p(ctx.element) === id;
    }
    if (typeof params.auth_asym_id !== 'undefined') {
        const p = Props.chain.auth_asym_id, id = '' + params.auth_asym_id;
        return ctx => p(ctx.element) === id;
    }
    return void 0;
}

function residueTest(params: any): QueryPredicate | undefined {
    const props: StructureElement.Property<any>[] = [], values: any[] = [];

    if (typeof params.label_seq_id !== 'undefined') {
        props.push(Props.residue.label_seq_id);
        values.push(+params.label_seq_id);
    }

    if (typeof params.auth_seq_id !== 'undefined') {
        props.push(Props.residue.auth_seq_id);
        values.push(+params.auth_seq_id);
    }

    if (typeof params.label_comp_id !== 'undefined') {
        props.push(Props.residue.label_comp_id);
        values.push(params.label_comp_id);
    }

    if (typeof params.auth_comp_id !== 'undefined') {
        props.push(Props.residue.auth_comp_id);
        values.push(params.auth_comp_id);
    }

    if (typeof params.pdbx_PDB_ins_code !== 'undefined') {
        props.push(Props.residue.pdbx_PDB_ins_code);
        values.push(params.pdbx_PDB_ins_code);
    }

    switch (props.length) {
        case 0: return void 0;
        case 1: return ctx => props[0](ctx.element) === values[0];
        case 2: return ctx => props[0](ctx.element) === values[0] && props[1](ctx.element) === values[1];
        case 3: return ctx => props[0](ctx.element) === values[0] && props[1](ctx.element) === values[1] && props[2](ctx.element) === values[2];
        default: {
            const len = props.length;
            return ctx => {
                for (let i = 0; i < len; i++) if (!props[i](ctx.element) !== values[i]) return false;
                return true;
            };
        }
    }
}

// function buildResiduesQuery(params: any): Query.Provider {
//     return Queries.generators.atoms({ entityTest: entityTest(params), chainTest: chainTest(params), residueTest: residueTest(params) });
// }

const QueryMap: { [id: string]: Partial<QueryDefinition> } = {
    'full': { niceName: 'Full Structure', query: () => Queries.generators.all, description: 'The full structure.' },
    'residueInteraction': {
        niceName: 'Residues Inside a Sphere',
        description: 'Identifies all residues within the given radius from the source residue.',
        query(p) {
            const center = Queries.generators.atoms({ entityTest: entityTest1_555(p), chainTest: chainTest(p), residueTest: residueTest(p) });
            return Queries.modifiers.includeSurroundings(center, { radius: p.radius, wholeResidues: true });
        },
        structureTransform(p, s) {
            return StructureSymmetry.builderSymmetryMates(s, p.radius).run();
        },
        params: [
            AtomSiteParameters.entity_id,
            AtomSiteParameters.label_asym_id,
            AtomSiteParameters.auth_asym_id,
            AtomSiteParameters.label_comp_id,
            AtomSiteParameters.auth_comp_id,
            AtomSiteParameters.pdbx_PDB_ins_code,
            AtomSiteParameters.label_seq_id,
            AtomSiteParameters.auth_seq_id,
            {
                name: 'radius',
                type: QueryParamType.Float,
                defaultValue: 5,
                exampleValue: '5',
                description: 'Value in Angstroms.',
                validation(v: any) {
                    if (v < 1 || v > 10) {
                        throw `Invalid radius for residue interaction query (must be a value between 1 and 10).`;
                    }
                }
            },
        ]
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
                case QueryParamType.String: ret[key] = value; break;
                case QueryParamType.Integer: ret[key] = parseInt(value); break;
                case QueryParamType.Float: ret[key] = parseFloat(value); break;
            }

            if (p.validation) p.validation(ret[key]);
        }
    }

    return ret;
}

export function normalizeQueryParams(query: QueryDefinition, params: any) {
    return _normalizeQueryParams(params, query.params);
}