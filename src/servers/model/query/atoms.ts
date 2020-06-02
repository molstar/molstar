/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { QueryPredicate, StructureElement, StructureProperties as Props } from '../../../mol-model/structure';
import { AtomsQueryParams } from '../../../mol-model/structure/query/queries/generators';
import { AtomSiteSchema, AtomSiteSchemaElement } from '../server/api';
import { ElementSymbol } from '../../../mol-model/structure/model/types';

export function getAtomsTests(params: AtomSiteSchema): Partial<AtomsQueryParams>[] {
    if (!params) return [{ }];
    if (Array.isArray(params)) {
        return params.map(p => atomsTest(p));
    } else {
        return [atomsTest(params)];
    }
}

function atomsTest(params: AtomSiteSchemaElement): Partial<AtomsQueryParams> {
    return {
        entityTest: entityTest(params),
        chainTest: chainTest(params),
        residueTest: residueTest(params),
        atomTest: atomTest(params)
    };
}

function entityTest(params: AtomSiteSchemaElement): QueryPredicate | undefined {
    if (!params || typeof params.label_entity_id === 'undefined') return void 0;
    const p = Props.entity.id, id = '' + params.label_entity_id;
    return ctx => p(ctx.element) === id;
}

function chainTest(params: AtomSiteSchemaElement): QueryPredicate | undefined {
    if (!params) return void 0;

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

function residueTest(params: AtomSiteSchemaElement): QueryPredicate | undefined {
    if (!params) return void 0;

    const props: StructureElement.Property<any>[] = [], values: any[] = [];

    if (typeof params.label_seq_id !== 'undefined') {
        props.push(Props.residue.label_seq_id);
        values.push(+params.label_seq_id);
    }

    if (typeof params.auth_seq_id !== 'undefined') {
        props.push(Props.residue.auth_seq_id);
        values.push(+params.auth_seq_id);
    }

    if (typeof params.pdbx_PDB_ins_code !== 'undefined') {
        props.push(Props.residue.pdbx_PDB_ins_code);
        values.push(params.pdbx_PDB_ins_code);
    }

    return andEqual(props, values);
}

function atomTest(params: AtomSiteSchemaElement): QueryPredicate | undefined {
    if (!params) return void 0;

    const props: StructureElement.Property<any>[] = [], values: any[] = [];

    if (typeof params.label_atom_id !== 'undefined') {
        props.push(Props.atom.label_atom_id);
        values.push(params.label_atom_id);
    }

    if (typeof params.auth_atom_id !== 'undefined') {
        props.push(Props.atom.auth_atom_id);
        values.push(params.auth_atom_id);
    }

    if (typeof params.type_symbol !== 'undefined') {
        props.push(Props.atom.type_symbol);
        values.push(ElementSymbol(params.type_symbol));
    }

    if (typeof params.label_comp_id !== 'undefined') {
        props.push(Props.atom.label_comp_id);
        values.push(params.label_comp_id);
    }

    if (typeof params.auth_comp_id !== 'undefined') {
        props.push(Props.atom.auth_comp_id);
        values.push(params.auth_comp_id);
    }

    return andEqual(props, values);
}

function andEqual(props: StructureElement.Property<any>[], values: any[]): QueryPredicate | undefined {
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