/**
 * Copyright (c) 2018 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Symbol, { Arguments, Argument } from '../../symbol'
import B from '../../builder'
import * as M from './macro'
import MolScript from '../../symbol-table'
import Type from '../../type'
import * as Struct from '../../symbol-table/structure-query'
import Expression from '../../expression'
import { UniqueArray } from 'mol-data/generic'

export type MolScriptSymbol =
    | { kind: 'alias', aliases: string[], symbol: Symbol }
    | { kind: 'macro', aliases: string[], symbol: Symbol, translate: (args: any) => Expression }

function Alias(symbol: Symbol<any>, ...aliases: string[]): MolScriptSymbol { return { kind: 'alias', aliases, symbol }; }
function Macro(symbol: Symbol<any>, translate: (args: any) => Expression, ...aliases: string[]): MolScriptSymbol {
    symbol.info.namespace = 'molscript-macro';
    symbol.id = `molscript-macro.${symbol.info.name}`;
    return { kind: 'macro', symbol, translate, aliases: [symbol.info.name, ...aliases] };
}

export function isMolScriptSymbol(x: any): x is MolScriptSymbol {
    return x.kind === 'alias' || x.kind === 'macro';
}

export const SymbolTable = [
    [
        'Core symbols',
        Alias(MolScript.core.type.bool, 'bool'),
        Alias(MolScript.core.type.num, 'num'),
        Alias(MolScript.core.type.str, 'str'),
        Alias(MolScript.core.type.regex, 'regex'),
        Alias(MolScript.core.type.list, 'list'),
        Alias(MolScript.core.type.set, 'set'),

        Alias(MolScript.core.type.compositeKey, 'composite-key'),
        Alias(MolScript.core.logic.not, 'not'),
        Alias(MolScript.core.logic.and, 'and'),
        Alias(MolScript.core.logic.or, 'or'),
        Alias(MolScript.core.ctrl.if, 'if'),
        Alias(MolScript.core.ctrl.fn, 'fn'),
        Alias(MolScript.core.ctrl.eval, 'eval'),
        Alias(MolScript.core.math.add, 'add', '+'),
        Alias(MolScript.core.math.sub, 'sub', '-'),
        Alias(MolScript.core.math.mult, 'mult', '*'),
        Alias(MolScript.core.math.div, 'div', '/'),
        Alias(MolScript.core.math.pow, 'pow', '**'),
        Alias(MolScript.core.math.mod, 'mod'),
        Alias(MolScript.core.math.min, 'min'),
        Alias(MolScript.core.math.max, 'max'),
        Alias(MolScript.core.math.floor, 'floor'),
        Alias(MolScript.core.math.ceil, 'ceil'),
        Alias(MolScript.core.math.roundInt, 'round'),
        Alias(MolScript.core.math.abs, 'abs'),
        Alias(MolScript.core.math.sqrt, 'sqrt'),
        Alias(MolScript.core.math.sin, 'sin'),
        Alias(MolScript.core.math.cos, 'cos'),
        Alias(MolScript.core.math.tan, 'tan'),
        Alias(MolScript.core.math.asin, 'asin'),
        Alias(MolScript.core.math.acos, 'acos'),
        Alias(MolScript.core.math.atan, 'atan'),
        Alias(MolScript.core.math.sinh, 'sinh'),
        Alias(MolScript.core.math.cosh, 'cosh'),
        Alias(MolScript.core.math.tanh, 'tanh'),
        Alias(MolScript.core.math.exp, 'exp'),
        Alias(MolScript.core.math.log, 'log'),
        Alias(MolScript.core.math.log10, 'log10'),
        Alias(MolScript.core.math.atan2, 'atan2'),
        Alias(MolScript.core.rel.eq, 'eq', '='),
        Alias(MolScript.core.rel.neq, 'neq', '!='),
        Alias(MolScript.core.rel.lt, 'lt', '<'),
        Alias(MolScript.core.rel.lte, 'lte', '<='),
        Alias(MolScript.core.rel.gr, 'gr', '>'),
        Alias(MolScript.core.rel.gre, 'gre', '>='),
        Alias(MolScript.core.rel.inRange, 'in-range'),
        Alias(MolScript.core.str.concat, 'concat'),
        Alias(MolScript.core.str.match, 'regex.match'),
        Alias(MolScript.core.list.getAt, 'list.get'),
        Alias(MolScript.core.set.has, 'set.has'),
        Alias(MolScript.core.set.isSubset, 'set.subset'),
    ],
    [
        'Structure',
        [
            'Types',
            Alias(MolScript.structureQuery.type.entityType, 'ent-type'),
            Alias(MolScript.structureQuery.type.authResidueId, 'auth-resid'),
            Alias(MolScript.structureQuery.type.labelResidueId, 'label-resid'),
            Alias(MolScript.structureQuery.type.ringFingerprint, 'ringfp'),
            Alias(MolScript.structureQuery.type.bondFlags, 'bond-flags'),
        ],
        [
            'Slots',
            Alias(MolScript.structureQuery.slot.elementSetReduce, 'atom.set.reduce.value'),
        ],
        [
            'Generators',
            Alias(MolScript.structureQuery.generator.atomGroups, 'sel.atom.atom-groups'),
            Alias(MolScript.structureQuery.generator.queryInSelection, 'sel.atom.query-in-selection'),
            Alias(MolScript.structureQuery.generator.rings, 'sel.atom.rings'),
            Alias(MolScript.structureQuery.generator.empty, 'sel.atom.empty'),

            Macro(Symbol('sel.atom.atoms', Arguments.Dictionary({
                0: Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test applied to each atom.' })
            }), Struct.Types.ElementSelection, 'A selection of singleton atom sets.'),
            args => B.struct.generator.atomGroups({ 'atom-test':  M.tryGetArg(args, 0, true) })),

            Macro(Symbol('sel.atom.res', Arguments.Dictionary({
                0: Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test applied to the 1st atom of each residue.' })
            }), Struct.Types.ElementSelection, 'A selection of atom sets grouped by residue.'),
            args => B.struct.generator.atomGroups({
                'residue-test':  M.tryGetArg(args, 0, true),
                'group-by': B.ammp('residueKey')
            })),

            Macro(Symbol('sel.atom.chains', Arguments.Dictionary({
                0: Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test applied to the 1st atom of each chain.' })
            }), Struct.Types.ElementSelection, 'A selection of atom sets grouped by chain.'),
            args => B.struct.generator.atomGroups({
                'chain-test': M.tryGetArg(args, 0, true),
                'group-by': B.ammp('chainKey')
            })),
        ],
        [
            'Modifiers',
            Alias(MolScript.structureQuery.modifier.queryEach, 'sel.atom.query-each'),
            Alias(MolScript.structureQuery.modifier.intersectBy, 'sel.atom.intersect-by'),
            Alias(MolScript.structureQuery.modifier.exceptBy, 'sel.atom.except-by'),
            Alias(MolScript.structureQuery.modifier.unionBy, 'sel.atom.union-by'),
            Alias(MolScript.structureQuery.modifier.union, 'sel.atom.union'),
            Alias(MolScript.structureQuery.modifier.cluster, 'sel.atom.cluster'),
            Alias(MolScript.structureQuery.modifier.includeSurroundings, 'sel.atom.include-surroundings'),
            Alias(MolScript.structureQuery.modifier.includeConnected, 'sel.atom.include-connected'),
            Alias(MolScript.structureQuery.modifier.expandProperty, 'sel.atom.expand-property'),

            Macro(Symbol('sel.atom.around', Arguments.Dictionary({
                0: Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test applied to the 1st atom of each chain.' })
            }), Struct.Types.ElementSelection, 'A selection of singleton atom sets with centers within "radius" of the center of any atom in the given selection.'),
            args => B.struct.modifier.exceptBy({
                '0': B.struct.filter.within({
                    '0': B.struct.generator.atomGroups(), target: M.tryGetArg(args, 0), 'max-radius': M.tryGetArg(args, 'radius')
                }),
                by: M.tryGetArg(args, 0)
            }))
        ],
        [
            'Filters',
            Alias(MolScript.structureQuery.filter.pick, 'sel.atom.pick'),
            Alias(MolScript.structureQuery.filter.withSameAtomProperties, 'sel.atom.with-same-atom-properties'),
            Alias(MolScript.structureQuery.filter.intersectedBy, 'sel.atom.intersected-by'),
            Alias(MolScript.structureQuery.filter.within, 'sel.atom.within'),
            Alias(MolScript.structureQuery.filter.isConnectedTo, 'sel.atom.is-connected-to'),
        ],
        [
            'Combinators',
            Alias(MolScript.structureQuery.combinator.intersect, 'sel.atom.intersect'),
            Alias(MolScript.structureQuery.combinator.merge, 'sel.atom.merge'),
            Alias(MolScript.structureQuery.combinator.distanceCluster, 'sel.atom.dist-cluster'),
        ],
        [
            'Atom Set Properties',
            Alias(MolScript.structureQuery.atomSet.atomCount, 'atom.set.atom-count'),
            Alias(MolScript.structureQuery.atomSet.countQuery, 'atom.set.count-query'),
            Alias(MolScript.structureQuery.atomSet.reduce, 'atom.set.reduce'),
            Alias(MolScript.structureQuery.atomSet.propertySet, 'atom.set.property'),

            Macro(Symbol('atom.set.max', Arguments.Dictionary({
                0: Argument(Type.Num, { description: 'Numeric atom property.'})
            }), Type.Num, 'Maximum of the given property in the current atom set.'),
            args => M.aggregate(M.tryGetArg(args, 0), B.core.math.max)),

            Macro(Symbol('atom.set.sum', Arguments.Dictionary({
                0: Argument(Type.Num, { description: 'Numeric atom property.'})
            }), Type.Num, 'Sum of the given property in the current atom set.'),
            args => M.aggregate(M.tryGetArg(args, 0), B.core.math.add, 0)),

            Macro(Symbol('atom.set.avg', Arguments.Dictionary({
                0: Argument(Type.Num, { description: 'Numeric atom property.'})
            }), Type.Num, 'Average of the given property in the current atom set.'),
            args => B.core.math.div([ M.aggregate(M.tryGetArg(args, 0), B.core.math.add, 0), B.struct.atomSet.atomCount() ])),

            Macro(Symbol('atom.set.min', Arguments.Dictionary({
                0: Argument(Type.Num, { description: 'Numeric atom property.'})
            }), Type.Num, 'Minimum of the given property in the current atom set.'),
            args => M.aggregate(M.tryGetArg(args, 0), B.core.math.min))
        ],
        [
            'Atom Properties',
            Alias(MolScript.structureQuery.atomProperty.core.elementSymbol, 'atom.el'),
            Alias(MolScript.structureQuery.atomProperty.core.vdw, 'atom.vdw'),
            Alias(MolScript.structureQuery.atomProperty.core.mass, 'atom.mass'),
            Alias(MolScript.structureQuery.atomProperty.core.atomicNumber, 'atom.atomic-number'),
            Alias(MolScript.structureQuery.atomProperty.core.x, 'atom.x'),
            Alias(MolScript.structureQuery.atomProperty.core.y, 'atom.y'),
            Alias(MolScript.structureQuery.atomProperty.core.z, 'atom.z'),
            Alias(MolScript.structureQuery.atomProperty.core.atomKey, 'atom.key'),
            Alias(MolScript.structureQuery.atomProperty.core.bondCount, 'atom.bond-count'),

            Alias(MolScript.structureQuery.atomProperty.topology.connectedComponentKey, 'atom.key.molecule'),

            Alias(MolScript.structureQuery.atomProperty.macromolecular.authResidueId, 'atom.auth-resid'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.labelResidueId, 'atom.label-resid'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.residueKey, 'atom.key.res'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.chainKey, 'atom.key.chain'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.entityKey, 'atom.key.entity'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.isHet, 'atom.is-het'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.id, 'atom.id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_atom_id, 'atom.label_atom_id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_alt_id, 'atom.label_alt_id', 'atom.altloc'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_comp_id, 'atom.label_comp_id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_asym_id, 'atom.label_asym_id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_entity_id, 'atom.label_entity_id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.label_seq_id, 'atom.label_seq_id'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.auth_atom_id, 'atom.auth_atom_id', 'atom.name'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.auth_comp_id, 'atom.auth_comp_id', 'atom.resname'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.auth_asym_id, 'atom.auth_asym_id', 'atom.chain'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.auth_seq_id, 'atom.auth_seq_id', 'atom.resno'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.pdbx_PDB_ins_code, 'atom.pdbx_PDB_ins_code', 'atom.inscode'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.pdbx_formal_charge, 'atom.pdbx_formal_charge'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.occupancy, 'atom.occupancy'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.B_iso_or_equiv, 'atom.B_iso_or_equiv', 'atom.bfactor'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.entityType, 'atom.entity-type'),

            Alias(MolScript.structureQuery.atomProperty.macromolecular.secondaryStructureKey, 'atom.key.sec-struct'),

            Alias(MolScript.structureQuery.atomProperty.macromolecular.isModified, 'atom.is-modified'),
            Alias(MolScript.structureQuery.atomProperty.macromolecular.modifiedParentName, 'atom.modified-parent'),

            Macro(Symbol('atom.sec-struct.is', Arguments.List(Struct.Types.SecondaryStructureFlag), Type.Bool,
                `Test if the current atom is part of an secondary structure. Optionally specify allowed sec. struct. types: ${Type.oneOfValues(Struct.Types.SecondaryStructureFlag).join(', ')}`),
            args => B.core.flags.hasAny([B.struct.atomProperty.macromolecular.secondaryStructureFlags(), B.struct.type.secondaryStructureFlags(args)])),
        ],
        [
            'Bond Properties',
            Alias(MolScript.structureQuery.bondProperty.order, 'bond.order'),
            Macro(Symbol('bond.is', Arguments.List(Struct.Types.BondFlag), Type.Bool,
                `Test if the current bond has at least one (or all if partial = false) of the specified flags: ${Type.oneOfValues(Struct.Types.BondFlag).join(', ')}`),
            args => B.core.flags.hasAny([B.struct.bondProperty.flags(), B.struct.type.bondFlags(M.getPositionalArgs(args))])),
        ]
    ]
];

const list: MolScriptSymbol[] = [];

function makeList(xs: any[]) {
    for (const x of xs) {
        if (isMolScriptSymbol(x)) list.push(x);
        else if (x instanceof Array) makeList(x);
    }
}

makeList(SymbolTable);

const normalized = (function () {
    const symbolList: [string, MolScriptSymbol][] = [];
    const symbolMap: { [id: string]: MolScriptSymbol | undefined } = Object.create(null);
    const namedArgs = UniqueArray.create<string, string>();
    const constants = UniqueArray.create<string, string>();

    for (const s of list) {
        for (const a of s.aliases) {
            symbolList.push([a, s]);
            if (symbolMap[a]) throw new Error(`Alias '${a}' already in use.`);
            symbolMap[a] = s;
        }
        const args = s.symbol.args;
        if (args.kind !== 'dictionary') {
            if (args.type.kind === 'oneof') {
                Type.oneOfValues(args.type).forEach(v => UniqueArray.add(constants, v, v));
            }
            continue;
        }
        for (const a of Object.keys(args.map)) {
            if (isNaN(a as any)) UniqueArray.add(namedArgs, a, a);
            const arg = ((args.map as any)[a]) as Argument;
            if (arg.type.kind === 'oneof') {
                Type.oneOfValues(arg.type).forEach(v => UniqueArray.add(constants, v, v));
            }
        }
    }

    return { symbolList, symbolMap, namedArgs: namedArgs.array, constants: constants.array }
})();

export const MolScriptSymbols = list;
export const Constants = normalized.constants;
export const NamedArgs = normalized.namedArgs;
export const SymbolMap = normalized.symbolMap;
export const SymbolList = normalized.symbolList;

const sortedSymbols = SymbolList.map(s => s[0]).sort((a, b) => {
    if (a.length === b.length) return (a < b) as any;
    return a.length - b.length;
});
export default [...sortedSymbols, ...NamedArgs.map(a => ':' + a), ...Constants];