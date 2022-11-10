/**
 * Copyright (c) 2017-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 *
 * Adapted from MolQL project
 */

import * as P from '../../../mol-util/monadic-parser';
import * as h from '../helper';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
const B = MolScriptBuilder;
import { OperatorList } from '../types';
import { Expression } from '../../language/expression';

export const operators: OperatorList = [
    {
        '@desc': 'Selects atoms that are not included in s1.',
        '@examples': [
            'NOT resn ALA',
            'not (resi 42 or chain A)',
            '!resi 42 or chain A',
        ],
        name: 'not',
        type: h.prefix,
        rule: P.MonadicParser.alt(
            P.MonadicParser.regexp(/NOT/i).skip(P.MonadicParser.whitespace),
            P.MonadicParser.string('!').skip(P.MonadicParser.optWhitespace)
        ),
        map: (op, selection) => h.invertExpr(selection),
    },
    {
        '@desc': 'Selects atoms included in both s1 and s2.',
        '@examples': ['chain A AND name CA'],
        name: 'and',
        type: h.binaryLeft,
        rule: h.infixOp(/AND|&/i),
        map: (op, selection, by) =>
            B.struct.modifier.intersectBy({ 0: selection, by }),
    },
    {
        '@desc': 'Selects atoms included in either s1 or s2.',
        '@examples': ['chain A OR chain B'],
        name: 'or',
        type: h.binaryLeft,
        rule: h.infixOp(/OR|\|/i),
        map: (op: string, s1: Expression, s2: Expression) => B.struct.combinator.merge([s1, s2]),
    },
    {
        '@desc':
            'Selects atoms in s1 whose identifiers name, resi, resn, chain and segi all match atoms in s2.',
        '@examples': ['chain A IN chain B'],
        name: 'in',
        type: h.binaryLeft,
        rule: h.infixOp(/IN/i),
        map: (op: string, selection: Expression, source: Expression) => {
            return B.struct.filter.withSameAtomProperties({
                0: selection,
                source,
                property: B.core.type.compositeKey([
                    B.ammp('label_atom_id'),
                    B.ammp('label_seq_id'),
                    B.ammp('label_comp_id'),
                    B.ammp('auth_asym_id'),
                    B.ammp('label_asym_id'),
                ]),
            });
        },
    },
    {
        '@desc':
            'Selects atoms in s1 whose identifiers name and resi match atoms in s2.',
        '@examples': ['chain A LIKE chain B'],
        name: 'like',
        type: h.binaryLeft,
        rule: h.infixOp(/LIKE|l\./i),
        map: (op: string, selection: Expression, source: Expression) => {
            return B.struct.filter.withSameAtomProperties({
                0: selection,
                source,
                property: B.core.type.compositeKey([
                    B.ammp('label_atom_id'),
                    B.ammp('label_seq_id'),
                ]),
            });
        },
    },
    {
        '@desc':
            'Selects all atoms whose van der Waals radii are separated from the van der Waals radii of s1 by a minimum of X Angstroms.',
        '@examples': ['solvent GAP 2'],
        name: 'gap',
        type: h.postfix,
        rule: h
            .postfixOp(/GAP\s+([-+]?[0-9]*\.?[0-9]+)/i, 1)
            .map((x: any) => parseFloat(x)),
        map: (distance: number, target: Expression) => {
            return B.struct.filter.within({
                '0': B.struct.generator.all(),
                target,
                'atom-radius': B.acp('vdw'),
                'max-radius': distance,
                invert: true,
            });
        },
    },
    {
        '@desc':
            'Selects atoms with centers within X Angstroms of the center of any atom in s1.',
        '@examples': ['resname LIG AROUND 1'],
        name: 'around',
        abbr: ['a.'],
        type: h.postfix,
        rule: h
            .postfixOp(/(AROUND|a\.)\s+([-+]?[0-9]*\.?[0-9]+)/i, 2)
            .map((x: any) => parseFloat(x)),
        map: (radius: number, target: Expression) => {
            return B.struct.modifier.exceptBy({
                '0': B.struct.filter.within({
                    '0': B.struct.generator.all(),
                    target,
                    'max-radius': radius,
                }),
                by: target,
            });
        },
    },
    {
        '@desc':
            'Expands s1 by all atoms within X Angstroms of the center of any atom in s1.',
        '@examples': ['chain A EXPAND 3'],
        name: 'expand',
        abbr: ['x.'],
        type: h.postfix,
        rule: h
            .postfixOp(/(EXPAND|x\.)\s+([-+]?[0-9]*\.?[0-9]+)/i, 2)
            .map((x: any) => parseFloat(x)),
        map: (radius: number, selection: Expression) => {
            return B.struct.modifier.includeSurroundings({ 0: selection, radius });
        },
    },
    {
        '@desc':
            'Selects atoms in s1 that are within X Angstroms of any atom in s2.',
        '@examples': ['chain A WITHIN 3 OF chain B'],
        name: 'within',
        abbr: ['w.'],
        type: h.binaryLeft,
        rule: h.ofOp('WITHIN', 'w.'),
        map: (radius: number, selection: Expression, target: Expression) => {
            return B.struct.filter.within({
                0: selection,
                target,
                'max-radius': radius,
            });
        },
    },
    {
        '@desc':
            'Same as within, but excludes s2 from the selection (and thus is identical to s1 and s2 around X).',
        '@examples': ['chain A NEAR_TO 3 OF chain B'],
        name: 'near_to',
        abbr: ['nto.'],
        type: h.binaryLeft,
        rule: h.ofOp('NEAR_TO', 'nto.'),
        map: (radius: number, selection: Expression, target: Expression) => {
            return B.struct.modifier.exceptBy({
                '0': B.struct.filter.within({
                    '0': selection,
                    target,
                    'max-radius': radius,
                }),
                by: target,
            });
        },
    },
    {
        '@desc': 'Selects atoms in s1 that are at least X Anstroms away from s2.',
        '@examples': ['solvent BEYOND 2 OF chain A'],
        name: 'beyond',
        abbr: ['be.'],
        type: h.binaryLeft,
        rule: h.ofOp('BEYOND', 'be.'),
        map: (radius: number, selection: Expression, target: Expression) => {
            return B.struct.modifier.exceptBy({
                '0': B.struct.filter.within({
                    '0': selection,
                    target,
                    'max-radius': radius,
                    invert: true,
                }),
                by: target,
            });
        },
    },
    {
        '@desc': 'Expands selection to complete residues.',
        '@examples': ['BYRESIDUE name N'],
        name: 'byresidue',
        abbr: ['byresi', 'byres', 'br.'],
        type: h.prefix,
        rule: h.prefixOp(/BYRESIDUE|byresi|byres|br\./i),
        map: (op: string, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.expandProperty({
                    '0': B.struct.modifier.union({ 0: selection }),
                    property: B.ammp('residueKey'),
                })
            );
        },
    },
    {
        '@desc':
            'Completely selects all alpha carbons in all residues covered by a selection.',
        '@examples': ['BYCALPHA chain A'],
        name: 'bycalpha',
        abbr: ['bca.'],
        type: h.prefix,
        rule: h.prefixOp(/BYCALPHA|bca\./i),
        map: (op: string, selection: Expression) => {
            return B.struct.generator.queryInSelection({
                '0': B.struct.modifier.expandProperty({
                    '0': B.struct.modifier.union({ 0: selection }),
                    property: B.ammp('residueKey'),
                }),
                query: B.struct.generator.atomGroups({
                    'atom-test': B.core.rel.eq([
                        B.atomName('CA'),
                        B.ammp('label_atom_id'),
                    ]),
                }),
            });
        },
    },
    {
        '@desc': 'Expands selection to complete molecules.',
        '@examples': ['BYMOLECULE resi 20-30'],
        name: 'bymolecule',
        isUnsupported: true, // structure-query.atom-property.topology.connected-component-key' is not implemented
        abbr: ['bymol', 'bm.'],
        type: h.prefix,
        rule: h.prefixOp(/BYMOLECULE|bymol|bm\./i),
        map: (op: string, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.expandProperty({
                    '0': B.struct.modifier.union({ 0: selection }),
                    property: B.atp('connectedComponentKey'),
                })
            );
        },
    },
    {
        '@desc': 'Expands selection to complete fragments.',
        '@examples': ['BYFRAGMENT resi 10'],
        name: 'byfragment',
        abbr: ['byfrag', 'bf.'],
        isUnsupported: true,
        type: h.prefix,
        rule: h.prefixOp(/BYFRAGMENT|byfrag|bf\./i),
        map: (op: string, selection: Expression) => [op, selection],
    },
    {
        '@desc': 'Expands selection to complete segments.',
        '@examples': ['BYSEGMENT resn CYS'],
        name: 'bysegment',
        abbr: ['bysegi', 'byseg', 'bs.'],
        type: h.prefix,
        rule: h.prefixOp(/BYSEGMENT|bysegi|byseg|bs\./i),
        map: (op: string, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.expandProperty({
                    '0': B.struct.modifier.union({ 0: selection }),
                    property: B.ammp('chainKey'),
                })
            );
        },
    },
    {
        '@desc': 'Expands selection to complete objects.',
        '@examples': ['BYOBJECT chain A'],
        name: 'byobject',
        abbr: ['byobj', 'bo.'],
        isUnsupported: true,
        type: h.prefix,
        rule: h.prefixOp(/BYOBJECT|byobj|bo\./i),
        map: (op: string, selection: Expression) => [op, selection],
    },
    {
        '@desc': 'Expands selection to unit cell.',
        '@examples': ['BYCELL chain A'],
        name: 'bycell',
        isUnsupported: true,
        type: h.prefix,
        rule: h.prefixOp(/BYCELL/i),
        map: (op: string, selection: Expression) => [op, selection],
    },
    {
        '@desc': 'All rings of size ≤ 7 which have at least one atom in s1.',
        '@examples': ['BYRING resn HEM'],
        name: 'byring',
        // isUnsupported: true, // structure-query.atom-set.atom-count' is not implemented.
        type: h.prefix,
        rule: h.prefixOp(/BYRING/i),
        map: (op: string, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.intersectBy({
                    '0': B.struct.filter.pick({
                        '0': B.struct.generator.rings(),
                        test: B.core.logic.and([
                            B.core.rel.lte([B.struct.atomSet.atomCount(), 7]),
                            B.core.rel.gr([B.struct.atomSet.countQuery([selection]), 1]),
                        ]),
                    }),
                    by: selection,
                })
            );
        },
    },
    {
        '@desc': 'Selects atoms directly bonded to s1, excludes s1.',
        '@examples': ['NEIGHBOR resn CYS'],
        name: 'neighbor',
        type: h.prefix,
        abbr: ['nbr.'],
        rule: h.prefixOp(/NEIGHBOR|nbr\./i),
        map: (op: string, selection: Expression) => {
            return B.struct.modifier.exceptBy({
                '0': h.asAtoms(
                    B.struct.modifier.includeConnected({
                        '0': B.struct.modifier.union({ 0: selection }),
                        'bond-test': true,
                    })
                ),
                by: selection,
            });
        },
    },
    {
        '@desc': 'Selects atoms directly bonded to s1, may include s1.',
        '@examples': ['BOUND_TO name CA'],
        name: 'bound_to',
        abbr: ['bto.'],
        type: h.prefix,
        rule: h.prefixOp(/BOUND_TO|bto\./i),
        map: (op: string, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.includeConnected({
                    '0': B.struct.modifier.union({ 0: selection }),
                })
            );
        },
    },
    {
        '@desc': 'Extends s1 by X bonds connected to atoms in s1.',
        '@examples': ['resname LIG EXTEND 3'],
        name: 'extend',
        abbr: ['xt.'],
        type: h.postfix,
        rule: h.postfixOp(/(EXTEND|xt\.)\s+([0-9]+)/i, 2).map((x: any) => parseInt(x)),
        map: (count: number, selection: Expression) => {
            return h.asAtoms(
                B.struct.modifier.includeConnected({
                    '0': B.struct.modifier.union({ 0: selection }),
                    'bond-test': true,
                    'layer-count': count,
                })
            );
        },
    },
];
