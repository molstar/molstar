/**
 * Copyright (c) 2018-2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Type from '../type';
import * as Core from './core';
import { Arguments, Argument } from '../symbol';
import { symbol } from '../helpers';

export namespace Types {
    export const ElementSymbol = Type.Value('Structure', 'ElementSymbol');
    export const AtomName = Type.Value('Structure', 'AtomName');

    export const BondFlag = Type.OneOf('Structure', 'BondFlag', Type.Str, ['covalent', 'metallic', 'ion', 'hydrogen', 'sulfide', 'computed', 'aromatic']);
    export const BondFlags = Core.Types.Flags(BondFlag, 'BondFlags');

    export const SecondaryStructureFlag = Type.OneOf('Structure', 'SecondaryStructureFlag', Type.Str, ['alpha', 'beta', '3-10', 'pi', 'sheet', 'strand', 'helix', 'turn', 'none']);
    export const SecondaryStructureFlags = Core.Types.Flags(SecondaryStructureFlag, 'SecondaryStructureFlag');

    export const RingFingerprint = Type.Value('Structure', 'RingFingerprint');
    export const EntityType = Type.OneOf('Structure', 'EntityType', Type.Str, ['polymer', 'non-polymer', 'water', 'branched']);
    export const EntitySubtype = Type.OneOf('Structure', 'EntitySubtype', Type.Str, ['other', 'polypeptide(D)', 'polypeptide(L)', 'polydeoxyribonucleotide', 'polyribonucleotide', 'polydeoxyribonucleotide/polyribonucleotide hybrid', 'cyclic-pseudo-peptide', 'peptide nucleic acid', 'oligosaccharide']);
    export const ObjectPrimitive = Type.OneOf('Structure', 'ObjectPrimitive', Type.Str, ['atomistic', 'sphere', 'gaussian', 'other']);
    export const ResidueId = Type.Value('Structure', 'ResidueId');

    export const ElementSet = Type.Value('Structure', 'ElementSet');
    export const ElementSelection = Type.Value('Structure', 'ElementSelection');
    export const ElementReference = Type.Value('Structure', 'ElementReference');

    export const ElementSelectionQuery = Core.Types.Fn(ElementSelection, 'ElementSelectionQuery');
}

const type = {
    '@header': 'Types',
    elementSymbol: symbol(
        Arguments.Dictionary({ 0: Argument(Type.Str) }),
        Types.ElementSymbol, 'Create element symbol representation from a string value.'),

    atomName: symbol(
        Arguments.Dictionary({ 0: Argument(Type.AnyValue) }), Types.AtomName, 'Convert a value to an atom name.'),

    entityType: symbol(
        Arguments.Dictionary({ 0: Argument(Types.EntityType) }),
        Types.EntityType,
        `Create normalized representation of entity type: ${Type.oneOfValues(Types.EntityType).join(', ')}.`),

    bondFlags: symbol(
        Arguments.List(Types.BondFlag),
        Types.BondFlags,
        `Create bond flags representation from a list of strings. Allowed flags: ${Type.oneOfValues(Types.BondFlag).join(', ')}.`),

    ringFingerprint: symbol(
        Arguments.List(Types.ElementSymbol, { nonEmpty: true }),
        Types.RingFingerprint,
        'Create ring fingerprint from the supplied atom element list.'),

    secondaryStructureFlags: symbol(
        Arguments.List(Types.SecondaryStructureFlag),
        Types.SecondaryStructureFlags,
        `Create secondary structure flags representation from a list of strings. Allowed flags: ${Type.oneOfValues(Types.SecondaryStructureFlag).join(', ')}.`),

    authResidueId: symbol(Arguments.Dictionary({
        0: Argument(Type.Str, { description: 'auth_asym_id' }),
        1: Argument(Type.Num, { description: 'auth_seq_id' }),
        2: Argument(Type.Str, { description: 'pdbx_PDB_ins_code', isOptional: true })
    }), Types.ResidueId, `Residue identifier based on "auth_" annotation.`),
    labelResidueId: symbol(Arguments.Dictionary({
        0: Argument(Type.Str, { description: 'label_entity_id' }),
        1: Argument(Type.Str, { description: 'label_asym_id' }),
        2: Argument(Type.Num, { description: 'label_seq_id' }),
        3: Argument(Type.Str, { description: 'pdbx_PDB_ins_code', isOptional: true })
    }), Types.ResidueId, `Residue identifier based on mmCIF's "label_" annotation.`)
};

const slot = {
    '@header': 'Iteration Slots',
    element: symbol(Arguments.None, Types.ElementReference, 'A reference to the current element.'),
    elementSetReduce: symbol(Arguments.None, Type.Variable('a', Type.AnyValue, true), 'Current value of the element set reducer.')
};

const generator = {
    '@header': 'Generators',
    all: symbol(Arguments.None, Types.ElementSelectionQuery, 'The entire structure.'),

    atomGroups: symbol(Arguments.Dictionary({
        'entity-test': Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test for the 1st atom of every entity' }),
        'chain-test': Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test for the 1st atom of every chain'  }),
        'residue-test': Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'Test for the 1st atom every residue'  }),
        'atom-test': Argument(Type.Bool, { isOptional: true, defaultValue: true }),
        'group-by': Argument(Type.Any, { isOptional: true, defaultValue: `atom-key`, description: 'Group atoms to sets based on this property. Default: each atom has its own set' }),
    }), Types.ElementSelectionQuery, 'Return all atoms for which the tests are satisfied, grouped into sets.'),

    bondedAtomicPairs: symbol(Arguments.Dictionary({
        0: Argument(Type.Bool, { isOptional: true, defaultValue: 'true for covalent bonds' as any, description: 'Test each bond with this predicate. Each bond is visited twice with swapped atom order.' }),
        // TODO: shoud we support this or just use queryEach to get similar behavior
        // 'group-by': Argument(Type.Any, { isOptional: true, defaultValue: ``, description: 'Group the bonds using the privided value' }),
    }), Types.ElementSelectionQuery, 'Return all pairs of atoms for which the test is satisfied.'),

    rings: symbol(Arguments.Dictionary({
        'fingerprint': Argument(Types.RingFingerprint, { isOptional: true }),
        'only-aromatic': Argument(Type.Bool, { isOptional: true, defaultValue: false }),
    }), Types.ElementSelectionQuery, 'Return all rings or those with the specified fingerprint and/or only aromatic rings.'),

    queryInSelection: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        query: Argument(Types.ElementSelectionQuery),
        'in-complement': Argument(Type.Bool, { isOptional: true, defaultValue: false })
    }), Types.ElementSelectionQuery, 'Executes query only on atoms that are in the source selection.'),

    empty: symbol(Arguments.None, Types.ElementSelectionQuery, 'Nada.'),
};

const modifier = {
    '@header': 'Selection Modifications',

    queryEach: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        query: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'Query every atom set in the input selection separately.'),

    intersectBy: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        by: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'Intersect each atom set from the first sequence from atoms in the second one.'),

    exceptBy: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        by: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, `Remove all atoms from 'selection' that occur in 'by'.`),

    unionBy: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        by: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'For each atom set A in the orginal sequence, combine all atoms sets in the target selection that intersect with A.'),

    union: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'Collects all atom sets in the sequence into a single atom set.'),

    cluster: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        'min-distance': Argument(Type.Num, { isOptional: true, defaultValue: 0 }),
        'max-distance': Argument(Type.Num),
        'min-size': Argument(Type.Num, { description: 'Minimal number of sets to merge, must be at least 2', isOptional: true, defaultValue: 2 }),
        'max-size': Argument(Type.Num, { description: 'Maximal number of sets to merge, if not set, no limit', isOptional: true }),
    }), Types.ElementSelectionQuery, 'Combines atom sets that have mutual distance in the interval [min-radius, max-radius]. Minimum/maximum size determines how many atom sets can be combined.'),

    includeSurroundings: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        radius: Argument(Type.Num),
        'atom-radius': Argument(Type.Num, { isOptional: true, defaultValue: 0, description: 'Value added to each atom before the distance check, for example VDW radius. Using this argument is computationally demanding.' }),
        'as-whole-residues': Argument(Type.Bool, { isOptional: true })
    }), Types.ElementSelectionQuery, 'For each atom set in the selection, include all surrouding atoms/residues that are within the specified radius.'),

    includeConnected: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        'bond-test': Argument(Type.Bool, { isOptional: true, defaultValue: 'true for covalent bonds' as any }),
        'layer-count': Argument(Type.Num, { isOptional: true, defaultValue: 1, description: 'Number of bonded layers to include.' }),
        'as-whole-residues': Argument(Type.Bool, { isOptional: true })
    }), Types.ElementSelectionQuery, 'Pick all atom sets that are connected to the target.'),

    wholeResidues: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
    }), Types.ElementSelectionQuery, 'Expand the selection to whole residues.'),

    expandProperty: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        property: Argument(Type.AnyValue)
    }), Types.ElementSelectionQuery, 'To each atom set in the selection, add all atoms that have the same property value that was already present in the set.')
};

const filter = {
    '@header': 'Selection Filters',
    pick: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        test: Argument(Type.Bool)
    }), Types.ElementSelectionQuery, 'Pick all atom sets that satisfy the test.'),

    first: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'Take the 1st atom set in the sequence.'),

    withSameAtomProperties: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        source: Argument(Types.ElementSelectionQuery),
        property: Argument(Type.Any)
    }), Types.ElementSelectionQuery, 'Pick all atom sets for which the set of given atom properties is a subset of the source properties.'),

    intersectedBy: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        by: Argument(Types.ElementSelectionQuery)
    }), Types.ElementSelectionQuery, 'Pick all atom sets that have non-zero intersection with the target.'),

    within: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        target: Argument(Types.ElementSelectionQuery),
        'min-radius': Argument(Type.Num, { isOptional: true, defaultValue: 0 }),
        'max-radius': Argument(Type.Num),
        'atom-radius': Argument(Type.Num, { isOptional: true, defaultValue: 0, description: 'Value added to each atom before the distance check, for example VDW radius. Using this argument is computationally demanding.' }),
        invert: Argument(Type.Bool, { isOptional: true, defaultValue: false, description: 'If true, pick only atom sets that are further than the specified radius.' }),
    }), Types.ElementSelectionQuery, 'Pick all atom sets from selection that have any atom within the radius of any atom from target.'),

    isConnectedTo: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery),
        target: Argument(Types.ElementSelectionQuery),
        'bond-test': Argument(Type.Bool, { isOptional: true, defaultValue: 'true for covalent bonds' as any }),
        disjunct: Argument(Type.Bool, { isOptional: true, defaultValue: true, description: 'If true, there must exist a bond to an atom that lies outside the given atom set to pass test.' }),
        invert: Argument(Type.Bool, { isOptional: true, defaultValue: false, description: 'If true, return atom sets that are not connected.' })
    }), Types.ElementSelectionQuery, 'Pick all atom sets that are connected to the target.'),
};

const combinator = {
    '@header': 'Selection Combinators',
    intersect: symbol(Arguments.List(Types.ElementSelectionQuery), Types.ElementSelectionQuery, 'Return all unique atom sets that appear in all of the source selections.'),
    merge: symbol(Arguments.List(Types.ElementSelectionQuery), Types.ElementSelectionQuery, 'Merges multiple selections into a single one. Only unique atom sets are kept.'),
    distanceCluster: symbol(Arguments.Dictionary({
        matrix: Argument(Core.Types.List(Core.Types.List(Type.Num)), { description: 'Distance matrix, represented as list of rows (num[][])). Lower triangle is min distance, upper triangle is max distance.' }),
        selections: Argument(Core.Types.List(Types.ElementSelectionQuery), { description: 'A list of held selections.' })
    }), Types.ElementSelectionQuery, 'Pick combinations of atom sets from the source sequences that are mutually within distances specified by a matrix.')
};

const atomSet = {
    '@header': 'Atom Sets',

    atomCount: symbol(Arguments.None, Type.Num),

    countQuery: symbol(Arguments.Dictionary({
        0: Argument(Types.ElementSelectionQuery)
    }), Type.Num, 'Counts the number of occurences of a specific query inside the current atom set.'),

    reduce: symbol(Arguments.Dictionary({
        initial: Argument(Type.Variable('a', Type.AnyValue, true), { description: 'Initial value assigned to slot.atom-set-reduce. Current atom is set to the 1st atom of the current set for this.' }),
        value: Argument(Type.Variable('a', Type.AnyValue, true), { description: 'Expression executed for each atom in the set' })
    }), Type.Variable('a', Type.AnyValue, true), 'Execute the value expression for each atom in the current atom set and return the result. Works the same way as Array.reduce in JavaScript (``result = value(value(...value(initial)))``)'),

    propertySet: symbol(Arguments.Dictionary({
        0: Argument(Core.Types.ConstrainedVar),
    }), Core.Types.Set(Core.Types.ConstrainedVar), 'Returns a set with all values of the given property in the current atom set.'),
};

const atomProperty = {
    '@header': 'Atom Properties',

    core: {
        '@header': 'Core Properties',

        elementSymbol: atomProp(Types.ElementSymbol),

        vdw: atomProp(Type.Num, 'Van der Waals radius'),
        mass: atomProp(Type.Num, 'Atomic weight'),
        atomicNumber: atomProp(Type.Num, 'Atomic number'),

        x: atomProp(Type.Num, 'Cartesian X coordinate'),
        y: atomProp(Type.Num, 'Cartesian Y coordinate'),
        z: atomProp(Type.Num, 'Cartesian Z coordinate'),

        atomKey: atomProp(Type.AnyValue, 'Unique value for each atom. Main use case is grouping of atoms.'),

        bondCount: symbol(Arguments.Dictionary({
            0: Argument(Types.ElementReference, { isOptional: true, defaultValue: 'slot.current-atom' }),
            flags: Argument(Types.BondFlags, { isOptional: true, defaultValue: 'covalent' as any }),
        }), Type.Num, 'Number of bonds (by default only covalent bonds are counted).'),

        sourceIndex: atomProp(Type.Num, 'Index of the atom/element in the input file.'),
        operatorName: atomProp(Type.Str, 'Name of the symmetry operator applied to this element.'),
        modelIndex: atomProp(Type.Num, 'Index of the model in the input file.'),
        modelLabel: atomProp(Type.Str, 'Label/header of the model in the input file.')
    },

    topology: {
        connectedComponentKey: atomProp(Type.AnyValue, 'Unique value for each connected component.')
    },

    macromolecular: {
        '@header': 'Macromolecular Properties (derived from the mmCIF format)',

        authResidueId: atomProp(Types.ResidueId, `type.auth-residue-id symbol executed on current atom's residue`),
        labelResidueId: atomProp(Types.ResidueId, `type.label-residue-id symbol executed on current atom's residue`),

        residueKey: atomProp(Type.AnyValue, 'Unique value for each tuple ``(label_entity_id,auth_asym_id, auth_seq_id, pdbx_PDB_ins_code)``, main use case is grouping of atoms'),
        chainKey: atomProp(Type.AnyValue, 'Unique value for each tuple ``(label_entity_id, auth_asym_id)``, main use case is grouping of atoms'),
        entityKey: atomProp(Type.AnyValue, 'Unique value for each tuple ``label_entity_id``, main use case is grouping of atoms'),

        isHet: atomProp(Type.Bool, 'Equivalent to atom_site.group_PDB !== ATOM'),

        id: atomProp(Type.Num, '_atom_site.id'),

        label_atom_id: atomProp(Types.AtomName),
        label_alt_id: atomProp(Type.Str),
        label_comp_id: atomProp(Type.Str),
        label_asym_id: atomProp(Type.Str),
        label_entity_id: atomProp(Type.Str),
        label_seq_id: atomProp(Type.Num),

        auth_atom_id: atomProp(Types.AtomName),
        auth_comp_id: atomProp(Type.Str),
        auth_asym_id: atomProp(Type.Str),
        auth_seq_id: atomProp(Type.Num),

        pdbx_PDB_ins_code: atomProp(Type.Str),
        pdbx_formal_charge: atomProp(Type.Num),

        occupancy: atomProp(Type.Num),
        B_iso_or_equiv: atomProp(Type.Num),

        entityType: atomProp(Types.EntityType, 'Type of the entity as defined in mmCIF (polymer, non-polymer, branched, water)'),
        entitySubtype: atomProp(Types.EntitySubtype, 'Subtype of the entity as defined in mmCIF _entity_poly.type and _pdbx_entity_branch.type (other, polypeptide(D), polypeptide(L), polydeoxyribonucleotide, polyribonucleotide, polydeoxyribonucleotide/polyribonucleotide hybrid, cyclic-pseudo-peptide, peptide nucleic acid, oligosaccharide)'),
        entityPrdId: atomProp(Type.Str, `The PRD ID of the entity.`),
        entityDescription: atomProp(Core.Types.List(Type.Str)),
        objectPrimitive: atomProp(Types.ObjectPrimitive, 'Type of the primitive object used to model this segment as defined in mmCIF/IHM (atomistic, sphere, gaussian, other)'),

        secondaryStructureKey: atomProp(Type.AnyValue, 'Unique value for each secondary structure element.'),
        secondaryStructureFlags: atomProp(Types.SecondaryStructureFlags),
        isModified: atomProp(Type.Bool, 'True if the atom belongs to modification of a standard residue.'),
        modifiedParentName: atomProp(Type.Str, `'3-letter' code of the modifed parent residue.`),
        isNonStandard: atomProp(Type.Bool, 'True if this is a non-standard residue.'),
        chemCompType: atomProp(Type.Str, `Type of the chemical component as defined in mmCIF.`),
    }
};

const bondProperty = {
    '@header': 'Bond Properties',

    flags: bondProp(Types.BondFlags),
    order: bondProp(Type.Num),
    length: bondProp(Type.Num),
    atomA: bondProp(Types.ElementReference),
    atomB: bondProp(Types.ElementReference)
};

function atomProp(type: Type, description?: string) {
    return symbol(Arguments.Dictionary({ 0: Argument(Types.ElementReference, { isOptional: true, defaultValue: 'slot.current-atom' }) }), type, description);
}

function bondProp(type: Type, description?: string) {
    return symbol(Arguments.None, type, description);
}

export default {
    '@header': 'Structure Queries',
    type,
    slot,
    generator,
    modifier,
    filter,
    combinator,
    atomSet,
    atomProperty,
    bondProperty: bondProperty
};