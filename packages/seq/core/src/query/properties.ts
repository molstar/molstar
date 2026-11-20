/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Property access (_spec/README.md "Queries", Properties; properties.md): a property is a typed symbol read per
// element of one level (molecule, unit, residue, atom), apart from the language that combines forms. The built-ins
// have a short name (printed) and a long one (MolViewSpec's, mmCIF's: read); an extension adds its own through the
// registry (`query-property`), namespaced. A name resolves to its definition here; how a value is read for an element
// is `structure.ts`'s.

import type { QueryModel, PropertyRead } from './model.js';

export type PropertyLevel = 'molecule' | 'unit' | 'residue' | 'atom';

/** `words` is a set: `=` on it is membership (a position's alternative components). */
export type PropertyType = 'word' | 'integer' | 'number' | 'boolean' | 'words';

/** Coarse to fine: a test reads its level's properties and coarser ones, never finer (`seq.query.property.level`). */
export const PropertyLevels: readonly PropertyLevel[] = ['molecule', 'unit', 'residue', 'atom'];

export const levelRank = (l: PropertyLevel) => PropertyLevels.indexOf(l);

/** A property's definition: its names, level and type, and for an extension's, how it is read. */
export interface PropertyDef {
  /** The short name, printed (`chain`, `af:plddt`). */
  readonly name: string;
  /** The long names read as it (`label_asym_id`). */
  readonly aliases: readonly string[];
  readonly level: PropertyLevel;
  readonly type: PropertyType;
  /** An extension's reader; the built-ins are read by `structure.ts`. */
  readonly read?:
    | ((model: QueryModel, system: string) => PropertyRead)
    | { readonly values: string; readonly field: string };
}

const def = (name: string, level: PropertyLevel, type: PropertyType, ...aliases: string[]): PropertyDef =>
  Object.freeze({ name, aliases: Object.freeze(aliases), level, type });

/** The built-in properties (properties.md, Built in), in the docs' order. */
export const BuiltInProperties: readonly PropertyDef[] = Object.freeze([
  def('entity', 'molecule', 'word', 'label_entity_id'),
  def('molecule-class', 'molecule', 'word'),
  def('identifier', 'molecule', 'words'),
  def('chain', 'unit', 'word', 'label_asym_id'),
  def('auth-chain', 'unit', 'word', 'auth_asym_id'),
  def('operator', 'unit', 'word', 'instance_id'),
  def('seq', 'residue', 'integer', 'label_seq_id'),
  def('auth-seq', 'residue', 'integer', 'auth_seq_id'),
  def('ins', 'residue', 'word', 'pdbx_PDB_ins_code'),
  def('index', 'residue', 'integer', 'residue_index'),
  def('component', 'residue', 'words', 'label_comp_id'),
  def('component-class', 'residue', 'words'),
  def('standard', 'residue', 'boolean'),
  def('modified', 'residue', 'boolean'),
  def('atom', 'atom', 'word', 'label_atom_id'),
  def('element', 'atom', 'word', 'type_symbol'),
  def('atom-index', 'atom', 'integer', 'atom_index'),
  def('b-factor', 'atom', 'number', 'B_iso_or_equiv'),
  def('occupancy', 'atom', 'number'),
]);

/** The properties a query reads: the built-ins and the registry's, by every name each has. */
export interface PropertyTable {
  readonly byName: ReadonlyMap<string, PropertyDef>;
  readonly defs: readonly PropertyDef[];
}

function table(defs: readonly PropertyDef[]): PropertyTable {
  const byName = new Map<string, PropertyDef>();
  for (const d of defs) for (const n of [d.name, ...d.aliases]) byName.set(n, d);
  return { byName, defs };
}

export const BuiltInPropertyTable: PropertyTable = table(BuiltInProperties);

/** The built-ins with an extension's properties beside them (checked at composition, `query-extensions.ts`). */
export function propertyTable(extra: readonly PropertyDef[]): PropertyTable {
  return extra.length ? table([...BuiltInProperties, ...extra]) : BuiltInPropertyTable;
}

/** A name's spelling for a near match: lowercase, `_` as `-`. */
const loose = (n: string) => n.toLowerCase().replace(/_/g, '-');

/**
 * A property by any of its names; else why not: a namespaced name nobody registered (`missing`), or a bare unknown
 * one with the names in a hint (the short name of a long one misspelled first).
 */
export function findProperty(
  t: PropertyTable,
  name: string,
): { readonly def: PropertyDef } | { readonly missing: string } | { readonly unknown: string } {
  const d = t.byName.get(name);
  if (d) return { def: d };
  const ns = name.indexOf(':');
  if (ns > 0) return { missing: name.slice(0, ns) };
  const near = t.defs.find((x) => [x.name, ...x.aliases].some((n) => loose(n) === loose(name)));
  const names = t.defs.map((x) => x.name).join(', ');
  return { unknown: near ? `did you mean ${near.name}? properties: ${names}` : `properties: ${names}` };
}
