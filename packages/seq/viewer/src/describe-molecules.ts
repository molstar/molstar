/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a molecule's place is, in words (_spec/README.md "HoverModel", Molecules; "Input", Announcements): a modified
// residue's component and name (`PSU (pseudouridine)`), its alternatives (`H2U | U`), a molecule's residue by its code,
// an atom by its name and element with its residue; and the molecule a place is on, with its identifiers (what a
// host's card draws a structure from). Read from the static model's molecule entries; nothing is held.

import type { StaticModel } from '@molstar/seqx-core/resolve';

/** A row of a card and the words the announcer says for it. */
export interface ChemistryRow {
  readonly label: string;
  readonly value: string;
  readonly words: string;
}

/** A component's code with its name when the spec or the host knows one: `PSU (pseudouridine)`. */
function named(model: StaticModel, code: string): string {
  const name = model.compound(code)?.name;
  return name ? `${code} (${name})` : code;
}

/**
 * The chemistry at a position of a system: on an atoms system the atom (name and element) and its residue; on a
 * molecule's residues (or a unit of one) the residue's component; on a sequence (or a unit of one) a modified
 * position's component, or its alternatives (microheterogeneity), none for a standard one.
 */
export function chemistryRows(model: StaticModel, system: string, position: number): ChemistryRow[] {
  const kind = model.systems.get(system)?.kind;
  if (kind === 'atoms') {
    const atoms = model.atoms?.(system);
    if (!atoms || position < 0 || position >= atoms.length) return [];
    const name = atoms.names.get(position) ?? '';
    const element = atoms.element.get(position) ?? '';
    const r = atoms.residue(position);
    const code = model.components?.(atoms.of)?.at(r)[0] ?? '';
    const atom = element ? `${name} (${element})` : name;
    const residue = `${code ? `${code} ` : ''}${r + 1}`;
    return [
      { label: 'Atom', value: atom, words: `atom ${atom}` },
      { label: 'Residue', value: residue, words: `residue ${residue}` },
    ];
  }
  const c = model.components?.(system);
  if (!c || position < 0 || position >= c.length) return [];
  const at = c.at(position);
  if (!at.length) return [];
  if (!c.polymer) {
    const value = named(model, at[0]);
    return [{ label: 'Component', value, words: `component ${value}` }];
  }
  // a polymer's position is said only where it is not its letter's standard component
  const std = c.standard?.(position);
  if (at.length === 1 && at[0] === std) return [];
  const value = at.length > 1 ? at.join(' | ') : named(model, at[0]);
  return [{ label: at.length > 1 ? 'Components' : 'Component', value, words: `component ${value}` }];
}

/** The molecule a place is on, as a host's card reads it: its address, name, class and identifiers. */
export interface PlaceMolecule {
  readonly address: string;
  readonly name?: string;
  readonly class: string;
  readonly identifiers: readonly { readonly namespace: string; readonly value: string }[];
  /** The unit it is shown as, when the place is on a unit of a complex. */
  readonly unit?: string;
}

/** The molecule (not a sequence) a system is of: a molecule's residues, its atoms, a unit of one or its atoms. */
export function moleculeOfSystem(model: StaticModel, system: string): PlaceMolecule | undefined {
  const m = model.molecular;
  const residues = m.atoms.get(system)?.of ?? system;
  const unit = m.units.get(residues);
  const e = m.molecules.get(unit ? unit.of : residues);
  if (!e) return undefined;
  return {
    address: e.address,
    ...(e.name !== undefined ? { name: e.name } : {}),
    class: e.class,
    identifiers: e.identifiers,
    ...(unit ? { unit: unit.address } : {}),
  };
}
