/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement, StructureProperties as P, Unit } from '@molstar/model/model/structure';

// A structure as PDB text: coordinate records only (ATOM, HETATM, TER, END), in the fixed columns of the PDB format
// (v3.3). Mol* has no PDB writer of its own; this one is for tools that read nothing else. What the format cannot hold
// (more than 99,999 atoms, chain ids longer than one character, residue names longer than four) is written as best it
// can and said in `warnings`: mmCIF holds it all.

export interface PdbText {
  readonly text: string;
  readonly warnings: readonly string[];
}

const pad = (s: string, n: number) => (s.length >= n ? s.slice(0, n) : s.padEnd(n));
const right = (s: string, n: number) => (s.length >= n ? s.slice(s.length - n) : s.padStart(n));

/** A number in a fixed column; one too wide for it (a coordinate beyond ±999.999) loses its leading characters: said. */
function fixed(x: number, digits: number, n: number, warnings: Set<string>): string {
  const s = x.toFixed(digits);
  if (s.length > n)
    warnings.add('Values too wide for their PDB columns (large coordinates, B-factors) were cut: use mmCIF.');
  return right(s, n);
}

/** Columns 13-16: a name of four fills them; a shorter one of a one-letter element starts at 14 (the convention). */
function atomName(name: string, element: string): string {
  if (name.length >= 4) return name.slice(0, 4);
  return element.length === 1 ? pad(` ${name}`, 4) : pad(name, 4);
}

function charge(c: number): string {
  if (!c) return '  ';
  return `${Math.min(9, Math.abs(c))}${c > 0 ? '+' : '-'}`;
}

/** The structure's atomic units as PDB records; coarse units (spheres, gaussians) are left out and said. */
export function to_pdb(structure: Structure): PdbText {
  const lines: string[] = [];
  const warnings = new Set<string>();
  const l = StructureElement.Location.create(structure);
  let serial = 0;
  // TER closes a polymer chain: after its last ATOM record, before the next chain or the end.
  let open: { resName: string; chain: string; resSeq: string; iCode: string } | undefined;
  let lastChain = -1;
  let lastUnit = -1;

  const ter = () => {
    if (!open) return;
    serial++;
    lines.push(
      `TER   ${right(String(serial % 100000), 5)}      ${right(open.resName, 3)} ${open.chain}${open.resSeq}${open.iCode}`,
    );
    open = undefined;
  };

  for (const unit of structure.units) {
    if (!Unit.isAtomic(unit)) {
      warnings.add('Coarse parts (spheres, gaussians) have no PDB records: they were left out.');
      continue;
    }
    l.unit = unit;
    const { elements } = unit;
    for (let i = 0, il = elements.length; i < il; i++) {
      l.element = elements[i];
      const chainIndex = unit.chainIndex[l.element];
      if (unit.id !== lastUnit || chainIndex !== lastChain) {
        ter();
        lastUnit = unit.id;
        lastChain = chainIndex;
      }

      const group = P.residue.group_PDB(l);
      const isAtom =
        group === 'ATOM' || (group !== 'HETATM' && P.entity.type(l) === 'polymer' && !P.residue.isNonStandard(l));
      const element = P.atom.type_symbol(l).toUpperCase();
      let chain = P.chain.auth_asym_id(l) || P.chain.label_asym_id(l);
      if (chain.length > 1) {
        warnings.add('Chain ids longer than one character were cut to their first.');
        chain = chain[0];
      }
      chain = chain || ' ';
      const altLoc = P.atom.label_alt_id(l);
      let resName = P.atom.auth_comp_id(l);
      let altCol = pad(altLoc, 1);
      if (resName.length === 4) {
        // A four-character name takes the alternate location's column too.
        if (altLoc) warnings.add('Four-character residue names replaced their atoms’ alternate locations.');
        altCol = resName[0];
        resName = resName.slice(1);
      } else if (resName.length > 4) {
        warnings.add('Residue names longer than four characters were cut.');
        resName = resName.slice(0, 3);
      }
      const seq = P.residue.auth_seq_id(l);
      if (seq > 9999 || seq < -999) warnings.add('Residue numbers beyond 4 digits were cut to their last ones.');
      const resSeq = right(String(seq), 4);
      const iCode = pad(P.residue.pdbx_PDB_ins_code(l), 1);

      // A ligand (water, ions) after its chain's polymer: the polymer ends first. A modified residue (MSE, TPO) is a
      // HETATM inside the polymer and does not end it.
      const polymer = isAtom || P.entity.type(l) === 'polymer';
      if (!polymer) ter();
      serial++;
      if (serial === 100000) warnings.add('More than 99,999 atoms: the serial numbers start again from 0.');
      lines.push(
        (isAtom ? 'ATOM  ' : 'HETATM') +
          right(String(serial % 100000), 5) +
          ' ' +
          atomName(P.atom.auth_atom_id(l), element) +
          altCol +
          right(resName, 3) +
          ' ' +
          chain +
          resSeq +
          iCode +
          '   ' +
          fixed(P.atom.x(l), 3, 8, warnings) +
          fixed(P.atom.y(l), 3, 8, warnings) +
          fixed(P.atom.z(l), 3, 8, warnings) +
          fixed(P.atom.occupancy(l), 2, 6, warnings) +
          fixed(P.atom.B_iso_or_equiv(l), 2, 6, warnings) +
          '          ' +
          right(element, 2) +
          charge(P.atom.pdbx_formal_charge(l)),
      );
      if (polymer) open = { resName, chain, resSeq, iCode };
    }
  }
  ter();
  lines.push('END', '');
  return { text: lines.join('\n'), warnings: [...warnings] };
}
