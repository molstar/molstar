/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Entry ids a load names (PyMOL's `fetch`, VMD's `mol pdbload`): the PDB's, as 4 characters or extended.

/** What a PDB id is, for a hint. */
export const PdbIdHint = 'PDB ids are 4 characters (1cbs) or extended: pdb_ and 8 (pdb_00001cbs).';

/**
 * The PDB id `text` names, lowercase, or undefined: `1cbs`, or the extended `pdb_00001cbs` (wwPDB's 12 characters).
 * An extended id of the four-character era (`pdb_0000…`) is the same entry as its last four: loaded by those, which
 * every mirror serves; a later one goes as written.
 */
export function pdbId(text: string): string | undefined {
  const id = text.toLowerCase();
  if (/^[0-9][a-z0-9]{3}$/.test(id)) return id;
  const extended = /^pdb_([0-9]{4})([0-9][a-z0-9]{3})$/.exec(id);
  if (!extended) return undefined;
  return extended[1] === '0000' ? extended[2] : id;
}
