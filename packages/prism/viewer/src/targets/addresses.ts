/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ComponentExpression } from './types.js';

// Address lists (targets.md §1, §5): `A`, `B:120`, `B:200-210`, `B:120A`, `C:122@OG1`, `A/120`, `A 120-130`,
// `id:1234`, `@#1234`, `index:1233`, `element:FE`, `C:122 element:S`, `A[ASM-X0-1]:120`; a list is separated by commas
// (or semicolons), a trailing one is fine. Author numbering throughout.

/**
 * Residue addresses (`A:315`, author numbering, an insertion code if any) folded into ranges: consecutive numbers of
 * one chain become `A:248-256`; with an insertion code a residue stays on its own. No Mol* here: stories and agents
 * use it where Mol* is not. (Ported from Stories Untold's `any-molql`.)
 */
export function compactResidues(addresses: readonly string[]): string[] {
  const out: string[] = [];
  let run: { chain: string; start: number; end: number } | undefined;
  const flush = () => {
    if (run) out.push(`${run.chain}:${run.start}${run.end !== run.start ? `-${run.end}` : ''}`);
    run = undefined;
  };
  for (const a of addresses) {
    const m = /^([^:]+):(-?\d+)$/.exec(a);
    if (!m) {
      flush();
      out.push(a);
      continue;
    }
    const chain = m[1] as string;
    const n = Number(m[2]);
    if (run && run.chain === chain && n === run.end + 1) run.end = n;
    else {
      flush();
      run = { chain, start: n, end: n };
    }
  }
  flush();
  return out;
}

/**
 * The residues of one chain (author numbers and insertion codes, in any order) as compact ranges: `['120-130', '135',
 * '140A']`. Consecutive numbers fold; a residue with an insertion code stays on its own.
 */
export function compactNumbers(residues: readonly { seq: number; ins: string }[]): string[] {
  const sorted = [...residues].sort((a, b) => a.seq - b.seq || (a.ins < b.ins ? -1 : a.ins > b.ins ? 1 : 0));
  const out: string[] = [];
  let start = 0;
  let end = 0;
  let open = false;
  const flush = () => {
    if (open) out.push(end !== start ? `${start}-${end}` : `${start}`);
    open = false;
  };
  let last: { seq: number; ins: string } | undefined;
  for (const r of sorted) {
    if (last && last.seq === r.seq && last.ins === r.ins) continue;
    last = r;
    if (r.ins) {
      flush();
      out.push(`${r.seq}${r.ins}`);
    } else if (open && r.seq === end + 1) end = r.seq;
    else if (open && r.seq === end) continue;
    else {
      flush();
      start = end = r.seq;
      open = true;
    }
  }
  flush();
  return out;
}

export type AddressResult =
  | { ok: true; items: ComponentExpression[] }
  | {
      ok: false;
      /** The part that did not read, and why. */
      message: string;
      hint: string;
      /** Whether the text looked like an address list (so its message, not a query's, is the one to show). */
      addressy: boolean;
    };

/** Words of selection languages that a chain id must not be (so `chain A` and `resi 100` are query text). */
const ReservedChains = new Set(
  (
    'chain resi resn resid resname resnum name elem element segi segid id index rank type alt altloc ss model object ' +
    'state none and or not of in by to within around near same like byres bychain byobject bymolecule bysegi bound ' +
    'protein nucleic polymer solvent organic metals hetatm backbone sidechain hydrogens atom residue helix sheet ' +
    'all not'
  ).split(' '),
);

const ChainPattern = /^[A-Za-z0-9]{1,4}$/;
const ResiduePattern = /^(-?\d+)([A-Za-z]?)(?:\s*(?:-|\.\.|–|—|to)\s*(-?\d+)([A-Za-z]?))?$/i;

class AddressError extends Error {
  constructor(
    message: string,
    readonly hint: string,
    readonly addressy: boolean,
  ) {
    super(message);
  }
}

const AddressHint =
  'Write chain and residues in author numbering: "A", "B:120", "B:200-210", "B:120A", "C:122@OG1", "A 120-130"; ' +
  'separate several with commas.';

/**
 * Reads an address list. Permissive as Stories reads small models' input: case, spaces, `A/120`, `A 120-130`, a
 * trailing comma. A residue-only item after an item with a chain belongs to that chain (`A:120, 130-140`). What does
 * not read fails with the part that did not.
 */
export function parseAddressList(text: string): AddressResult {
  const parts = text
    .split(/[,;\n]+/)
    .map((p) => p.trim())
    .filter(Boolean);
  if (parts.length === 0) {
    return { ok: false, message: 'The address list is empty.', hint: AddressHint, addressy: false };
  }
  const items: ComponentExpression[] = [];
  let chain: string | undefined;
  try {
    for (const part of parts) {
      const item = parseItem(part, chain);
      if (item.chain !== undefined) chain = item.chain;
      items.push(item.expression);
    }
  } catch (e) {
    if (e instanceof AddressError) return { ok: false, message: e.message, hint: e.hint, addressy: e.addressy };
    throw e;
  }
  return { ok: true, items };
}

interface Item {
  expression: ComponentExpression;
  /** The chain this item names (the next residue-only item belongs to it). */
  chain?: string;
}

function parseItem(part: string, stickyChain: string | undefined): Item {
  // `C:122@OG1` is `C:122 @OG1`; key:value pairs lose the space after the colon.
  const normalized = part
    .replace(/\s*@\s*/g, ' @')
    .replace(/\b(id|index|element|atom)\s*:\s*/gi, '$1:')
    .replace(/\s+/g, ' ')
    .trim();
  const expression: ComponentExpression = {};
  const rest: string[] = [];
  const once = <K extends keyof ComponentExpression>(field: K, value: NonNullable<ComponentExpression[K]>) => {
    if (expression[field] !== undefined) {
      throw new AddressError(`"${part}" gives ${field} twice.`, 'Put one atom, id, index or element in an item.', true);
    }
    expression[field] = value;
  };
  for (const token of normalized.split(' ')) {
    let m: RegExpExecArray | null;
    if ((m = /^@#(.*)$/.exec(token))) once('atom_id', count(part, 'atom id', m[1] as string));
    else if ((m = /^@(.+)$/.exec(token))) once('auth_atom_id', (m[1] as string).toUpperCase());
    else if (token === '@') throw new AddressError(`"${part}" has "@" without an atom name.`, AddressHint, true);
    else if ((m = /^id:(.*)$/i.exec(token))) once('atom_id', count(part, 'atom id', m[1] as string));
    else if ((m = /^index:(.*)$/i.exec(token))) once('atom_index', count(part, 'atom index', m[1] as string));
    else if ((m = /^element:(.*)$/i.exec(token))) {
      const symbol = (m[1] as string).toUpperCase();
      if (!/^[A-Z][A-Z]?[A-Z]?$/.test(symbol)) {
        throw new AddressError(
          `"${m[1]}" is not an element symbol (in "${part}").`,
          'Write element:FE, element:S.',
          true,
        );
      }
      once('type_symbol', symbol);
    } else if ((m = /^atom:(.+)$/i.exec(token))) once('auth_atom_id', (m[1] as string).toUpperCase());
    else rest.push(token);
  }

  let named: string | undefined;
  let residueOnly = false;
  const address = rest.join(' ');
  if (address) {
    const residues = readResidues(part, address, expression, stickyChain);
    named = residues.chain;
    residueOnly = residues.residueOnly;
  }
  if (residueOnly && expression.auth_asym_id === undefined && stickyChain !== undefined) {
    expression.auth_asym_id = stickyChain;
  }
  if (Object.keys(expression).length === 0) {
    throw new AddressError(`Nothing to select in "${part}".`, AddressHint, false);
  }
  return { expression, chain: named };
}

/** A whole number written in an address ("1234"), or an error saying which part is not one. */
function count(part: string, what: string, text: string): number {
  if (!/^\d+$/.test(text)) {
    throw new AddressError(
      `"${text}" is not an ${what} (in "${part}").`,
      `An ${what} is a whole number: ${what === 'atom id' ? 'id:1234 or @#1234' : 'index:1233'}. List several as separate items.`,
      true,
    );
  }
  return Number(text);
}

/** Reads `A`, `A:120`, `A/120`, `A 120-130`, `A[ASM-X0-1]:120`, `120-130`, `A120` into `expression`. */
function readResidues(
  part: string,
  address: string,
  expression: ComponentExpression,
  stickyChain: string | undefined,
): { chain?: string; residueOnly: boolean } {
  const addressy = /[:/[\]@]|\d/.test(part);
  // Residues only: `120`, `-5`, `120-130`, `120A`.
  if (/^-?\d/.test(address) && !/[:/]/.test(address)) {
    setResidues(part, address, expression);
    return { residueOnly: true, chain: stickyChain };
  }
  let chain: string;
  let instance: string | undefined;
  let residues = '';
  // "A120", "A120-130": one letter and the residues, no separator.
  const glued = /^([A-Za-z])(-?\d+[A-Za-z]?(?:-\d+[A-Za-z]?)?)$/.exec(address);
  let m = /^([A-Za-z0-9]+)\[([^\]]+)\]\s*(?:[:/]|\s)?\s*(.*)$/.exec(address);
  if (glued) {
    chain = glued[1] as string;
    residues = glued[2] as string;
  } else if (m) {
    chain = m[1] as string;
    instance = m[2] as string;
    residues = m[3] as string;
  } else {
    m = /^([A-Za-z0-9]+)\s*(?:([:/])\s*(.*)|\s+(.*))?$/.exec(address);
    if (!m) {
      throw new AddressError(
        `Could not read "${address}" as chain and residues (in "${part}").`,
        AddressHint,
        addressy,
      );
    }
    chain = m[1] as string;
    residues = (m[3] ?? m[4] ?? '') as string;
  }
  const bare = !instance && !residues && !part.includes(':') && !part.includes('/');
  if (ReservedChains.has(chain.toLowerCase()) || (bare && !ChainPattern.test(chain))) {
    throw new AddressError(
      `"${chain}" is not a chain id (in "${part}").`,
      `${AddressHint} Chain ids are author chain ids of up to 4 characters.`,
      false,
    );
  }
  if (!ChainPattern.test(chain)) {
    throw new AddressError(
      `"${chain}" is not a chain id (in "${part}").`,
      `${AddressHint} Chain ids are author chain ids of up to 4 characters.`,
      addressy,
    );
  }
  expression.auth_asym_id = chain;
  if (instance) expression.instance_id = instance;
  if (residues) setResidues(part, residues, expression);
  return { chain, residueOnly: false };
}

function setResidues(part: string, text: string, expression: ComponentExpression) {
  const m = ResiduePattern.exec(text.trim());
  if (!m) {
    throw new AddressError(
      `"${text}" is not a residue number or range (in "${part}").`,
      'Write residues as 120, 120A (insertion code) or 200-210.',
      true,
    );
  }
  const begSeq = Number(m[1]);
  const begIns = (m[2] ?? '').toUpperCase();
  if (m[3] === undefined) {
    expression.auth_seq_id = begSeq;
    if (begIns) expression.pdbx_PDB_ins_code = begIns;
    return;
  }
  const endSeq = Number(m[3]);
  const endIns = (m[4] ?? '').toUpperCase();
  if (begIns || endIns) {
    // Insertion codes make it a range in the order of the file.
    expression.beg_auth_seq_id = begSeq;
    expression.end_auth_seq_id = endSeq;
    if (begIns) expression.beg_auth_ins_code = begIns;
    if (endIns) expression.end_auth_ins_code = endIns;
  } else {
    expression.beg_auth_seq_id = Math.min(begSeq, endSeq);
    expression.end_auth_seq_id = Math.max(begSeq, endSeq);
  }
}
