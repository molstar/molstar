/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MVSAnnotationRow } from './schemas.js';
import { isDefined } from './utils.js';

/** Longest label made of the descriptions of individual rows (the labels are shown in narrow places, e.g. buttons) */
const MaxListedLength = 32;
/** Most rows whose descriptions are listed one by one */
const MaxListedRows = 3;
/** Most chains listed by name when each row is just a chain */
const MaxListedChains = 6;

const AuthSuffix = ' (auth)';

/** A piece of a description. `auth` is `true` for `auth_*` criteria, `false` for `label_*` ones, and `undefined` for criteria which do not have the two variants. */
interface Token {
  text: string;
  auth?: boolean;
}

function describeRange(beg: number | undefined | null, end: number | undefined | null): string | undefined {
  if (isDefined(beg) && isDefined(end)) return beg === end ? `${beg}` : `${beg}-${end}`;
  if (isDefined(beg)) return `${beg}-`;
  if (isDefined(end)) return `-${end}`;
  return undefined;
}

/** The pieces of the description of a row, in the order in which they are shown. */
function tokensOfRow(row: MVSAnnotationRow): Token[] {
  const tokens: Token[] = [];
  if (isDefined(row.label_entity_id)) tokens.push({ text: `Entity ${row.label_entity_id}` });

  const chains: Token[] = [];
  if (isDefined(row.label_asym_id)) chains.push({ text: row.label_asym_id, auth: false });
  if (isDefined(row.auth_asym_id)) chains.push({ text: row.auth_asym_id, auth: true });

  // residues: `single` ones can be written right after the chain name (A432), ranges follow it after a colon (A:10-50)
  const residues: { token: Token; single: boolean }[] = [];
  if (isDefined(row.label_seq_id)) residues.push({ token: { text: `${row.label_seq_id}`, auth: false }, single: true });
  if (isDefined(row.auth_seq_id)) {
    residues.push({ token: { text: `${row.auth_seq_id}${row.pdbx_PDB_ins_code ?? ''}`, auth: true }, single: true });
  } else if (isDefined(row.pdbx_PDB_ins_code)) {
    const text = row.pdbx_PDB_ins_code === '' ? 'no ins code' : `ins ${row.pdbx_PDB_ins_code}`;
    residues.push({ token: { text }, single: false });
  }
  const labelRange = describeRange(row.beg_label_seq_id, row.end_label_seq_id);
  if (labelRange) residues.push({ token: { text: labelRange, auth: false }, single: false });
  const authRange = describeRange(row.beg_auth_seq_id, row.end_auth_seq_id);
  if (authRange) residues.push({ token: { text: authRange, auth: true }, single: false });
  if (isDefined(row.residue_index)) residues.push({ token: { text: `index ${row.residue_index}` }, single: false });

  if (chains.length === 1 && residues.length === 1 && chains[0].auth === residues[0].token.auth) {
    // the common cases: chain + residue (A432) and chain + range (A:10-50)
    const { token, single } = residues[0];
    tokens.push({ text: `${chains[0].text}${single ? '' : ':'}${token.text}`, auth: token.auth });
  } else {
    tokens.push(...chains);
    for (const { token } of residues) {
      tokens.push({
        ...token,
        text: chains.length === 0 && token.auth !== undefined ? `res ${token.text}` : token.text,
      });
    }
  }

  if (isDefined(row.label_comp_id)) tokens.push({ text: row.label_comp_id, auth: false });
  if (isDefined(row.auth_comp_id) && row.auth_comp_id !== row.label_comp_id) {
    tokens.push({ text: row.auth_comp_id, auth: true });
  }
  if (isDefined(row.label_atom_id)) tokens.push({ text: row.label_atom_id, auth: false });
  if (isDefined(row.auth_atom_id)) tokens.push({ text: row.auth_atom_id, auth: true });
  if (isDefined(row.type_symbol)) tokens.push({ text: `${row.type_symbol.toUpperCase()} atoms` });
  if (isDefined(row.atom_id)) tokens.push({ text: `atom ${row.atom_id}` });
  if (isDefined(row.atom_index)) tokens.push({ text: `atom index ${row.atom_index}` });
  if (isDefined(row.instance_id)) tokens.push({ text: `instance ${row.instance_id}` });
  if (isDefined(row.operator_name)) tokens.push({ text: `operator ${row.operator_name}` });
  return tokens;
}

/** The description of a row without the `(auth)` marker, and whether the row is about `auth_*` criteria only
 * (in that case the marker applies to all of it and is written once, otherwise each `auth_*` piece is marked). */
function describeRowParts(row: MVSAnnotationRow): { text: string; auth: boolean } {
  const tokens = tokensOfRow(row);
  if (tokens.length === 0) return { text: 'All', auth: false };
  const variants = tokens.filter((t) => t.auth !== undefined);
  const auth = variants.length > 0 && variants.every((t) => t.auth);
  const text = tokens.map((t) => (t.auth && !auth ? `${t.text}${AuthSuffix}` : t.text)).join(' ');
  return { text, auth };
}

/** A short description of the criteria in one row: `A` (a chain), `A432` (a residue), `A:10-50` (residues), `A432 CA` (an atom),
 * `HEM` (a component), with ` (auth)` after criteria which use the `auth_*` numbering. A row without criteria selects everything (`All`). */
export function describeRow(row: MVSAnnotationRow): string {
  const { text, auth } = describeRowParts(row);
  return auth ? `${text}${AuthSuffix}` : text;
}

type RowKind = 'chain' | 'residue' | 'other';

/** Is the row just a chain, just a single residue in a chain, or something else? */
function kindOf(row: MVSAnnotationRow): RowKind {
  const defined = Object.keys(row).filter((key) => key !== 'group_id' && isDefined((row as any)[key]));
  const set = new Set(defined);
  const chain = set.has('label_asym_id') ? 'label_asym_id' : set.has('auth_asym_id') ? 'auth_asym_id' : undefined;
  if (!chain) return 'other';
  if (defined.length === 1) return 'chain';
  const seq = chain === 'label_asym_id' ? 'label_seq_id' : 'auth_seq_id';
  const allowed = new Set([chain, seq, ...(seq === 'auth_seq_id' ? ['pdbx_PDB_ins_code'] : [])]);
  return set.has(seq) && defined.every((key) => allowed.has(key)) ? 'residue' : 'other';
}

/** A short name for a selection made of `rows` (an element is selected if it satisfies any of them), for labels.
 * A few rows are described one by one (`A`, `A:10-50`, `A, B, C`, `A432 + A500`),
 * more of them are only counted (`120 residues`, `A: 120 residues`, `30 chains`, `120 selections`). */
export function describeRows(rows: readonly MVSAnnotationRow[]): string {
  if (rows.length === 0) return 'Empty Selection';
  if (rows.length === 1) return describeRow(rows[0]);

  const parts = rows.map(describeRowParts);
  const allAuth = parts.every((p) => p.auth);
  const kinds = rows.map(kindOf);
  const allChains = kinds.every((k) => k === 'chain');
  // when every row is about `auth_*` criteria, the marker is written once at the end
  const marker = allAuth ? AuthSuffix : '';
  const texts = parts.map((p) => (allAuth ? p.text : p.auth ? `${p.text}${AuthSuffix}` : p.text));

  const maxListed = allChains ? MaxListedChains : MaxListedRows;
  if (rows.length <= maxListed) {
    const listed = `${texts.join(allChains ? ', ' : ' + ')}${marker}`;
    if (listed.length <= MaxListedLength) return listed;
  }

  // counted: chains, residues (of one chain or of several), or anything
  const noun = allChains ? 'chains' : kinds.every((k) => k === 'residue') ? 'residues' : 'selections';
  const chain = rows[0].label_asym_id ?? rows[0].auth_asym_id;
  const sameChain =
    isDefined(chain) &&
    !allChains &&
    rows.every(
      (r) =>
        (r.label_asym_id ?? r.auth_asym_id) === chain &&
        isDefined(r.label_asym_id) === isDefined(rows[0].label_asym_id),
    );
  const prefix = sameChain ? `${chain}: ` : '';
  return `${prefix}${rows.length} ${noun}${marker}`;
}
