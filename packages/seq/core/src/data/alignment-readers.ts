/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The alignment readers (_spec/README.md "Data", Alignment readers): aligned FASTA, Stockholm and Newick, in their own
// entry (a reader of a spec never pays for them). Each reads text into plain data and reports what it cannot read, by
// line, never throwing (EC-40): a row of another width is left out and said, a tree that fails is said and the
// alignment kept (EC-43).

import { issueOf } from '@molstar/core/repr/tree';
import { firstNumberOf, openAlignment } from './alignment-source.js';
import { SourceIssues } from './issues.js';
import { type FormatRuntime, isAbortError } from './protocol.js';
import { SourceOpenError } from './seq-runtime.js';

/** What a reader could not read: a line (1-based, 0 for the whole text) and why. */
export interface ReadIssue {
  readonly line: number;
  readonly message: string;
}

/** A row of an alignment: its label, its aligned letters (gaps `-` and `.`), its first and last residue numbers. */
export interface AlignmentRow {
  readonly label: string;
  readonly letters: string;
  /** The first residue's number from a `NAME/27-137` label (EC-39); 1 when the label says none. */
  readonly first: number;
  /** Its own fields (Stockholm's `#=GS <name> <tag> <value>`: `OS`, `DE`, `AC`), by tag; absent when none. */
  readonly fields?: Readonly<Record<string, string>>;
}

/** An alignment as read: its rows, its width, its insert columns and the pairs a consensus structure names. */
export interface AlignmentRead {
  readonly rows: readonly AlignmentRow[];
  readonly columns: number;
  /**
   * Insert columns (EC-38): `#=GC RF`'s `.` and `-` where it is given, else a column whose letters are all lowercase;
   * absent when none is.
   */
  readonly inserts?: Uint8Array;
  /** Base pairs of `#=GC SS_cons` as column pairs `[i, j]` with `i < j` (nested brackets and pseudoknot letters). */
  readonly pairs?: readonly (readonly [number, number])[];
  readonly issues: readonly ReadIssue[];
}

export { firstNumberOf };

/**
 * The rows of the alignment's width (the most common one, on a tie the first to reach it: one stray row never sets it), the
 * rest said and left out; a label twice said (its rows kept: a row is addressed by its index too); insert columns from
 * lowercase.
 */
function finish(
  read: { label: string; letters: string; line: number }[],
  issues: ReadIssue[],
  rf?: string,
  ss?: string,
  gs?: ReadonlyMap<string, Readonly<Record<string, string>>>,
): AlignmentRead {
  const widths = new Map<number, number>();
  let columns = 0;
  let most = 0;
  for (const r of read) {
    const n = (widths.get(r.letters.length) ?? 0) + 1;
    widths.set(r.letters.length, n);
    if (n > most) {
      most = n;
      columns = r.letters.length;
    }
  }
  const rows: AlignmentRow[] = [];
  const labels = new Set<string>();
  for (const r of read) {
    if (labels.has(r.label)) issues.push({ line: r.line, message: `the label ${r.label} names more than one row` });
    labels.add(r.label);
    if (r.letters.length !== columns) {
      issues.push({
        line: r.line,
        message: `the row ${r.label} is ${r.letters.length} columns wide, not ${columns}: left out`,
      });
      continue;
    }
    if (/[^A-Za-z*.\-~]/.test(r.letters)) {
      issues.push({ line: r.line, message: `the row ${r.label} holds a character that is no letter or gap: left out` });
      continue;
    }
    const fields = gs?.get(r.label);
    rows.push({
      label: r.label,
      letters: r.letters.replace(/~/g, '-'),
      first: firstNumberOf(r.label),
      ...(fields ? { fields } : {}),
    });
  }
  let inserts: Uint8Array | undefined;
  if (rf !== undefined && rf.length === columns) {
    inserts = new Uint8Array(columns);
    for (let j = 0; j < columns; j++) inserts[j] = rf[j] === '.' || rf[j] === '-' ? 1 : 0;
  } else {
    if (rf !== undefined) issues.push({ line: 0, message: 'the #=GC RF line is not the alignment’s width: not read' });
    // a column whose letters are all lowercase (gaps aside) is an insert (A3M's lowercase is on no column)
    for (let j = 0; j < columns && rows.length; j++) {
      let lower = false;
      let upper = false;
      for (const r of rows) {
        const c = r.letters.charCodeAt(j);
        if (c >= 97 && c <= 122) lower = true;
        else if (c >= 65 && c <= 90) upper = true;
      }
      if (lower && !upper) (inserts ??= new Uint8Array(columns))[j] = 1;
    }
  }
  let pairs: [number, number][] | undefined;
  if (ss !== undefined) {
    if (ss.length !== columns) issues.push({ line: 0, message: 'the #=GC SS_cons line is not the alignment’s width' });
    else {
      const p = structurePairs(ss);
      pairs = p.pairs;
      if (p.unpaired) issues.push({ line: 0, message: `SS_cons leaves ${p.unpaired} brackets unpaired` });
    }
  }
  return { rows, columns, ...(inserts ? { inserts } : {}), ...(pairs ? { pairs } : {}), issues };
}

const Opens: Readonly<Record<string, string>> = { '<': '>', '(': ')', '[': ']', '{': '}' };

/**
 * The pairs of a dot-bracket structure: `<>`, `()`, `[]`, `{}` each a nesting of its own, and a pseudoknot as an
 * uppercase letter paired with its lowercase; what is left open or closed is counted.
 */
export function structurePairs(ss: string): { pairs: [number, number][]; unpaired: number } {
  const stacks = new Map<string, number[]>();
  const pairs: [number, number][] = [];
  let unpaired = 0;
  for (let i = 0; i < ss.length; i++) {
    const c = ss[i];
    let open: string | undefined;
    let close: string | undefined;
    if (c in Opens) open = c;
    else if (c >= 'A' && c <= 'Z') open = c;
    else if (c >= 'a' && c <= 'z') close = c.toUpperCase();
    else
      for (const o in Opens)
        if (Opens[o] === c) {
          close = o;
          break;
        }
    if (open !== undefined) {
      let s = stacks.get(open);
      if (!s) stacks.set(open, (s = []));
      s.push(i);
    } else if (close !== undefined) {
      const j = stacks.get(close)?.pop();
      if (j === undefined) unpaired++;
      else pairs.push([j, i]);
    }
  }
  for (const s of stacks.values()) unpaired += s.length;
  pairs.sort((a, b) => a[0] - b[0]);
  return { pairs, unpaired };
}

/** Aligned FASTA (and A3M's uppercase columns): `>label` lines, each followed by its aligned letters over lines. */
export function readAlignedFasta(text: string): AlignmentRead {
  const issues: ReadIssue[] = [];
  const read: { label: string; letters: string; line: number }[] = [];
  const lines = typeof text === 'string' ? text.split(/\r?\n/) : [];
  if (typeof text !== 'string') issues.push({ line: 0, message: 'the text is not a string' });
  let cur: { label: string; parts: string[]; line: number } | undefined;
  const flush = () => {
    if (cur) read.push({ label: cur.label, letters: cur.parts.join(''), line: cur.line });
  };
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line || line.startsWith(';')) continue;
    if (line.startsWith('>')) {
      flush();
      const label = line.slice(1).split(/\s/)[0];
      if (!label) issues.push({ line: i + 1, message: 'a record without a label' });
      cur = { label: label || `row ${read.length + 1}`, parts: [], line: i + 1 };
    } else if (!cur) issues.push({ line: i + 1, message: 'letters before any record: left out' });
    else cur.parts.push(line.replace(/\s+/g, ''));
  }
  flush();
  if (!read.length && !issues.length) issues.push({ line: 0, message: 'no records' });
  return finish(read, issues);
}

/**
 * Stockholm 1.0: `# STOCKHOLM 1.0`, `#=GS <name> <tag> <value>` (a row's fields, a tag given twice joined by a
 * space), `#=GF`/`#=GR` annotation (read past), `#=GC RF` (insert columns) and `#=GC SS_cons` (pairs), `name letters`
 * lines in blocks (a name's letters joined across them), `//` its end. The first alignment of the file.
 */
export function readStockholm(text: string): AlignmentRead {
  const issues: ReadIssue[] = [];
  if (typeof text !== 'string')
    return { rows: [], columns: 0, issues: [{ line: 0, message: 'the text is not a string' }] };
  const lines = text.split(/\r?\n/);
  let i = 0;
  while (i < lines.length && !lines[i].trim()) i++;
  if (!/^# STOCKHOLM 1\.\d/.test(lines[i] ?? '')) issues.push({ line: i + 1, message: 'no `# STOCKHOLM 1.0` header' });
  else i++;
  const order: string[] = [];
  const parts = new Map<string, { parts: string[]; line: number }>();
  const gc = new Map<string, string[]>();
  const gs = new Map<string, Record<string, string>>();
  let ended = false;
  for (; i < lines.length; i++) {
    const line = lines[i];
    const t = line.trim();
    if (!t) continue;
    if (t === '//') {
      ended = true;
      break;
    }
    if (t.startsWith('#=GC')) {
      const m = /^#=GC\s+(\S+)\s+(\S+)$/.exec(t);
      if (!m) issues.push({ line: i + 1, message: 'a #=GC line without a feature and its text' });
      else {
        const g = gc.get(m[1]);
        if (g) g.push(m[2]);
        else gc.set(m[1], [m[2]]);
      }
      continue;
    }
    if (t.startsWith('#=GS')) {
      const m = /^#=GS\s+(\S+)\s+(\S+)\s+(.*)$/.exec(t);
      if (!m) issues.push({ line: i + 1, message: 'a #=GS line without a name, a tag and its text' });
      else {
        let f = gs.get(m[1]);
        if (!f) gs.set(m[1], (f = Object.create(null) as Record<string, string>));
        const v = m[3].trim();
        f[m[2]] = Object.hasOwn(f, m[2]) ? `${f[m[2]]} ${v}` : v;
      }
      continue;
    }
    if (t.startsWith('#')) continue;
    const m = /^(\S+)\s+(\S+)$/.exec(t);
    if (!m) {
      issues.push({ line: i + 1, message: 'a line that is no name and letters: left out' });
      continue;
    }
    const p = parts.get(m[1]);
    if (p) p.parts.push(m[2]);
    else {
      order.push(m[1]);
      parts.set(m[1], { parts: [m[2]], line: i + 1 });
    }
  }
  if (!ended) issues.push({ line: lines.length, message: 'the alignment does not end with `//`' });
  const read = order.map((label) => {
    const p = parts.get(label)!;
    return { label, letters: p.parts.join(''), line: p.line };
  });
  if (!read.length) issues.push({ line: 0, message: 'no rows' });
  return finish(read, issues, gc.get('RF')?.join(''), gc.get('SS_cons')?.join(''), gs);
}

/** A tree as read: nodes in pre-order (the root 0), each its parent, branch length, label and children. */
export interface TreeRead {
  /** The parent of each node, -1 for the root. */
  readonly parent: Int32Array;
  /** Each node's branch length to its parent (0 when the text gives none). */
  readonly length: Float64Array;
  /** Each node's label (an internal node's too: its id), `''` when none. */
  readonly label: readonly string[];
  readonly children: readonly (readonly number[])[];
  /** The leaves in the order the text gives them (tree order). */
  readonly leaves: readonly number[];
  readonly issues: readonly ReadIssue[];
}

/** The most nodes a tree is read with: a tree is a figure's, never a genome's (refused past it, said). */
export const MaxTreeNodes = 1_000_000;

/**
 * Newick: `(A:0.1,(B:0.2,C)n1:0.4)root;`, labels bare or quoted (`'a b'`, `''` an escaped quote), `[comments]` read
 * past. Iterative (a deep tree costs no stack); what does not parse is said and the tree is empty.
 */
export function readNewick(text: string): TreeRead {
  const parent: number[] = [];
  const length: number[] = [];
  const label: string[] = [];
  const children: number[][] = [];
  const issues: ReadIssue[] = [];
  const empty = (message: string): TreeRead => ({
    parent: new Int32Array(0),
    length: new Float64Array(0),
    label: [],
    children: [],
    leaves: [],
    issues: [{ line: 0, message }],
  });
  if (typeof text !== 'string') return empty('the text is not a string');
  const add = (p: number) => {
    const n = parent.length;
    if (n >= MaxTreeNodes) return -1;
    parent.push(p);
    length.push(0);
    label.push('');
    children.push([]);
    if (p >= 0) children[p].push(n);
    return n;
  };
  let i = 0;
  const skip = () => {
    for (;;) {
      while (i < text.length && /\s/.test(text[i])) i++;
      if (text[i] !== '[') return;
      const end = text.indexOf(']', i);
      i = end < 0 ? text.length : end + 1;
    }
  };
  const name = (): string | undefined => {
    skip();
    if (text[i] === "'") {
      let out = '';
      i++;
      for (;;) {
        if (i >= text.length) return undefined;
        if (text[i] === "'") {
          if (text[i + 1] === "'") {
            out += "'";
            i += 2;
            continue;
          }
          i++;
          return out;
        }
        out += text[i++];
      }
    }
    const s = i;
    while (i < text.length && !/[\s(),:;[\]']/.test(text[i])) i++;
    return text.slice(s, i);
  };
  const tail = (n: number): string | undefined => {
    const l = name();
    if (l === undefined) return 'an unclosed quote';
    label[n] = l;
    skip();
    if (text[i] === ':') {
      i++;
      skip();
      const s = i;
      while (i < text.length && /[0-9eE+\-.]/.test(text[i])) i++;
      const v = Number(text.slice(s, i));
      if (i === s || !Number.isFinite(v)) return 'a branch length that is no number';
      length[n] = v;
    }
    return undefined;
  };
  let cur = add(-1);
  skip();
  let fail: string | undefined;
  for (;;) {
    skip();
    const c = text[i];
    if (c === '(') {
      i++;
      const n = add(cur);
      if (n < 0) {
        fail = `more than ${MaxTreeNodes} nodes`;
        break;
      }
      cur = n;
      continue;
    }
    // a node ends: its label and length, then `,` a sibling, `)` its parent's end, `;` the tree's
    fail = tail(cur);
    if (fail) break;
    skip();
    const d = text[i];
    if (d === ',') {
      i++;
      const p = parent[cur];
      if (p < 0) {
        fail = 'a `,` outside any parentheses';
        break;
      }
      const n = add(p);
      if (n < 0) {
        fail = `more than ${MaxTreeNodes} nodes`;
        break;
      }
      cur = n;
    } else if (d === ')') {
      i++;
      const p = parent[cur];
      if (p < 0) {
        fail = 'a `)` closes nothing';
        break;
      }
      cur = p;
    } else if (d === ';' || d === undefined) {
      if (parent[cur] >= 0) fail = 'the tree ends inside parentheses';
      else if (d === undefined) issues.push({ line: 0, message: 'the tree does not end with `;`' });
      else i++;
      break;
    } else {
      fail = `an unexpected character ${JSON.stringify(d)}`;
      break;
    }
  }
  if (fail) return empty(`the tree cannot be read at ${i}: ${fail}`);
  skip();
  if (i < text.length) issues.push({ line: 0, message: 'text after the tree’s `;`: read past' });
  const leaves: number[] = [];
  for (let n = 0; n < children.length; n++) if (!children[n].length) leaves.push(n);
  return {
    parent: Int32Array.from(parent),
    length: Float64Array.from(length),
    label,
    children,
    leaves,
    issues,
  };
}

/**
 * The rows of an alignment each leaf names (EC-43): by its label, else by the label before its last `/` (a tree of
 * names against rows of `NAME/27-137`); -1 for a leaf that names none. A row named twice is the first leaf's.
 */
export function matchLeaves(tree: TreeRead, rows: readonly { readonly label: string }[]): Int32Array {
  const byLabel = new Map<string, number>();
  const byName = new Map<string, number>();
  rows.forEach((r, k) => {
    if (!byLabel.has(r.label)) byLabel.set(r.label, k);
    const s = r.label.lastIndexOf('/');
    const n = s < 0 ? r.label : r.label.slice(0, s);
    if (!byName.has(n)) byName.set(n, k);
  });
  const out = new Int32Array(tree.leaves.length).fill(-1);
  const taken = new Set<number>();
  tree.leaves.forEach((leaf, k) => {
    const l = tree.label[leaf];
    const s = l.lastIndexOf('/');
    const row = byLabel.get(l) ?? byName.get(l) ?? (s < 0 ? undefined : byName.get(l.slice(0, s)));
    if (row !== undefined && !taken.has(row)) {
      out[k] = row;
      taken.add(row);
    }
  });
  return out;
}

/** The rows' own fields as columns, a tag a column (a row without the tag missing there); none when no row has one. */
function rowFields(rows: readonly AlignmentRow[]): { fields?: Record<string, (string | null)[]> } {
  const tags: string[] = [];
  for (const r of rows) if (r.fields) for (const t of Object.keys(r.fields)) if (!tags.includes(t)) tags.push(t);
  if (!tags.length) return {};
  const fields: Record<string, (string | null)[]> = Object.create(null);
  for (const t of tags) fields[t] = rows.map((r) => (r.fields && Object.hasOwn(r.fields, t) ? r.fields[t] : null));
  return { fields };
}

/**
 * The runtime of alignment files (`aligned-fasta`, `stockholm`): the file read whole through the assets port (or its
 * text inline as `value`) and opened as an alignment source (`openAlignment`); a file that leaves no row is refused
 * with the reader's first issue.
 */
export const AlignmentFileRuntime: FormatRuntime = {
  name: 'alignment-file',
  formats: ['aligned-fasta', 'stockholm'],
  async open(source, ctx) {
    const fail = (code: 'seq.source.unreadable' | 'seq.source.no_data' | 'seq.source.malformed', message: string) =>
      new SourceOpenError([issueOf(SourceIssues, code, message, { address: source.address })]);
    if (source.uri === undefined && typeof source.value !== 'string')
      throw fail('seq.source.no_data', 'An alignment file needs a uri or its text inline');
    let text: string;
    if (typeof source.value === 'string') text = source.value;
    else
      try {
        text = await ctx.assets.readText(source.uri!, ctx.signal);
      } catch (e) {
        if (isAbortError(e) || ctx.signal.aborted) throw e;
        throw fail('seq.source.unreadable', `The file ${source.uri} could not be read`);
      }
    const read = source.format === 'stockholm' ? readStockholm(text) : readAlignedFasta(text);
    if (!read.rows.length)
      throw fail('seq.source.malformed', `The alignment has no row: ${read.issues[0]?.message ?? 'empty'}`);
    const rows = read.rows;
    return openAlignment(
      {
        label: rows.map((r) => r.label),
        letters: rows.map((r) => r.letters),
        first: rows.map((r) => r.first),
        ...rowFields(rows),
      },
      source,
    );
  },
};
