/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's text formats (seqspec:sources#data-by-reference): FASTA (aligned FASTA and A3M are FASTA records),
// Stockholm with its tags, WUSS pairs, Clustal and Newick. Parsers run at the edges: they allocate, report issues and
// never throw. Codes mirror the document's (`model.source.format` is `svs.source.format`, and so on).

import type { AlignmentInput, ParsedTree } from "./alignment.js";
import { Issue } from "./common.js";

export interface ParseOptions {
  /** A size bound in UTF-8 bytes; absent: none (a source's `value` is bound at 1 MiB). */
  maxBytes?: number;
  /** The JSON pointer issues are reported at. */
  path?: string;
}

/** The record id is the header's first word; the rest is its description. Sequences have no whitespace. */
export interface FastaRecord {
  id: string;
  description: string;
  sequence: string;
}

/** One alignment of a Stockholm file; the tags are kept as read (per-column ones with their blocks joined). */
export interface StockholmRecord {
  /** Row names in first-seen order ("HBA_HUMAN/1-141"), and their gapped rows. */
  names: string[];
  rows: string[];
  /** `#=GF feature text`: lines of one feature joined with "\n". */
  gf: Map<string, string>;
  /** `#=GC feature data` per column (`RF`, `SS_cons`, …). */
  gc: Map<string, string>;
  /** `#=GS name feature text`, by row name. */
  gs: Map<string, Map<string, string>>;
  /** `#=GR name feature data` per column, by row name. */
  gr: Map<string, Map<string, string>>;
}

/** A Clustal alignment (`.aln`): its header line, and its rows by name with their blocks joined in order. */
export interface ClustalRecord {
  header: string;
  names: string[];
  rows: string[];
}

const FORMAT = "model.source.format";
const EMPTY_RECORD = "model.source.empty_record";
/** A record's id missing or repeated: the records are read, and the id may not name one of them. */
const ID = "model.source.id";
/** A file whose last record or tree lacks its terminator: read as if it had one. */
const UNTERMINATED = "model.source.unterminated";
/** `;` (Pearson) and `#` lines are comments, anywhere in a FASTA file. */
const COMMENT = /^\s*[;#]/;
/**
 * HH-suite's annotation records in an A3M: secondary structure (`ss_pred`, `ss_conf`, `ss_dssp`), DSSP's solvent
 * accessibility (`sa_dssp`, which `addss.pl` writes too) and the consensus.
 */
const HHSUITE = /^s[sa]_|^Consensus$/;

/** UTF-8 bytes of `text` within the bound; else an error. Counts exactly only when the length alone cannot tell. */
function withinBound(text: string, maxBytes: number | undefined, issues: Issue[], path: string): boolean {
  if (maxBytes == null || text.length * 3 <= maxBytes) return true;
  let bytes = text.length > maxBytes ? text.length : 0;
  for (let i = 0; i < text.length && bytes <= maxBytes; i++) {
    const c = text.charCodeAt(i);
    const surrogate = c >= 0xd800 && c < 0xdc00;
    bytes += c < 0x80 ? 1 : c < 0x800 ? 2 : surrogate ? 4 : 3;
    if (surrogate) i++;
  }
  if (bytes <= maxBytes) return true;
  issues.push(Issue.error("model.source.size_bound", path, `The source is over its bound of ${maxBytes} bytes.`));
  return false;
}

/** The first few of a repeated cause's examples, and how many more (GV-07: one issue per cause, not one per record). */
function firstFew(items: readonly string[], max = 5): string {
  const shown = items.slice(0, max).join(", ");
  return items.length > max ? `${shown} and ${items.length - max} more` : shown;
}

/** LF, CRLF and CR line ends. */
const lines = (text: string) => text.split(/\r\n|\r|\n/);

function parseFasta(text: string, options: ParseOptions = {}): { records: FastaRecord[]; issues: Issue[] } {
  const issues: Issue[] = [];
  const path = options.path ?? "";
  const records: FastaRecord[] = [];
  if (!withinBound(text, options.maxBytes, issues, path)) return { records, issues };
  let parts: string[] | null = null;
  let id = "";
  let description = "";
  let headerLine = 0;
  // Each cause said once for the file, with its count and first examples (GV-07: a hit list repeats ids by the
  // thousand).
  const empty: string[] = [];
  const noId: string[] = [];
  const before: string[] = [];
  const flush = () => {
    if (!parts) return;
    const sequence = parts.join("");
    records.push({ id, description, sequence });
    // Kept, so that a record index does not shift; a node that reads it is `model.molecule.content`.
    if (!sequence) empty.push(`"${id}" (line ${headerLine})`);
  };
  lines(text).forEach((line, k) => {
    if (line.startsWith(">")) {
      flush();
      const header = line.slice(1).trim();
      const space = header.search(/\s/);
      id = space < 0 ? header : header.slice(0, space);
      description = space < 0 ? "" : header.slice(space + 1).trim();
      if (!id) noId.push(String(k + 1));
      headerLine = k + 1;
      parts = [];
    } else if (!COMMENT.test(line)) {
      if (parts) parts.push(line.replace(/\s+/g, ""));
      else if (line.trim()) before.push(String(k + 1));
    }
  });
  flush();
  if (before.length > 0)
    issues.push(Issue.error(FORMAT, path, `FASTA: a sequence before the first header (line ${firstFew(before)}).`));
  if (noId.length > 0)
    issues.push(Issue.warning(ID, path, `FASTA: ${noId.length} header(s) with no id (line ${firstFew(noId)}).`));
  if (empty.length > 0)
    issues.push(
      Issue.info(EMPTY_RECORD, path, `FASTA: ${empty.length} record(s) with no sequence: ${firstFew(empty)}.`),
    );
  if (!records.length && !issues.length) issues.push(Issue.error(FORMAT, path, "FASTA: no records."));
  const ids = new Set<string>();
  const repeated = new Set<string>();
  for (const r of records) {
    if (ids.has(r.id)) repeated.add(r.id);
    ids.add(r.id);
  }
  if (repeated.size > 0) {
    const named = firstFew([...repeated].map((x) => `"${x}"`));
    issues.push(Issue.warning(ID, path, `FASTA: ${repeated.size} id(s) name several records: ${named}.`));
  }
  return { records, issues };
}

/**
 * Aligned FASTA or A3M records as an alignment's input; labels are the record ids. In an A3M, HH-suite's `ss_*`, `sa_*`
 * and `Consensus` records are not rows: one `model.source.skipped` note says how many were left out.
 */
function fastaAlignment(
  records: readonly FastaRecord[],
  key: string,
  flavor: "aligned" | "a3m",
  options: ParseOptions = {},
): { input: AlignmentInput; issues: Issue[] } {
  const issues: Issue[] = [];
  const rows = flavor === "a3m" ? records.filter((r) => !HHSUITE.test(r.id)) : records;
  if (rows.length < records.length) {
    const m = `${records.length - rows.length} HH-suite record(s) (ss_*, sa_*, Consensus) are not alignment rows.`;
    issues.push(Issue.info("model.source.skipped", options.path ?? "", `FASTA: ${m}`));
  }
  return { input: { key, flavor, gapped: rows.map((r) => r.sequence), labels: rows.map((r) => r.id) }, issues };
}

function parseStockholm(text: string, options: ParseOptions = {}): { records: StockholmRecord[]; issues: Issue[] } {
  const issues: Issue[] = [];
  const path = options.path ?? "";
  const records: StockholmRecord[] = [];
  if (!withinBound(text, options.maxBytes, issues, path)) return { records, issues };
  let rec: StockholmRecord | null = null;
  let rowOf = new Map<string, number>();
  const bad = (where: string, m: string) => issues.push(Issue.error(FORMAT, path, `Stockholm ${where}: ${m}`));
  const append = (m: Map<string, string>, key: string, s: string, sep: string) => {
    const prev = m.get(key);
    m.set(key, prev === undefined ? s : prev + sep + s);
  };
  const nested = (m: Map<string, Map<string, string>>, name: string) => {
    const inner = m.get(name) ?? new Map<string, string>();
    m.set(name, inner);
    return inner;
  };
  /** Per-column tags must span the columns; others are dropped with an error. */
  const checkColumns = (tags: Map<string, string>, what: string, n: number, where: string) => {
    for (const [f, v] of tags) {
      if (v.length === n) continue;
      tags.delete(f);
      bad(where, `${what} ${f} has ${v.length} columns, not ${n}.`);
    }
  };
  const close = (r: StockholmRecord) => {
    const where = `record ${records.length}`;
    const n = r.rows[0]?.length ?? 0;
    if (!r.rows.length) bad(where, "it has no rows.");
    checkColumns(r.gc, "#=GC", n, where);
    for (const [name, tags] of r.gr) {
      if (!rowOf.has(name)) bad(where, `#=GR for "${name}", which has no row.`);
      checkColumns(tags, `#=GR ${name}`, n, where);
    }
    records.push(r);
    rec = null;
    rowOf = new Map();
  };
  const empty = (): StockholmRecord => ({
    names: [],
    rows: [],
    gf: new Map(),
    gc: new Map(),
    gs: new Map(),
    gr: new Map(),
  });
  const all = lines(text);
  for (let k = 0; k < all.length; k++) {
    const line = all[k].trim();
    if (!line) continue;
    if (line.startsWith("# STOCKHOLM")) {
      if (rec) close(rec);
      rec = empty();
      continue;
    }
    if (line === "//") {
      if (rec) close(rec);
      else bad(`line ${k + 1}`, "`//` with no record open.");
      continue;
    }
    rec ??= empty();
    const tag = /^#=(GF|GC|GS|GR)\s+(\S+)\s*(.*)$/.exec(line);
    if (tag) {
      const [, kind, name, rest] = tag;
      if (kind === "GF") append(rec.gf, name, rest, "\n");
      else if (kind === "GC") append(rec.gc, name, rest.replace(/\s+/g, ""), "");
      else {
        const [, feature, value] = /^(\S*)\s*(.*)$/.exec(rest) as RegExpExecArray;
        const gr = kind === "GR";
        append(nested(gr ? rec.gr : rec.gs, name), feature, gr ? value.replace(/\s+/g, "") : value, gr ? "" : "\n");
      }
      continue;
    }
    if (line.startsWith("#")) continue;
    const m = /^(\S+)\s+(\S+)$/.exec(line);
    if (!m) {
      bad(`line ${k + 1}`, "expected a row name and its residues.");
      continue;
    }
    const at = rowOf.get(m[1]);
    if (at === undefined) {
      rowOf.set(m[1], rec.rows.length);
      rec.names.push(m[1]);
      rec.rows.push(m[2]);
    } else rec.rows[at] += m[2];
  }
  if (rec) {
    issues.push(Issue.info(UNTERMINATED, path, "Stockholm: the last record does not end with `//`."));
    close(rec);
  }
  if (!records.length && !issues.length) issues.push(Issue.error(FORMAT, path, "Stockholm: no records."));
  return { records, issues };
}

/** Insert columns from `#=GC RF`: 1 where it marks the column with `.` or `-` (seqspec:model#alignments-and-trees). */
function insertColumns(rf: string): Uint8Array {
  return Uint8Array.from(rf, (c) => (c === "." || c === "-" ? 1 : 0));
}

/** A Stockholm record as an alignment's input: aligned, labels the row names, insert columns from `RF`. */
function stockholmAlignment(record: StockholmRecord, key: string): AlignmentInput {
  const rf = record.gc.get("RF");
  const input: AlignmentInput = { key, flavor: "aligned", gapped: record.rows, labels: record.names };
  if (rf !== undefined) input.insertColumns = insertColumns(rf);
  return input;
}

/** The header lines of the Clustal writers (CLUSTAL W, Clustal Omega, MUSCLE, PROBCONS): their first word. */
const CLUSTAL_HEADER = /^(CLUSTAL|MUSCLE|PROBCONS)\b/;
/** A conservation line: indented, `*`, `:` and `.` under the columns. */
const CONSERVATION = /^\s+[*:.][*:.\s]*$/;

/**
 * A Clustal alignment (seqspec:sources#formats-and-services): the header, then blocks of `name residues [count]`
 * lines, each block listing every row in the first block's order, its rows of one length, a conservation line and a
 * blank line after it; a row's blocks joined in order. Residues may be split by spaces; the count is not read. Any
 * other shape is an error that says what to write, and nothing is read.
 */
function parseClustal(text: string, options: ParseOptions = {}): { record: ClustalRecord | null; issues: Issue[] } {
  const issues: Issue[] = [];
  const path = options.path ?? "";
  if (!withinBound(text, options.maxBytes, issues, path)) return { record: null, issues };
  const fail = (m: string) => {
    issues.push(Issue.error(FORMAT, path, `Clustal: ${m}`));
    return { record: null, issues };
  };
  const all = lines(text);
  let k = 0;
  while (k < all.length && !all[k].trim()) k++;
  const header = (all[k] ?? "").trim();
  if (!CLUSTAL_HEADER.test(header)) {
    return fail(
      `line ${k + 1} is not a Clustal header; the file starts with "CLUSTAL" (or "MUSCLE", "PROBCONS"), e.g. "CLUSTAL W (2.1) multiple sequence alignment".`,
    );
  }
  const names: string[] = [];
  const parts: string[][] = [];
  /** The block being read: its first line, its rows so far, its width. */
  let [start, row, width] = [-1, 0, -1];
  let first = true;
  const close = () => {
    if (start < 0) return null;
    if (!first && row !== names.length) {
      return `the block at line ${start} has ${row} row(s); the first block has ${names.length}: every block lists every row.`;
    }
    [start, row, width, first] = [-1, 0, -1, false];
    return null;
  };
  for (k++; k < all.length; k++) {
    const line = all[k];
    if (!line.trim() || CONSERVATION.test(line)) {
      const bad = close();
      if (bad) return fail(bad);
      continue;
    }
    if (/^\s/.test(line)) {
      return fail(
        `line ${k + 1} starts with a space; a row starts with its name, and a conservation line holds only "*", ":" and ".".`,
      );
    }
    const words = line.trim().split(/\s+/);
    if (words.length > 2 && /^\d+$/.test(words[words.length - 1])) words.pop();
    if (words.length < 2) return fail(`line ${k + 1} has a name and no residues; a row is "name residues".`);
    const name = words[0];
    const residues = words.slice(1).join("");
    if (start < 0) start = k + 1;
    if (width >= 0 && residues.length !== width) {
      const prev = first ? names[names.length - 1] : names[row - 1];
      return fail(
        `the block at line ${start} is ragged: "${prev}" has ${width} columns there, "${name}" ${residues.length}; pad rows with "-".`,
      );
    }
    width = residues.length;
    if (first) {
      names.push(name);
      parts.push([residues]);
    } else {
      if (row >= names.length) {
        return fail(
          `the block at line ${start} has more rows than the first block (${names.length}): every block lists every row once.`,
        );
      }
      if (names[row] !== name) {
        return fail(
          `line ${k + 1} names "${name}" where the first block has "${names[row]}": every block lists the rows in the first block's order.`,
        );
      }
      parts[row].push(residues);
    }
    row++;
  }
  const bad = close();
  if (bad) return fail(bad);
  if (!names.length) return fail("no rows after the header.");
  const repeated = names.filter((n, i) => names.indexOf(n) !== i);
  if (repeated.length > 0) {
    const named = firstFew([...new Set(repeated)].map((x) => `"${x}"`));
    issues.push(Issue.warning(ID, path, `Clustal: ${repeated.length} name(s) repeat another row's: ${named}.`));
  }
  return { record: { header, names, rows: parts.map((p) => p.join("")) }, issues };
}

/** A Clustal alignment as an alignment's input: aligned, labels the row names. */
function clustalAlignment(record: ClustalRecord, key: string): AlignmentInput {
  return { key, flavor: "aligned", gapped: record.rows, labels: record.names };
}

const OPEN = "([{<";
const CLOSE = ")]}>";

/**
 * WUSS pairs: per column its partner, -1 unpaired. `()`, `<>`, `[]`, `{}` and `Aa` letter pairs pair; anything else is
 * unpaired. Unbalanced input is null, with an error.
 */
function wussPairs(ss: string, issues: Issue[], path: string): Int32Array | null {
  const partner = new Int32Array(ss.length).fill(-1);
  const stacks: number[][] = Array.from({ length: 30 }, () => []);
  for (let i = 0; i < ss.length; i++) {
    const c = ss[i];
    const code = ss.charCodeAt(i);
    let open = OPEN.indexOf(c);
    let close = CLOSE.indexOf(c);
    if (code >= 65 && code <= 90) open = 4 + code - 65;
    else if (code >= 97 && code <= 122) close = 4 + code - 97;
    if (open >= 0) stacks[open].push(i);
    else if (close >= 0) {
      const j = stacks[close].pop();
      if (j === undefined) {
        issues.push(Issue.error(FORMAT, path, `WUSS: "${c}" at column ${i} closes nothing.`));
        return null;
      }
      partner[i] = j;
      partner[j] = i;
    }
  }
  const left = stacks.findIndex((s) => s.length > 0);
  if (left >= 0) {
    issues.push(Issue.error(FORMAT, path, `WUSS: the pair opened at column ${stacks[left][0]} is not closed.`));
    return null;
  }
  return partner;
}

const NUMBER = /[-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?/y;
const UNQUOTED = /[^\s()[\]',:;]+/y;

/** One Newick tree, nodes in preorder (the root is 0). Labels are kept as written (`_` stays), quotes undone. */
function parseNewick(text: string, options: ParseOptions = {}): { tree: ParsedTree | null; issues: Issue[] } {
  const issues: Issue[] = [];
  const path = options.path ?? "";
  if (!withinBound(text, options.maxBytes, issues, path)) return { tree: null, issues };
  const parent: number[] = [-1];
  const length: number[] = [Number.NaN];
  const label: (string | null)[] = [null];
  const fail = (i: number, m: string) => {
    issues.push(Issue.error(FORMAT, path, `Newick at ${i}: ${m}`));
    return { tree: null, issues };
  };
  const node = (p: number) => {
    parent.push(p);
    length.push(Number.NaN);
    label.push(null);
    return parent.length - 1;
  };
  let cur = 0;
  let depth = 0;
  let i = 0;
  let done = false;
  while (i < text.length && !done) {
    const c = text[i];
    if (/\s/.test(c)) i++;
    else if (c === "[") {
      const end = text.indexOf("]", i);
      if (end < 0) return fail(i, "a comment is not closed.");
      i = end + 1;
    } else if (c === "(") {
      // A group opens a node that has nothing yet: not after a group, a label or a length.
      if (cur + 1 < parent.length || label[cur] !== null || !Number.isNaN(length[cur])) {
        return fail(i, "a group follows a group, a label or a length.");
      }
      cur = node(cur);
      depth++;
      i++;
    } else if (c === "," || c === ")") {
      if (depth === 0) return fail(i, `"${c}" outside parentheses.`);
      if (c === ",") cur = node(parent[cur]);
      else {
        cur = parent[cur];
        depth--;
      }
      i++;
    } else if (c === ";") {
      done = true;
      i++;
    } else if (c === ":") {
      NUMBER.lastIndex = i + 1;
      const m = NUMBER.exec(text);
      if (!m) return fail(i, "a branch length is not a number.");
      if (!Number.isNaN(length[cur])) return fail(i, "a node has two lengths.");
      length[cur] = Number(m[0]);
      i = NUMBER.lastIndex;
    } else {
      if (label[cur] !== null) return fail(i, "a node has two labels.");
      if (c === "'") {
        let s = "";
        let j = i + 1;
        for (;;) {
          const q = text.indexOf("'", j);
          if (q < 0) return fail(i, "a quoted label is not closed.");
          s += text.slice(j, q);
          if (text[q + 1] !== "'") {
            i = q + 1;
            break;
          }
          s += "'";
          j = q + 2;
        }
        label[cur] = s;
      } else {
        UNQUOTED.lastIndex = i;
        const m = UNQUOTED.exec(text);
        if (!m) return fail(i, `unexpected "${c}".`);
        label[cur] = m[0];
        i = UNQUOTED.lastIndex;
      }
    }
  }
  if (!done && parent.length === 1 && label[0] === null && Number.isNaN(length[0])) return fail(i, "there is no tree.");
  if (depth !== 0) return fail(i, "parentheses are not balanced.");
  if (
    text
      .slice(i)
      .replace(/\[[^\]]*\]/g, "")
      .trim()
  )
    return fail(i, "text after the tree's `;`.");
  if (!done) issues.push(Issue.info(UNTERMINATED, path, "Newick: the tree does not end with `;`."));
  const tree = { parent: Int32Array.from(parent), branchLength: Float64Array.from(length), label };
  return { tree, issues };
}

export const Fasta = { parse: parseFasta, alignment: fastaAlignment };
export const Stockholm = { parse: parseStockholm, alignment: stockholmAlignment, insertColumns, pairs: wussPairs };
export const Newick = { parse: parseNewick };
export const Clustal = { parse: parseClustal, alignment: clustalAlignment };
