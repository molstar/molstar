/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's formats (seqspec:sources#data-by-reference): FASTA, aligned FASTA and A3M, Stockholm with its tags,
// Clustal, Newick, TSV, CSV, JSON and MessagePack tables; and the `uniprot` service. Read from bytes, or from a
// source's `value` encoded as UTF-8.

import {
  Clustal,
  columnNumbers,
  decodeMsgPack,
  decodeTypedArrays,
  Fasta,
  type Geometry,
  type Issue,
  Newick,
  Stockholm,
  type StockholmRecord,
} from "@molstar/seq-core/model";
import type { FormatResult, StageContext } from "../contract.js";
import type { FormatEntry, ServiceEntry } from "../registry.js";
import { isList } from "../tables/table.js";

const decoder = new TextDecoder();
const text = (bytes: Uint8Array) => decoder.decode(bytes);
type Cells = Record<string, unknown[]>;
/** An object keyed by data (column names): no prototype, so `constructor` or `__proto__` is a key like any other. */
const cells = (): Cells => Object.create(null);
const table = (columns: Record<string, ArrayLike<unknown>>, length: number, text = false): FormatResult =>
  text ? { kind: "table", columns, length, text } : { kind: "table", columns, length };
const fail = (ctx: StageContext, m: string): FormatResult => {
  ctx.report({ code: "svs.source.format", severity: "error", path: "", message: m });
  return table({}, 0);
};

/**
 * The last FASTA file's parse, by its bytes' identity: a FASTA alignment is read for its source's default (molecules)
 * and for the alignment, two requests, one text (GV-02: parsed twice), and those two reads meet one after the other.
 * One slot, not a cache by the bytes: the source's `read` closure holds the bytes for as long as the source is open,
 * so anything reachable from them lived as long and every open FASTA source kept its parse (LA-38). A new file's parse
 * evicts the old one; the memory is one parse and the last file's bytes, and a read of an earlier file parses again.
 * Its issues are reported on every read, and `open` says each once.
 */
let fastaSlot: { bytes: Uint8Array; parsed: ReturnType<typeof Fasta.parse> } | null = null;

function fasta(flavor: "aligned" | "a3m" | null): FormatEntry["read"] {
  return (bytes, request, ctx) => {
    if (fastaSlot?.bytes !== bytes) fastaSlot = { bytes, parsed: Fasta.parse(text(bytes)) };
    const { records, issues } = fastaSlot.parsed;
    for (const i of issues) ctx.report(i);
    if (request.tag === "header") {
      return table(
        { row: records.map((_, k) => k), id: records.map((r) => r.id), header: records.map((r) => r.description) },
        records.length,
      );
    }
    if (flavor || request.want === "alignments") {
      const a = Fasta.alignment(records, "", flavor ?? request.flavor);
      for (const i of a.issues) ctx.report(i);
      return { kind: "alignments", alignments: [a.input] };
    }
    const molecules = records.map((r) => ({ key: r.id, residues: r.sequence, label: r.id }));
    return { kind: "molecules", molecules };
  };
}

const isGap = (code: number) => code === 45 || code === 46;

/** A pass's columns; `value` has items for the values geometry, `start` and `stop` for the others, `label` for ranges. */
interface Out {
  start: Int32Array;
  stop: Int32Array;
  value: Float64Array;
  label: string[];
  target: string[] | null;
}

/**
 * A Stockholm tag (`GC:SS_cons`, `GR:SS`, `GS:OS`) as a table in the reading annotation's geometry: the items counted
 * first, then the typed columns filled. Anything it cannot read is `svs.source.format`.
 */
function tagTable(rec: StockholmRecord, tag: string, geometry: Geometry | null, ctx: StageContext): FormatResult {
  const [kind, feature] = [tag.slice(0, 2), tag.slice(3)];
  if (kind === "GS") {
    if (geometry != null && geometry !== "rows") return fail(ctx, `${tag} is read as rows, not ${geometry}`);
    const row: number[] = [];
    const values: string[] = [];
    rec.names.forEach((name, k) => {
      const v = rec.gs.get(name)?.get(feature);
      if (v === undefined) return;
      row.push(k);
      values.push(v);
    });
    return table({ row: Int32Array.from(row), [feature]: values }, row.length);
  }
  if (kind !== "GC" && kind !== "GR") return fail(ctx, `"${tag}" is not GC:<tag>, GR:<tag> or GS:<tag>`);
  const geo = geometry ?? "ranges";
  if (geo !== "ranges" && geo !== "values" && geo !== "pairs") {
    return fail(ctx, `${tag} is read as ranges, values or pairs, not ${geo}`);
  }
  /** The tag's strings: the columns' one, or each row's, its `target` one string. */
  const strings: { s: string; row: number; target: string | null }[] = [];
  if (kind === "GC") {
    const s = rec.gc.get(feature);
    if (s === undefined) return fail(ctx, `the alignment has no #=GC ${feature}`);
    strings.push({ s, row: -1, target: null });
  } else {
    const rowOf = new Map(rec.names.map((name, row) => [name, row]));
    for (const [name, tags] of rec.gr) {
      const [s, row] = [tags.get(feature), rowOf.get(name) ?? -1];
      if (s !== undefined && row >= 0) strings.push({ s, row, target: String(row) });
    }
    if (!strings.length) return fail(ctx, `no row has #=GR ${feature}`);
  }
  const pos = new Int32Array(strings.reduce((w, x) => Math.max(w, x.s.length), 0));
  const bad: Issue[] = [];
  let numbers = 0;
  /** One pass over the items: their count (-1: unbalanced WUSS, reported in `bad`), and with `out` their columns. */
  const pass = (out: Out | null): number => {
    let n = 0;
    for (const x of strings) {
      const [n0, gapped] = [n, x.row < 0 ? null : rec.rows[x.row]];
      // A column's index among the kept ones: all of them, or a row's residues with its gap columns dropped.
      for (let c = 0, k = 0; c < x.s.length; c++) pos[c] = gapped && isGap(gapped.charCodeAt(c)) ? -1 : k++;
      if (geo === "pairs") {
        const partner = Stockholm.pairs(x.s, bad, "");
        if (!partner) return -1;
        for (let c = 0; c < x.s.length; c++) {
          const d = partner[c];
          if (d <= c || pos[c] < 0 || pos[d] < 0) continue;
          if (out) {
            out.start[n] = pos[c];
            out.stop[n] = pos[d];
          }
          n++;
        }
      } else if (geo === "values") {
        for (let c = 0; c < x.s.length; c++) {
          if (pos[c] < 0) continue;
          const ch = x.s.charCodeAt(c);
          const v = ch >= 48 && ch <= 57 ? ch - 48 : ch === 42 ? 10 : Number.NaN;
          if (!Number.isNaN(v)) numbers++;
          if (out) out.value[n] = v;
          n++;
        }
      } else {
        for (let c = 0; c < x.s.length; ) {
          let e = c + 1;
          while (e < x.s.length && x.s.charCodeAt(e) === x.s.charCodeAt(c)) e++;
          let first = -1;
          let last = -1;
          for (let i = c; i < e; i++) {
            if (pos[i] < 0) continue;
            if (first < 0) first = pos[i];
            last = pos[i];
          }
          if (!isGap(x.s.charCodeAt(c)) && first >= 0) {
            if (out) {
              out.start[n] = first;
              out.stop[n] = last + 1;
              out.label[n] = x.s[c];
            }
            n++;
          }
          c = e;
        }
      }
      // A row's one string for all its items.
      if (out?.target && x.target !== null) out.target.fill(x.target, n0, n);
    }
    return n;
  };
  const n = pass(null);
  if (n < 0) {
    for (const i of bad) ctx.report({ ...i, message: `${tag}: ${i.message}` });
    return table({}, 0);
  }
  if (geo === "values" && numbers === 0) return fail(ctx, `${tag} has no digit or * to read as a value`);
  const out: Out = {
    start: new Int32Array(geo === "values" ? 0 : n),
    stop: new Int32Array(geo === "values" ? 0 : n),
    value: new Float64Array(geo === "values" ? n : 0),
    label: new Array<string>(geo === "ranges" ? n : 0),
    target: kind === "GR" ? new Array<string>(n) : null,
  };
  pass(out);
  const cols: Record<string, ArrayLike<unknown>> = out.target ? { target: out.target } : {};
  if (geo === "values") return table({ ...cols, value: out.value }, n);
  if (geo === "pairs") {
    return table({ ...cols, start: out.start, b_start: out.stop, category: Array(n).fill("base_pair") }, n);
  }
  return table({ ...cols, start: out.start, stop: out.stop, label: out.label }, n);
}

const stockholm: FormatEntry["read"] = (bytes, request, ctx) => {
  const { records, issues } = Stockholm.parse(text(bytes));
  for (const i of issues) ctx.report(i);
  if (request.tag == null) return { kind: "alignments", alignments: records.map((r) => Stockholm.alignment(r, "")) };
  const r = records[typeof request.record === "number" ? request.record : 0];
  return r ? tagTable(r, request.tag, request.geometry, ctx) : fail(ctx, `no record ${request.record ?? 0}`);
};

/** Clustal (`.aln`): one alignment; a file it cannot read reported once, nothing read. */
const clustal: FormatEntry["read"] = (bytes, _, ctx) => {
  const { record, issues } = Clustal.parse(text(bytes));
  for (const i of issues) ctx.report(i);
  return record ? { kind: "alignments", alignments: [Clustal.alignment(record, "")] } : table({}, 0);
};

const newick: FormatEntry["read"] = (bytes, _, ctx) => {
  const { tree, issues } = Newick.parse(text(bytes));
  for (const i of issues) ctx.report(i);
  // No tree: the parse said why, once.
  return tree ? { kind: "tree", tree } : issues.length ? table({}, 0) : fail(ctx, "no tree");
};

/**
 * Rows of a delimited file, each passed to `onRow` as the first `count` cells of one reused array (false stops): quotes as RFC 4180 (a doubled quote
 * is one; a cell is quoted when its first character other than a space or tab is `"`, and what follows the closing quote
 * joins it), unquoted cells trimmed, blank lines skipped. A cell is a slice of the text, not built a character at a time,
 * and no array is made per row, so a large file costs its cells.
 */
function delimitedRows(s: string, sep: string, onRow: (cells: readonly string[], count: number) => boolean): void {
  const [S, n] = [sep.charCodeAt(0), s.length];
  // The row's cells are its first `count`: the array is never shortened, which would drop its store each row.
  const row: string[] = [];
  let count = 0;
  /** The end of the cell's text from `i`: the separator, a line end or the end of the text. */
  const until = (i: number) => {
    for (let c = s.charCodeAt(i); i < n && c !== S && c !== 10 && c !== 13; c = s.charCodeAt(++i));
    return i;
  };
  let i = 0;
  for (;;) {
    let j = i;
    for (let c = s.charCodeAt(j); j < n && c !== S && (c === 32 || c === 9); c = s.charCodeAt(++j));
    let cell: string;
    if (j < n && s.charCodeAt(j) === 34) {
      cell = "";
      for (let from = j + 1; ; ) {
        const q = s.indexOf('"', from);
        cell += s.slice(from, q < 0 ? n : q);
        i = q < 0 ? n : q + 1;
        if (q < 0 || s.charCodeAt(q + 1) !== 34) break;
        cell += '"';
        from = q + 2;
      }
      const e = until(i);
      if (e > i) cell += s.slice(i, e);
      i = e;
    } else {
      const e = until(i);
      cell = s.slice(i, e).trim();
      i = e;
    }
    row[count++] = cell;
    const c = i < n ? s.charCodeAt(i) : -1;
    if (c === S) {
      i++;
      continue;
    }
    if ((count > 1 || row[0] !== "") && !onRow(row, count)) return;
    if (c < 0) return;
    count = 0;
    i += c === 13 && s.charCodeAt(i + 1) === 10 ? 2 : 1;
  }
}

/**
 * A table by its header, its cells text (empty ones null), typed after `fields` (seqspec:format#annotation-tables); a
 * repeated name is an error.
 */
function delimited(sep: string): FormatEntry["read"] {
  return (bytes, _, ctx) => {
    const s = text(bytes);
    // Columns are sized once, by an upper bound on the rows, and cut to them after: grown by a push a cell, each
    // column's store was copied at every growth (23 bytes a cell of garbage, load.heap.ts). The bound is the lines, and
    // at most the text's length over the header's cells, since a row of `k` cells takes at least `k` characters (its
    // `k - 1` separators and its line end, the last row's aside): the slots are then at most the text's length, whatever
    // the header (GI-04: sized by the lines alone, 20,000 names over 200,000 empty lines asked for 4 billion).
    let lines = 1;
    for (let at = s.indexOf("\n"); at >= 0; at = s.indexOf("\n", at + 1)) lines++;
    let names: string[] | null = null;
    let cols: (string | null)[][] = [];
    let [length, error] = [0, ""];
    delimitedRows(s, sep, (cells, count) => {
      if (!names) {
        names = cells.slice(0, count).map((name) => name.trim());
        const rows = Math.min(lines, Math.floor(s.length / names.length) + 1);
        cols = names.map(() => new Array(rows));
        const seen = new Set<string>();
        const twice = names.find((name) => seen.size === seen.add(name).size);
        if (twice !== undefined) error = `the header names "${twice}" twice`;
      } else if (count !== names.length) {
        error = `row ${length + 2} has ${count} cells; the header has ${names.length}`;
      } else {
        for (let k = 0; k < count; k++) cols[k][length] = cells[k] === "" ? null : cells[k];
        length++;
      }
      return !error;
    });
    if (!names) return fail(ctx, "no header row");
    if (error) return fail(ctx, error);
    const out = cells();
    (names as string[]).forEach((name, k) => {
      cols[k].length = length;
      out[name] = cols[k];
    });
    return table(out, length, true);
  };
}

const isObject = (v: unknown): v is Record<string, unknown> => v !== null && typeof v === "object" && !Array.isArray(v);

/** A `Table` (seqspec:format#positions) as columns, dictionaries decoded; null when the value is neither form. */
export function tableColumns(v: unknown): { columns: Cells; length: number } | null {
  if (Array.isArray(v)) {
    if (!v.every(isObject)) return null;
    const columns = cells();
    for (const row of v)
      for (const k of Object.keys(row)) columns[k] ??= v.map((r) => (Object.hasOwn(r, k) ? r[k] : null) ?? null);
    return { columns, length: v.length };
  }
  if (!isObject(v)) return null;
  const columns = cells();
  for (const [k, c] of Object.entries(v)) {
    if (Array.isArray(c)) columns[k] = Array.from(c);
    else if (ArrayBuffer.isView(c)) columns[k] = columnNumbers(c as Float64Array, null);
    else if (isObject(c) && Array.isArray(c.json)) columns[k] = c.json;
    else if (isObject(c) && Array.isArray(c.dict) && isList(c.codes)) {
      const dict = c.dict;
      columns[k] = Array.from(c.codes, (i) => (i == null ? null : (dict[i as number] ?? null)));
    } else return null;
  }
  return { columns, length: Object.values(columns)[0]?.length ?? 0 };
}

/** A `Table` in either form, else an object, from a parse of the bytes (JSON or MessagePack, seqspec:format#values). */
function tableOrObject(parse: (bytes: Uint8Array) => unknown, what: string): FormatEntry["read"] {
  return (bytes, _, ctx) => {
    let v: unknown;
    try {
      v = parse(bytes);
    } catch (e) {
      return fail(ctx, `not ${what}: ${e instanceof Error ? e.message : e}`);
    }
    const t = tableColumns(v);
    return t ? table(t.columns, t.length) : { kind: "object", value: v };
  };
}

export const KernelFormats: readonly FormatEntry[] = [
  {
    name: "fasta",
    extensions: ["fasta", "fa", "faa", "fna"],
    gives: ["molecules", "alignments", "table"],
    reads: ["want", "flavor", "tag"],
    read: fasta(null),
  },
  { name: "a3m", extensions: ["a3m"], gives: ["alignments"], reads: ["want", "tag"], read: fasta("a3m") },
  {
    name: "stockholm",
    extensions: ["sto", "stk"],
    gives: ["alignments", "table"],
    reads: ["want", "flavor", "tag", "geometry", "record"],
    read: stockholm,
  },
  { name: "clustal", extensions: ["aln", "clustal", "clw"], gives: ["alignments"], reads: [], read: clustal },
  { name: "newick", extensions: ["nwk", "newick", "tree"], gives: ["tree"], reads: [], read: newick },
  { name: "tsv", extensions: ["tsv", "tab"], gives: ["table"], reads: [], read: delimited("\t") },
  { name: "csv", extensions: ["csv"], gives: ["table"], reads: [], read: delimited(",") },
  {
    name: "json",
    extensions: ["json"],
    gives: ["table", "object"],
    reads: [],
    read: tableOrObject((b) => JSON.parse(text(b)), "JSON"),
  },
  {
    name: "msgpack",
    extensions: ["msgpack", "mpk"],
    gives: ["table", "object"],
    reads: [],
    read: tableOrObject((b) => decodeMsgPack(b, decodeTypedArrays), "MessagePack"),
  },
];

type Feature = {
  type?: string;
  description?: string;
  featureId?: string;
  location?: { start?: { value?: number }; end?: { value?: number } };
  begin?: string;
  end?: string;
  /** The REST entry's object, or the Proteins API's letters (with `wildType`). */
  alternativeSequence?: string | { originalSequence?: string; alternativeSequences?: string[] };
  wildType?: string;
  consequenceType?: string;
};

/**
 * UniProt's features (or the Proteins API's variants) as ranges in indices, one row per located feature and per
 * alternative sequence (`alt`, "" a deletion; `original` the residues replaced), the rows of one feature sharing its `id`.
 * Where the response gives the protein's sequence (the Proteins API's a string, the entry's `{ value, length }`), a
 * feature starting past its end (a stop-lost extension, a frameshift's new tail) is left out and one running past it is
 * clipped, one `svs.source.skipped` note counting both (seqspec:sources#formats-and-services).
 */
function features(value: unknown, ctx: StageContext | undefined, what: string): FormatResult {
  const list = (isObject(value) && Array.isArray(value.features) ? value.features : []) as Feature[];
  const seq = isObject(value) ? value.sequence : undefined;
  const length =
    typeof seq === "string" ? seq.length : isObject(seq) && Number.isInteger(seq.length) ? (seq.length as number) : -1;
  const cols: Cells = { start: [], stop: [], category: [], label: [], id: [], alt: [], original: [] };
  let [past, clipped] = [0, 0];
  for (const f of list) {
    const s = Number(f.location?.start?.value ?? f.begin);
    let e = Number(f.location?.end?.value ?? f.end);
    if (!(Number.isInteger(s) && Number.isInteger(e) && s >= 1 && e >= s)) continue;
    if (length >= 0 && s > length) {
      past++;
      continue;
    }
    if (length >= 0 && e > length) {
      clipped++;
      e = length;
    }
    const a = f.alternativeSequence;
    const original = typeof a === "string" ? (f.wildType ?? null) : (a?.originalSequence ?? null);
    const given = typeof a === "string" ? [a] : (a?.alternativeSequences ?? []);
    for (const alt of given.length ? given : [original != null ? "" : null]) {
      cols.start.push(s - 1);
      cols.stop.push(e);
      cols.category.push(f.type ?? null);
      cols.label.push(f.description || f.consequenceType || null);
      cols.id.push(f.featureId ?? null);
      cols.alt.push(alt);
      cols.original.push(original);
    }
  }
  for (const k of ["id", "alt", "original"]) if (cols[k].every((v) => v === null)) delete cols[k];
  if (ctx && (past || clipped)) {
    const said = [
      past ? `${past} start past the sequence's end (${length}) and are left out` : "",
      clipped ? `${clipped} run past ${past ? "it" : `the sequence's end (${length})`} and are clipped` : "",
    ];
    const message = `${what}: ${said.filter(Boolean).join("; ")}`;
    ctx.report({ code: "svs.source.skipped", severity: "info", path: "", message });
  }
  return table(cols, cols.start.length);
}

/** `uniprot`: an accession or isoform; a sequence it reads gets the numbering `uniprot:<id>`. */
export const UniprotService: ServiceEntry = {
  name: "uniprot",
  resources: ["sequence", "features", "variants"],
  resolve(id, resource) {
    const base = "https://rest.uniprot.org/uniprotkb";
    if (resource === "sequence") return { uri: `${base}/${id}.fasta`, format: "fasta" };
    if (resource === "features") return { uri: `${base}/${id}.json`, format: "json" };
    if (resource === "variants") return { uri: `https://www.ebi.ac.uk/proteins/api/variation/${id}`, format: "json" };
    return null;
  },
  after(result, id, resource, ctx) {
    if (result.kind === "object") return features(result.value, ctx, `uniprot ${id} ${resource}`);
    if (result.kind !== "molecules") return result;
    const molecules = result.molecules.map((m) => ({
      ...m,
      label: id,
      identifiers: [{ namespace: "uniprot", value: id }],
      numberings: [...(m.numberings ?? []), { name: `uniprot:${id}`, method: "offset" as const, first: 1 }],
    }));
    return { kind: "molecules", molecules };
  },
};
