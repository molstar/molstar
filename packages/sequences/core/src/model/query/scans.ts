/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QR2b's scans (seqspec:seqql#performance): `letters` through a 256-entry mask, `positions` with its test per position,
// `search` with one matcher of @molstar/seq-core/algorithms/search. Each reads its `:in` (default: every unit and row; for
// `search`, every unit) or, as an `and`'s direct operand with no `:in` of its own, the `and`'s selection so far
//. A unit's letters are its molecule's; a row's its own residues, inserts included; the
// columns have none. Letters and positions give one group per run, `search` one per hit, scored by its edits.

import {
  type ApproximateEnds,
  approximateMatcher,
  classMatcher,
  compileMotif,
  exactMatcher,
  type Matcher,
  regexMatcher,
  toCodes,
} from "@molstar/seq-core/algorithms/search";
import { CoordinateSystemKind as K } from "../coordinates.js";
import type { SeqSet } from "../set.js";
import type { CompileContext, QueryFn } from "./compile.js";
import type { QueryContext } from "./context.js";
import { eachRun, type SelectionBuilder, sameSystem } from "./kernels.js";
import { rt, unionOf } from "./runtime.js";
import type { SeqSelection } from "./selection.js";
import { domainOf } from "./systems.js";

const dom = new Int32Array(2);
const count = (n: number, what: string) => `${n.toLocaleString("en-US")} ${what}${n === 1 ? "" : "s"}`;

/** Where a scan reads, one canonical group: its `:in`; else the `and`'s selection; else every unit (and held row). */
function scope(ctx: QueryContext, within: QueryFn<SeqSelection> | undefined, rows: boolean): SeqSelection {
  const s = ctx.scope;
  ctx.scope = null;
  if (within) return unionOf(ctx, within(ctx));
  if (s) return s;
  const set = ctx.set;
  const b = ctx.builder();
  b.begin();
  for (let u = 0; u < set.units.count; u++) if (domainOf(set, K.unit, u, -1, dom)) b.add(K.unit, u, -1, dom[0], dom[1]);
  for (let a = 0; rows && a < set.alignments.length; a++) {
    for (let r = 0; r < (set.alignments[a]?.rowCount ?? 0); r++)
      if (domainOf(set, K.row, a, r, dom)) b.add(K.row, a, r, dom[0], dom[1]);
  }
  b.end();
  return ctx.finish(b);
}

/** The letters run k's system reads: `code[off + index]`. */
const text = { code: new Uint8Array(0) as Uint8Array, off: 0, cyclic: false };
/** Run k's letters into `text`: 1 held; 0 none (the columns); -1 not loaded (a lazy molecule, a row not held). */
function lettersOf(set: SeqSet, u: SeqSelection, k: number): number {
  const kind = u.kind[k];
  if (kind === K.unit || kind === K.molecule) {
    const m = set.molecules[kind === K.unit ? set.units.molecule[u.id[k]] : u.id[k]];
    if (!m || m.lazy) return -1;
    text.code = m.code;
    text.off = 0;
    text.cyclic = m.cyclic;
    return 1;
  }
  if (kind !== K.row) return 0;
  const rows = set.alignments[u.id[k]]?.rows;
  const s = rows ? rows.slotOf(u.row[k] as never) : -1;
  if (!rows || s < 0) return -1;
  text.code = rows.residues;
  text.off = rows.residueStart[s];
  text.cyclic = false;
  return 1;
}

/**
 * The positions a scan of `u` visits, counted against the limit once, before it starts (so the message names them
 * all); with `letters`, the systems without letters are skipped, each kind with one warning.
 */
function visits(ctx: QueryContext, u: SeqSelection, path: string, letters: boolean): number {
  let [n, seqs, rows, cols, lazy, held] = [0, 0, 0, 0, 0, 1];
  for (let k = 0; k < u.runCount; k++) {
    if (k === 0 || !sameSystem(u, k - 1, k)) {
      held = letters ? lettersOf(ctx.set, u, k) : 1;
      if (held < 0) lazy++;
      else if (held === 0) cols++;
      else if (u.kind[k] === K.row) rows++;
      else seqs++;
    }
    if (held > 0) n += u.stop[k] - u.start[k];
  }
  if (letters && cols > 0) {
    const say = "alignment columns have no letters and were skipped; name the rows: (and (rows :of A) (letters C))";
    ctx.warn("seqql.letters.columns", path, say);
  }
  if (lazy > 0)
    ctx.warn("seqql.lazy", path, `the letters of ${count(lazy, "sequence")} are not loaded and were not read`);
  const what = [seqs && count(seqs, "sequence"), rows && count(rows, "row"), cols && count(cols, "column range")];
  ctx.visit(n, path, what.filter((x) => x).join(" and ") || "nothing");
  return n;
}

/** The runs of a mask's letters in [s, t) of `text`, appended to `b`'s open group. */
function scanMask(b: SelectionBuilder, mask: Uint8Array, kind: number, id: number, row: number, s: number, t: number) {
  const code = text.code;
  const off = text.off;
  for (let i = s; i < t; ) {
    while (i < t && mask[code[off + i]] === 0) i++;
    const start = i;
    while (i < t && mask[code[off + i]] !== 0) i++;
    if (i > start) b.add(kind, id, row, start, i);
  }
}

/** `(letters L :in Q)`: positions whose letter is one of L, case-insensitive; one group per run. */
function letters(c: CompileContext): QueryFn<SeqSelection> {
  const [L, within] = [String(c.raw("0") ?? ""), c.query("in")];
  const mask = new Uint8Array(256);
  for (const ch of L.toUpperCase()) if (ch.charCodeAt(0) < 256) mask[ch.charCodeAt(0)] = 1;
  return (ctx) => {
    const u = scope(ctx, within, true);
    const n = visits(ctx, u, c.path, true);
    const set = ctx.set;
    const b = ctx.builder();
    b.begin();
    for (let k = 0, held = 0; k < u.runCount; k++) {
      if (k === 0 || !sameSystem(u, k - 1, k)) {
        ctx.tick();
        held = lettersOf(set, u, k);
      }
      if (held > 0) scanMask(b, mask, u.kind[k], u.id[k], u.row[k], u.start[k], u.stop[k]);
    }
    b.end();
    const out = ctx.finish(b);
    if (out.runCount > ctx.limits.groups) ctx.over("groups", out.runCount);
    if (out.runCount === 0) ctx.because(c.path, `no ${[...L.toUpperCase()].join(", ")} in ${count(n, "position")}`);
    return eachRun(out);
  };
}

/** Run k's positions passing `where` (the slot's index set per position), appended to `b`'s open group. */
function testRun(ctx: QueryContext, b: SelectionBuilder, where: QueryFn<unknown>, u: SeqSelection, k: number) {
  const slot = ctx.slot;
  const t = u.stop[k];
  let open = -1;
  for (let x = u.start[k]; x < t; x++) {
    slot.index = x;
    if (where(ctx) === true) {
      if (open < 0) open = x;
    } else if (open >= 0) {
      b.add(u.kind[k], u.id[k], u.row[k], open, x);
      open = -1;
    }
  }
  if (open >= 0) b.add(u.kind[k], u.id[k], u.row[k], open, t);
}

/** `(positions :in Q :where T)`: the general case, a test per position over the reused slot; one group per run. */
function positions(c: CompileContext): QueryFn<SeqSelection> {
  const [within, where] = [c.query("in"), c.value("where")];
  return (ctx) => {
    const u = scope(ctx, within, true);
    const n = visits(ctx, u, c.path, false);
    const [b, slot, missed] = [ctx.builder(), ctx.slot, ctx.missed];
    b.begin();
    for (let k = 0; k < u.runCount; k++) {
      if (k === 0 || !sameSystem(u, k - 1, k)) {
        ctx.tick();
        slot.kind = u.kind[k];
        slot.id = u.id[k];
        slot.row = u.row[k];
      }
      if (where) testRun(ctx, b, where, u, k);
      else b.add(u.kind[k], u.id[k], u.row[k], u.start[k], u.stop[k]);
    }
    b.end();
    const out = ctx.finish(b);
    if (out.runCount > ctx.limits.groups) ctx.over("groups", out.runCount);
    if (out.runCount === 0) ctx.tested(c.path, n, "position", missed);
    return eachRun(out);
  };
}

/** The hit visitor's target: the builder, the system searched, and its letters' offset in the searched array. */
const H = { b: null as SelectionBuilder | null, kind: 0, id: 0, row: 0, off: 0 };
function hit(start: number, stop: number, score: number): void {
  const b = H.b as SelectionBuilder;
  b.begin(-1, -1, score);
  // A hit across the origin of a cyclic molecule (start > stop): its two runs, one group.
  if (start < stop) b.add(H.kind, H.id, H.row, start - H.off, stop - H.off);
  else {
    b.add(H.kind, H.id, H.row, 0, stop);
    b.add(H.kind, H.id, H.row, start, text.code.length);
  }
  b.end();
}

interface Pattern {
  readonly text: string;
  readonly edits: number;
  readonly classes: string | undefined;
  readonly regex: boolean;
  readonly ends: ApproximateEnds;
}
/** `:regex` a regular expression; `:classes` a motif (with `:edits`, its classes as equal pairs); `:edits`
 * the approximate matcher; else the exact one. Letters are upper-cased, as the codes are. */
function matcherOf(ctx: QueryContext, c: CompileContext, p: Pattern): Matcher {
  const table = p.classes === undefined ? null : c.env.classes.get(p.classes);
  if (p.classes !== undefined && !table) {
    const have = [...c.env.classes.keys()];
    const say = have.length ? `its tables: ${have.join(", ")}` : "this host registers none";
    ctx.fail("seqql.classes.unknown", `${c.path}/args/classes`, `no letter-class table "${p.classes}" here; ${say}`);
  }
  const upper = p.text.toUpperCase();
  try {
    if (p.regex) return regexMatcher(new RegExp(p.text, "i"));
    if (!table) return p.edits > 0 ? approximateMatcher(toCodes(upper), p.edits, p) : exactMatcher(toCodes(upper));
    if (p.edits === 0) return classMatcher(compileMotif(upper, table));
    const equal: number[] = [];
    for (const [k, v] of Object.entries(table)) for (const x of toCodes(v)) equal.push(k.charCodeAt(0), x);
    return approximateMatcher(toCodes(upper), p.edits, { ends: p.ends, equal: Uint8Array.from(equal) });
  } catch (e) {
    if (!(e instanceof RangeError || e instanceof SyntaxError)) throw e;
    const why = e.message.replace(/^\w+: /, "");
    return ctx.fail("seqql.type", `${c.path}/args/0`, `search cannot use the pattern "${p.text}": ${why}`);
  }
}

/** `(search P :in Q :edits n :classes C :regex b :ends e)`: one group per hit, in the order of their ends per run. */
function search(c: CompileContext): QueryFn<SeqSelection> {
  const within = c.query("in");
  const p: Pattern = {
    text: String(c.raw("0") ?? ""),
    edits: c.raw<number>("edits") ?? 0,
    classes: c.raw<string>("classes"),
    regex: c.raw<boolean>("regex") === true,
    ends: c.raw<ApproximateEnds>("ends") ?? "local_best_last",
  };
  let m: Matcher | null = null;
  return (ctx) => {
    m ??= matcherOf(ctx, c, p);
    const u = scope(ctx, within, false);
    const n = visits(ctx, u, c.path, true);
    const b = ctx.builder();
    H.b = b;
    for (let k = 0, held = 0; k < u.runCount; k++) {
      if (k === 0 || !sameSystem(u, k - 1, k)) {
        ctx.tick();
        held = lettersOf(ctx.set, u, k);
        H.kind = u.kind[k];
        H.id = u.id[k];
        H.row = u.row[k];
        H.off = text.off;
      }
      // A cyclic molecule read whole is searched as a circle: a hit may cross its origin.
      const circle = text.cyclic && u.start[k] === 0 && u.stop[k] === text.code.length;
      if (held > 0) m.find(text.code, text.off + u.start[k], text.off + u.stop[k], circle, hit);
    }
    H.b = null;
    const out = ctx.finish(b);
    const more = p.edits > 0 || p.regex ? "" : "; allow edits with :edits 1";
    if (out.groupCount === 0) ctx.because(c.path, `no hit in ${count(n, "position")}${more}`);
    return out;
  };
}

export const SCAN_RUNTIMES = [
  rt("seq.generator.letters", letters),
  rt("seq.generator.positions", positions),
  rt("seq.generator.search", search),
];
