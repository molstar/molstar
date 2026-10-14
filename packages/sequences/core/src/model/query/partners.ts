/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `partners`: one side of a pair by its partner. Q is unioned once; its runs on units ask the pair index which pairs
// touch each position (`Pairs.touching`, allocation-free), each end visited once (an epoch mark per row and side); ends
// on a molecule with no unit or on an alignment are found by a scan of the pair rows. For every pair passing `:where`
// (its item) with an end in Q, the other end is returned, when `:side` allows it.

import { type Annotation, Pairs } from "../annotation.js";
import type { SeqIndex, UnitId } from "../common.js";
import { CoordinateSystemKind as K } from "../coordinates.js";
import type { CompileContext, QueryFn } from "./compile.js";
import type { QueryContext } from "./context.js";
import { eachRun, firstAfter, type SelectionBuilder } from "./kernels.js";
import { addTarget, annotationOf, rt, unionOf, Wants, wrongKind } from "./runtime.js";
import type { SeqSelection } from "./selection.js";
import { unitsOf } from "./systems.js";

const WHOLE = 0x7fffffff;
/** The run in progress: per (row, side) the epoch it was touched in; per row the epoch its test ran in, and its result. */
const P = {
  ctx: null as QueryContext | null,
  b: null as SelectionBuilder | null,
  a: null as Annotation | null,
  where: undefined as QueryFn<unknown> | undefined,
  /** Ends that may be returned: 1 a, 2 b, 3 both. */
  sides: 3,
  epoch: 0,
  /** Pair ends Q touched in this run (seqspec:seqql#empty-results). */
  reached: 0,
  touched: new Int32Array(64),
  tested: new Int32Array(64),
  passed: new Uint8Array(64),
};

/** Row `r`'s end on `side` touches Q: its other end is returned, once, when the pair passes and `:side` allows it. */
function touch(r: number, side: number, p: typeof P): void {
  const at = 2 * r + side;
  if (p.touched[at] === p.epoch) return;
  p.touched[at] = p.epoch;
  p.reached++;
  const other = 1 - side;
  if ((p.sides & (other ? 2 : 1)) === 0) return;
  const a = p.a as Annotation;
  const ctx = p.ctx as QueryContext;
  if (p.where) {
    if (p.tested[r] !== p.epoch) {
      p.tested[r] = p.epoch;
      ctx.slot.annotation = a.id;
      ctx.slot.item = a.item[r];
      p.passed[r] = p.where(ctx) === true ? 1 : 0;
    }
    if (!p.passed[r]) return;
  }
  const b = p.b as SelectionBuilder;
  const pc = a.pairs;
  if (other === 0) {
    // Side a always has a position: one without is `model.position.missing` at build (only side b may be whole).
    addTarget(b, ctx.set, a.slots, a.slot[r], a.start[r], Math.max(a.stop[r], a.start[r] + 1));
  } else if (pc) {
    const s = pc.bStart[r];
    addTarget(b, ctx.set, pc.bTarget, r, s < 0 ? 0 : s, s < 0 ? WHOLE : Math.max(pc.bStop[r], s + 1));
  }
}

/** An end not in the pair index (a molecule with no unit, a row, the columns) against Q's union, by its target. */
function scanEnd(u: SeqSelection, a: Annotation, r: number, side: number): void {
  const set = (P.ctx as QueryContext).set;
  const tg = side ? (a.pairs as NonNullable<Annotation["pairs"]>).bTarget : a.slots;
  const i = side ? r : a.slot[r];
  const mol = tg.molecule[i];
  const units = unitsOf(set);
  if (tg.unit[i] >= 0 || (mol >= 0 && units.start[mol + 1] > units.start[mol])) return;
  const kind = mol >= 0 ? K.molecule : tg.row[i] >= 0 ? K.row : K.columns;
  const id = mol >= 0 ? mol : tg.alignment[i];
  const row = kind === K.row ? tg.row[i] : -1;
  const s0 = side ? (a.pairs as NonNullable<Annotation["pairs"]>).bStart[r] : a.start[r];
  const s = s0 < 0 ? 0 : s0;
  const e = s0 < 0 ? WHOLE : Math.max(side ? (a.pairs as NonNullable<Annotation["pairs"]>).bStop[r] : a.stop[r], s + 1);
  const j = firstAfter(u, kind, id, row, s);
  if (j < u.runCount && u.kind[j] === kind && u.id[j] === id && u.row[j] === row && u.start[j] < e) touch(r, side, P);
}

const SIDES: Record<string, number> = { any: 3, a: 1, b: 2 };

/** `(partners Q :via A :side S :where T)`: the ends whose other end is in Q, one group per run, untagged. */
function partners(c: CompileContext): QueryFn<SeqSelection> {
  const [q, where] = [c.query("0") as QueryFn<SeqSelection>, c.value("where")];
  const [via, side] = [c.raw<string>("via") as string, SIDES[c.raw<string>("side") ?? "any"]];
  const path = `${c.path}/args/via`;
  return (ctx) => {
    const a = ctx.set.annotations[annotationOf(ctx, via, path, Wants.annotation)] as Annotation;
    if (a.geometry !== "pairs") wrongKind(ctx, via, path, "", `${a.geometry} annotation`);
    const u = unionOf(ctx, q(ctx));
    const b = ctx.builder();
    if (P.touched.length < 2 * a.count) {
      P.touched = new Int32Array(4 * a.count);
      P.tested = new Int32Array(2 * a.count);
      P.passed = new Uint8Array(2 * a.count);
    }
    // Saved and restored: a :where may hold another partners (a selection inside a test).
    const saved = { ctx: P.ctx, b: P.b, a: P.a, where: P.where, sides: P.sides, reached: P.reached };
    Object.assign(P, { ctx, b, a, where, sides: side, epoch: P.epoch + 1, reached: 0 });
    b.begin();
    let other = false;
    for (let k = 0; k < u.runCount; k++) {
      if (u.kind[k] !== K.unit) {
        other = true;
        continue;
      }
      ctx.visit(u.stop[k] - u.start[k], path, via);
      ctx.tick();
      for (let x = u.start[k]; x < u.stop[k]; x++) Pairs.touching(a, u.id[k] as UnitId, x as SeqIndex, touch, P);
    }
    for (let r = 0; other && r < a.count; r++) {
      scanEnd(u, a, r, 0);
      scanEnd(u, a, r, 1);
    }
    b.end();
    const reached = P.reached;
    Object.assign(P, saved);
    const out = ctx.finish(b);
    if (out.runCount > ctx.limits.groups) ctx.over("groups", out.runCount);
    if (out.runCount === 0 && u.runCount > 0) {
      const ends = `${reached} pair end${reached === 1 ? ` of ${via} touches` : `s of ${via} touch`} it`;
      ctx.because(c.path, reached ? `${ends}; none passed :where or :side` : `no pair of ${via} touches it`);
    }
    return eachRun(out);
  };
}

export const PARTNERS_RUNTIMES = [rt("seq.modifier.partners", partners)];
