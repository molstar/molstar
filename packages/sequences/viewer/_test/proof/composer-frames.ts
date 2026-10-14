/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's proof, for a change meant to draw exactly as before (a move, a split): the node tier run with this
// setup (`vitest run --config vitest.proof.config.ts`, `SEQ_PROOF_OUT=<folder>`) writes a line for every composer
// frame any test draws: the test, the frame's number in it, its flags, the hash of what it drew (`hashRecording` of
// the backend's recording; a tee's two), and the composer's outputs that are not drawn: its stats, `animating`,
// `again` and `animStep`, the issues it reported during the frame, the ring labels it placed for the hit
// (`placedRingLabels`), and the section's hit grid (a `hitTest` every 16 CSS px). Run it before and after the change,
// then `diff -r` the two folders: equal folders, the same frames. `SEQ_PROOF_COVER=1` also writes each
// `SectionComposer` method's call count per test file (the corpus reaches every method). A frame whose test reads
// the clock may differ between two runs of the same code: run the before twice to know them; the files known to vary
// are named in `VARIES` with why, their drawing, stats and hit grid written as `varies`.

import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, join, relative } from "node:path";
import { afterAll, beforeAll, expect } from "vitest";
import type { ComposerStats, RenderBackend } from "../../src/draw/index.js";
import type { SeqViewerContext } from "../../src/index.js";

/**
 * The viewer's modules, imported in `beforeAll`: after the test file's own `vi.mock`s are registered, so a test that
 * mocks a module the composer reads draws mocked under the proof too (LC-44).
 */
let m: {
  draw: typeof import("../../src/draw/index.js");
  rings: typeof import("../../src/draw/ring-labels.js");
  root: typeof import("../../src/index.js");
  testing: typeof import("../../src/testing/index.js");
  svg: typeof import("../../src/svg/backend.js");
};

const OUT = process.env.SEQ_PROOF_OUT;
const COVER = process.env.SEQ_PROOF_COVER === "1";
const GRID = 16;
const ROOT = join(import.meta.dirname, "../..");

/** The composer's fields the proof reads (some private). */
interface Composer {
  readonly ctx: SeqViewerContext;
  readonly section: string;
  readonly stats: ComposerStats;
  readonly animating: boolean;
  readonly again: boolean;
  readonly animStep: Float64Array;
}

const lines: string[] = [];
const calls: Record<string, number> = {};
const counts = new Map<string, number>();
let file = "";

/** A string's FNV-1a hash (8 hex digits): the outputs that are long. */
function digest(s: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < s.length; i++) h = Math.imul(h ^ s.charCodeAt(i), 0x01000193);
  return (h >>> 0).toString(16).padStart(8, "0");
}
/** An int array's FNV-1a hash, element by element. */
function digestInts(a: ArrayLike<number>): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < a.length; i++) h = Math.imul(h ^ (a[i] | 0), 0x01000193);
  return (h >>> 0).toString(16).padStart(8, "0");
}
/**
 * What a backend drew: a recording's `hashRecording`, a tee's two sides, the SVG backend's body and stats, a pixel
 * backend's pixels (`px`, as the arrivals' `Stacks`); another, its name alone (LC-38: SVG and pixel frames were
 * hashed by name).
 */
function drawn(be: RenderBackend): string {
  const { RecordingBackend, TeeBackend, hashRecording } = m.testing;
  if (be instanceof RecordingBackend) return hashRecording(be);
  if (be instanceof TeeBackend) return `tee(${drawn(be.a)},${drawn(be.b)})`;
  if (be instanceof m.svg.SvgBackend) return `svg(${digest(be.out.join(""))},${digest(JSON.stringify(be.stats))})`;
  const px = (be as { px?: unknown }).px;
  if (ArrayBuffer.isView(px)) {
    const { w, h } = be as { w?: number; h?: number };
    return `${be.constructor.name}(${w}x${h},${digestInts(px as unknown as ArrayLike<number>)})`;
  }
  return be.constructor.name;
}
function rings(c: Composer): string {
  const l = c.ctx.layout.get(c.section)?.current;
  const r = l ? m.rings.placedRingLabels.get(l) : undefined;
  if (!r) return "-";
  return digest(
    JSON.stringify([r.n, r.text, [...r.side], [...r.pts], [...r.item], [...r.ux], r.at, r.hidden], (_, v) =>
      typeof v === "number" ? v.toString() : v,
    ),
  );
}
let HIT_KEYS: string[] = [];
function grid(c: Composer): string {
  const l = c.ctx.layout.get(c.section)?.current;
  if (!l) return "-";
  const hit = m.root.Hit.create() as unknown as Record<string, number> & ReturnType<typeof m.root.Hit.create>;
  const out: string[] = [];
  for (let y = GRID / 2; y < l.height; y += GRID)
    for (let x = GRID / 2; x < l.width; x += GRID) {
      let s: string;
      try {
        s = c.ctx.layout.hitTest(c.section, x, y, hit)
          ? `${HIT_KEYS.map((k) => hit[k]).join(",")},${hit.span[0]},${hit.span[1]}`
          : "";
      } catch (e) {
        s = `throw ${e instanceof Error ? e.message : String(e)}`;
      }
      out.push(s);
    }
  return digest(out.join(";"));
}

/**
 * Test files whose frames vary from run to run of the same code, by what: their drawing, stats and hit grid are written
 * as `varies` (LC-46). The jbrowse fixtures' bbi blocks are inflated through `DecompressionStream`, and their arrival
 * order, with `MaxInFlight` 6, sometimes leaves `gc_src`'s two edge regions at level 4 unasked (15 requests, not 17),
 * so those ends are drawn coarse; a longer settle did not remove it.
 */
const VARIES = ["_test/jbrowse-zoom.test.ts", "_test/jbrowse-collapse.test.ts", "_test/jbrowse-nav.test.ts"];

type Frame = (this: unknown, be: RenderBackend, flags: number) => void;
let frame: Frame;
beforeAll(async () => {
  m = {
    draw: await import("../../src/draw/index.js"),
    rings: await import("../../src/draw/ring-labels.js"),
    root: await import("../../src/index.js"),
    testing: await import("../../src/testing/index.js"),
    svg: await import("../../src/svg/backend.js"),
  };
  HIT_KEYS = Object.entries(m.root.Hit.create())
    .filter(([k, v]) => typeof v === "number" && k !== "x" && k !== "y")
    .map(([k]) => k);
  const proto = m.draw.SectionComposer.prototype as unknown as Record<string, (...a: unknown[]) => unknown>;
  if (Object.hasOwn(proto.frame, "proof")) return;
  frame = proto.frame as unknown as Frame;
  const framedOf = function (this: unknown, be: RenderBackend, flags: number) {
    return framed(this as Composer, be, flags);
  };
  proto.frame = Object.assign(framedOf, { proof: true }) as unknown as (...a: unknown[]) => unknown;
  if (COVER)
    for (const k of Object.getOwnPropertyNames(proto)) {
      const d = Object.getOwnPropertyDescriptor(proto, k);
      if (k === "constructor" || !d || typeof d.value !== "function") continue;
      const f = d.value as (...a: unknown[]) => unknown;
      calls[k] = 0;
      proto[k] = function (this: unknown, ...a: unknown[]) {
        calls[k]++;
        return f.apply(this, a);
      };
    }
});
function framed(c: Composer, be: RenderBackend, flags: number) {
  const ctx = c.ctx as { report: SeqViewerContext["report"] };
  const reports: string[] = [];
  const own = Object.hasOwn(ctx, "report");
  const report = ctx.report;
  ctx.report = (code: string, message: string) => {
    reports.push(`${code} ${message}`);
    report.call(ctx, code, message);
  };
  try {
    frame.call(c, be, flags);
  } finally {
    if (own) ctx.report = report;
    else delete (ctx as Partial<typeof ctx>).report;
  }
  if (!OUT) return;
  const st = expect.getState();
  file ||= st.testPath ?? "";
  const test = st.currentTestName ?? "(outside a test)";
  const n = counts.get(test) ?? 0;
  counts.set(test, n + 1);
  const varies = VARIES.some((f) => file.endsWith(f));
  const fields = [
    test,
    n,
    flags,
    varies ? "varies" : drawn(be),
    varies ? "varies" : JSON.stringify(c.stats),
    +c.animating,
    +c.again,
    c.animStep[0],
    digest(reports.join("\n")),
    reports.length,
    rings(c),
    varies ? "varies" : grid(c),
  ];
  lines.push(fields.join("\t"));
}

afterAll(() => {
  if (!OUT) return;
  const at = join(OUT, relative(ROOT, file || (expect.getState().testPath ?? "unknown")));
  mkdirSync(dirname(at), { recursive: true });
  if (lines.length > 0) writeFileSync(`${at}.tsv`, `${lines.join("\n")}\n`);
  if (COVER) writeFileSync(`${at}.cover.json`, JSON.stringify(calls));
});
