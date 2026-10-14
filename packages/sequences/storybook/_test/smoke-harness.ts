/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headless viewer the Seq* story smoke test (smoke.test.ts, smoke-more.test.ts) opens each story document on: the stories' own
// spec (the default registry, the fixture structures, the interactions extension) with a manual clock and the fixture
// text measurer, no canvas. It settles the graph, shows every view, applies every snapshot and every patch a story's
// buttons make, and after each step reads what a person would see go wrong: the issues, and a visible section with
// nothing in it.

import type { Issue } from "@molstar/seq-core/model";
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  manualScheduler,
  type SeqViewerContext,
  type SeqViewerRegistryEntry,
} from "@molstar/seq-viewer";
import { InteractionsViewerRegistry } from "@molstar/seq-viewer/interactions/viewer";
import { fixtureMeasurer, streamsOf } from "@molstar/seq-viewer/testing";
import { fixtureLoaders, MolstarFixture } from "../../core/_test/format/load-fixtures.ts";
import type { AllowedWarning, SmokeCase } from "../src/_stories/support/seqstar-smoke.ts";

export interface Finding {
  readonly step: string;
  readonly problem: string;
}

const SIZE: readonly [number, number] = [1000, 400];
const flush = () => new Promise((ok) => setTimeout(ok, 0));

/** Settles the stream manager: passes, clock ticks and microtasks until a few quiet rounds. */
export async function settle(ctx: SeqViewerContext, scheduler: ReturnType<typeof manualScheduler>) {
  for (let k = 0; k < 4; k++) {
    streamsOf(ctx).pass();
    scheduler.advance(16);
    await flush();
  }
  await ctx.graph.update();
}

function size(ctx: SeqViewerContext) {
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(SIZE[0], SIZE[1], 1);
}

const text = (i: Issue) => `${i.severity} ${i.code} ${i.path}: ${i.message}`;

const allowed = (i: Issue, allow: readonly AllowedWarning[]) =>
  allow.some((a) => a.code === i.code && (a.includes === undefined || i.message.includes(a.includes)));

/** The issues the context shows now, those that are not allowed: errors, and warnings not on the story's list. */
function badIssues(ctx: SeqViewerContext, allow: readonly AllowedWarning[]): string[] {
  return ctx.issues.current.value
    .filter((i) => i.severity === "error" || (i.severity === "warning" && !allowed(i, allow)))
    .map(text);
}

/** The current view's visible sections, each with a block and something laid out. */
function emptySections(ctx: SeqViewerContext): string[] {
  const view = ctx.graph.assembled.value?.view?.resolved;
  const set = ctx.graph.assembled.value?.set;
  if (!view || !set) return ["no view assembled"];
  const out: string[] = [];
  for (let k = 0; k < view.sections.count; k++) {
    const ref = set.strings.get(view.sections.ref[k]);
    if (view.sections.hidden[k]) continue;
    if (view.sections.blockCount[k] === 0) out.push(`section ${ref}: no block`);
    else if ((ctx.layout.get(ref)?.current.entryCount ?? 0) === 0) out.push(`section ${ref}: nothing laid out`);
  }
  return out;
}

/** What a step gave back: `false` (an action refused, an expectation not met) or issues, as findings unless allowed. */
function stepProblems(result: unknown, allow: readonly AllowedWarning[]): string[] {
  if (result === false) return ["the step failed (an action refused, or an expectation that did not hold)"];
  if (!Array.isArray(result)) return [];
  return (result as Issue[])
    .filter((i) => i.severity === "error" || (i.severity === "warning" && !allowed(i, allow)))
    .map(text);
}

/** Opens the case and runs every step; the findings, each once, empty when the story is sound. */
export async function smoke(c: SmokeCase): Promise<Finding[]> {
  const findings: Finding[] = [];
  const seen = new Set<string>();
  const note = (step: string, problem: string) => {
    if (seen.has(problem)) return;
    seen.add(problem);
    findings.push({ step, problem });
  };
  const allow = c.allow ?? [];
  const scheduler = manualScheduler();
  const spec = DefaultSeqViewerSpec();
  const registry: SeqViewerRegistryEntry[] = [
    DefaultRegistry,
    ...(c.fixtureStructures === false ? [] : [{ svs: MolstarFixture }]),
    InteractionsViewerRegistry,
    ...(c.registry ?? []),
  ];
  const ctx = createSeqViewer(
    { ...spec, registry },
    { scheduler, loaders: c.loaders ?? fixtureLoaders(), text: fixtureMeasurer() },
  );
  const check = async (step: string) => {
    size(ctx);
    await settle(ctx, scheduler);
    for (const p of badIssues(ctx, allow)) note(step, p);
    for (const p of emptySections(ctx)) note(step, p);
  };
  try {
    const done = ctx.setDocument(c.doc);
    size(ctx);
    const r = await done;
    if (!r.ok) note("open", `the document did not open: ${r.issues.map(text).join(" | ")}`);
    if (c.view) ctx.views.setCurrent(c.view);
    await check("open");
    const first = ctx.views.current.value;
    for (const v of ctx.views.list.value.map((x) => x.ref)) {
      if (v === first) continue;
      ctx.views.setCurrent(v);
      await check(`view ${v}`);
    }
    if (first) ctx.views.setCurrent(first);
    await check("back to the first view");
    for (const s of ctx.snapshots.list.value) {
      await ctx.snapshots.apply(s.ref, { animate: false });
      await check(`snapshot ${s.ref}`);
    }
    await ctx.snapshots.apply(null, { animate: false });
    await check("back to the base");
    for (const s of c.steps ?? []) {
      const step = `step "${s.label}"`;
      try {
        for (const p of stepProblems(await s.run(ctx), allow)) note(step, p);
      } catch (e) {
        note(step, `threw: ${e instanceof Error ? e.message : String(e)}`);
      }
      await check(step);
    }
    if (c.streams) {
      const st = ctx.streams.stats.current;
      if (st.requests === 0) note("end", "the story streams, but no region was asked for");
      if (st.errors > 0) note("end", `${st.errors} region requests failed`);
    }
  } finally {
    ctx.dispose();
  }
  return findings;
}
