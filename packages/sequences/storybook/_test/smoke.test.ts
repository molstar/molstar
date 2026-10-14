/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every Seq* story's document opened on a headless viewer (smoke-harness.ts): the story registry
// (lib/seqstar/storybook/src/_stories/support/seqstar-smoke.ts) lists what each story opens; each document is set, its graph settled,
// every view shown, every snapshot and every patch a story's buttons make applied, and after each step there must be
// no issue of severity error, no warning the story does not show on purpose, and no visible section without a block
// or a laid-out track. A story file with an export the registry does not name fails here, so a new story brings its
// entry. Several stories broke silently before this: a lazy alignment sharing a ref with its view, buttons that
// patched a track instead of its layer, a filter that kept rows no track draws, a query naming one that had not run.
// The stories' documents are opened half here and half in smoke-more.test.ts (smoke-stories.ts), so the halves run at
// once; the four with the largest documents are the slow tier's, every other runs in the fast tier.

import { readdirSync, readFileSync } from "node:fs";
import { join } from "node:path";
import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { describe, expect, it } from "vitest";
import { CASES } from "../../core/_test/format/cases.js";
import { keyedFamilyDoc } from "../../viewer/_test/keyed-family.js";
import { BAD, styled } from "../src/_stories/support/seqstar-docs-views.js";
import { storyPatch } from "../src/_stories/support/seqstar-patch.js";
import { SMOKE } from "../src/_stories/support/seqstar-smoke.js";
import { smoke } from "./smoke-harness.js";
import { root, SLOW_STORIES, storyDocuments } from "./smoke-stories.js";

const stories = join(root, "lib/seqstar/storybook/src/_stories");

/** The exports of every story file, as `file#Export`. */
function storyKeys(): string[] {
  return readdirSync(stories)
    .filter((f) => f.endsWith(".stories.tsx"))
    .flatMap((f) =>
      [...readFileSync(join(stories, f), "utf8").matchAll(/^export const (\w+)\b/gm)].map((m) => `${f}#${m[1]}`),
    );
}

describe("the story registry", () => {
  const keys = new Set(SMOKE.map((s) => s.key));
  it("names every story export, and only story exports", () => {
    const real = new Set(storyKeys());
    expect([...real].filter((k) => !keys.has(k))).toEqual([]);
    // An entry that is skipped may name a story another branch has not committed yet.
    expect(SMOKE.filter((s) => !s.skip && !real.has(s.key)).map((s) => s.key)).toEqual([]);
  });
  it("gives each story documents, the key of one that has them, or a reason", () => {
    for (const s of SMOKE) {
      expect(
        [s.cases, s.same, s.skip].filter((x) => x !== undefined),
        s.key,
      ).toHaveLength(1);
      if (s.same) expect(SMOKE.find((x) => x.key === s.same)?.cases, `${s.key} -> ${s.same}`).toBeTypeOf("function");
    }
  });
  it("names the slow tier's stories among those with documents", () => {
    expect([...SLOW_STORIES].filter((k) => !SMOKE.find((s) => s.key === k)?.cases)).toEqual([]);
  });
});

storyDocuments(0, 2);

describe("the smoke test catches what broke stories", () => {
  const problems = async (c: Parameters<typeof smoke>[0]) => (await smoke(c)).map((f) => `${f.step}: ${f.problem}`);
  it("a warning the story does not allow", async () => {
    const doc = styled(CASES.case3(), "family.main", BAD);
    const found = await problems({ label: "bad", doc });
    expect(found.some((p) => p.includes("svs.section.style"))).toBe(true);
  });
  it("a button's patch naming a ref that is not there", async () => {
    const found = await problems({
      label: "bad patch",
      doc: CASES.case1(),
      steps: [
        {
          label: "wrong ref",
          run: async (ctx) =>
            (
              await storyPatch(ctx, {
                format: "seq_view_spec_patch",
                ops: [{ op: "set", ref: "no-such-node", params: {} }],
              })
            ).issues,
        },
      ],
    });
    expect(found.some((p) => p.startsWith('step "wrong ref"'))).toBe(true);
  });
  it("a keyed-set button that patches the track instead of its rows layer", async () => {
    const family = [
      { label: "A", gapped: "ACDE" },
      { label: "B", gapped: "ACDF" },
    ];
    const found = await problems({
      label: "track patch",
      doc: keyedFamilyDoc(family),
      steps: [
        {
          label: "by identity",
          run: async (ctx) =>
            (
              await storyPatch(ctx, {
                format: "seq_view_spec_patch",
                ops: [{ op: "set", ref: "members", params: { order: "data" } }],
              })
            ).issues,
        },
      ],
    });
    expect(found.some((p) => p.startsWith('step "by identity"'))).toBe(true);
  });
  it("a refused action", async () => {
    const found = await problems({
      label: "refused",
      doc: CASES.case1(),
      steps: [{ label: "no such action", run: (ctx) => ctx.actions.run("no.such.action", null) }],
    });
    expect(found.some((p) => p.startsWith('step "no such action": the step failed'))).toBe(true);
  });
  it("a view with no block", async () => {
    const b = createSeqViewSpec({ title: "empty" });
    const v = b.view({ ref: "v", label: "V" });
    v.track({ label: "" }).ruler();
    const found = await problems({ label: "no block", doc: b.getDocument() });
    expect(found.length).toBeGreaterThan(0);
  });
});
