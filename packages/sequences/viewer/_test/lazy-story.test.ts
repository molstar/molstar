/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The lazy alignment story's document and loaders (lazy-story-fixtures.ts) opened headlessly: the document validates,
// the alignment resolves into the view's section, and the first pass asks for its tiles; answered, the rows are held.

import { validate } from "@molstar/seq-core/format";
import { expect, it } from "vitest";
import { registries } from "../../core/_test/format/load-fixtures.js";
import { KernelRenderers, type LayerRenderer } from "../src/index.js";
import type { Streams } from "../src/streams.js";
import { opened } from "./layout-fixtures.js";
import { LAZY_ROWS, lazyStoryDoc, lazyStoryLoaders } from "./lazy-story-fixtures.js";

const flush = () => new Promise((ok) => setTimeout(ok, 0));

it("the lazy alignment story's document resolves its alignment and asks for its first tiles", async () => {
  const doc = lazyStoryDoc();
  expect(validate(doc, registries()).filter((i) => i.severity === "error")).toEqual([]);
  const t = await opened(doc, { renderers: [...KernelRenderers] as LayerRenderer[], loaders: lazyStoryLoaders() });
  const m = t.ctx.graph.assembled.value?.view;
  expect(m?.set.alignments.some((a) => a?.lazy)).toBe(true);
  const section = t.ctx.viewports.sections.value[0];
  expect(section).toBeTruthy();
  expect(t.ctx.layout.get(section)?.current.rows.count).toBe(LAZY_ROWS);
  (t.ctx.streams as Streams).pass();
  expect(t.ctx.streams.stats.current.requests).toBeGreaterThan(0);
  for (let k = 0; k < 4; k++) {
    t.scheduler.advance(16);
    await flush();
  }
  expect(t.ctx.streams.stats.current.coarse + t.ctx.streams.stats.current.fine).toBeGreaterThan(0);
  // Nothing fails with the failure knob at 0: the reference row and the labels are served too.
  expect(t.ctx.streams.stats.current.errors).toBe(0);
  expect(t.ctx.issues.current.value.map((i) => `${i.code}: ${i.message}`)).toEqual([]);
  t.ctx.dispose();
});

it("the story's failure knob fails tiles only: labels and rows still load", async () => {
  const t = await opened(lazyStoryDoc(), {
    renderers: [...KernelRenderers] as LayerRenderer[],
    loaders: lazyStoryLoaders(() => Promise.reject(new Error("503 from the story's knob"))),
  });
  (t.ctx.streams as Streams).pass();
  for (let k = 0; k < 4; k++) {
    t.scheduler.advance(16);
    await flush();
  }
  const issues = t.ctx.issues.current.value.map((i) => i.message);
  expect(issues.some((m) => m.includes("503"))).toBe(true);
  expect(issues.some((m) => /rows\/|labels|404/.test(m))).toBe(false);
  t.ctx.dispose();
});
