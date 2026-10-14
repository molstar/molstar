/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A streamed alignment's row labels are held as its regions are (seqspec:sources#the-viewers-side; VT-07): a labels
// chunk counts its strings' bytes under `Stream.CacheBytes` and is evicted as a region is, its rows' labels going with
// it; the chunks in view are kept, and a chunk scrolled back into view is asked again.

import { type StreamedData, StreamProp } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { case9 } from "../../core/_test/format/cases.ts";
import { SeqViewerConfig } from "../src/index.ts";
import type { Streams } from "../src/streams.ts";
import { opened, preparedRenderers } from "./layout-fixtures.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));

describe("a streamed alignment's labels under the cache bound (VT-07)", () => {
  /** Case 9's 18,000 rows scrolled through `rows` at a time, under a 64 KB cache, and back to the top. */
  async function tour(rows: number) {
    const t = await opened(case9(), { renderers: preparedRenderers(), size: { "msa.main": [1000, 500] } });
    t.ctx.config.set(SeqViewerConfig.Stream.CacheBytes, 64 * 1024);
    const s = t.ctx.streams as Streams;
    const vp = t.vp("msa.main");
    const scroll = async (y: number) => {
      vp.scrollTo(0, y, "api");
      t.layout("msa.main");
      for (let k = 0; k < 4; k++) {
        s.pass();
        await flush();
      }
    };
    const labels = () => {
      let n = 0;
      let bytes = 0;
      let first: string | undefined;
      for (const a of t.ctx.graph.assembled.value?.set.alignments ?? []) {
        const d = a?.custom.props.get(StreamProp) as StreamedData | undefined;
        if (!d) continue;
        n += d.labels.size;
        for (const x of d.labels.values()) bytes += 2 * x.length;
        first ??= d.labels.get(0);
      }
      return { n, bytes, first };
    };
    await scroll(0);
    expect(labels().first).toBeTruthy();
    for (let y = 0; y < 18_000 * 18; y += rows * 18) await scroll(y);
    // Scrolled through: the chunks left behind went, the ones in view are held and counted.
    expect(labels().n).toBeLessThan(5_000);
    expect(labels().first).toBeUndefined();
    expect(s.stats.current.bytes).toBeGreaterThanOrEqual(labels().bytes);
    await scroll(0);
    expect(labels().first).toBeTruthy();
    t.ctx.dispose();
  }

  it("are counted and evicted with their chunks, and asked again in view: 4,000 rows a scroll", () => tour(4000));

  it("are counted and evicted with their chunks, and asked again in view", { timeout: 30_000, tags: ["slow"] }, () =>
    tour(400),
  );
});
