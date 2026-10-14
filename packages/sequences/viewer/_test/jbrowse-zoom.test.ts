/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST4b: the JBrowse 2 demo's sources (jbrowse-story-fixtures.ts, @molstar/seq-core/stream/gmod) at a whole chromosome read only
// what zoom levels need: no 2bit sequence (the letters gate), no VCF data (the counts are the tabix index's estimate),
// and for BigBed and BigWig nothing from their data sections (the counts and summaries are their zoom records). The
// reads replay recorded slices of the public files (core/_test/stream/gmod/fixtures/whole/; the .tbi is the main fixtures');
// `SEQSTAR_RECORD=1` records them again, by hand, over curl (never in a test run).

import { afterAll, expect, it } from "vitest";
import type { RangeReader } from "../../core/src/stream/gmod/index.ts";
import { Dirty } from "../src/draw/composer.ts";
import type { Streams } from "../src/streams.ts";
import { composed, recording } from "./composer-fixtures.ts";
import { demoSources, RECORD } from "./jbrowse-sources.ts";
import { JBROWSE_SECTION, jbrowseDoc } from "./jbrowse-story-fixtures.ts";

const CHR21 = 46_709_983;
const src = demoSources("whole");
afterAll(src.save);

/** A bbi file's data section (raw items or values): [fullDataOffset, the first zoom level's data). */
async function dataSection(read: RangeReader): Promise<[number, number]> {
  const head = await read(0, 64 + 24 * 10);
  const dv = new DataView(head.buffer, head.byteOffset, head.byteLength);
  const zooms = dv.getUint16(6, true);
  let first = Number.POSITIVE_INFINITY;
  for (let k = 0; k < zooms; k++) first = Math.min(first, Number(dv.getBigUint64(64 + 24 * k + 8, true)));
  return [Number(dv.getBigUint64(16, true)), first];
}

it("a whole chromosome reads zoom levels only: no 2bit sequence, no VCF data, no BigBed or BigWig data blocks", async () => {
  const { io, asked, loaders } = src;
  const t = await composed(jbrowseDoc("chr21", 0, CHR21), {
    width: 1200,
    height: 400,
    dpr: 1,
    loaders,
  });
  const opened = Object.fromEntries(Object.entries(io).map(([k, x]) => [k, x.reads.length]));
  for (let k = 0; k < 40 && (k < 3 || t.ctx.streams.stats.current.inFlight > 0); k++) {
    (t.ctx.streams as Streams).pass();
    await new Promise((ok) => setTimeout(ok, RECORD ? 200 : 0));
  }
  expect(t.ctx.issues.current.value.map((i) => i.message)).toEqual([]);
  expect(t.ctx.streams.stats.current.inFlight).toBe(0);
  // Every items and values request is a level's; nothing of the letters is asked.
  const data = asked.filter((r) => r.system !== null && r.level < 0);
  expect(data).toEqual([]);
  expect(asked.some((r) => r.system !== null && r.source === "genes_src")).toBe(true);
  // The transports: nothing of the 2bit or the VCF after their manifests; the bbi files' reads outside their data.
  expect(io.twobit.reads.length).toBe(opened.twobit);
  expect(io.vcf.reads.length).toBe(opened.vcf);
  for (const key of ["bigwig", "bigbed"] as const) {
    const [a, e] = await dataSection(io[key].read);
    const after = io[key].reads.slice(opened[key]);
    expect(after.length).toBeGreaterThan(0);
    expect(
      after.filter(([o, l]) => o < e && o + l > a),
      key,
    ).toEqual([]);
  }
  // Drawn: GC and genes from zoom records, ClinVar's estimate with its hint, the letters' hatch with theirs.
  const be = t.frame(recording(), Dirty.tracks);
  const texts = be.textItems.text.slice(0, be.textItems.n);
  expect(texts).toContain("zoom in to see variants");
  expect(texts).toContain("zoom in to see letters");
  expect(be.rectItems.n).toBeGreaterThan(200);
  expect(JBROWSE_SECTION).toBe(t.section);
  t.ctx.dispose();
}, 600_000);
