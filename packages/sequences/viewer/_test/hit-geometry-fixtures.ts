/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Sized shapes for the hits by drawn geometry (HR2), shared by the test and the heap guard: a track per sized type over
// 400 residues: dots of radius 2 to 8 by score, bars of 0.2 to 1 of their lane by confidence, boundaries' glyphs of
// radius 2 to 8 by weight (one at each end of the sequence too, also drawn as points and as variants), pair ends of
// radius 2 to 8 and arcs 1 to 16 px wide by count.

import { createSeqViewSpec } from "@molstar/seq-core/builder";

const SEQ =
  "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ".repeat(4);
export const S = "v.main";
const scale = (field: string, range: [number, number]) => ({ theme: "field_scale", params: { field, range } });

/** One track per sized type: dots radius 2–8 by score, bars 0.2–1 of their lane by confidence, pairs' ends and stroke. */
export function hitGeometryDoc() {
  const b = createSeqViewSpec();
  const p = b.data.sequence({ ref: "p", value: SEQ });
  const sites = b.data.annotation({
    ref: "sites",
    geometry: "ranges",
    target: p,
    value: { start: [20, 60], score: [0, 100], category: ["site", "site"] },
  } as never);
  const bars = b.data.annotation({
    ref: "bars",
    geometry: "ranges",
    target: p,
    value: { start: [10, 50], stop: [30, 70], confidence: [0.2, 1] },
  } as never);
  const ticks = b.data.annotation({
    ref: "ticks",
    geometry: "ranges",
    target: p,
    // Inner boundaries, and one at each end of the sequence (a terminal insertion's).
    value: { boundary: [0, 20, 60, SEQ.length], weight: [5, 0, 10, 5] },
  } as never);
  const arcs = b.data.annotation({
    ref: "arcs",
    geometry: "pairs",
    target: p,
    value: { start: [10, 150], b_target: ["p", "p"], b_start: [40, 180], count: [1, 10] },
  } as never);
  const pairs = b.data.annotation({
    ref: "pairs",
    geometry: "pairs",
    target: p,
    value: { start: [10, 40], b_target: ["p", "p"], b_start: [30, 80], count: [1, 10] },
  } as never);
  const sec = b.view({ ref: "v", label: "V" }).section({ ref: S }).block({ sequence: p });
  sec
    .track({ ref: "t-dots", label: "Dots", height: 2 })
    .points({ data: sites, radius: scale("score", [2, 8]) } as never);
  sec
    .track({ ref: "t-bars", label: "Bars", height: 2 })
    .ranges({ data: bars, thickness: scale("confidence", [0.2, 1]) } as never);
  sec
    .track({ ref: "t-ticks", label: "Ticks", height: 2 })
    .ranges({ data: ticks, radius: scale("weight", [2, 8]) } as never);
  // The same boundaries as points and as stacked variants: glyphs on the edges, hit there.
  sec.track({ ref: "t-tick-points", label: "Tick points", height: 2 }).points({ data: ticks });
  sec.track({ ref: "t-tick-variants", label: "Tick variants", height: 2 }).variants({ data: ticks });
  sec.track({ ref: "t-pairs", label: "Pairs", height: 4 }).pairs({
    data: pairs,
    display: "ends",
    radius: scale("count", [2, 8]),
  } as never);
  sec.track({ ref: "t-arcs", label: "Arcs", height: 4 }).pairs({
    data: arcs,
    display: "arcs",
    width: scale("count", [1, 16]),
  } as never);
  return b.getDocument();
}
