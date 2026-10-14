/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Two sequences on one axis, side by side (the strip navigator's two blocks) or stacked: a ruler, letters, ranges, values
// and points, in `tracks` rows. By value, columnar annotations, for the browser suite.

import { createSeqViewSpec } from "@molstar/seq-core/builder";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";
/** `n` residues from a fixed protein, rotated per seed. */
export const residues = (n: number, seed = 0) =>
  Array.from({ length: n }, (_, i) => PROTEIN[(i * 7 + seed * 13) % PROTEIN.length]).join("");

export function twoChains(
  opts: { stacked?: boolean; length?: number; tracks?: number; layers?: readonly string[] } = {},
) {
  const layers = new Set(opts.layers ?? ["letters", "ranges", "values", "points"]);
  const n = opts.length ?? 400;
  const m = Math.round(n * 0.7);
  const b = createSeqViewSpec({ title: "Two chains" });
  const a = b.data.sequence({ ref: "a", label: "Chain A", value: residues(n, 1) });
  const c = b.data.sequence({ ref: "c", label: "Chain C", value: residues(m, 2) });
  const ranges = (ref: string, target: unknown, len: number, seed: number) => {
    const start = Array.from({ length: 12 }, (_, k) => (k * 37 + seed * 11) % (len - 40));
    return b.data.annotation({
      ref,
      geometry: "ranges",
      target: target as never,
      value: {
        start,
        stop: start.map((s, k) => s + 6 + ((k * 5) % 30)),
        label: start.map((_, k) => `F${k}`),
        category: start.map((_, k) => `c${k % 4}`),
      },
    } as never);
  };
  const values = (ref: string, target: unknown, len: number, seed: number) =>
    b.data.annotation({
      ref,
      geometry: "values",
      target: target as never,
      value: { value: Array.from({ length: len }, (_, i) => ((i * 31 + seed * 7) % 100) / 100) },
    } as never);
  const ra = ranges("r1", a, n, 1);
  const rc = ranges("r2", c, m, 2);
  const va = values("v1", a, n, 1);
  const vc = values("v2", c, m, 2);
  const v = b.view({
    ref: "v",
    label: "Chains",
    arrangement: opts.stacked ? "stacked" : "side",
    height: "box",
  } as never);
  v.block({ sequence: a }).block({ sequence: c });
  v.track({ label: "" }).ruler();
  for (let k = 0; k < (opts.tracks ?? 3); k++) {
    if (layers.has("letters")) v.track({ label: `Letters ${k}` }).letters();
    if (layers.has("ranges")) v.track({ label: `Features ${k}` }).ranges({ data: [ra, rc] } as never);
    if (layers.has("values")) v.track({ label: `Values ${k}` }).values({ data: [va, vc] } as never);
    if (layers.has("points")) v.track({ label: `Sites ${k}` }).points({ data: [ra, rc] } as never);
  }
  return b.getDocument();
}
