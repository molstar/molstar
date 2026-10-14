/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// GenomeScale's adapter (CS2; the Storybook story `Seq*/Performance/Large documents › Genome-scale, 3M nucleotides` has
// its own copy): the synthetic genome of core/_test/format/genome-fixture.ts on a context whose renderers are the kernel's
// with each `prepare` and `draw` timed. The old story wrapped `linear.draw` of the canvas renderers; the next's
// renderers have one `draw` (the composer's `DrawView` is the same object for every pass), so the wrapper is on the
// entry itself. Sums in a typed array, read back as `__genome` in the console: per renderer the draws' calls, total ms
// and slowest call, and the prepares' ms by type.

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { genomeDoc, synthGenome } from "../../../core/_test/format/genome-fixture.js";
import type { DrawView } from "../../src/draw/view.js";
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  KernelRenderers,
  type LayerInput,
  type LayerRenderer,
  type Scheduler,
  type SeqViewerContext,
} from "../../src/index.js";
export const SECTION = "genome.main";

/** Per renderer: draw ms summed, calls, the slowest call; prepare ms by type. */
export const PROBE = {
  names: [] as string[],
  draw: new Float64Array(3 * 64),
  prepare: {} as Record<string, number>,
  /** The document last built (the Storybook story sets it, to find a coding gene's exons from the console). */
  doc: null as SeqViewSpecDoc | null,
  /** The draws per renderer since the last reset: [name, calls, total ms, max ms]. */
  report() {
    return PROBE.names
      .map((n, i) => [n, PROBE.draw[3 * i + 1], +PROBE.draw[3 * i].toFixed(1), +PROBE.draw[3 * i + 2].toFixed(2)])
      .filter((r) => (r[1] as number) > 0);
  },
  reset() {
    PROBE.draw.fill(0);
  },
};

type Draw = (v: DrawView, s: never, p: never) => void;
function timedDraw(i: number, fn: Draw | undefined): Draw | undefined {
  if (!fn) return fn;
  return (v, st, p) => {
    const t = performance.now();
    fn(v, st, p);
    const ms = performance.now() - t;
    PROBE.draw[3 * i] += ms;
    PROBE.draw[3 * i + 1]++;
    if (ms > PROBE.draw[3 * i + 2]) PROBE.draw[3 * i + 2] = ms;
  };
}

/** The kernel's renderers with their prepares, draws and overviews timed. */
export const TIMED: LayerRenderer[] = (KernelRenderers as LayerRenderer[]).map((r, i) => {
  PROBE.names[i] = r.name;
  return {
    ...r,
    prepare: (input: LayerInput, o?: unknown) => {
      const t = performance.now();
      const out = (r.prepare as (a: LayerInput, b?: unknown) => unknown)(input, o);
      PROBE.prepare[r.name] = (PROBE.prepare[r.name] ?? 0) + performance.now() - t;
      return out;
    },
    draw: timedDraw(i, r.draw as Draw | undefined) as LayerRenderer["draw"],
  };
});

/** Everything is by value: a load is a bug, never a fetch. */
const NO_LOADS = { load: () => Promise.reject(new Error("the genome story loads nothing")) };

/** The genome document (`genomeDoc` over `synthGenome`), shared by the stack's story and tests. */
export function genomeDocument(o: Parameters<typeof synthGenome>[0] = {}): SeqViewSpecDoc {
  return genomeDoc(synthGenome(o));
}

/** A context showing the genome with the timed renderers. */
export function genomeContext(
  scheduler: Scheduler,
  doc: SeqViewSpecDoc,
): { ctx: SeqViewerContext; ready: Promise<void> } {
  const spec = DefaultSeqViewerSpec();
  const [first, ...rest] = spec.registry ?? [];
  const ctx = createSeqViewer(
    { ...spec, registry: [{ ...first, renderers: TIMED }, ...rest] },
    { scheduler, loaders: NO_LOADS },
  );
  (globalThis as { __genome?: typeof PROBE }).__genome = PROBE;
  const ready = ctx.setDocument(doc).then((r) => {
    if (!r.ok) throw new Error(JSON.stringify(r.issues));
  });
  return { ctx, ready };
}
