/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LC-02 (seqspec:viewer#extension-points): every type the renderer contract's members take is `/draw`'s, so a host's
// renderer or navigator names them from the entry (`ringLabels(s, p, out: RingLabelSink)`, `stream?(s): StreamLayer`,
// a navigator's `draw(n: NavView)`, the sinks, pyramids and views the helpers take). The typecheck is the guard; the
// test holds a renderer written against them.

import type {
  FieldTable,
  ItemOrder,
  LayerRenderer,
  LocalCells,
  MarkSet,
  NavView,
  RingLabelSink,
  SizeView,
  StreamLayer,
  TreeRows,
  ValuePyramid,
} from "@molstar/seq-viewer/draw";
import { expect, it } from "vitest";

it("a renderer's ring labels and stream are typed by `/draw`'s names", () => {
  const stream = (_s: unknown): StreamLayer | null => null;
  const ringLabels = (_s: unknown, _p: unknown, _out: RingLabelSink) => {};
  const parts: Pick<LayerRenderer, "stream" | "ringLabels"> = { stream, ringLabels };
  expect(parts.stream?.(null)).toBeNull();
});

it("the helpers' and a navigator's types are named from `/draw`", () => {
  type Named = [FieldTable, ItemOrder, LocalCells, MarkSet, NavView, SizeView, TreeRows, ValuePyramid];
  const none: Named | null = null;
  expect(none).toBeNull();
});
