/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The lazy alignment of `Seq*/Streaming/Streamed data › A lazy alignment by view` (lib/seqstar/storybook, ST3b): its document and its
// sampled `seqstar-tiles/1` transport, shared by the story and lazy-story.test.ts, which opens them headlessly so the
// story cannot silently load nothing again. Tiles of letters (and their share) from a hash of the cell, labels by
// index; nothing kept (the atlas story's transport).

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import { type Loaders, type SeqViewSpecDoc, TILE } from "@molstar/seq-core/format";
import { type TileTransport, tilesStream } from "@molstar/seq-core/stream/tiles";

export const LAZY_URI = "https://example.org/streamed-lazy";
export const LAZY_ROWS = 20_000;
export const LAZY_COLUMNS = 1_200;
const AMINO = "ACDEFGHIKLMNPQRSTVWY";
const encoder = new TextEncoder();

export const sampledTiles: TileTransport = async (path) => {
  const [kind, ...rest] = path.split("/");
  const a = rest.map(Number);
  if (kind === "tiles") {
    const out = new Uint8Array(2 * TILE * TILE);
    for (let i = 0; i < TILE * TILE; i++) {
      const h = Math.imul(i * 2654435761 + a[2] * 97 + a[3] * 31 + a[0], 40503) >>> 0;
      out[i] = h % 7 === 0 ? 45 : AMINO.charCodeAt(h % 20);
      out[TILE * TILE + i] = 110 + (h % 140);
    }
    return out;
  }
  if (kind === "labels")
    return encoder.encode(JSON.stringify(Array.from({ length: a[1] - a[0] }, (_, i) => `SEQ${a[0] + i}`)));
  // A row whole (its label and gapped string): the viewer asks for the reference row.
  if (kind === "rows") {
    let row = "";
    for (let c = 0; c < LAZY_COLUMNS; c++) {
      const h = Math.imul((a[0] + 1) * 73856093 + c * 19349663, 2654435761) >>> 0;
      row += h % 7 === 0 ? "-" : AMINO[(h >>> 8) % 20];
    }
    return encoder.encode(JSON.stringify({ label: `SEQ${a[0]}`, row }));
  }
  throw new Error(`404 ${path}`);
};

/**
 * The loaders over the transport, each tile request through `wrap` (the story's latency and failure knobs: regions
 * only, so a failure is an error band with Retry); labels and rows are answered at once.
 */
export function lazyStoryLoaders(
  wrap: (fetch: () => Promise<Uint8Array>) => Promise<Uint8Array> = (f) => f(),
): Loaders {
  return {
    load: () => Promise.reject(new Error("no files")),
    stream: tilesStream((u) =>
      u === LAZY_URI
        ? (path) => (path.startsWith("tiles/") ? wrap(() => sampledTiles(path)) : sampledTiles(path))
        : null,
    ),
  };
}

/** The story's document: a ruler and the rows of a lazy alignment, `LAZY_ROWS` × `LAZY_COLUMNS`. */
export function lazyStoryDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "A lazy alignment" });
  const aln = b.data.alignment({
    ref: "lazy",
    row_count: LAZY_ROWS,
    column_count: LAZY_COLUMNS,
    source: { uri: LAZY_URI, protocol: "seqstar-tiles/1" },
  } as never);
  // Its own ref: the view named "lazy" like the alignment was `svs.ref.duplicate`, and the story showed nothing (ST3b).
  const v = b.view({ ref: "rows", label: "Lazy", height: "box" });
  v.block({ alignment: aln });
  v.track({ label: "" }).ruler();
  v.track({ label: "Rows" }).rows();
  return b.getDocument();
}
