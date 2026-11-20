/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The 2OCJ view and a counting side for the structure link's performance cases (its time and its allocation).

import { readFileSync } from 'node:fs';
import { P53Views } from '@molstar/seqx-examples';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { EntityStart, Observed } from '@molstar/seqx-examples/p53/structure';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import type { StructureTarget } from '@molstar/seqx-viewer/link-structure';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';

/** The 2OCJ view, its entity's file read from the examples' fixtures. */
export async function open2ocj(): Promise<SeqViewer> {
  const read = (uri: string) => readFileSync(`${PackageFolder}${uri}`);
  const v = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    ports: {
      assets: () => ({
        reader: {
          readText: async (uri) => read(uri).toString('utf8'),
          readBytes: async (uri) => new Uint8Array(read(uri)),
        },
      }),
    },
  });
  v.resize(1440);
  await v.load(P53Views.find((x) => x.id === 'structure')!.build());
  await v.settle();
  return v;
}

/** A side that keeps the last highlight and counts them. */
export const countingSide = () => {
  const side = {
    told: 0,
    last: null as StructureTarget | null,
    highlight(t: StructureTarget | null) {
      side.told++;
      side.last = t;
    },
    select() {},
    focus() {},
  };
  return side;
};

/** Hover `k`'s position on P04637: modeled residues of the entity, a new one each move. */
export const hoverAt = (k: number) => EntityStart + Observed[0] + 1 + (k % 150);
