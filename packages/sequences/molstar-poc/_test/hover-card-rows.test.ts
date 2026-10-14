/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Every coverage row gets the hover card's "Show in 3D" (seqspec:libraries#the-mol-adapter): the PDB coverage opened,
// each structure opened, the pointer resting on each chain row's segment, through the real layout, hit test and card
// controller; coming from outside, and moving straight from one row to the next (a card open on the last one). With
// 1TUP chain D shown, the button says what it would do, by the author's chain IDs: "Shown in 3D" on 1TUP B
// (label_asym_id D), "Bring chain A to front" on 1TUP A (label_asym_id C), "Show in 3D" on the others.

import type { SeqSet } from '@molstar/seq-core/model';
import {
  createSeqViewer,
  DefaultRegistry,
  DefaultSeqViewerSpec,
  Hit,
  HoverKind,
  manualScheduler,
  SeqViewerConfig,
} from '@molstar/seq-viewer';
import { expect, it } from 'vitest';
import { helpers } from '../../viewer/_test/layout-fixtures.js';
import { createHoverCard } from '../../viewer/src/hover-card.js';
import { fixtureMeasurer } from '../../viewer/src/testing/index.js';
import { coverageRows, MolstarSvsEntry } from '../src/index.js';
import { type Shown, showInStructureAction } from '../src/viewer.js';
import { demoLoaders, P53 } from './fixtures.js';
import { P53_STRUCTURES, P53_TRACKS, p53CoverageOptions, p53DemoDoc } from './p53-demo.js';

const SECTION = 'main.main';

it("every chain row's segment opens the card with Show in 3D", async () => {
  const scheduler = manualScheduler();
  const shown: Shown = { structure: { source: { service: 'pdbe', id: '1tup' }, chains: ['D'] }, colorBy: null };
  const action = showInStructureAction(
    () => {},
    () => shown,
  );
  const ctx = createSeqViewer(
    {
      ...DefaultSeqViewerSpec(),
      registry: [DefaultRegistry, { svs: MolstarSvsEntry, actions: [action] }],
      config: [[SeqViewerConfig.UI.ItemActions, [action.name]]],
    },
    { scheduler, loaders: demoLoaders(), text: fixtureMeasurer() },
  );
  const done = ctx.setDocument(p53DemoDoc(P53.sequence));
  for (const ref of ctx.viewports.sections.value) ctx.viewports.get(ref)?.setSize(1200, 2000, 1);
  expect((await done).ok).toBe(true);
  ctx.ui.toggleGroup('g-pdb');
  for (const s of P53_STRUCTURES) ctx.ui.toggleGroup(`${P53_TRACKS.coverage}-${s.complex}`);
  await ctx.graph.update();
  scheduler.advance(50);
  const card = createHoverCard(ctx);
  const t = helpers(ctx);
  const set = ctx.graph.assembled.value?.set as SeqSet;
  const a = set.annotations.find((x) => x && set.strings.get(x.key) === 'coverage');
  if (!a) throw new Error('no coverage');
  const entry = a.fields.find((f) => set.strings.get(f.name) === 'entry');
  const chain = a.fields.find((f) => set.strings.get(f.name) === 'chain');
  const p53 = coverageRows(set, p53CoverageOptions);
  const missed: string[] = [];
  const titles = new Map<string, string>();
  const hit = Hit.create();
  for (const straight of [false, true]) {
    for (let r = 0; r < a.count; r++) {
      const s = P53_STRUCTURES.find(
        (x) => x.label === set.strings.get((entry?.data as Int32Array | undefined)?.[r] ?? -1),
      );
      const c = set.strings.get((chain?.data as Int32Array | undefined)?.[r] ?? -1);
      const track = `${P53_TRACKS.coverage}-${s?.complex}`;
      const view = ctx.graph.assembled.value?.view;
      const k = view?.resolved.tracks.split[t.trackIndex(track)] ?? -1;
      const labels = (view?.prepared[k]?.state as { split?: { labels: string[] } } | undefined)?.split?.labels ?? [];
      // Its middle on P04637 (the rows sit on the chains' units; the provider's own rows are in UniProt numbers).
      const mid = (p53.beg[r] - 1 + p53.end[r]) >> 1;
      if (!straight) {
        ctx.hover.setFromHit(Hit.create(), 'pointer');
        scheduler.advance(1000);
      }
      ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, mid), t.yOf(SECTION, track, labels.indexOf(c)), hit);
      ctx.hover.setFromHit(hit, 'pointer');
      scheduler.advance(1000);
      const open = card.current.value;
      const label = `${s?.label} ${c}${straight ? ', straight from the last row' : ''} (hit ${hit.kind === HoverKind.item ? 'item' : hit.kind})`;
      if (!open || open.actions[0]?.name !== action.name || open.target.item !== a.item[r]) missed.push(label);
      else titles.set(`${s?.label} ${c}`, open.actions[0].title);
    }
  }
  expect(missed).toEqual([]);
  // By the author's chain IDs: 1TUP's B is label_asym_id D, its A is C.
  expect([titles.get('1TUP B'), titles.get('1TUP A'), titles.get('2OCJ A'), titles.get('AlphaFold A')]).toEqual([
    'Shown in 3D',
    'Bring chain A to front',
    'Show in 3D',
    'Show in 3D',
  ]);
  card.dispose();
  ctx.dispose();
});
