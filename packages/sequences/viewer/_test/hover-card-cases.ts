/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The case of the hover card's guard (seqspec:interaction#the-hover-card; hover-card.heap.ts) and its quick twin
// (hover-card.quick.ts): a card open on an item, and a step that bumps the hover's version with the pointer on that
// item (what each move within the item does to the handler). The heap guard and the quick twin share it.

import { case13 } from '../../core/_test/format/cases.js';
import { createHoverCard } from '../src/hover-card.js';
import { type ActionEntry, Hit, type ItemTarget, SeqViewerConfig } from '../src/index.js';
import { state } from '../src/state.js';
import { opened } from './layout-fixtures.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 20_000;
const SECTION = 'features.main';
const Show: ActionEntry<ItemTarget> = { name: 'test:show', title: 'Show', when: () => true, run: () => {} };

/** The card open on the first domain, and the step: one version bump. */
export async function openCard() {
  const t = await opened(case13(), { size: { [SECTION]: [1000, 900] }, registry: [{ actions: [Show] }] });
  t.ctx.config.set(SeqViewerConfig.UI.ItemActions, ['test:show']);
  // The card's own view of the hover: the context's hover record, a version only it listens to (the viewer's other
  // hover subscribers, guarded elsewhere, stay out of this measurement).
  const version = state(0);
  const card = createHoverCard(
    Object.assign(Object.create(t.ctx), { hover: { live: { current: t.ctx.hover.live.current, version } } }),
  );
  const a = Hit.create();
  t.ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, 150), t.yOf(SECTION, 'features.domain'), a);
  t.ctx.hover.setFromHit(a, 'pointer');
  version.set(1);
  t.scheduler.advance(500);
  return { card, step: () => version.set(version.value + 1) };
}
