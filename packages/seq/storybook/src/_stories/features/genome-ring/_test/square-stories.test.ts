/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens } from '@molstar/seqx-core/render';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { partnersOf, rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import {
  cellPoint,
  partnerKeys,
  RingSection,
  squareClicked,
  squareDragged,
  squareOf,
  squareResized,
  squareScattered,
  stemMenuOpen,
} from '../square-gestures.js';
import { ribbonsOf, ringStoryView, strongestStem } from '../views.js';

// The square's and the stem menu's stories (Features/Genome Ring), headless twins on a fake clock: each story's opening
// replayed as the story runs it, what it leaves checked, the view's SVG with its navigators light and dark.

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open() {
  viewer = await openExample(ringStoryView('ring'), {
    delayMs: 0,
    wait: 3000,
    fixtures: true,
    extensions: [rondoViewer],
  });
  viewer.draw();
  return viewer;
}

/** An opening run as the story runs it, its patches landed and the view drawn again. */
async function replay(v: SeqViewer, opening: (v: SeqViewer) => unknown) {
  await opening(v);
  for (let i = 0; i < 10; i++) {
    await vi.advanceTimersByTimeAsync(100);
    v.draw();
  }
  squareOf(v);
}

const include = (v: SeqViewer) => {
  const views = v.spec.root.value.children!.find((c) => c.kind === 'views')!;
  return views.children![0].children!.find((c) => c.id === 'ring')!.params?.include;
};

/** The SVG with the navigators, light and dark: settled, no error, the two differ. */
async function exportedBoth(v: SeqViewer) {
  const svg = async (dark: boolean) => {
    const r = v.renderSvg({ idPrefix: 's-', tokens: dark ? DarkTokens : LightTokens, navigators: 'include' });
    await vi.advanceTimersByTimeAsync(400);
    return r;
  };
  const light = await svg(false);
  const dark = await svg(true);
  expect(light.settled).toBe(true);
  expect(light.issues.filter((i) => i.severity === 'error')).toEqual([]);
  expect(dark.svg).not.toBe(light.svg);
}

describe('The square beside the ring, its gestures replayed', () => {
  it('dragged: chr4 to chr6 moved three chromosomes on', async () => {
    const v = await open();
    await replay(v, squareDragged);
    expect(include(v)).toEqual(['chr7', 'chr8', 'chr9']);
    await exportedBoth(v);
  });

  it('resized: chr4 and chr5 grown by the corner to chr9', async () => {
    const v = await open();
    await replay(v, squareResized);
    expect(include(v)).toEqual(['chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']);
    await exportedBoth(v);
  });

  it('clicked: chr4 × chr7 at once', async () => {
    const v = await open();
    await replay(v, squareClicked);
    expect(include(v)).toEqual(['chr4', 'chr7']);
    await exportedBoth(v);
  });

  it('a scattered set: chr4, then Cmd+click chr9', async () => {
    const v = await open();
    await replay(v, squareScattered);
    expect(include(v)).toEqual(['chr4', 'chr9']);
    expect(v.view.resolved.value!.sections.find((s) => s.address === RingSection)!.blocks!.length).toBe(2);
    await exportedBoth(v);
  });
});

describe('The square’s moves are steps of the view’s history', () => {
  /** The gesture run, its patches landed, the view drawn. */
  async function settled(v: SeqViewer, gesture: (v: SeqViewer) => unknown) {
    await replay(v, gesture);
    for (let i = 0; i < 10; i++) {
      await vi.advanceTimersByTimeAsync(100);
      v.draw();
    }
  }

  const each = [
    ['a click', squareClicked, ['chr4', 'chr7']],
    ['a drag', squareDragged, ['chr7', 'chr8', 'chr9']],
    ['a resize', squareResized, ['chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']],
  ] as const;

  for (const [name, gesture, shown] of each)
    it(`${name}: Back is possible, shows the blocks it left, Forward again`, async () => {
      const v = await open();
      // the opening's own include (a drag and a resize begin on a run), one step already
      const steps = v.view.history.value.back;
      await settled(v, gesture);
      expect(include(v)).toEqual(shown);
      expect(v.view.canBack).toBe(true);
      // one step for the whole gesture: a drag's frames add none
      const made = v.view.history.value.back - steps;
      expect(made).toBeGreaterThanOrEqual(1);
      expect(made).toBeLessThanOrEqual(2);
      expect(v.view.back()).toBe(true);
      await settled(v, () => undefined);
      expect(include(v)).not.toEqual(shown);
      expect(v.view.canForward).toBe(true);
      expect(v.view.forward()).toBe(true);
      await settled(v, () => undefined);
      expect(include(v)).toEqual(shown);
    });

  it('a drag’s frames push nothing: the step is the release’s', async () => {
    const v = await open();
    await v.execute('section.include', { address: RingSection, include: ['chr4', 'chr5', 'chr6'] });
    await settled(v, () => undefined);
    const before = v.view.history.value.back;
    const n = squareOf(v);
    const from = cellPoint(v, n, 'chr5', 'chr5');
    const to = cellPoint(v, n, 'chr8', 'chr8');
    n.down(from[0], from[1]);
    for (let s = 1; s <= 6; s++) {
      n.move(from[0] + ((to[0] - from[0]) * s) / 6, from[1] + ((to[1] - from[1]) * s) / 6);
      expect(v.view.history.value.back).toBe(before);
    }
    n.up();
    await settled(v, () => undefined);
    expect(v.view.history.value.back).toBe(before + 1);
    expect(include(v)).toEqual(['chr7', 'chr8', 'chr9']);
    expect(v.view.back()).toBe(true);
    await settled(v, () => undefined);
    expect(include(v)).toEqual(['chr4', 'chr5', 'chr6']);
  });
});

describe('The stem menu and the partner keys, replayed', () => {
  it('the strongest stem’s menu open: its partners, select, show only', async () => {
    const v = await open();
    expect(await stemMenuOpen(v)).toBe(true);
    const menu = v.behavior<HeadersBehavior>('headers')!.menu.value!.menu;
    expect(menu.groups[0].map((i) => i.name)).toEqual([
      'rondo:show-partner',
      'rondo:select-stem-partners',
      'rondo:show-chromosome',
    ]);
    await exportedBoth(v);
  });

  it('`.` twice: the cursor on the strongest stem’s second partner', async () => {
    const v = await open();
    const l = ribbonsOf(v)!;
    const { stem } = strongestStem(l);
    partnerKeys(v);
    const second = partnersOf(l, stem)[1];
    const c = v.selection.cursor.value!;
    expect(c.section).toBe(l.section);
    expect(c.system).toBe(l.systems[l.stemSystem[second.stem]]);
    expect(c.position).toBe(l.stemStart[second.stem]);
    await exportedBoth(v);
  });
});
