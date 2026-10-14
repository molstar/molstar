/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureElement } from '@molstar/model/model/structure';
import type { StructureView3D } from '@molstar/seq-molstar-poc/viewer';
import { filter, firstValueFrom } from 'rxjs';
import { describe, expect, it } from 'vitest';
import { linkSequence } from '../link.js';
import { headlessViewers, hostStructure } from './fixtures.js';

/** A 3D view that records what the link asks of it and what it is subscribed to. */
function fakeView() {
  const calls: string[] = [];
  const hover = new Set<(l: StructureElement.Loci | null) => void>();
  const select = new Set<(l: StructureElement.Loci | null) => void>();
  const view: StructureView3D = {
    highlight: (l) => calls.push(`highlight(${l ? 'loci' : 'null'})`),
    select: (l) => calls.push(`select(${l ? 'loci' : 'null'})`),
    focus: () => calls.push('focus'),
    onHover: (fn) => {
      hover.add(fn);
      return () => hover.delete(fn);
    },
    onSelect: (fn) => {
      select.add(fn);
      return () => select.delete(fn);
    },
  };
  return { view, calls, hover, select };
}

const shown = async (h: ReturnType<typeof headlessViewers>) => {
  await firstValueFrom(h.viewer.structure.pipe(filter((s) => s !== undefined)));
  await h.viewer.seq.graph.update();
  h.scheduler.advance(50);
};

describe('linking the sequence to a 3D view', () => {
  it('subscribes when there is a structure, sends a Seq* selection to 3D and focuses it, and drops on cleanup', async () => {
    const h = headlessViewers();
    const fake = fakeView();
    const stop = linkSequence(h.viewer, () => fake.view);
    // Nothing shown yet: nothing subscribed.
    expect(fake.hover.size + fake.select.size).toBe(0);
    h.structure.set(await hostStructure());
    await shown(h);
    expect(fake.hover.size).toBe(1);
    expect(fake.select.size).toBe(1);
    // The replay at link time sends no clear to the 3D view.
    expect(fake.calls).toEqual([]);
    const seq = h.viewer.seq;
    const set = seq.graph.assembled.value?.set;
    const ref = set?.strings.get(set.units.key[0]) ?? '';
    const loci = seq.refs.toLoci(ref, { start: 0, stop: 5 });
    expect(loci).not.toBeNull();
    seq.selection.set(loci, 'pointer');
    expect(fake.calls.slice(0, 2)).toEqual(['select(loci)', 'focus']);
    stop();
    expect(fake.hover.size + fake.select.size).toBe(0);
    h.viewers.dispose();
  });

  it('is remade for another structure, and not made without a 3D view of it', async () => {
    const h = headlessViewers();
    const fake = fakeView();
    let views = 0;
    const stop = linkSequence(h.viewer, () => {
      views++;
      return fake.view;
    });
    h.structure.set(await hostStructure('1tup'));
    await shown(h);
    h.structure.set(await hostStructure('1c26'));
    await shown(h);
    expect(views).toBe(2);
    expect(fake.hover.size).toBe(1);
    stop();
    expect(fake.hover.size).toBe(0);

    const none = linkSequence(h.viewer, () => undefined);
    h.structure.set(await hostStructure('2ocj'));
    await shown(h);
    none();
    h.viewers.dispose();
  });
});
