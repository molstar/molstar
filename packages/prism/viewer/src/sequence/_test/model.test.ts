/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { inkForColor } from '@molstar/prism-ui';
import { filter, firstValueFrom } from 'rxjs';
import { describe, expect, it } from 'vitest';
import { skinFor } from '../skin.js';
import { headlessViewers, hostStructure } from './fixtures.js';

/** Waits for the document to be in (the viewer says what it shows), then lets the frames run. */
const settled = async (h: ReturnType<typeof headlessViewers>) => {
  await Promise.race([
    firstValueFrom(h.viewer.structure.pipe(filter((s) => s !== undefined))),
    new Promise((_, reject) => setTimeout(() => reject(new Error(JSON.stringify(h.viewer.seq.issues))), 5000)),
  ]);
  await h.viewer.seq.graph.update();
  h.scheduler.advance(50);
};

describe('the sequence model', () => {
  it('creates the default viewer, empty until the host shows a structure', async () => {
    const h = headlessViewers();
    expect(h.viewers.ids.value).toEqual(['default']);
    expect(h.viewer.id).toBe('default');
    expect(h.viewer.structure.value).toBeUndefined();
    expect(h.viewer.style.value).toBe('compact');
    expect(() => h.viewers.get('other')).toThrow();
    h.viewers.dispose();
    expect(h.viewers.ids.value).toEqual([]);
  });

  it('builds one view with one section over the structure, and says what it shows', async () => {
    const h = headlessViewers();
    const host = await hostStructure();
    h.structure.set(host);
    await settled(h);
    const shown = h.viewer.structure.value;
    expect(shown?.title).toBe('1TUP');
    expect(shown?.structure).toBe(host.structure);
    const set = h.viewer.seq.graph.assembled.value?.set;
    expect(set?.units.rangeStart.length).toBeGreaterThan(0);
    expect(h.viewer.seq.viewports.sections.value).toEqual(['main']);
    h.viewers.dispose();
  });

  it("reads the host's parse and refuses the network: no issue for the structure, a loaded file would fail", async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    expect(h.viewer.seq.issues.current.value.filter((i) => i.severity === 'error')).toEqual([]);
    const set = h.viewer.seq.graph.assembled.value?.set;
    // 1TUP: p53 chains and DNA chains are units of the complex.
    expect(set?.units.rangeStart.length).toBeGreaterThanOrEqual(3);
    h.viewers.dispose();
  });

  it('switches the style by UI state alone, and keeps it across a new document', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const doc = h.viewer.seq.document.shown.value;
    h.viewer.setStyle('tracks');
    expect(h.viewer.style.value).toBe('tracks');
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('tracks');
    expect(h.viewer.seq.document.shown.value).toBe(doc);
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('tracks');
    h.viewers.dispose();
  });

  it('writes the theme into the document and keeps the set (so the link stays valid)', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const set = h.viewer.seq.graph.assembled.value?.set;
    expect(h.viewer.theme.value.letters).toBe('zappo');
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('zappo');
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    expect(h.viewer.theme.value.letters).toBe('taylor');
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('taylor');
    expect(h.viewer.seq.graph.assembled.value?.set).toBe(set);
    expect(h.viewer.legend.value.title).toBe('Taylor');
    h.viewers.dispose();
  });

  it('follows the ink: the skin is the ink, no CSS read', async () => {
    const h = headlessViewers('#ffffff');
    expect(h.viewer.seq.skin.value).toEqual(skinFor(inkForColor('#ffffff')));
    expect(h.viewer.seq.skin.value.dark).toBe(false);
    h.ink.set(inkForColor('#111318'));
    expect(h.viewer.seq.skin.value).toEqual(skinFor(inkForColor('#111318')));
    expect(h.viewer.seq.skin.value.dark).toBe(true);
    expect(h.viewer.seq.skin.value.surface).toBe('rgb(17 19 24)');
    h.viewer.setFonts({ mono: 'M', ui: 'U' });
    expect(h.viewer.seq.skin.value.fontMono).toBe('M');
    h.viewers.dispose();
  });

  it('shows a second structure in place of the first', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure('1tup'));
    await settled(h);
    const other = await hostStructure('1c26');
    h.structure.set(other);
    await settled(h);
    expect(h.viewer.structure.value?.title).toBe('1C26');
    expect(h.viewer.structure.value?.structure).toBe(other.structure);
    h.structure.set(undefined);
    expect(h.viewer.structure.value).toBeUndefined();
    h.viewers.dispose();
  });

  it('dispose leaves nothing running', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    h.viewers.dispose();
    expect(h.viewer.seq.disposed).toBe(true);
    expect(h.scheduler.pendingFrames).toBe(0);
    expect(h.scheduler.pendingTimers).toBe(0);
    // Later changes of what it followed reach nothing.
    h.ink.set(inkForColor('#000000'));
    h.structure.set(undefined);
    expect(h.viewer.seq.skin.value.surface).not.toBe('rgb(0 0 0)');
  });
});
