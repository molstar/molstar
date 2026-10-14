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
import { ComplexRef, sequenceDocument } from '../document.js';
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
    // The tracks are the default; the text is a click away.
    expect(h.viewer.style.value).toBe('tracks');
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
    // The theme applies (a scheme that did not fit the first block would be reported and drawn as the letters' own).
    expect(h.viewer.seq.issues.current.value).toEqual([]);
    const set = h.viewer.seq.graph.assembled.value?.set;
    // 1TUP: p53 chains and DNA chains are units of the complex.
    expect(set?.units.rangeStart.length).toBeGreaterThanOrEqual(3);
    h.viewers.dispose();
  });

  it('switches the style, the text with only the letters (no lane for the tracks), and keeps it across a new document', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const doc = h.viewer.seq.document.shown.value;
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('prism_tracks');
    h.viewer.setStyle('compact');
    expect(h.viewer.style.value).toBe('compact');
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('compact');
    await settled(h);
    const text = JSON.stringify(h.viewer.seq.document.shown.value);
    expect(text).not.toContain('molstar:secondary_structure');
    expect(h.viewer.seq.document.shown.value).not.toBe(doc);
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('compact');
    h.viewer.setStyle('tracks');
    await settled(h);
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('molstar:secondary_structure');
    h.viewers.dispose();
  });

  it('writes the theme into the document and keeps the set (so the link stays valid)', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const set = h.viewer.seq.graph.assembled.value?.set;
    expect(h.viewer.theme.value.letters).toBe('prism:residues');
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('prism:residues');
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    expect(h.viewer.theme.value.letters).toBe('taylor');
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('taylor');
    expect(h.viewer.seq.graph.assembled.value?.set).toBe(set);
    expect(h.viewer.legend.value.title).toBe('Taylor');
    h.viewers.dispose();
  });

  it("the tracks version has no header column; its rail is the host's choice", async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    h.viewer.seq.viewports.get('main')?.setSize(1000, 200, 1);
    const layout = () => h.viewer.seq.layout.get('main')?.current;
    const props = () => h.viewer.seq.layout.props('main');
    // Tracks is the default.
    h.scheduler.advance(50);
    expect(h.viewer.seq.layout.presetOf(0)?.geometry).toBe('linear');
    expect(props().headerWidth).toBe(0);
    expect(layout()?.navHeight).toBe(0);
    h.viewer.setRail({ navigator: 'strip', height: 30 });
    h.scheduler.advance(50);
    expect(h.viewer.rail.value).toEqual({ navigator: 'strip', height: 30 });
    expect(layout()?.navHeight).toBe(30);
    h.viewer.setRail({ navigator: 'none' });
    h.scheduler.advance(50);
    expect(layout()?.navHeight).toBe(0);
    // The text version: wrapped, centered by the model's render props.
    h.viewer.setStyle('compact');
    h.scheduler.advance(50);
    expect(h.viewer.seq.layout.presetOf(0)?.geometry).toBe('wrapped');
    expect(props().wrappedAlign).toBe('center');
    h.viewers.dispose();
  });

  it("Shift adds or removes each stretch (`shiftPositions: 'toggle'`), not Seq*'s own `extend`", () => {
    const h = headlessViewers();
    expect(h.viewer.seq.behaviors.get<{ shiftPositions: string }>('select')?.params.value.shiftPositions).toBe(
      'toggle',
    );
    h.viewers.dispose();
  });

  it('a predicted model gets a pLDDT track in the AlphaFold bands, and the legend lists the four bands', async () => {
    const h = headlessViewers();
    const host = await hostStructure();
    expect(h.viewer.legend.value.sections).toBeUndefined();
    h.structure.set({ ...host, predicted: true });
    await settled(h);
    // The recorded entry is no prediction, but its B-factors are read as the confidence: the track is there, with no issue.
    expect(h.viewer.seq.issues.current.value).toEqual([]);
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('molstar:plddt');
    const sections = h.viewer.legend.value.sections;
    expect(sections?.map((x) => x.title)).toEqual(['Secondary structure', 'Model confidence']);
    expect(sections?.[0].items.map((i) => i.label)).toEqual(['Helix', 'Strand', 'Turn']);
    const bands = sections?.[1];
    expect(bands?.items.map((i) => i.label)).toEqual([
      'Very high (pLDDT > 90)',
      'Confident (70–90)',
      'Low (50–70)',
      'Very low (< 50)',
    ]);
    // An experimental structure after it: the secondary structure, no bands.
    h.structure.set(host);
    await settled(h);
    expect(h.viewer.legend.value.sections?.map((x) => x.title)).toEqual(['Secondary structure']);
    h.viewers.dispose();
  });

  it('the tracks have their menus: the model says what each offers, they run on the 3D view', async () => {
    const calls: string[] = [];
    const h = headlessViewers({ actions: { colorStructure: (c) => void calls.push(c) } });
    const host = await hostStructure();
    h.structure.set({ ...host, predicted: true });
    await settled(h);
    const menu = (track: string) => h.viewer.trackActions(track).map((a) => a.label);
    expect(menu('letters')).toEqual(['Apply default coloring']);
    expect(menu('ss-track')).toEqual(['Color by secondary structure']);
    expect(menu('plddt-track')).toEqual(['Color by pLDDT']);
    expect(menu('ruler')).toEqual([]);
    for (const t of ['letters', 'ss-track', 'plddt-track']) h.viewer.trackActions(t)[0].run();
    expect(calls).toEqual(['default', 'secondary', 'plddt']);
    // The tracks exist in the document under those refs.
    const refs = JSON.stringify(h.viewer.seq.document.shown.value);
    for (const t of ['letters', 'ss-track', 'plddt-track']) expect(refs).toContain(`"ref":"${t}"`);
    h.viewers.dispose();
    // Without a host's actions no track has a menu.
    const bare = headlessViewers();
    expect(bare.viewer.trackActions('letters')).toEqual([]);
    bare.viewers.dispose();
  });

  it('follows the ink: the skin is the ink, no CSS read', async () => {
    const h = headlessViewers();
    expect(h.viewer.seq.skin.value).toEqual(skinFor(inkForColor('#ffffff')));
    expect(h.viewer.seq.skin.value.dark).toBe(false);
    h.ink.set(inkForColor('#111318'));
    expect(h.viewer.seq.skin.value).toEqual(skinFor(inkForColor('#111318')));
    expect(h.viewer.seq.skin.value.dark).toBe(true);
    expect(h.viewer.seq.skin.value.surface).toBe('rgb(17 19 24)');
    h.viewer.setFonts({ mono: 'M', ui: 'U', size: 13 });
    expect(h.viewer.seq.skin.value.fontMono).toBe('M');
    // The text version's letters: at most the panel text size, from the host's tokens; the line the whole width.
    const props = h.viewer.seq.renderProps.value;
    expect([
      props.wrappedFontPx,
      props.maxLetterPx,
      props.wrappedMaxWidth,
      props.wrappedColors,
      props.unobservedGray,
    ]).toEqual([13, 13, 0, true, true]);
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

  it('the document lists only the ligands the 3D view draws (the same list), and every one when none is given', async () => {
    const host = await hostStructure();
    const base = { title: 'T', key: host.key, format: host.source.format, theme: { letters: null } };
    const options = (ligands?: string[]) => {
      const text = JSON.stringify(sequenceDocument({ ...base, ...(ligands ? { ligands } : {}) }));
      const node = (JSON.parse(text) as { root: unknown }).root;
      const found: { ref?: string; params?: { options?: { ligand_comps?: string[] } } }[] = [];
      const walk = (n: { ref?: string; children?: unknown[] }) => {
        if (n.ref === ComplexRef) found.push(n);
        for (const c of n.children ?? []) walk(c as never);
      };
      walk(node as never);
      return found[0]?.params?.options;
    };
    expect(options()).toBeUndefined();
    expect(options(['HEM', 'ERY'])?.ligand_comps).toEqual(['HEM', 'ERY']);
  });
});
