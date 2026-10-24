/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { SecondaryStructureProvider as ComputedSecondaryStructure } from '@molstar/model/props/computed/secondary-structure';
import { inkForColor } from '@molstar/prism-ui';
import { filter, firstValueFrom } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { ComplexRef, sequenceDocument } from '../document.js';
import { DefaultTracksRail } from '../presets.js';
import { SecondaryTrackProvider } from '../tracks.js';
import { applySaved } from '../runtime.js';
import { skinFor } from '../skin.js';
import type { SequenceBehaviorContext } from '../behaviors.js';
import { headlessViewers, hostStructure, ValuesTrack, valuesBehavior } from './fixtures.js';

async function until(check: () => boolean, ms = 3000): Promise<void> {
  const t0 = Date.now();
  while (!check()) {
    if (Date.now() - t0 > ms) throw new Error(`timed out after ${ms} ms`);
    await new Promise((r) => setTimeout(r, 5));
  }
}

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

  it('switches the style on one document (the same object, secondary structure in both), and keeps it across a new theme', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const doc = h.viewer.seq.document.shown.value;
    const set = h.viewer.seq.graph.assembled.value?.set;
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('prism_tracks');
    expect(JSON.stringify(doc)).toContain('molstar:secondary_structure');
    h.viewer.setStyle('compact');
    expect(h.viewer.style.value).toBe('compact');
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('compact');
    await settled(h);
    expect(h.viewer.seq.document.shown.value).toBe(doc);
    expect(h.viewer.seq.graph.assembled.value?.set).toBe(set);
    expect(h.viewer.seq.layout.presetOf(0)?.geometry).toBe('wrapped');
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('compact');
    h.viewer.setStyle('tracks');
    await settled(h);
    expect(h.viewer.seq.ui.style.value.get('main')).toBe('prism_tracks');
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('molstar:secondary_structure');
    h.viewers.dispose();
  });

  it("the text style's render props are the model's: centered, colored, the letters' size", async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    h.viewer.setStyle('compact');
    h.scheduler.advance(50);
    const props = h.viewer.seq.layout.props('main');
    expect(props.wrappedAlign).toBe('center');
    expect(props.wrappedColors).toBe(true);
    expect(props.unobservedGray).toBe(true);
    h.viewers.dispose();
  });

  it('draws no colors as plain letters: cells in the surface color, the text in the ink', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    h.viewer.setStyle('compact');
    h.viewer.setTheme({ letters: null });
    await settled(h);
    h.scheduler.advance(50);
    const shown = JSON.stringify(h.viewer.seq.document.shown.value);
    expect(shown).toContain('"theme":"uniform"');
    expect(shown).not.toContain('prism:residues');
    expect(h.viewer.seq.layout.props('main').wrappedColors).toBe(false);
    // Colors again: the text is tinted again.
    h.viewer.setTheme({ letters: 'taylor' });
    await settled(h);
    h.scheduler.advance(50);
    expect(h.viewer.seq.layout.props('main').wrappedColors).toBe(true);
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

  it("Shift adds or removes each stretch (Seq*'s default `shiftPositions: 'toggle'`)", () => {
    const h = headlessViewers();
    expect(h.viewer.seq.behaviors.get<{ shiftPositions: string }>('select')?.params.value.shiftPositions).toBe(
      'toggle',
    );
    h.viewers.dispose();
  });

  it("a behavior adds its track to each document by patch, its legend block and its track's menu", async () => {
    const patched: boolean[] = [];
    const h = headlessViewers({ behaviors: [valuesBehavior(patched)] });
    const host = await hostStructure();
    h.structure.set(host);
    await settled(h);
    await until(() => patched.length === 1);
    expect(patched).toEqual([true]);
    await h.viewer.seq.graph.update();
    // The recorded entry is no prediction, but its B-factors are read as the confidence: the track is there, with no issue.
    expect(h.viewer.seq.issues.current.value).toEqual([]);
    const doc = JSON.stringify(h.viewer.seq.document.shown.value);
    expect(doc).toContain(`"ref":"${ValuesTrack}"`);
    expect(doc.indexOf('"ref":"ss-track"')).toBeLessThan(doc.indexOf(`"ref":"${ValuesTrack}"`));
    const sections = h.viewer.legend.value.sections;
    expect(sections?.map((x) => x.title)).toEqual(['Secondary structure', 'Values']);
    expect(sections?.[0].items.map((i) => i.label)).toEqual(['Helix', 'Strand', 'Turn']);
    // Another structure: a new document, and the behavior adds its track to it again.
    h.structure.set(await hostStructure('2ocj'));
    await settled(h);
    await until(() => patched.length === 2);
    expect(patched).toEqual([true, true]);
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain(`"ref":"${ValuesTrack}"`);
    h.viewers.dispose();
  });

  it('drops a patch made for a document that was replaced, and a removed behavior withdraws its offer', async () => {
    const h = headlessViewers();
    const host = await hostStructure();
    let ctx: SequenceBehaviorContext | undefined;
    const remove = h.viewer.addBehavior({
      id: 'test.sequence.probe',
      title: 'Probe',
      attach: (c) => ((ctx = c), () => {}),
    });
    expect(() => h.viewer.addBehavior({ id: 'test.sequence.probe', title: 'Again', attach: () => () => {} })).toThrow();
    h.structure.set(host);
    await settled(h);
    const first = ctx!.shown.value!;
    ctx!.offer({ legend: [{ title: 'Probe', items: [] }] });
    expect(h.viewer.legend.value.sections?.map((x) => x.title)).toEqual(['Secondary structure', 'Probe']);
    h.structure.set(await hostStructure('2ocj'));
    await settled(h);
    expect(ctx!.shown.value).not.toBe(first);
    expect(await ctx!.patch(first, [{ op: 'remove', ref: 'letters' }])).toBe(false);
    expect(JSON.stringify(h.viewer.seq.document.shown.value)).toContain('"ref":"letters"');
    remove();
    expect(h.viewer.legend.value.sections?.map((x) => x.title)).toEqual(['Secondary structure']);
    h.viewers.dispose();
  });

  it('the tracks have their menus: the model says what each offers, they run on the 3D view', async () => {
    const calls: string[] = [];
    const h = headlessViewers({
      actions: { colorStructure: (c) => void calls.push(c) },
      behaviors: [valuesBehavior()],
    });
    const host = await hostStructure();
    h.structure.set(host);
    await settled(h);
    // Asked through Seq*'s `interaction.menu` (sequence-control.tsx): no behavior of Prism's in the spec for it.
    expect(h.viewer.seq.behaviors.get('contextMenu')).toBeNull();
    const menu = (track: string) => h.viewer.trackActions(track).map((a) => a.label);
    expect(menu('letters')).toEqual(['Apply default coloring']);
    expect(menu('ss-track')).toEqual(['Color by secondary structure']);
    expect(menu(ValuesTrack)).toEqual(['Color by values']);
    expect(menu('ruler')).toEqual([]);
    for (const t of ['letters', 'ss-track', ValuesTrack]) h.viewer.trackActions(t)[0].run();
    expect(calls).toEqual(['default', 'secondary-structure', 'plddt-confidence']);
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

  it('announces every document it builds, though the structure and the complex ref are the same (the 3D link is remade)', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await settled(h);
    const before = h.viewer.structure.value;
    let emitted = 0;
    const sub = h.viewer.structure.subscribe(() => emitted++);
    emitted = 0;
    // A new document over the same structure (the text version's colors).
    h.viewer.setTheme({ letters: null });
    await until(() => h.viewer.structure.value !== before);
    expect(emitted).toBeGreaterThan(0);
    expect(h.viewer.structure.value?.structure).toBe(before?.structure);
    expect(h.viewer.structure.value?.complex).toBe(before?.complex);
    sub.unsubscribe();
    h.viewers.dispose();
  });

  it('restores the default rail with the default state (a snapshot without the sequence, Undo to the start)', () => {
    const h = headlessViewers();
    h.viewer.setRail({ navigator: 'strip', height: 30 });
    applySaved(h.viewer, () => {}, undefined, undefined);
    expect(h.viewer.rail.value).toEqual(DefaultTracksRail);
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
    const base = { title: 'T', key: host.key, format: host.source.format, tracks: [] };
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

  it("follows the 3D view's secondary structure settings: the track's data names them, and a change recomputes it", async () => {
    const h = headlessViewers();
    const host = await hostStructure('2ocj');
    host.source.drawn = host.structure;
    h.structure.set(host);
    await settled(h);
    const rows = () => {
      const set = h.viewer.seq.graph.assembled.value?.set;
      const a = set?.annotations.find((x) => x && set.strings.get(x.key).endsWith('ss'));
      return a?.count ?? -1;
    };
    // Automatic: 2OCJ's own assignment, bends and all.
    const own = rows();
    expect(own).toBeGreaterThan(0);
    // The 3D view set to DSSP (its custom property's settings), then the host tells the same structure with them.
    await ComputedSecondaryStructure.attach(
      { runtime: SyncRuntimeContext, assetManager: new DefaultAssetManager() },
      host.structure.root,
      { type: { name: 'dssp', params: {} } } as never,
    );
    h.structure.set({ ...host, assignment: 'dssp' });
    await until(() => JSON.stringify(h.viewer.seq.document.shown.value).includes('dssp'));
    await settled(h);
    const dssp = rows();
    expect(dssp).toBeGreaterThan(0);
    expect(dssp).not.toBe(own);
    // Told again with the settings as they are: no patch (the document stays the same object).
    const doc = h.viewer.seq.document.shown.value;
    h.structure.set({ ...host, assignment: 'dssp, again' });
    await settled(h);
    expect(h.viewer.seq.document.shown.value).toBe(doc);
    // A frame (the same key, the same settings): no sync at all, so no provider builds (a frame comes at animation
    // rate).
    const build = vi.spyOn(SecondaryTrackProvider, 'build');
    h.structure.set({ ...host, assignment: 'dssp, again', title: 'a frame' });
    await settled(h);
    expect(build).not.toHaveBeenCalled();
    h.viewers.dispose();
  });
});
