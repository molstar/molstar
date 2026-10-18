/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import { describe, expect, it, vi } from 'vitest';
import { cbs, crambin, headlessViewer, until } from '../../_test/harness.js';
import { memoryStore } from '../../persist/stores.js';
import type { PrismViewer } from '../../viewer.js';
import { type LociPreview, lociPreviewSink } from '../../sources/select-queries.js';
import { countLine, watchSelectCommands } from '../select.js';

// The Select group (commands.md §5): generated selection-query commands with the group's shared mode, previews and
// counts for the highlighted row; and the keymap's layers in the viewer (commands.md §6, viewer-spec.md §2.5).

const selected = (v: PrismViewer) => {
  let n = 0;
  v.plugin.managers.structure.selection.entries.forEach((e) => {
    n += StructureElement.Loci.size(e.selection);
  });
  return n;
};

describe('the Select group', () => {
  it('lists the registry’s queries by category, and what the loaded structures hold', async () => {
    const v = await headlessViewer();
    const commands = v.managers.commands;
    expect(commands.groupPath('select.query.type').map((g) => g.title)).toEqual(['Select', 'By query', 'Type']);
    expect(commands.get('select.query.type.protein')?.group).toBe('select.query.type');
    // Nothing loaded: the structure's own rows are not there, the registry's say why they cannot run.
    expect(commands.get('select.chain.A')).toBeUndefined();
    expect(commands.blocked(commands.get('select.query.type.protein')!)).toBe('Nothing is loaded.');
    // `Select…` (by text) sits in the group, so its form starts from the shared mode.
    expect(commands.get('select.set')?.group).toBe('select');

    expect((await v.api.load.files({ files: [cbs()] })).ok).toBe(true);
    await until(() => !!commands.get('select.chain.A'));
    expect(commands.get('select.residue.REA')?.group).toBe('select.structure.residues');
    expect(commands.get('select.element.C')?.group).toBe('select.structure.elements');
    expect(commands.children('select.structure.entities').commands.length).toBeGreaterThan(0);
    // Gone with the structure.
    expect((await v.api.session.clear()).ok).toBe(true);
    await until(() => !commands.get('select.chain.A'));
  });

  it('runs with the shared mode, and counts and previews what it would select', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const commands = v.managers.commands;
    await until(() => !!commands.get('select.residue.REA'));
    const ligand = commands.get('select.residue.REA')!;
    const ctx = commands.context.value;
    const count = await ligand.count!({}, ctx, new AbortController().signal);
    expect(count).toMatch(/^1 residue · \d+ atoms$/);
    const preview = ligand.preview!({}, ctx) as LociPreview;
    expect(preview.kind).toBe('loci');
    expect((await preview.resolve()).length).toBe(1);

    const run = commands.open('select.residue.REA');
    expect(run.kind).toBe('ran');
    if (run.kind === 'ran') await run.result;
    const one = selected(v);
    expect(one).toBeGreaterThan(0);
    // Add: the protein joins the ligand; Remove: the ligand goes again.
    commands.setShared('select', 'add');
    const protein = commands.open('select.query.type.protein');
    if (protein.kind === 'ran') expect((await protein.result)?.message).toMatch(/^Added|^Selected/);
    const both = selected(v);
    expect(both).toBeGreaterThan(one);
    commands.setShared('select', 'remove');
    const removed = commands.open('select.residue.REA');
    if (removed.kind === 'ran') await removed.result;
    expect(selected(v)).toBe(both - one);
  });

  it('lights a preview with the hover highlight, and drops one that a newer preview or a clear replaced', async () => {
    const v = await headlessViewer({ files: [crambin()] });
    const highlights = v.plugin.managers.interactivity.lociHighlights;
    const only = vi.spyOn(highlights, 'highlightOnly');
    const clear = vi.spyOn(highlights, 'clearHighlights');
    const sink = lociPreviewSink(v.plugin);
    const structure = v.plugin.managers.structure.hierarchy.current.structures[0]!.cell.obj!.data;
    const loci = StructureElement.Loci.all(structure);
    let release!: () => void;
    const slow: LociPreview = {
      kind: 'loci',
      resolve: () => new Promise((r) => (release = () => r([loci]))),
    };
    sink.set(slow);
    sink.set({ kind: 'loci', resolve: async () => [loci] });
    await until(() => only.mock.calls.length === 1);
    release();
    await new Promise((r) => setTimeout(r, 0));
    expect(only).toHaveBeenCalledTimes(1);
    sink.clear();
    expect(clear).toHaveBeenCalled();
  });

  it('resolves a row again once the structures’ data changed, though the summary did not', async () => {
    const resolve = vi.fn(async () => ({ ok: true, parts: [] }));
    const defs = new Map<string, any>();
    let data: object[] = [{}];
    const entry = { id: 'select.chain.A', title: 'Chain A', group: 'g', icon: 'x', expression: {} };
    const drop = watchSelectCommands({
      queries: { list: () => ({ groups: [], entries: [entry] }), ensure: async () => {} } as any,
      api: { describe: () => ({}) } as any,
      targets: { resolve } as any,
      structures: { subscribe: (f: () => void) => (f(), { unsubscribe() {} }) } as any,
      data: () => data,
      registerLazy: (source) => {
        for (const d of source.make()?.commands ?? []) defs.set(d.id, d);
        return { stale() {}, dispose() {} };
      },
      registerGroups: () => () => {},
    });
    const count = () => defs.get('select.chain.A').count({}, {}, new AbortController().signal);
    await count();
    await count();
    expect(resolve).toHaveBeenCalledTimes(1);
    data = [{}];
    await count();
    expect(resolve).toHaveBeenCalledTimes(2);
    drop();
  });

  it('says counts in a few words', () => {
    expect(countLine(0, 0)).toBe('nothing');
    expect(countLine(1, 20)).toBe('1 residue · 20 atoms');
    expect(countLine(3, 3)).toBe('3 atoms');
    expect(countLine(1200, 9800)).toBe('1,200 residues · 9,800 atoms');
  });
});

describe('keymap in the viewer', () => {
  it('takes the spec’s keys as the app layer and keeps the person’s with the preferences', async () => {
    const preferences = memoryStore({ version: '1', budget: 64_000 });
    const spec = {
      persist: { key: 'keymap-test', preferences },
      ui: { keybindings: [{ command: 'camera.reset', key: { key: 'r', alt: true } }] },
    };
    const v = await headlessViewer({ spec });
    const commands = v.managers.commands;
    const reset = commands.get('camera.reset')!;
    expect(commands.shortcut(reset)).toEqual({ key: 'r', alt: true });
    commands.keymap.remap('camera.reset', { key: 'h', alt: true });
    expect(preferences.get('keybindings')).toEqual({ 'camera.reset': { key: 'h', alt: true } });
    await v.dispose();

    const again = await headlessViewer({ spec });
    expect(again.managers.commands.shortcut(again.managers.commands.get('camera.reset')!)).toEqual({
      key: 'h',
      alt: true,
    });
    again.managers.commands.keymap.resetAll();
    expect(preferences.get('keybindings')).toBeUndefined();
  });
});
