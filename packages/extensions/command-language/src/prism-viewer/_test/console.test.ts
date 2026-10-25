/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { PrismStore } from '@molstar/prism-viewer/spec';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { CommandLanguageExtension } from '../../prism-viewer.js';
import { completeLine } from '../console.js';
import { ConsoleModel, HistoryPreference } from '../console-model.js';

// The console's model on a headless viewer (_spec §9): submissions and what they said, the history, script files with
// the files they load, `undo` alone. The view is a thin binding of this.

const cif = readFileSync(new URL('../../../../../../data/examples/1cbs_updated.cif', import.meta.url));
const structureFile = () => new File([cif], '1cbs.cif');

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer() {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), CommandLanguageExtension] },
  });
  viewers.push(v);
  return v;
}

function store(): PrismStore & { readonly data: Map<string, unknown> } {
  const data = new Map<string, unknown>();
  return { data, get: (k) => data.get(k), set: (k, v) => void data.set(k, v), remove: (k) => void data.delete(k) };
}

describe('the console', () => {
  it('runs submissions, says what they did, and keeps the history', async () => {
    const v = await viewer();
    const prefs = store();
    const c = new ConsoleModel(v.api, prefs);
    await c.openFiles([
      new File(['load 1cbs.cif\nshow sticks, resi 10-20\ncount_atoms resn REA'], 'setup.pml'),
      structureFile(),
    ]);
    const [loaded] = c.entries.value;
    expect(loaded, JSON.stringify(loaded)).toMatchObject({ input: 'setup.pml', state: 'done', dialect: 'pymol' });
    expect(loaded.lines.map((l) => l.text)).toContain('22');
    await c.submit('color red, resn REA');
    expect(c.entries.value[1].state).toBe('done');
    await c.submit('shw sticks\ncolor notacolor');
    const failed = c.entries.value[2];
    expect(failed.state).toBe('failed');
    expect(failed.lines.map((l) => l.text)).toEqual([
      expect.stringContaining("line 1: PyMOL has no command 'shw'"),
      expect.stringContaining("line 2: 'notacolor' is not a PyMOL color"),
    ]);
    expect(prefs.data.get(HistoryPreference)).toEqual(['color red, resn REA', 'shw sticks\ncolor notacolor']);
  });

  it('walks the history with ↑ and ↓, back to the draft', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api, store());
    await c.submit('get_names');
    await c.submit('help');
    c.setDraft('half typed');
    c.previous();
    expect(c.draft.value).toBe('help');
    c.previous();
    expect(c.draft.value).toBe('get_names');
    c.previous();
    expect(c.draft.value).toBe('get_names');
    c.next();
    c.next();
    expect(c.draft.value).toBe('half typed');
  });

  it('runs a line in the dialect its first word names, and undo alone as the viewer’s', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api, store());
    await c.openFiles([new File(['load 1cbs.cif'], 'a.pml'), structureFile()]);
    const steps = v.api.history.steps.value.length;
    await c.submit('vmd mol color Chain');
    expect(c.entries.value[1]).toMatchObject({ dialect: 'vmd', state: 'done' });
    await c.submit('undo');
    expect(c.entries.value[2].state).toBe('done');
    expect(v.api.history.steps.value.length).toBeLessThanOrEqual(steps + 1);
  });

  it('loads files dropped without a script, as the line that loads them', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    await c.openFiles([structureFile()]);
    expect(c.entries.value[0]).toMatchObject({ input: 'load 1cbs.cif', state: 'done' });
    expect(v.plugin.managers.structure.hierarchy.current.structures).toHaveLength(1);
    c.setDialect('vmd');
    await c.openFiles([structureFile()]);
    expect(c.entries.value[1]).toMatchObject({ input: 'mol new 1cbs.cif', state: 'done' });
  });

  it('completes the line with Tab’s first completion', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    expect(completeLine(c, 'show stic')).toBe('show sticks, ');
    expect(completeLine(c, 'vmd mol addr')).toBe('vmd mol addrep ');
  });
});

describe('the console, from a second review', () => {
  it('runs the dropped scripts no other one includes, once', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    await c.openFiles([
      new File(['load 1cbs.cif\n@helper.pml'], 'main.pml'),
      new File(['show sticks, resn REA'], 'helper.pml'),
      structureFile(),
    ]);
    expect(c.entries.value.map((e) => [e.input, e.state])).toEqual([['main.pml', 'done']]);
  });

  it('names the included file and its line when a command there fails as it runs', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    await c.openFiles([
      new File(['load 1cbs.cif\n\n@helper.pml'], 'main.pml'),
      new File(['show sticks, resn REA\n\ncolor red, resn XYZ'], 'helper.pml'),
      structureFile(),
    ]);
    expect(c.entries.value[0].lines[0].text).toContain('helper.pml line 3');
  });

  it('runs a submission made while one runs after it, on the scene it leaves', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    const first = c.submit('load 1cbs.cif', [structureFile()]);
    const second = c.submit('color red, 1cbs');
    await Promise.all([first, second]);
    expect(c.entries.value.map((e) => e.state)).toEqual(['done', 'done']);
  });
});

describe('editing an entry', () => {
  it('puts its line back in its own dialect; a script file is not edited', async () => {
    const v = await viewer();
    const c = new ConsoleModel(v.api);
    await c.openFiles([new File(['load 1cbs.cif'], 'a.pml'), structureFile()]);
    c.setDialect('vmd');
    await c.submit('mol list');
    c.setDialect('pymol');
    const [file, line] = c.entries.value;
    expect(c.editable(file)).toBe(false);
    c.edit(line);
    expect(c.dialect.value).toBe('vmd');
    expect(c.draft.value).toBe('mol list');
  });
});
