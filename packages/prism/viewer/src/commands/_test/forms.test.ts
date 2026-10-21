/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';

// Commands over operations, as the command menu runs them (api.md §9): the target is asked for;
// fields that default (a structure left empty, a language, a mode) do not stop the run.

describe('operation commands', () => {
  it('run with only the target given', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    for (const id of ['focus.set', 'select.set']) {
      const opened = v.managers.commands.open(id, 'ligand');
      expect(opened.kind).toBe('form');
      if (opened.kind !== 'form') continue;
      expect(opened.invocation.state.value.missing).toEqual([]);
      const result = await opened.invocation.confirm();
      expect(result?.message).toMatch(/ligand|residue|atom/i);
    }
    expect(v.plugin.managers.structure.selection.elementCount()).toBeGreaterThan(0);
  });

  it('show the query language only while the target is query text, optional', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const opened = v.managers.commands.open('focus.set', 'ligand');
    if (opened.kind !== 'form') throw new Error(opened.kind);
    const fields = () => opened.invocation.state.value.fields.map((f) => f.id);
    expect(fields()).not.toContain('language');
    opened.invocation.set('target', 'resn REA');
    expect(fields()).toContain('language');
    expect(opened.invocation.state.value.missing).toEqual([]);
    opened.invocation.set('target', 'A:10-20');
    expect(fields()).not.toContain('language');
  });

  it('read Show as… and Color by theme… as sentences: the subject, its structure, then the choice', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const read = (id: string) => {
      const opened = v.managers.commands.open(id);
      if (opened.kind !== 'form') throw new Error(opened.kind);
      const s = opened.invocation.state.value;
      v.managers.commands.close();
      return s.fields.map((f) => (f.before ? `${f.before} ${f.id}` : f.id));
    };
    expect(read('representation.set')).toEqual(['component', 'in structure', 'as type']);
    expect(read('representation.color-theme')).toEqual(['target', 'in structure', 'by theme']);
  });

  it('ask for the target when it is missing', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const opened = v.managers.commands.open('focus.set');
    if (opened.kind !== 'form') throw new Error(opened.kind);
    expect(opened.invocation.state.value.missing).toEqual(['target']);
  });
});
