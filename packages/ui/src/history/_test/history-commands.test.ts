/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { historyCommands } from '../history-commands.js';
import { HistoryModel, type HistoryPort } from '../history-model.js';

// Undo and Redo as commands, in plain Node (history §1): `when` says why not, a run says what it did.

function setup() {
  let a = 0;
  const port: HistoryPort<number> = {
    checkpoint: () => a,
    restore: async (c) => {
      a = c;
    },
  };
  const history = new HistoryModel(port);
  const [undo, redo] = historyCommands(history);
  const signal = new AbortController().signal;
  return { history, undo, redo, get: () => a, set: (v: number) => (a = v), signal };
}

describe('historyCommands', () => {
  it('is greyed with nothing to undo or redo; undoes and redoes, naming the step', async () => {
    const h = setup();
    expect(h.undo.when!(undefined)).toBe('Nothing to undo');
    expect(h.redo.when!(undefined)).toBe('Nothing to redo');
    await h.history.run('Set a to 1', { kind: 'person' }, async () => void h.set(1));
    expect(h.undo.when!(undefined)).toBe(true);
    expect(await h.undo.run({}, undefined, h.signal)).toEqual({ message: 'Undone: Set a to 1' });
    expect(h.get()).toBe(0);
    expect(h.redo.when!(undefined)).toBe(true);
    expect(await h.redo.run({}, undefined, h.signal)).toEqual({ message: 'Redone: Set a to 1' });
    expect(h.get()).toBe(1);
  });

  it('a refused undo fails with its reason', async () => {
    const h = setup();
    await expect(h.undo.run({}, undefined, h.signal)).rejects.toThrow('Nothing to undo');
  });

  it('greys both while a run is in flight, Redo too', async () => {
    const h = setup();
    await h.history.run('Set a to 1', { kind: 'person' }, async () => void h.set(1));
    await h.history.undo();
    let open!: () => void;
    const running = h.history.run('Slow', { kind: 'person' }, () => new Promise<void>((r) => (open = r)));
    expect(h.undo.when!(undefined)).toBe('A change is running');
    expect(h.redo.when!(undefined)).toBe('A change is running');
    await new Promise((r) => setTimeout(r, 0));
    open();
    await running;
    expect(h.redo.when!(undefined)).toBe('Nothing to redo');
  });

  it('says what it did in the app’s words', async () => {
    const h = setup();
    const [undo] = historyCommands(h.history, { labels: { undo: 'Annuler', undone: (s) => `Annulé : ${s}` } });
    expect(undo.title).toBe('Annuler');
    await h.history.run('Set a to 1', { kind: 'person' }, async () => void h.set(1));
    expect(await undo.run({}, undefined, h.signal)).toEqual({ message: 'Annulé : Set a to 1' });
  });
});
