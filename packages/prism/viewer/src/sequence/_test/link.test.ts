/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '@molstar/model/model/structure';
import { observableState } from '@molstar/prism-ui/core';
import { filter, firstValueFrom } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import type { Mark, MarkHover, MarksSource } from '../../sources/marks.js';
import { linkSequence, type SequenceView3D } from '../link.js';
import { headlessViewers, hostStructure } from './fixtures.js';

/** A 3D view that records what the link asks of it and what it is subscribed to. */
function fakeView() {
  const calls: string[] = [];
  const hover = new Set<(l: StructureElement.Loci | null) => void>();
  const select = new Set<(l: StructureElement.Loci | null) => void>();
  const view: SequenceView3D = {
    highlight: (l) => calls.push(`highlight(${l ? 'loci' : 'null'})`),
    select: (l) => calls.push(`select(${l ? 'loci' : 'null'})`),
    focus: () => calls.push('focus'),
    resetCamera: () => calls.push('reset'),
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
    fake.calls.length = 0;
    stop();
    expect(fake.hover.size + fake.select.size).toBe(0);
    // What the sequence lit in 3D goes with the link (and the host's hover state with it).
    expect(fake.calls).toEqual(['highlight(null)']);
    h.viewers.dispose();
  });

  it('follows a frame step without clearing the 3D hover or sending the selection again (perf-audit-2 10/B6)', async () => {
    const h = headlessViewers();
    const fake = fakeView();
    const views: Structure[] = [];
    const stop = linkSequence(h.viewer, (s) => {
      views.push(s);
      return fake.view;
    });
    const host = await hostStructure();
    h.structure.set(host);
    await shown(h);
    const seq = h.viewer.seq;
    const set = seq.graph.assembled.value?.set;
    const ref = set?.strings.get(set.units.key[0]) ?? '';
    seq.selection.set(seq.refs.toLoci(ref, { start: 0, stop: 5 }), 'pointer');
    expect(fake.calls.slice(0, 2)).toEqual(['select(loci)', 'focus']);
    fake.calls.length = 0;

    // a frame: a new Structure of the same atoms, under the same key
    for (let i = 0; i < 3; i++) {
      const frame = Structure.create(host.structure.units);
      h.structure.set({ ...host, structure: frame });
      expect(h.viewer.frame.value).toBe(frame);
    }
    expect(views.length).toBe(4);
    expect(fake.calls).toEqual([]);
    expect(fake.select.size).toBe(1);

    // the views of the new frame still send what the person does
    seq.selection.set(null, 'pointer');
    expect(fake.calls[0]).toBe('select(null)');
    stop();
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

  it('resets the camera once when the person empties the selection in the sequence, and not for an emptying from 3D', async () => {
    const h = headlessViewers();
    const fake = fakeView();
    const stop = linkSequence(h.viewer, () => fake.view);
    h.structure.set(await hostStructure());
    await shown(h);
    const seq = h.viewer.seq;
    const set = seq.graph.assembled.value?.set;
    const ref = set?.strings.get(set.units.key[0]) ?? '';
    const loci = () => seq.refs.toLoci(ref, { start: 0, stop: 5 });
    const resets = () => fake.calls.filter((c) => c === 'reset').length;

    // Selecting, and emptying an empty selection: no reset.
    seq.selection.clear('pointer');
    seq.selection.set(loci(), 'pointer');
    expect(resets()).toBe(0);
    // Emptied by a click in the sequence: one reset, after the 3D selection is cleared.
    seq.selection.set(null, 'pointer');
    expect(resets()).toBe(1);
    expect(fake.calls.slice(-2)).toEqual(['select(null)', 'reset']);
    // Empty again: nothing more.
    seq.selection.clear('pointer');
    expect(resets()).toBe(1);

    // Emptied from 3D (the selection bar, a click in 3D): the camera stays.
    seq.selection.set(loci(), 'pointer');
    seq.selection.set(null, 'host:3d');
    expect(resets()).toBe(1);
    // A data edit's carry or the API: the same.
    seq.selection.set(loci(), 'pointer');
    seq.selection.set(null, 'api');
    expect(resets()).toBe(1);
    stop();
    h.viewers.dispose();
  });

  it('reveals what a link channel lights after the dwell (the revealHighlight behavior), and not a channel it does not follow', async () => {
    const h = headlessViewers();
    h.structure.set(await hostStructure());
    await shown(h);
    const seq = h.viewer.seq;
    const reveal = vi.spyOn(seq, 'reveal');
    const set = seq.graph.assembled.value?.set;
    const ref = set?.strings.get(set.units.key[0]) ?? '';
    const loci = seq.refs.toLoci(ref, { start: 0, stop: 5 });
    expect(loci).not.toBeNull();
    seq.highlight.set('other:x', loci);
    h.scheduler.advance(500);
    expect(reveal).not.toHaveBeenCalled();
    seq.highlight.set('links:3d', loci);
    h.scheduler.advance(100);
    expect(reveal).not.toHaveBeenCalled();
    h.scheduler.advance(200);
    expect(reveal).toHaveBeenCalledTimes(1);
    h.viewers.dispose();
  });

  it("shows a view's marks as highlight channels in their tones and its hover as the 3D view's, and clears them on cleanup", async () => {
    const h = headlessViewers();
    const fake = fakeView();
    const list = observableState<readonly Mark[]>([]);
    const hover = observableState<MarkHover | undefined>(undefined);
    const marks = {
      marks: list.asReadonly,
      hover: hover.asReadonly,
    } as MarksSource;
    const stop = linkSequence(h.viewer, () => fake.view, { marks });
    const host = await hostStructure();
    h.structure.set(host);
    await shown(h);
    const seq = h.viewer.seq;
    const setHover = vi.spyOn(seq.hover, 'setLoci');
    const all = Structure.toStructureElementLoci(host.structure);
    const mark = (id: string, tone: Mark['tone']): Mark => ({
      owner: 'plot',
      id,
      label: id,
      tone,
      structure: 's',
      loci: all,
    });

    list.set([mark('a', 'highlight'), mark('b', 'accent')]);
    const channels = () => [...seq.highlight.channels.value].filter(([c]) => c.startsWith('marks:'));
    expect(channels().map(([c, e]) => [c, e.style])).toEqual([
      ['marks:plot:a', 'highlight'],
      ['marks:plot:b', 'accent'],
    ]);
    list.set([mark('b', 'accent')]);
    expect(channels().map(([c]) => c)).toEqual(['marks:plot:b']);

    hover.set({ owner: 'plot', structure: 's', loci: all });
    expect(setHover).toHaveBeenLastCalledWith(expect.anything(), '3d');
    // The hover is not sent back to the 3D view, which lit it.
    expect(fake.calls).toEqual([]);
    hover.set(undefined);
    expect(setHover).toHaveBeenLastCalledWith(null, '3d');

    // Another model's loci are not this one's residues (they would remap by unit id and atom index).
    const other = Structure.toStructureElementLoci((await hostStructure('1c26')).structure);
    list.set([{ ...mark('c', 'accent'), loci: other }]);
    expect(channels()).toEqual([]);
    hover.set({ owner: 'plot', structure: 's', loci: other });
    expect(setHover).toHaveBeenLastCalledWith(null, '3d');
    hover.set(undefined);

    hover.set({ owner: 'plot', structure: 's', loci: all });
    stop();
    expect(channels()).toEqual([]);
    expect(setHover).toHaveBeenLastCalledWith(null, '3d');
    h.viewers.dispose();
  });

  it("focuses the sequence on an owner's new marks together when they ask for it (a box's two ranges)", async () => {
    const h = headlessViewers();
    const fake = fakeView();
    const host = await hostStructure();
    const all = Structure.toStructureElementLoci(host.structure);
    const part = (k: number) => StructureElement.Loci(all.structure, [all.elements[k]]);
    const mark = (id: string, k: number, focusSequence?: boolean): Mark => ({
      owner: 'plot',
      id,
      label: id,
      tone: 'accent',
      structure: 's',
      loci: part(k),
      ...(focusSequence ? { focusSequence } : {}),
    });
    // Marks there when the link is made are shown, not focused.
    const list = observableState<readonly Mark[]>([mark('old', 0, true)]);
    const marks = { marks: list.asReadonly, hover: observableState(undefined).asReadonly } as MarksSource;
    const stop = linkSequence(h.viewer, () => fake.view, { marks });
    h.structure.set(host);
    await shown(h);
    const focus = vi.spyOn(h.viewer.seq.focus, 'fromLoci');
    expect(focus).not.toHaveBeenCalled();
    // Two ranges in two chains: one focus over both.
    const box = [mark('aligned', 0, true), mark('scored', 1, true)];
    list.set(box);
    expect(focus).toHaveBeenCalledTimes(1);
    const target = focus.mock.calls[0][0] as { systems: string[]; runs: Int32Array };
    expect(target.systems.length).toBe(2);
    // The same marks set again (new objects, as each highlight.mark makes), or marks that do not ask: no move.
    list.set(box.map((m) => ({ ...m })));
    list.set([...box, mark('other', 0)]);
    expect(focus).toHaveBeenCalledTimes(1);
    // Another range (its label says it): a move.
    list.set([{ ...box[0], label: 'aligned, moved' }, box[1]]);
    expect(focus).toHaveBeenCalledTimes(2);
    stop();
    h.viewers.dispose();
  });

  it("adds the hovered track and Seq*'s lines to Mol*'s label, escaped, and labels again when they change", async () => {
    const h = headlessViewers();
    const fake = fakeView();
    let read: (() => string | undefined) | undefined;
    const refresh = vi.fn();
    const dispose = vi.fn();
    fake.view.annotateHover = (r) => {
      read = r;
      return { refresh, dispose };
    };
    const stop = linkSequence(h.viewer, () => fake.view);
    h.structure.set(await hostStructure());
    await shown(h);
    expect(read?.()).toBeUndefined();
    const labels = h.viewer.seq.labels.current;
    labels.set({
      title: 'A ALA 12',
      lines: [
        { key: 'copies', value: '× 2 (A, B)' },
        { key: 'pLDDT', value: '91.5' },
      ],
      from: ['copies', 'value'],
      track: 'Confidence <b>',
    });
    // The residue's lines are Mol*'s to say: only the track's own.
    expect(read?.()).toBe('<b>Confidence &lt;b&gt;</b><br/>pLDDT: 91.5');
    expect(refresh).toHaveBeenCalledTimes(1);
    // The same lines again: nothing to label again.
    labels.set({ title: 'A ALA 13', lines: [{ key: 'pLDDT', value: '91.5' }], track: 'Confidence <b>' });
    expect(refresh).toHaveBeenCalledTimes(1);
    // A value keyed by the track's name: no track line above it.
    labels.set({
      title: 'A ALA 13',
      lines: [{ key: 'Cα spread', value: '0.037' }],
      from: ['value'],
      track: 'Cα spread',
    });
    expect(read?.()).toBe('Cα spread: 0.037');
    // On one line each (a line break would make the label Markdown); the copies differing here is kept.
    labels.set({
      title: 'A ALA 13',
      lines: [
        { key: 'copies', value: 'varies' },
        { key: 'note', value: 'two\nlines' },
      ],
      from: ['copies', 'item'],
    });
    expect(read?.()).toBe('copies: varies<br/>note: two lines');
    // Only the residue's lines: nothing added, not even the track's name.
    labels.set({ title: 'A ALA 13', lines: [{ key: 'auth', value: '97' }], from: ['position'], track: 'Cα spread' });
    expect(read?.()).toBeUndefined();
    // Only a title (Mol* names the residue): nothing added.
    labels.set({ title: 'A ALA 13', lines: [] });
    expect(read?.()).toBeUndefined();
    labels.set(null);
    expect(refresh).toHaveBeenCalledTimes(4);
    stop();
    expect(dispose).toHaveBeenCalledTimes(1);
    h.viewers.dispose();
  });
});
