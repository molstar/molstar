/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  AddedTrack,
  addTrack,
  cancelDrag,
  composeEditing,
  dragHeight,
  type EditingSession,
  endDrag,
  historyRows,
  proteinSection,
  setDomainsColor,
  toggleTrack,
  trackHeight,
  trackOrder,
} from '../compose.js';

// The Editing story's additions, headless on a fake clock, each composed as the story composes it: the host's stack,
// the viewer and its SVG export (an added track, a hidden one, a theme, Undo and Redo; a drag of a height).

const tracks = P53Views.find((v) => v.id === 'tracks')!;

let session: EditingSession | undefined;

// a frame tick asked for by an earlier file on timers since dropped would never run: this file's own host function
const raf = DomHost.requestAnimationFrame;
beforeEach(() => {
  DomHost.requestAnimationFrame = (fn) => {
    setTimeout(() => fn(performance.now()), 16);
  };
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});

afterEach(() => {
  DomHost.requestAnimationFrame = raf;
  session?.dispose();
  session = undefined;
  vi.useRealTimers();
});

/** The viewer's SVG export, its asynchronous work let go on the fake clock. */
async function svgOf(v: SeqViewer, idPrefix: string): Promise<string> {
  const p = v.renderSvg({ idPrefix });
  await vi.advanceTimersByTimeAsync(300);
  return (await p).svg;
}

describe('Editing: adding a track and a height drag, headless', () => {
  async function make() {
    const s = composeEditing(tracks);
    session = s;
    s.viewer.resize(720);
    await s.viewer.load(tracks.build());
    await vi.advanceTimersByTimeAsync(400);
    return s;
  }

  it('add a track, hide one, set a theme: three steps, undone and redone, the SVG following', async () => {
    const s = await make();
    const svg0 = await svgOf(s.viewer, 'e-');
    expect(trackOrder(s)).toEqual(['letters', 'domains', 'conservation']);

    expect(await addTrack(s)).toBe('Added');
    expect(await addTrack(s)).toBe('Already added');
    expect(await toggleTrack(s, 'conservation')).toBe('Hidden');
    expect(await setDomainsColor(s, '#d62728')).toBe('Colored');
    expect(trackOrder(s)).toEqual(['letters', 'domains', AddedTrack, 'conservation']);
    expect(historyRows(s).done).toEqual(['Color a layer', 'Hide a track', 'Add a track']);
    await vi.advanceTimersByTimeAsync(400);
    const svg1 = await svgOf(s.viewer, 'e-');
    expect(svg1).not.toBe(svg0);
    expect(svg1).toContain('Domains (added)');
    expect(svg1.toLowerCase()).toContain('#d62728');

    for (let i = 0; i < 3; i++) await s.stack.undo();
    await vi.advanceTimersByTimeAsync(400);
    expect(trackOrder(s)).toEqual(['letters', 'domains', 'conservation']);
    expect(await svgOf(s.viewer, 'e-')).toBe(svg0);
    await s.stack.redo();
    expect(trackOrder(s)).toContain(AddedTrack);
  });

  it('a drag previews a height in the spec shown, Esc ends it in nothing, a release in one step', async () => {
    const s = await make();
    const heightOf = (shown: boolean) =>
      proteinSection(s, shown)?.children?.find((c) => c.id === 'domains')?.params?.height;
    const original = heightOf(true);
    expect(original).toBeUndefined();
    // the slider follows the spec: the laid-out height while nothing is set
    expect(trackHeight(s, 'domains')).toBeGreaterThan(0);
    expect(dragHeight(s, 'domains', 40)).toBe(true);
    expect(dragHeight(s, 'domains', 52)).toBe(true);
    expect(heightOf(true)).toBe(52);
    expect(trackHeight(s, 'domains')).toBe(52);
    expect(heightOf(false)).toBe(original);
    expect(cancelDrag(s)).toBe('Cancelled');
    expect(heightOf(true)).toBe(original);
    expect(historyRows(s).done).toEqual([]);

    dragHeight(s, 'domains', 52);
    await endDrag(s);
    expect(heightOf(false)).toBe(52);
    expect(heightOf(true)).toBe(52);
    expect(historyRows(s).done).toEqual(['Resize a track']);
    await s.stack.undo();
    expect(heightOf(true)).toBe(original);
  });
});
