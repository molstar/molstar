/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's stories, headless on a fake clock, composed as the stories compose them: p53's tracks with a host's
// own card adapter in place of `@molstar/ui`'s get the same data (the viewer's description, the announcer's text) and
// the same place; the keys' data the panel and the host's list draw lists the toolbar's command with its key.

import { P53Views } from '@molstar/seqx-examples';
import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import type { HeadersBehavior } from '@molstar/seqx-viewer/behaviors/headers';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { describeKeys } from '@molstar/seqx-viewer/keys';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { hostCardLines, recordingCard, withCardPort } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

async function open(port: ReturnType<typeof recordingCard>['port']) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  viewer = createSeqViewer(withCardPort({ ...DefaultSeqViewerSpec(), ui: undefined }, port));
  viewer.resize(720);
  void viewer.load(P53Views.find((v) => v.id === 'tracks')!.build());
  await vi.advanceTimersByTimeAsync(2000);
  // the headers' lazy chunks (a dynamic import: real I/O, never on the fake clock, slow in a cold worker under load):
  // one landing after the card shows relays the headers out, and a frame that re-laid out takes a hover's card away
  // until it rests again (`framed(relaid)`), so the card is measured once they have landed
  await viewer.behavior<HeadersBehavior>('headers')!.loaded();
  await vi.advanceTimersByTimeAsync(200);
  viewer.draw();
  return viewer;
}

describe('the hover card’s stories, headless', () => {
  it('a host’s adapter is told the viewer’s card and place on rest', async () => {
    const rec = recordingCard();
    const v = await open(rec.port);
    const s = v.frame!.layout.sections[0];
    const t = s.tracks[0];
    v.input.move((s.axis!.x(9) + s.axis!.x(10)) / 2, t.box.y + t.box.height / 2);
    expect(rec.calls).toHaveLength(0);
    // the move itself asks for a chunk (the headers' icons, on the track's first hover): landed before the clock runs
    await v.behavior<HeadersBehavior>('headers')!.loaded();
    await vi.advanceTimersByTimeAsync(600);
    expect(rec.calls).toHaveLength(1);
    const { card, place } = rec.calls[0];
    expect(card).toEqual(v.behavior(HoverCardProvider)!.shown.value!.card);
    expect(hostCardLines(card)[0]).toBe(card.title);
    expect(card.text.length).toBeGreaterThan(0);
    expect(place.view.width).toBe(v.frame!.layout.width);
  });

  it('the keys’ data lists the keys panel’s command with its key', async () => {
    const v = await open(recordingCard().port);
    const view = describeKeys(v, false).find((g) => g.id === 'view')!;
    expect(view.entries.find((e) => e.action === 'help.keys')).toMatchObject({
      name: 'Keys and gestures',
      keys: ['?'],
    });
  });
});
