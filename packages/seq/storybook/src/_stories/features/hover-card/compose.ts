/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HoverCardProvider } from '@molstar/seqx-viewer/behaviors';
import type { HoverCard, HoverCardPlace, HoverCardPortFactory } from '@molstar/seqx-viewer/hover-card';
import type { SeqViewerSpec } from '@molstar/seqx-viewer';

// The hover card's stories, DOM-free (the stories and their headless twin call these): a viewer's spec with a host's
// own card adapter in place of `@molstar/ui`'s (the hover card port: the same data, the host's components and styles),
// and the lines a host's card writes for a card.

/** The default spec with the hover card drawn by `port` (a host's adapter). */
export function withCardPort(spec: SeqViewerSpec, port: HoverCardPortFactory): SeqViewerSpec {
  return {
    ...spec,
    behaviors: spec.behaviors.map((b) =>
      b.provider === HoverCardProvider ? { provider: HoverCardProvider, params: { port } } : b,
    ),
  };
}

/** A host's card as lines: its title, then `label = value` a row (what the story's terminal-like card prints). */
export function hostCardLines(card: HoverCard): string[] {
  return [card.title, ...card.rows.map((r) => `${r.label} = ${r.value}`)];
}

/** A host adapter that records what it is told (the story's twin reads it). */
export function recordingCard(): { calls: { card: HoverCard; place: HoverCardPlace }[]; port: HoverCardPortFactory } {
  const calls: { card: HoverCard; place: HoverCardPlace }[] = [];
  return { calls, port: () => ({ show: (card, place) => void calls.push({ card, place }), hide() {} }) };
}
