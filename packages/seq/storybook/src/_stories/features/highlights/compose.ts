/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { HighlightStyle } from '@molstar/seqx-core/render';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The highlights story's composition, DOM-free: the story's buttons and its headless twin call these, so the twin acts
// on the viewer as the buttons do (through the viewer's commands, never around them).

/** A channel the story sets: its name, what it lights (an address on p53) and the style it starts with. */
export interface StoryChannel {
  readonly name: string;
  readonly title: string;
  readonly address: string;
  readonly style: HighlightStyle;
}

/** Three channels that overlap on p53: a domain, a motif search result, the host's own. */
export const Channels: readonly StoryChannel[] = [
  { name: 'domain', title: 'DNA-binding domain', address: '/data/p53/(range :from 95 :to 292)', style: 'tint' },
  { name: 'motif', title: 'Motif search result', address: '/data/p53/(range :from 150 :to 240)', style: 'band' },
  { name: 'host', title: 'Host’s own', address: '/data/p53/(range :from 200 :to 330)', style: 'outline' },
];

/** The channels' colors (a token each): without them every named channel takes the one default color. */
export const ChannelColors = { domain: 'token:cat-1', motif: 'token:cat-3', host: 'token:cat-6' };

/** The viewer config the story composes with: the channels' colors. */
export const HighlightConfig = [[viewerItem('highlight.colors'), ChannelColors]] as const;

export const Styles: readonly HighlightStyle[] = ['tint', 'band', 'outline'];

/** The next style in the cycle. */
export const nextStyle = (s: HighlightStyle): HighlightStyle => Styles[(Styles.indexOf(s) + 1) % Styles.length];

/** Sets a channel (the story's button), in the style it has now, else its own. */
export async function setChannel(viewer: SeqViewer, channel: StoryChannel, style?: HighlightStyle): Promise<void> {
  await viewer.execute('address.highlight', {
    address: channel.address,
    channel: channel.name,
    style: style ?? viewer.highlight.get(channel.name)?.style ?? channel.style,
  });
}

/** Clears a channel (the button, and the Legend part's clear). */
export async function clearChannel(viewer: SeqViewer, channel: StoryChannel): Promise<void> {
  await viewer.execute('highlight.clear', { channel: channel.name });
}

/** Sets or clears a channel. */
export async function toggleChannel(viewer: SeqViewer, channel: StoryChannel): Promise<void> {
  if (viewer.highlight.get(channel.name)) await clearChannel(viewer, channel);
  else await setChannel(viewer, channel);
}

/** Gives a set channel its next style. */
export async function cycleStyle(viewer: SeqViewer, channel: StoryChannel): Promise<void> {
  const now = viewer.highlight.get(channel.name)?.style ?? channel.style;
  await setChannel(viewer, channel, nextStyle(now));
}

/** What the status line says: the channels set, in order, with their styles. */
export function describeChannels(viewer: SeqViewer): string {
  const names = viewer.highlight.names.value;
  return names.length
    ? `channels: ${names.map((n) => `${n} (${viewer.highlight.get(n)?.style ?? 'tint'})`).join(', ')}`
    : 'channels: none';
}
