/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ThemingBehavior } from '@molstar/seqx-viewer/behaviors';
import type { EditingSession } from '../../features/editing/compose.js';

// The themes editing story's composition, DOM-free: the story and its headless twin call these. A scheme is previewed
// while hovered in the picker and picked in one patch; the shared scale's domain is set for both signals in one step.
// Every edit goes through the viewer's commands or its legend behavior's preview.

export const Letters = '/views/main/protein/letters-classes/letters';
export const Signal = '/views/signal';
export const SignalLayers = ['/views/main/protein/conservation/values', '/views/main/protein/accessibility/values'];

const theming = (s: EditingSession) => s.viewer.behaviors.get('legend') as ThemingBehavior;

/** The schemes the picker offers for the classes track's letters. */
export function pickerSchemes(s: EditingSession): string[] {
  return theming(s)
    .fittingSchemes(Letters)
    .map((e) => e.name);
}

/** A scheme hovered in the picker: previewed, not a step. */
export function previewScheme(s: EditingSession, name: string): boolean {
  return theming(s).preview(Letters, 'color', { kind: 'scheme', params: { name } });
}

/** The hover left the picker (or Esc): the preview ends in nothing. */
export function cancelScheme(s: EditingSession): void {
  theming(s).cancelPreview();
}

/** A scheme clicked: one patch, one step of the host's stack. */
export function pickScheme(s: EditingSession, name: string): Promise<string> {
  return Promise.resolve(
    s.viewer.execute('layer.channel', {
      address: Letters,
      channel: 'color',
      value: { kind: 'scheme', params: { name } },
    }),
  );
}

/** The shared scale's domain set: both signals follow, one step. */
export function setSignalDomain(s: EditingSession, domain: readonly [number | null, number | null]): Promise<string> {
  return Promise.resolve(s.viewer.execute('scale.set', { address: Signal, scale: { domain: [...domain] } }));
}

/** Each signal's y domain as laid out. */
export function signalDomains(s: EditingSession): (readonly [number, number] | undefined)[] {
  const layers = s.viewer.frame?.layout.sections.flatMap((x) => x.tracks.flatMap((t) => t.layers)) ?? [];
  return SignalLayers.map((a) => {
    const l = layers.find((x) => x.key === a);
    return l ? l.type.domain?.(l.lay) : undefined;
  });
}
