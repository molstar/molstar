/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { legendRows } from '@molstar/seqx-viewer/legends';
import { type LegendPartOptions, mountLegend } from '@molstar/seqx-viewer/mount-legend';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The legends story's composition, DOM-free but for the elements it is given: the story and its headless twin call it.
// One legend part placed three ways by the host (a track header's line, a side panel, a band under the canvas) and the
// export's; the host decides where each legend sits, the layer's `legend: none` what a figure keys at all.

/** A track's address from a layer's. */
export const trackOf = (layer: string) => layer.slice(0, layer.lastIndexOf('/'));

export interface Placements {
  /** The legend beside a track's header: compact, the one track's. */
  readonly header: LegendPartOptions;
  /** The side panel: every legend but the header's, stacked, full. */
  readonly panel: LegendPartOptions;
  /** The band under the canvas: every legend, in a row, compact. */
  readonly band: LegendPartOptions;
}

/** The options of each placement for a viewer whose legends are derived: the header takes the first legend's track. */
export function placementsOf(viewer: SeqViewer): Placements | undefined {
  const first = legendRows(viewer)[0];
  if (!first) return undefined;
  const track = trackOf(first.layer);
  return {
    header: { of: track, form: 'compact', channels: 'collapse' },
    panel: { omit: [track], form: 'full', orient: 'vertical', channels: 'collapse' },
    band: { form: 'compact', orient: 'horizontal', channels: 'collapse' },
  };
}

export interface Hosts {
  readonly header: HTMLElement;
  readonly panel: HTMLElement;
  readonly band: HTMLElement;
}

/** Mounts the three placements into their hosts; returns what takes them away. */
export function mountPlacements(viewer: SeqViewer, hosts: Hosts): (() => void) | undefined {
  const p = placementsOf(viewer);
  if (!p) return undefined;
  const offs = [
    mountLegend(hosts.header, viewer, p.header),
    mountLegend(hosts.panel, viewer, p.panel),
    mountLegend(hosts.band, viewer, p.band),
  ];
  return () => offs.forEach((off) => off());
}
