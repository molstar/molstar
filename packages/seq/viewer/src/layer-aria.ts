/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewLayout } from '@molstar/seqx-core/render';

// A layer's ARIA words (_spec/README.md "Input", Announcements, A layer's summary): each layer whose type says a
// `summary` (core's `LayerType`: an extension's ribbons, "612 ribbons, 23 chromosomes") a group beside the canvas,
// labelled by its track's label (else its key), set when a frame settles. DOM-free words; the mount owns the nodes.

/** Each summarized layer of a layout: its label and its type's words. */
export function layerAria(layout: ViewLayout | undefined): { label: string; description: string }[] {
  const out: { label: string; description: string }[] = [];
  for (const s of layout?.sections ?? []) {
    for (const t of s.tracks)
      for (const l of t.layers) {
        const d = l.type.summary?.(l.lay);
        if (d) out.push({ label: t.track.label || l.key, description: d });
      }
    for (const l of s.sectionLays) {
      const d = l.type.summary?.(l.lay);
      if (d) out.push({ label: l.key, description: d });
    }
  }
  return out;
}

const Hosts = new WeakMap<HTMLElement, HTMLElement>();

/** The groups in a mounted viewer's container: visually hidden, a `role="group"` a layer, replaced when a frame settles. */
export function setLayerAria(container: HTMLElement, layout: ViewLayout | undefined): void {
  const words = layerAria(layout);
  let host = Hosts.get(container);
  if (!words.length) {
    host?.replaceChildren();
    return;
  }
  const doc = container.ownerDocument;
  if (!host) {
    host = doc.createElement('div');
    host.style.cssText = 'position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap';
    container.append(host);
    Hosts.set(container, host);
  }
  host.replaceChildren(
    ...words.map((w) => {
      const g = doc.createElement('div');
      g.setAttribute('role', 'group');
      g.setAttribute('aria-label', w.label);
      g.textContent = w.description;
      return g;
    }),
  );
}
