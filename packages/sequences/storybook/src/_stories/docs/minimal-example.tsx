/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Introduction's example, shown on the page as written (imported raw) and run beside it: one document, mounted by
// the React UI, by the DOM UI, rendered by Nightingale, and drawn headless as SVG.
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { compile, createRegistries, KernelRegistry, register, resolveView } from '@molstar/seq-core/format';
import { materialize, NightingaleSupported } from '@molstar/seq-nightingale';
import { mount } from '@molstar/seq-nightingale/mount';
import { defineSeqViewer } from '@molstar/seq-viewer/dom';
import { SeqViewer } from '@molstar/seq-viewer/react';
import { renderSvg } from '@molstar/seq-viewer/svg';
import '@molstar/seq-viewer/skin.css';

// A document: data (a sequence and an annotation over it) and a view of tracks. Synthetic data.
const b = createSeqViewSpec({ title: 'A synthetic protein' });
const protein = b.data.sequence({
  ref: 'protein',
  label: 'Protein',
  value: 'MSDNELKQLAEEGRKLVTALSPEDAIRWLHQNGFDVSKAIELLKGGDPSTVYYAVKIGNTEEAKQFLRDHGF',
});
const domains = b.data.annotation({
  ref: 'domains',
  geometry: 'ranges',
  target: protein,
  value: { start: [4, 40], stop: [30, 68], label: ['Helical bundle', 'Binding loop'] },
});
const view = b.view({ ref: 'tracks', label: 'Tracks' }).block({ sequence: protein });
view.track({ label: '' }).ruler();
view.track({ label: 'Residues' }).letters();
view.track({ label: 'Domains' }).ranges({ data: domains, labels: true });
export const doc = b.getDocument();

// React: the default UI as a component, `clean` by default (the canvas, hover and the tooltip); `ui={{ chrome: "full" }}`
// adds the tabs, section headers, track menus and the legend.
export const ReactExample = () => <SeqViewer document={doc} />;

// DOM: the <seq-viewer> element, no framework.
export function mountDomExample(host: HTMLElement) {
  defineSeqViewer();
  const element = document.createElement('seq-viewer');
  element.document = doc;
  host.append(element);
  return () => element.remove();
}

// Nightingale, the format's second renderer: the same document compiled (all by value: nothing to load), its view
// resolved for what Nightingale draws, materialized as Nightingale's components and mounted. `issues`: what it cannot draw.
const registries = createRegistries();
register(registries, [KernelRegistry]);
export async function mountNightingaleExample(host: HTMLElement) {
  const { set } = await compile(doc, registries, { load: () => Promise.reject(new Error('all by value')) });
  const resolved = resolveView(doc, set, 'tracks', registries, NightingaleSupported);
  const { root, issues } = materialize(resolved.view, set);
  const mounted = mount(root, host, { view: resolved.view, set });
  return { mounted, issues: [...resolved.issues, ...issues] };
}

// Headless: the same document as an SVG string, on its skin's surface (Node or the browser; no DOM).
export const svgExample = () => renderSvg(doc, { width: 640, background: true });
