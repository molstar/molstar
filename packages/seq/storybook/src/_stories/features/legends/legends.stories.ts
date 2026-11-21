/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { followSvg } from '../../follow-svg.js';
import { composeEditing, type EditingSession } from '../editing/compose.js';
import { mountPlacements } from './legends-compose.js';

// The legend part placed four ways by its host (the legends design): a track header's line (compact, the one track's), a
// side panel (the rest, stacked), a band under the canvas (compact, in a row) and the export's band below its figure.
// The same entries everywhere: hover or click a swatch to highlight or select its items.

const view = FeatureViews.find((v) => v.id === 'themes-editing')!;

const Page =
  'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

function box(label: string, css = ''): HTMLElement {
  const el = document.createElement('section');
  el.setAttribute('aria-label', label);
  el.style.cssText = `border: 1px solid var(--msui-border, currentColor); border-radius: 6px; padding: var(--msui-space-3); min-width: 0; ${css}`;
  const h = document.createElement('h3');
  h.textContent = label;
  h.style.cssText = 'margin: 0 0 var(--msui-space-2); font-size: var(--msui-text-xs); color: var(--msui-text-muted);';
  el.append(h);
  return el;
}

export class LegendPlacementsRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-legend-placements-root';

  session: EditingSession | undefined;
  loaded: Promise<void> | undefined;
  hosts: { header: HTMLElement; panel: HTMLElement; band: HTMLElement; svg: HTMLElement } | undefined;

  protected render(model: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText = Page;
    const note = document.createElement('p');
    note.style.cssText = Muted;
    note.textContent =
      'One legend part, placed by the host: the first track’s legend compact beside its header line, the others in a side panel, all of them compact in a band under the canvas, and the export’s band below its figure. A layer’s `legend: none` would take its key out of every one.';
    const header = box('Track header line', 'display: flex; gap: var(--msui-space-3); align-items: center;');
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const panel = box('Side panel');
    const main = document.createElement('div');
    main.style.cssText = 'display: grid; grid-template-columns: minmax(0, 1fr) 260px; gap: var(--msui-space-4);';
    const left = document.createElement('div');
    const band = box('Band under the canvas');
    left.append(canvasHost, band);
    main.append(left, panel);
    const svg = box('Export');
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'overflow-x: auto;';
    svg.append(svgHost);
    root.append(note, header, main, svg);

    const s = composeEditing(model, { target: canvasHost, legends: false });
    this.session = s;
    this.hosts = { header, panel, band, svg: svgHost };
    const follow = followSvg(s.viewer, svgHost, { prefix: 'lg', svg: { legend: 'below' } });
    let unmount: (() => void) | undefined;
    onCleanup(() => {
      unmount?.();
      follow.dispose();
      s.dispose();
      this.session = undefined;
    });
    s.viewer.resize(960);
    this.loaded = s.viewer.load(model.build()).then(async () => {
      await s.viewer.settle();
      unmount = mountPlacements(s.viewer, { header, panel, band });
    });
    return root;
  }
}

defineElements([LegendPlacementsRoot]);

const meta: Meta = { title: 'Features/Legends' };
export default meta;

/** The legend part placed four ways: a track header's line (compact), a side panel, a band under the canvas and the export's. */
export const LegendPlacements: StoryObj = {
  render: () => {
    const el = document.createElement(LegendPlacementsRoot.tag) as LegendPlacementsRoot;
    el.model = view;
    return el;
  },
};
