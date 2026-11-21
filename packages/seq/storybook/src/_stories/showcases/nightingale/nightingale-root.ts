/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import '@nightingale-elements/nightingale-linegraph-track';
import '@nightingale-elements/nightingale-sequence';
import '@nightingale-elements/nightingale-track';
import type { ExampleView } from '@molstar/seqx-examples';
import { mountLegend } from '@molstar/seqx-viewer/mount-legend';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import { followSvg } from '../../follow-svg.js';
import { composeViewer } from '../../features/selection/compose.js';
import { highlightOf, nightingaleData, type NightingaleData, Window } from './compose.js';
import type { ThemeName } from './nightingale-adapter.js';

// The Nightingale showcase's root (composition.md: other viewers read the resolved view). Seq* opens p53's tracks view;
// the adapter reads its letters, its domains and its conservation from the resolved view and the store, and what it
// reads is set on Nightingale's own web components (`nightingale-sequence`, `nightingale-track`,
// `nightingale-linegraph-track`), in the colors of Seq*'s theme tokens. The `link` story adds Seq*'s canvas, whose
// hover and selection are the elements' `highlight`.

export type NightingaleStoryKind = 'features' | 'values' | 'link';

export interface NightingaleStory {
  readonly view: ExampleView;
  readonly kind: NightingaleStoryKind;
  /** What is read from where, a sentence or two. */
  readonly reads: string;
}

const NgWidth = 640;

const Box = 'display: grid; gap: var(--msui-space-2); min-width: 0; align-content: start;';
const Title =
  'margin: 0; font-size: var(--msui-text-sm); font-weight: var(--msui-weight-strong); color: var(--msui-text-muted);';

// Nightingale draws its text with SVG attributes: the theme's text color is set by the stylesheet, which wins over them.
const NgStyle = `
  nightingale-sequence text, nightingale-track text, nightingale-linegraph-track text { fill: var(--msui-text); }
  nightingale-sequence .highlighted rect, nightingale-track .highlighted rect { opacity: 0.35; }
  nightingale-sequence rect.base_bg { fill: var(--ng-track, #eee); }
  nightingale-sequence, nightingale-track, nightingale-linegraph-track { display: block; }
`;

function section(title: string, body: HTMLElement): HTMLElement {
  const el = document.createElement('div');
  el.style.cssText = Box;
  const h = document.createElement('h4');
  h.style.cssText = Title;
  h.textContent = title;
  el.append(h, body);
  return el;
}

type NgElement = HTMLElement & { data?: unknown };

/** Nightingale's element for a tag, sized to the shared width and the window it shows (1-based, inclusive). */
function ngElement(tag: string, length: number, height: number, window?: { start: number; end: number }): NgElement {
  const el = document.createElement(tag) as NgElement;
  const attrs: Record<string, string | number> = {
    length,
    width: NgWidth,
    height,
    'margin-left': 0,
    'margin-right': 0,
    'margin-top': 0,
    'margin-bottom': 0,
    'display-start': window?.start ?? 1,
    'display-end': window?.end ?? length,
  };
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
  return el;
}

const themeOf = (el: HTMLElement): ThemeName =>
  (el.closest('[data-msui-theme]')?.getAttribute('data-msui-theme') ||
    el.ownerDocument.documentElement.getAttribute('data-msui-theme')) === 'dark'
    ? 'dark'
    : 'light';

class NightingaleRoot extends UiElement<NightingaleStory> {
  static readonly tag = 'seqx-nightingale-root';

  protected render({ view, kind, reads }: NightingaleStory): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const style = document.createElement('style');
    style.textContent = NgStyle;
    const hint = document.createElement('p');
    hint.style.cssText = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';
    hint.textContent = reads;
    const status = document.createElement('p');
    status.style.cssText =
      'margin: 0; font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono); color: var(--msui-text-muted);';

    const ngHost = document.createElement('div');
    ngHost.style.cssText = `width: ${NgWidth}px; max-width: 100%; display: grid; gap: var(--msui-space-1); overflow-x: auto;`;
    const seqHost = document.createElement('div');
    seqHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    const legendHost = document.createElement('div');
    const pair = document.createElement('div');
    pair.style.cssText =
      'display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5);';

    const viewer = composeViewer(view, { target: kind === 'link' ? seqHost : undefined, noRails: kind !== 'link' });
    onCleanup(() => viewer.retire());
    // what is made after the load: `onCleanup` is for the render scope only
    const late: (() => void)[] = [];
    onCleanup(() => {
      for (const d of late) d();
    });

    // the three elements, made once the data is known
    let data: NightingaleData | undefined;
    const els: NgElement[] = [];
    const paint = async () => {
      const theme = themeOf(root);
      data = await nightingaleData(viewer, theme);
      const { read, features, lines, colors } = data;
      root.style.setProperty('--ng-track', colors.track);
      for (const e of els) e.remove();
      els.length = 0;
      const window = kind === 'link' ? undefined : kind === 'features' ? Window : undefined;
      const seq = ngElement('nightingale-sequence', read.length, 36, window);
      seq.setAttribute('sequence', read.letters);
      const track = ngElement('nightingale-track', read.length, 32, window);
      track.setAttribute('layout', 'non-overlapping');
      track.data = features;
      const line = ngElement('nightingale-linegraph-track', read.length, 80, window);
      line.data = lines.map((l) => ({ ...l, range: [...l.range] }));
      for (const e of [seq, track, line]) {
        e.setAttribute('highlight-color', colors.highlight);
        e.setAttribute('margin-color', 'transparent');
      }
      const shown = kind === 'features' || kind === 'link' ? [seq, track] : [seq, line];
      els.push(...shown);
      ngHost.append(...shown);
      status.textContent = `read from the resolved view and the store: ${read.length} letters, ${read.ranges.length} ranges, ${read.values.length} values (${theme})`;
      showLink();
    };
    const showLink = () => {
      if (kind !== 'link' || !data) return;
      const highlight = highlightOf(viewer, data.read.system);
      for (const e of els) e.setAttribute('highlight', highlight);
    };

    if (kind === 'link') {
      onCleanup(viewer.hover.version.listen(showLink));
      onCleanup(viewer.selection.value.listen(showLink));
      onCleanup(mountLegend(legendHost, viewer));
    }
    let redrawSvg: (() => Promise<void>) | undefined;
    const mutations = new MutationObserver(() => {
      void paint();
      void redrawSvg?.();
    });
    mutations.observe(document.documentElement, { attributes: true, attributeFilter: ['data-msui-theme'] });
    onCleanup(() => mutations.disconnect());

    pair.append(section('Nightingale', ngHost));
    if (kind === 'link') pair.append(section('Seq*', seqHost));
    else pair.append(section('Seq*, the same positions', seqHost));
    root.append(style, hint, pair, status);

    void viewer
      .load(view.build())
      .then(() => viewer.settle())
      .then(async () => {
        await paint();
        if (kind === 'link') return;
        const section0 = viewer.view.resolved.value?.sections[0];
        if (!section0) return;
        // the same positions on both sides: Seq*'s export of the window Nightingale shows (0-based, half-open)
        const window = kind === 'features' ? ([Window.start - 1, Window.end] as const) : undefined;
        const svg = followSvg(viewer, seqHost, {
          prefix: 'ng',
          svg: {
            size: { width: NgWidth + 96, height: 200 },
            // Seq*'s export in the theme shown, read at each draw
            get theme() {
              return themeOf(root);
            },
            ...(window ? { viewport: { positions: { [section0.address]: window } } } : {}),
          },
        });
        late.push(() => svg.dispose());
        redrawSvg = () => svg.draw();
        await svg.draw();
      })
      .catch((e: unknown) => {
        status.textContent = `The view did not load: ${e instanceof Error ? e.message : String(e)}`;
      });
    return root;
  }
}

defineElements([NightingaleRoot]);

/** A Storybook `render` result. */
export function nightingaleStory(story: NightingaleStory): HTMLElement {
  const el = document.createElement(NightingaleRoot.tag) as NightingaleRoot;
  el.model = story;
  return el;
}
