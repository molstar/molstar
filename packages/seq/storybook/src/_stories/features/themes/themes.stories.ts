/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { Color } from '@molstar/core/util/color/color';
import { type ColorVision, ColorVisions, simulateColorVision } from '@molstar/core/util/color/cvd';
import { DarkTokens, LightTokens, type SeqTokens } from '@molstar/seqx-core/render';
import {
  BuiltInSchemes,
  darkVariant,
  resolveColor,
  type SchemeEntry,
  schemeColors,
  toRgb,
} from '@molstar/seqx-core/theme/index';
import { type ExampleView, FeatureViews, P53Views } from '@molstar/seqx-examples';
import { mountLegend } from '@molstar/seqx-viewer/mount-legend';
import { Button } from '@molstar/ui/components/button/button';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { type EditingSession, composeEditing } from '../editing/compose.js';
import { followSvg } from '../../follow-svg.js';
import { exampleStory } from '../../examples/example-root.js';
import { cancelScheme, pickerSchemes, pickScheme, previewScheme, setSignalDomain } from './themes-compose.js';

// p53's themes: the view with every scheme that fits its letters, its legends and their export; the view edited (a scheme
// previewed while hovered, picked in one step; the shared scale's domain set for both signals in one step; the reader's
// color-blind substitution); and the gallery of the built-in schemes in the light and dark tables, each seen through the
// three simulated color-vision deficiencies. Strings from the spec go in as text, never markup.

const view = P53Views.find((v) => v.id === 'themes')!;
const editingView = FeatureViews.find((v) => v.id === 'themes-editing')!;

const Page =
  'display: grid; align-content: start; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
const Row = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
const Muted = 'margin: 0; color: var(--msui-text-muted); max-width: 72ch;';

function text(tag: string, content: string, css = ''): HTMLElement {
  const el = document.createElement(tag);
  el.style.cssText = css;
  el.textContent = content;
  return el;
}

// Each root renders once it has a model (`UiElement` renders nothing without one): the view, the schemes.
export class ThemesEditingRoot extends UiElement<ExampleView> {
  static readonly tag = 'seqx-p53-themes-editing-root';

  /** The composition the story drives (its headless twin drives the same one). */
  session: EditingSession | undefined;
  /** The view's load, then the picker filled. */
  loaded: Promise<void> | undefined;

  protected render(view: ExampleView): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText = Page;
    root.append(
      text(
        'p',
        'Hover a scheme to preview it on the classes track (no step); click to pick it (one step). Set the shared scale’s top: both signals follow in one step. Undo takes either back whole. The legends below are the viewer’s Legend part: click a value to select its items.',
        Muted,
      ),
    );
    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    const s: EditingSession = composeEditing(view, { target: canvasHost });
    this.session = s;
    onCleanup(() => {
      s.dispose();
      this.session = undefined;
    });
    const status = text('p', '', Muted);
    const say = (p: Promise<string> | string) =>
      void Promise.resolve(p).then((m) => {
        status.textContent = m;
      });

    // the picker: the schemes that fit the letters, hovered to preview
    const picker = document.createElement('div');
    picker.style.cssText = Row;
    picker.setAttribute('role', 'listbox');
    picker.setAttribute('aria-label', 'Letter schemes');
    picker.addEventListener('mouseleave', () => cancelScheme(s));
    picker.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') cancelScheme(s);
    });

    const scale = document.createElement('label');
    scale.style.cssText = Row;
    const top = document.createElement('input');
    top.type = 'number';
    top.step = '0.1';
    top.placeholder = 'in view';
    top.style.cssText = 'width: 8ch;';
    top.addEventListener('change', () => say(setSignalDomain(s, [0, top.value === '' ? null : Number(top.value)])));
    scale.append('Signal scale top (empty: in view)', top);

    const tools = document.createElement('div');
    tools.style.cssText = Row;
    tools.append(
      Button({
        size: 'sm',
        children: 'Undo',
        disabled: derive(s.stack.steps, (x) => x.length === 0),
        onClick: () => void s.stack.undo(),
      }) as HTMLElement,
      Button({
        size: 'sm',
        children: 'Redo',
        disabled: derive(s.stack.redoable, (x) => x.length === 0),
        onClick: () => void s.stack.redo(),
      }) as HTMLElement,
      Button({
        size: 'sm',
        children: 'Color-blind safe',
        onClick: () => say(s.viewer.execute('theme.cvd-safe', {})),
      }) as HTMLElement,
    );

    const legend = document.createElement('div');
    onCleanup(mountLegend(legend, s.viewer));
    const svgHost = document.createElement('div');
    svgHost.style.cssText = 'min-width: 0; overflow-x: auto;';
    root.append(canvasHost, picker, scale, tools, status, legend, svgHost);
    const svg = followSvg(s.viewer, svgHost, { prefix: 'p53t' });
    onCleanup(() => svg.dispose());

    s.viewer.resize(960);
    this.loaded = s.viewer.load(view.build()).then(() => {
      for (const name of pickerSchemes(s)) {
        const b = Button({ size: 'sm', children: name, onClick: () => say(pickScheme(s, name)) }) as HTMLElement;
        b.setAttribute('role', 'option');
        b.addEventListener('mouseenter', () => previewScheme(s, name));
        b.addEventListener('focus', () => previewScheme(s, name));
        picker.append(b);
      }
    });
    return root;
  }
}

// ----- the schemes gallery -----

const Visions: readonly ColorVision[] = ColorVisions;

/** A scheme's literal and token colors in one table, as RGB. */
function colorsIn(entry: SchemeEntry, tokens: SeqTokens, dark: boolean): Color[] {
  const out: Color[] = [];
  for (const c of schemeColors(entry)) {
    const spelled = dark && entry.dark === 'derive' && c.startsWith('#') ? darkVariant(c) : c;
    const v = resolveColor(spelled, tokens);
    if (v !== undefined) out.push(toRgb(v));
  }
  return out;
}

function swatches(colors: readonly Color[], vision: ColorVision, label: string): HTMLElement {
  const row = document.createElement('div');
  row.style.cssText = 'display: flex; flex-wrap: wrap; gap: 1px;';
  row.setAttribute('role', 'img');
  row.setAttribute('aria-label', label);
  for (const c of colors) {
    const sw = document.createElement('span');
    sw.style.cssText = `width: 14px; height: 14px; background: ${Color.toHexStyle(simulateColorVision(c, vision))};`;
    row.append(sw);
  }
  return row;
}

export class SchemesGalleryRoot extends UiElement<readonly SchemeEntry[]> {
  static readonly tag = 'seqx-schemes-gallery-root';

  protected render(schemes: readonly SchemeEntry[]): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText = Page;
    root.append(
      text(
        'p',
        'Every built-in scheme in the light and the dark table, each as normal vision and simulated protanopia, deuteranopia and tritanopia (Machado 2009). A scheme marked safe keeps its colors apart under all three; one that is not names the scheme the reader’s substitution draws instead.',
        Muted,
      ),
    );
    const table = document.createElement('table');
    table.style.cssText = 'border-collapse: collapse; font-size: var(--msui-text-xs);';
    const head = document.createElement('tr');
    for (const h of ['Scheme', 'Table', ...Visions]) head.append(text('th', h, 'text-align: left; padding: 4px 8px;'));
    table.append(head);
    for (const entry of schemes) {
      for (const [name, tokens, dark] of [
        ['light', LightTokens, false],
        ['dark', DarkTokens, true],
      ] as const) {
        const colors = colorsIn(entry, tokens, dark);
        if (!colors.length) continue;
        const tr = document.createElement('tr');
        const bg = Color.toHexStyle(toRgb(tokens.colors['msui-bg']));
        tr.style.cssText = `background: ${bg};`;
        const safe = entry.cvdSafe ? 'safe' : `not safe${entry.substitute ? ` (substitute: ${entry.substitute})` : ''}`;
        tr.append(
          text('td', `${entry.label} · ${safe}`, `padding: 4px 8px; color: ${dark ? '#eee' : '#222'};`),
          text('td', name, `padding: 4px 8px; color: ${dark ? '#eee' : '#222'};`),
        );
        for (const v of Visions) {
          const td = document.createElement('td');
          td.style.cssText = 'padding: 4px 8px;';
          td.append(swatches(colors, v, `${entry.label}, ${name} table, ${v}`));
          tr.append(td);
        }
        table.append(tr);
      }
    }
    root.append(table);
    return root;
  }
}

defineElements([ThemesEditingRoot, SchemesGalleryRoot]);

const meta: Meta = { title: 'Features/Themes' };
export default meta;

/** p53 under every theme form: schemes on the letters, categories with their legend, pLDDT bands and a ramp, a shared scale, a log-sized variant track; the export draws the legends below. */
export const Themes: StoryObj = { render: () => exampleStory(view) };

/** The themes edited: a scheme previewed on hover and picked in one step, the shared scale set for both signals, the reader's substitution. */
export const EditingThemes: StoryObj = {
  render: () => {
    const el = document.createElement(ThemesEditingRoot.tag) as ThemesEditingRoot;
    el.model = editingView;
    return el;
  },
};

/** The built-in schemes in both tables, each under the three simulated color-vision deficiencies. */
export const SchemesGallery: StoryObj = {
  render: () => {
    const el = document.createElement(SchemesGalleryRoot.tag) as SchemesGalleryRoot;
    el.model = BuiltInSchemes;
    return el;
  },
};
