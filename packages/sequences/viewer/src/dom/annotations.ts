/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Annotations` part (seqspec:interaction#the-default-ui), DOM: a section's bookmarks set as annotations the person
// adds from the viewer's selection, with a label and a color, lists, reveals and deletes. A thin client of the patch
// API (src/annotations.ts: `bookmark.add`, `bookmark.remove`, `ctx.reveal`); the set is drawn by the document's own
// ranges track over it. The label and the color are leaf state, in the part's two fields; so is the place chosen when
// the panel has several sets.

import type { LineStyle } from '@molstar/seq-core/format';
import {
  ANNOTATION_COLOR,
  ANNOTATION_GLYPH,
  ANNOTATION_GLYPHS,
  ANNOTATION_STYLE,
  ANNOTATION_STYLES,
  type AnnotationEntry,
  type AnnotationPlace,
  addAnnotation,
  annotationEntries,
  annotationPlaces,
  annotationStates,
  annotationTarget,
  deleteName,
  placeGlyphs,
  removeAnnotation,
  revealAnnotation,
  showName,
} from '../annotations.js';
import { sameData } from '../same.js';
import { button, type DomPart, derived, h, off, type SectionProps, setAttr, setText, spot, watch } from './core.js';
import { CROSS } from './shared.js';

export interface AnnotationsPanelProps extends SectionProps {
  /** Another set than `<section>.bookmarks` (`bookmark.add`'s `data`). */
  readonly data?: string;
  /**
   * Several sets to choose between (each `bookmark.add`'s `data`, with the name a person sees): the add row offers the
   * choice, the list shows the sets one after another with their names, a delete acts on its row's set. Overrides `data`.
   */
  readonly places?: readonly AnnotationPlace[];
}

const SVG = 'http://www.w3.org/2000/svg';
/** The delete button's ✕ (`CROSS`), decorative: the button says "Delete …". */
function crossIcon(): SVGSVGElement {
  const svg = document.createElementNS(SVG, 'svg');
  svg.setAttribute('class', 'seqv-pin');
  svg.setAttribute('aria-hidden', 'true');
  svg.setAttribute('viewBox', '0 0 16 16');
  const path = document.createElementNS(SVG, 'path');
  for (const [k, v] of Object.entries({
    d: CROSS,
    fill: 'none',
    stroke: 'currentColor',
    'stroke-width': '1.4',
    'stroke-linecap': 'round',
  }))
    path.setAttribute(k, v);
  svg.append(path);
  return svg;
}

export const AnnotationsPanel: DomPart<AnnotationsPanelProps> = (host, ctx, props) => {
  let p = props;
  const at = spot(host);
  const label = h('input', {
    type: 'text',
    class: 'seqv-annotations-label',
    'aria-label': 'Annotation label',
    placeholder: 'Label',
    'data-seqv-key': 'label',
  });
  const color = h('input', {
    type: 'color',
    class: 'seqv-annotations-color',
    'aria-label': 'Annotation color',
    'data-seqv-key': 'color',
  });
  color.value = ANNOTATION_COLOR;
  const style = h(
    'select',
    { class: 'seqv-annotations-style', 'aria-label': 'Annotation line style', 'data-seqv-key': 'style' },
    ...ANNOTATION_STYLES.map((x) => h('option', { value: x }, x)),
  );
  style.value = ANNOTATION_STYLE;
  const lineStyle = () => style.value as LineStyle;
  style.addEventListener('change', () => refresh());
  // A place of markers (`kind` `points`) takes a glyph in the line style's place.
  const shape = h(
    'select',
    { class: 'seqv-annotations-style', 'aria-label': 'Annotation glyph', 'data-seqv-key': 'glyph' },
    ...ANNOTATION_GLYPHS.map((x) => h('option', { value: x }, x)),
  );
  shape.value = ANNOTATION_GLYPH;
  shape.addEventListener('change', () => refresh());
  const add = button(
    { class: 'seqv-switch-trigger seqv-annotations-add', 'data-seqv-key': 'add' },
    () => run(),
    'Add from selection',
  );
  const where = h('fieldset', { class: 'seqv-switch seqv-annotations-place', 'aria-label': 'Annotation placement' });
  const status = h('span', { class: 'seqv-annotations-status seqv-quiet', 'aria-live': 'polite' });
  const body = h('div', { class: 'seqv-annotations-body' });
  const addRow = h('div', { class: 'seqv-annotations-add-row' }, color, style, label, add);
  const root = h('div', { class: 'seqv-annotations', 'data-seqv-part': 'annotations' }, addRow, status, body);
  const places = () => annotationPlaces(p.section, p.data, p.places);
  /** The set the next add goes to: the chosen place while it is one of the places, else the first. */
  let chosen = '';
  const place = () => {
    const all = places();
    return all.find((x) => x.data === chosen) ?? all[0];
  };
  const target = () => place().data;
  const glyph = () => (placeGlyphs(place()) ? shape.value : undefined);
  /** The line style's picker, or the glyph's on a place of markers. */
  const picker = () => {
    const want = glyph() === undefined ? style : shape;
    if (want.isConnected) return;
    (want === style ? shape : style).replaceWith(want);
  };
  const marks = () => {
    const now = target();
    for (const b of where.querySelectorAll('button')) setAttr(b, 'aria-pressed', b.dataset.seqvPlace === now);
  };
  let drawn = '';
  const choices = () => {
    const all = places();
    const sig = all.length > 1 ? JSON.stringify(all) : '';
    if (sig !== drawn) {
      drawn = sig;
      where.replaceChildren(
        ...all.map((x, i) =>
          button(
            {
              class: 'seqv-choice',
              'data-seqv-key': `place:${i}`,
              'data-seqv-group': 'annotation-place',
              'data-seqv-place': x.data,
            },
            () => {
              chosen = x.data;
              marks();
              picker();
              refresh();
            },
            x.label,
          ),
        ),
      );
      if (sig) addRow.prepend(where);
      else where.remove();
    }
    marks();
  };
  const run = () => {
    if (add.getAttribute('aria-disabled') === 'true') return;
    if (addAnnotation(ctx, p.section, label.value, color.value, target(), lineStyle(), glyph())) label.value = '';
  };
  label.addEventListener('keydown', (e) => {
    if (e.key !== 'Enter') return;
    e.preventDefault();
    run();
  });
  const list = spot(body);

  const row = (a: AnnotationEntry) => {
    const swatch = h('span', { class: 'seqv-swatch', 'aria-hidden': 'true' });
    swatch.style.background = a.color;
    const go = button(
      {
        class: 'seqv-annotation-go',
        'data-seqv-key': `go:${a.key}`,
        'data-seqv-group': 'annotation-go',
        'aria-label': showName(a),
      },
      () => revealAnnotation(ctx, p.section, a.id, a.data),
      swatch,
      h('span', { class: 'seqv-annotation-label' }, a.label),
      a.glyph
        ? h('span', { class: 'seqv-quiet seqv-annotation-style' }, a.glyph)
        : a.lineStyle !== ANNOTATION_STYLE
          ? h('span', { class: 'seqv-quiet seqv-annotation-style' }, a.lineStyle)
          : null,
      a.place ? h('span', { class: 'seqv-quiet seqv-annotation-place' }, a.place) : null,
      h('span', { class: 'seqv-quiet' }, a.span),
    );
    const del = button(
      {
        class: 'seqv-icon',
        'data-seqv-key': `del:${a.key}`,
        'data-seqv-group': 'annotation-delete',
        'aria-label': deleteName(a),
      },
      () => removeAnnotation(ctx, p.section, a.id, a.data),
      crossIcon(),
    );
    return h('li', { class: 'seqv-annotation' }, go, del);
  };
  const show = (entries: readonly AnnotationEntry[]) => {
    setText(
      status,
      entries.length === 0 ? 'No annotations' : `${entries.length} annotation${entries.length > 1 ? 's' : ''}`,
    );
    list.set(
      entries.length === 0
        ? null
        : h('ul', { class: 'seqv-annotation-list', 'aria-label': 'Annotations' }, ...entries.map(row)),
    );
  };
  const refresh = () => {
    setAttr(
      add,
      'aria-disabled',
      annotationTarget(ctx, p.section, label.value, color.value, target(), lineStyle(), glyph()) === null,
    );
  };
  at.set(root);
  choices();
  picker();
  const subs = [
    ...derived([ctx.document.shown], () => annotationEntries(ctx, places()), show, sameData),
    ...watch(annotationStates(ctx), refresh),
  ];
  return {
    update(next) {
      p = next;
      choices();
      picker();
      show(annotationEntries(ctx, places()));
      refresh();
    },
    dispose() {
      off([...subs, list, at]);
    },
  };
};
