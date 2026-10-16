/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `Annotations` part (seqspec:interaction#the-default-ui), React: the DOM part's twin (src/dom/annotations.ts) from
// the same helpers (src/annotations.ts), with its structure, class names and text. A thin client of the patch API; the
// label and the color are leaf state of the part.

import type { LineStyle } from '@molstar/seq-core/format';
import { useMemo, useRef, useState } from 'react';
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
import { CROSS } from '../dom/shared.js';
import { sameData } from '../same.js';
import { type SectionProps, useDerived, useLanding } from './hooks.js';

export interface AnnotationsPanelProps extends SectionProps {
  /** Another set than `<section>.bookmarks` (`bookmark.add`'s `data`). */
  readonly data?: string;
  /**
   * Several sets to choose between (each `bookmark.add`'s `data`, with the name a person sees): the add row offers the
   * choice, the list shows the sets one after another with their names, a delete acts on its row's set. Overrides `data`.
   */
  readonly places?: readonly AnnotationPlace[];
}

/** A row: its range (reveals it) and a delete; focus a removed row held lands on the next one (`useLanding`). */
function Annotation({ a, onGo, onDelete }: { a: AnnotationEntry; onGo(): void; onDelete(): void }) {
  const go = useRef<HTMLButtonElement>(null);
  const del = useRef<HTMLButtonElement>(null);
  useLanding(go);
  useLanding(del);
  return (
    <li className="seqv-annotation">
      <button
        ref={go}
        type="button"
        className="seqv-annotation-go"
        data-seqv-key={`go:${a.key}`}
        data-seqv-group="annotation-go"
        aria-label={showName(a)}
        onClick={onGo}
      >
        <span className="seqv-swatch" aria-hidden="true" style={{ background: a.color }} />
        <span className="seqv-annotation-label">{a.label}</span>
        {a.glyph ? (
          <span className="seqv-quiet seqv-annotation-style">{a.glyph}</span>
        ) : (
          a.lineStyle !== ANNOTATION_STYLE && <span className="seqv-quiet seqv-annotation-style">{a.lineStyle}</span>
        )}
        {a.place && <span className="seqv-quiet seqv-annotation-place">{a.place}</span>}
        <span className="seqv-quiet">{a.span}</span>
      </button>
      <button
        ref={del}
        type="button"
        className="seqv-icon"
        data-seqv-key={`del:${a.key}`}
        data-seqv-group="annotation-delete"
        aria-label={deleteName(a)}
        onClick={onDelete}
      >
        <svg className="seqv-pin" aria-hidden="true" viewBox="0 0 16 16">
          <path d={CROSS} fill="none" stroke="currentColor" strokeWidth={1.4} strokeLinecap="round" />
        </svg>
      </button>
    </li>
  );
}

export function AnnotationsPanel({ context, section, data, places: given }: AnnotationsPanelProps) {
  const [label, setLabel] = useState('');
  const [color, setColor] = useState(ANNOTATION_COLOR);
  const [lineStyle, setLineStyle] = useState<LineStyle>(ANNOTATION_STYLE);
  const [glyphChosen, setGlyph] = useState(ANNOTATION_GLYPH);
  const [chosen, setChosen] = useState('');
  const key = given?.map((x) => `${x.data}\n${x.label}\n${x.kind ?? ''}`).join('\n');
  // The places as given, kept while their content is the same (a host may pass a new array each render).
  const places = useMemo(() => annotationPlaces(section, data, given), [section, data, key]);
  const place = places.find((x) => x.data === chosen) ?? places[0];
  const to = place.data;
  // A place of markers takes a glyph, the others a line style.
  const glyph = placeGlyphs(place) ? glyphChosen : undefined;
  const entries = useDerived([context.document.shown], () => annotationEntries(context, places), {
    compare: sameData,
    deps: [context, places],
  });
  const idle = useDerived(
    annotationStates(context),
    () => annotationTarget(context, section, label, color, to, lineStyle, glyph) === null,
    { deps: [context, section, to, label, color, lineStyle, glyph] },
  );
  const run = () => {
    if (!idle && addAnnotation(context, section, label, color, to, lineStyle, glyph)) setLabel('');
  };
  return (
    <div className="seqv-annotations" data-seqv-part="annotations">
      <div className="seqv-annotations-add-row">
        {places.length > 1 && (
          <fieldset aria-label="Annotation placement" className="seqv-switch seqv-annotations-place">
            {places.map((x, i) => (
              <button
                key={x.data}
                type="button"
                className="seqv-choice"
                aria-pressed={x.data === to}
                data-seqv-key={`place:${i}`}
                data-seqv-group="annotation-place"
                data-seqv-place={x.data}
                onClick={() => setChosen(x.data)}
              >
                {x.label}
              </button>
            ))}
          </fieldset>
        )}
        <input
          type="color"
          className="seqv-annotations-color"
          aria-label="Annotation color"
          data-seqv-key="color"
          value={color}
          onChange={(e) => setColor(e.target.value)}
        />
        {glyph === undefined ? (
          <select
            className="seqv-annotations-style"
            aria-label="Annotation line style"
            data-seqv-key="style"
            value={lineStyle}
            onChange={(e) => setLineStyle(e.target.value as LineStyle)}
          >
            {ANNOTATION_STYLES.map((x) => (
              <option key={x} value={x}>
                {x}
              </option>
            ))}
          </select>
        ) : (
          <select
            className="seqv-annotations-style"
            aria-label="Annotation glyph"
            data-seqv-key="glyph"
            value={glyph}
            onChange={(e) => setGlyph(e.target.value)}
          >
            {ANNOTATION_GLYPHS.map((x) => (
              <option key={x} value={x}>
                {x}
              </option>
            ))}
          </select>
        )}
        <input
          type="text"
          className="seqv-annotations-label"
          aria-label="Annotation label"
          placeholder="Label"
          data-seqv-key="label"
          value={label}
          onChange={(e) => setLabel(e.target.value)}
          onKeyDown={(e) => {
            if (e.key !== 'Enter') return;
            e.preventDefault();
            run();
          }}
        />
        <button
          type="button"
          className="seqv-switch-trigger seqv-annotations-add"
          data-seqv-key="add"
          aria-disabled={idle}
          onClick={run}
        >
          Add from selection
        </button>
      </div>
      <span className="seqv-annotations-status seqv-quiet" aria-live="polite">
        {entries.length === 0 ? 'No annotations' : `${entries.length} annotation${entries.length > 1 ? 's' : ''}`}
      </span>
      <div className="seqv-annotations-body">
        {entries.length > 0 && (
          <ul className="seqv-annotation-list" aria-label="Annotations">
            {entries.map((a) => (
              <Annotation
                key={a.key}
                a={a}
                onGo={() => revealAnnotation(context, section, a.id, a.data)}
                onDelete={() => removeAnnotation(context, section, a.id, a.data)}
              />
            ))}
          </ul>
        )}
      </div>
    </div>
  );
}
