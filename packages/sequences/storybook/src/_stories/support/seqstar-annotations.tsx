/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The controls of the Seq*/Text view stories (text-view.stories.tsx: Publication tracks, Alignment figure), over their
// live viewer: Text or Tracks (the one document drawn by two presets of its section, `ctx.ui.setStyle`: the document is
// not rebuilt, so what a person adds shows in both), and the `Annotations` part (`@molstar/seq-viewer/react`) in a
// column beside the viewer (`ExportCard`'s `aside`): select a range, choose a place (Above or Below, a Frame, a Box on a
// row, Markers), add it with a label, a color and a line style or a glyph, and it is a row of that place's bookmark set,
// drawn by the set's track.

import { type AnnotationPlace, AnnotationsPanel } from '@molstar/seq-viewer/react';
import { useEffect } from 'react';
import { Choices } from '../ui/index.js';
import { PrismCanvas, prism } from './seqstar-prism.js';
import {
  PUBLICATION_PLACES,
  PUBLICATION_PRESETS,
  PUBLICATION_SECTION,
  type PublicationView,
} from './seqstar-publication.js';
import f from './seqstar-format.module.css';
import type { Mounted } from './seqstar-viewer.js';

export function PublicationControls({
  m,
  view,
  onView,
  section = PUBLICATION_SECTION,
  presets = PUBLICATION_PRESETS,
}: {
  m: Mounted | null;
  view: PublicationView;
  onView(view: PublicationView): void;
  /** Another story's section and its two presets (default: the publication's). */
  section?: string;
  presets?: Readonly<Record<PublicationView, string | null>>;
}) {
  // UI state only: the document, its patches and its undo history stay as they are. A new context takes it again.
  useEffect(() => m?.ctx.ui.setStyle(section, presets[view]), [m, view, section, presets]);
  return (
    <div className={prism.narrow}>
      <Choices
        label="View"
        choices={[
          { id: 'text', label: 'Text', title: 'The wrapped text, tracks under every line' },
          { id: 'tracks', label: 'Tracks', title: 'The same tracks on one line, zoom and pan' },
        ]}
        active={view}
        onPick={(id) => onView(id as PublicationView)}
      />
    </div>
  );
}

/** The column beside the viewer: the Annotations panel, with the way to use it (each add, delete or Undo is a patch). */
export function PublicationAnnotations({
  m,
  section = PUBLICATION_SECTION,
  places = PUBLICATION_PLACES,
}: {
  m: Mounted;
  section?: string;
  places?: readonly AnnotationPlace[];
}) {
  return (
    <>
      <div className={f.subheading}>Annotations: select a range in the viewer, then add it</div>
      <PrismCanvas>
        <div style={{ paddingInline: 'var(--msp-inset)' }}>
          <AnnotationsPanel context={m.ctx} section={section} places={places} />
        </div>
      </PrismCanvas>
    </>
  );
}
