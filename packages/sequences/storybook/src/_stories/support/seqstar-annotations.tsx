/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The controls of the Seq*/Text view › Publication tracks story (text-view.stories.tsx), over its live viewer: Text or
// Tracks (the one document drawn by two presets of its section, `ctx.ui.setStyle`: the document is not rebuilt, so what
// a person adds shows in both), and the `Annotations` part (`@molstar/seq-viewer/react`) in a column beside the viewer
// (`ExportCard`'s `aside`): select a range, choose Above or Below, add it with a label and a color, and it is a row of
// one of the section's two bookmarks sets, drawn by that set's ranges track (before or after the letters).

import { AnnotationsPanel } from '@molstar/seq-viewer/react';
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
}: {
  m: Mounted | null;
  view: PublicationView;
  onView(view: PublicationView): void;
}) {
  // UI state only: the document, its patches and its undo history stay as they are. A new context takes it again.
  useEffect(() => m?.ctx.ui.setStyle(PUBLICATION_SECTION, PUBLICATION_PRESETS[view]), [m, view]);
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
export function PublicationAnnotations({ m }: { m: Mounted }) {
  return (
    <>
      <div className={f.subheading}>Annotations: select a range in the viewer, then add it</div>
      <PrismCanvas>
        <div style={{ paddingInline: 'var(--msp-inset)' }}>
          <AnnotationsPanel context={m.ctx} section={PUBLICATION_SECTION} places={PUBLICATION_PLACES} />
        </div>
      </PrismCanvas>
    </>
  );
}
