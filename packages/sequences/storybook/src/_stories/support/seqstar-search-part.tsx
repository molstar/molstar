/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Case 18 in the Viewer cases story: p53's Tracks view with the viewer's own `Search` part, in the bar (the UI spec
// gives `search(text, signal)`) or placed by the app above the viewer (the part's own `search` prop). The search is the
// story's, not the viewer's: the names of the items of p53's annotations through a `TextIndex` (anywhere in the name),
// and the letters of the sequence through an exact matcher, or a regular expression when the text is one (seqstar-story-search.ts). Its answer goes to `ctx.results`;
// the part steps through it and keeps it as a location set.

import { Search, SeqViewerLayout, SeqViewerProvider } from '@molstar/seq-viewer/react';
import { type RefObject, use, useMemo, useRef, useState } from 'react';
import { Choices } from '../ui/index.js';
import { compiledCase, described } from './seqstar-cases.js';
import f from './seqstar-format.module.css';
import { PrismCanvas, prism } from './seqstar-prism.js';
import q from './seqstar-search.module.css';
import { storySearch } from './seqstar-story-search.js';
import { useCaseContext, useStoryUi, type ViewerOptions } from './seqstar-viewer.js';

/** What the story suggests: each finds something the Tracks view shows (names of items, then letters). */
const TRY = [
  { text: 'zinc', hint: 'the four zinc ligands' },
  { text: 'DNA', hint: 'the DNA-binding domain and the DNA-contact residues' },
  { text: 'RHSVV', hint: 'letters: where p53 reads RHSVV' },
  { text: 'tetra', hint: 'the tetramerization domain' },
  { text: '[KR]{2}.{10,12}[KR]{3,}', hint: 'a regular expression: the bipartite nuclear localization signal' },
];

/** Types `text` in the first search box under `root`, as a person would: the part's own handlers take it from there. */
function fill(root: HTMLElement | null, text: string) {
  const input = root?.querySelector<HTMLInputElement>('.seqv-search-input');
  if (!input) return;
  Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')?.set?.call(input, text);
  input.dispatchEvent(new Event('input', { bubbles: true }));
  input.focus();
}

/** The suggested searches, each a chip that fills the box under `root` and searches. */
function Try({ root }: { root: RefObject<HTMLElement | null> }) {
  return (
    <div className={q.chips}>
      <span className={q.caption}>Try</span>
      {TRY.map((t) => (
        <button key={t.text} type="button" className={q.chip} onClick={() => fill(root.current, t.text)}>
          {t.text}
          <span className={q.at}>{t.hint}</span>
        </button>
      ))}
    </div>
  );
}

const PLACES = [
  { id: 'bar', label: "In the viewer's bar" },
  { id: 'app', label: 'Placed by the app, above the viewer' },
];

/** Case 18's viewer: the Search part where the choice puts it, the suggested searches, the legend (kept results). */
export function Case18Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase('case18'));
  const doc = useMemo(() => described('case18', c.doc), [c]);
  const [place, setPlace] = useState('bar');
  // The legend lists the kept results under Highlights (left out by `clean`).
  const m = useCaseContext(doc, { ...options, view: 'tracks', navigator: 'strip', parts: ['Legend'] });
  const search = useMemo(() => (m ? storySearch(m.ctx) : undefined), [m]);
  const extra = useMemo(() => (place === 'bar' ? { search } : {}), [place, search]);
  const ui = useStoryUi(m, extra);
  const card = useRef<HTMLElement>(null);
  return (
    <section className={f.card} ref={card}>
      <h2 className={f.heading}>p53, its tracks and location sets</h2>
      <div className={prism.controls}>
        <div className={prism.wide}>
          <Choices label="The Search part" choices={PLACES} active={place} onPick={setPlace} />
        </div>
      </div>
      <Try root={card} />
      {m && search ? (
        <SeqViewerProvider context={m.ctx} ui={ui}>
          {place === 'app' ? (
            <div className={q.outside}>
              <Search context={m.ctx} search={search} />
            </div>
          ) : null}
          <PrismCanvas>
            <SeqViewerLayout />
          </PrismCanvas>
        </SeqViewerProvider>
      ) : (
        <div className={f.status}>Loading…</div>
      )}
    </section>
  );
}
