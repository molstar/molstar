/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, State } from '@molstar/core/reactive';
import type { SearchStrands } from '@molstar/seqx-core/search';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { ToolbarMount } from '@molstar/seqx-viewer/viewer-spec';
import { Button } from '@molstar/ui/components/button/button';
import { Segmented } from '@molstar/ui/components/segmented/segmented';
import { TextInput } from '@molstar/ui/components/text-input/text-input';
import { DomHost } from '@molstar/ui/dom/host';
import { render } from '@molstar/ui/dom/render';
import { toolbarWithHiddenTracks } from './toolbar.js';

// The stories' search box (the host composition, as the hidden tracks' list): the library ships the search's model and
// commands, never a box, so this one is built on the public surface alone (a test holds its imports to public exports).
// An input, a matcher choice (IUPAC, exact, regex), a strand switch on nucleotides, the count ("2 of 14", "14+" while
// running, "No hits"), previous and next, Save and Clear. The text runs `search.run` once it rests 200 ms (a run in
// progress cancelled by the next); Enter and Shift+Enter step; Esc empties the box and runs `search.clear`.

/** How long the text rests before it is searched, ms. */
export const SearchRestMs = 200;

const Matchers = [
  { value: 'degenerate', label: 'IUPAC' },
  { value: 'exact', label: 'Exact' },
  { value: 'regex', label: 'Regex' },
];
const Strands = [
  { value: 'both' as SearchStrands, label: 'Both strands' },
  { value: 'forward' as SearchStrands, label: 'Forward' },
];

/** The count as the box shows it: "k of n", "n hits", "n+" while running, "No hits", or nothing before a search. */
export function searchCount(viewer: SeqViewer): string {
  const s = viewer.search;
  if (!s?.query.value) return '';
  const n = s.count.value;
  const fmt = (x: number) => x.toLocaleString('en-US');
  if (s.status.value === 'running') return n ? `${fmt(n)}+` : 'Searching…';
  if (!n) return 'No hits';
  const k = s.current.value;
  const more = s.complete.value ? '' : '+';
  return k >= 0 ? `${fmt(k + 1)} of ${fmt(n)}${more}` : `${fmt(n)}${more} hit${n === 1 ? '' : 's'}`;
}

/** Whether the view draws a nucleotide sequence (the strand switch's): read from the model the view resolves over. */
function hasNucleotides(viewer: SeqViewer): boolean {
  for (const s of viewer.view.resolved.value?.sections ?? []) {
    const a = viewer.spec.model.value.systems.get(s.axis.system)?.alphabet;
    if (a && a !== 'protein') return true;
  }
  return false;
}

export interface SearchBoxOptions {
  /** A pattern searched once the viewer has its spec. */
  readonly initial?: string;
  readonly matcher?: string;
}

/** Mounts the box into `target` for `viewer`; the function it returns takes it away. */
export function mountSearchBox(viewer: SeqViewer, target: HTMLElement, options: SearchBoxOptions = {}): () => void {
  const doc = target.ownerDocument;
  const host = doc.createElement('div');
  host.className = 'seqx-search-box';
  host.style.cssText = 'display: flex; flex-wrap: wrap; align-items: center; gap: var(--msui-space-2, 8px);';
  target.append(host);
  const search = viewer.search;
  if (!search) {
    host.textContent = 'This viewer does not list the search behavior.';
    return () => host.remove();
  }
  const text = new State(options.initial ?? '');
  const matcher = new State(options.matcher ?? viewer.config.search.matcher);
  const strands = new State<SearchStrands>(viewer.config.search.strands);
  const nucleotides = new State(hasNucleotides(viewer));
  // every state the count reads: one derived text
  const count = derive([search.query, search.status, search.count, search.current, search.complete], () =>
    searchCount(viewer),
  );
  const noHits = derive(search.hits, (r) => r.hits.length === 0);

  let rest: unknown;
  const offs: (() => void)[] = [];
  const runNow = () => {
    DomHost.clearTimeout(rest);
    rest = undefined;
    const pattern = text.value.trim();
    if (!pattern) {
      void viewer.execute('search.clear', {});
      return;
    }
    void viewer.execute('search.run', { pattern, matcher: matcher.value, strands: strands.value });
  };
  // the text searched once it rests (the previous run cancelled by the next)
  const soon = () => {
    DomHost.clearTimeout(rest);
    rest = DomHost.setTimeout(runNow, SearchRestMs);
  };
  const step = (back: boolean) => {
    // a pattern typed and not yet searched is searched first
    if (rest !== undefined) runNow();
    else void viewer.execute(back ? 'search.previous' : 'search.next', {});
  };

  const dispose = render(() => {
    const input = TextInput({
      value: text,
      search: true,
      placeholder: 'Search a motif',
      'aria-label': 'Search the sequences',
      onInput: (v) => {
        text.set(v);
        if (!v) runNow();
        else soon();
      },
    });
    // Enter steps on, Shift+Enter back (the input's own Enter is its submit)
    input.addEventListener('keydown', (e) => {
      if (e.key !== 'Enter' || e.isComposing) return;
      e.preventDefault();
      step(e.shiftKey);
    });
    const kinds = Segmented<string>({
      options: Matchers,
      value: matcher,
      label: 'Matcher',
      size: 'sm',
      onChange: () => text.value && runNow(),
    } as never);
    const strand = Segmented<SearchStrands>({
      options: Strands,
      value: strands,
      label: 'Strands',
      size: 'sm',
      onChange: () => text.value && runNow(),
    } as never);
    const strandWrap = doc.createElement('span');
    strandWrap.append(strand as Node);
    const showStrand = (on: boolean) => {
      strandWrap.style.display = on ? '' : 'none';
    };
    showStrand(nucleotides.value);
    offs.push(nucleotides.listen(showStrand));
    const counter = doc.createElement('span');
    counter.className = 'seqx-search-count';
    counter.setAttribute('aria-live', 'off');
    counter.style.cssText =
      'min-width: 9ch; font: var(--msui-text-sm, 13px) var(--msui-font-sans, sans-serif); color: var(--msui-text-muted, inherit); font-variant-numeric: tabular-nums;';
    counter.textContent = count.value;
    offs.push(
      count.listen((t) => {
        counter.textContent = t;
      }),
    );
    const previous = Button({
      size: 'sm',
      children: 'Previous',
      title: 'Previous hit (Shift+Enter)',
      disabled: noHits,
      onClick: () => step(true),
    });
    const next = Button({
      size: 'sm',
      children: 'Next',
      title: 'Next hit (Enter)',
      disabled: noHits,
      onClick: () => step(false),
    });
    const save = Button({
      size: 'sm',
      variant: 'ghost',
      children: 'Save',
      title: 'Save the hits as a set',
      disabled: noHits,
      onClick: () => void viewer.execute('search.save', {}),
    });
    const clear = Button({
      size: 'sm',
      variant: 'ghost',
      children: 'Clear',
      title: 'Clear the search',
      onClick: () => {
        text.set('');
        runNow();
      },
    });
    const wrap = doc.createElement('div');
    wrap.style.cssText = 'display: contents;';
    wrap.append(
      input as Node,
      kinds as Node,
      strandWrap,
      counter,
      previous as Node,
      next as Node,
      save as Node,
      clear as Node,
    );
    return wrap;
  }, host);

  // the strand switch follows the view (a protein-only view has none)
  const offView = viewer.view.resolved.listen(() => nucleotides.set(hasNucleotides(viewer)));
  if (options.initial) {
    void viewer.settle().then(() => {
      if (text.value === options.initial) runNow();
    });
  }
  return () => {
    DomHost.clearTimeout(rest);
    offView();
    for (const off of offs) off();
    dispose();
    host.remove();
  };
}

/**
 * `ui.toolbar` for a story with the box: the stories' toolbar (`toolbarWithHiddenTracks`), the search box under it.
 */
export function toolbarWithSearch(options: SearchBoxOptions = {}): ToolbarMount {
  return (viewer, target) => {
    const offBar = toolbarWithHiddenTracks(viewer, target);
    const slot = target.ownerDocument.createElement('div');
    slot.style.cssText = 'padding: var(--msui-space-2, 8px) 0;';
    target.append(slot);
    const offBox = mountSearchBox(viewer, slot, options);
    return () => {
      offBox();
      slot.remove();
      offBar();
    };
  };
}
