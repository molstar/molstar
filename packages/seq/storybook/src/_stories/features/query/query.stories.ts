/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
import { Button } from '@molstar/ui/components/button/button';
import { Segmented } from '@molstar/ui/components/segmented/segmented';
import { TextInput } from '@molstar/ui/components/text-input/text-input';
import { defineElements, onCleanup } from '@molstar/ui/dom';
import { UiElement } from '@molstar/ui/elements/ui-element';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import {
  clearQuery,
  composeViewer,
  type Emphasis,
  type QueryExample,
  QueryExamples,
  QueryRunner,
  type QueryReport,
  readableRanges,
  rerunWhileWaiting,
  runQuery,
  statusLine,
} from './compose.js';

// Query (Features): the address language over p53's data and the globins family. Whatever is typed is parsed, checked and
// resolved as typed by the library's own paths (`parseAddress`, `printAddress`, the viewer's address model: no parser of
// the story's), and its result is selected or highlighted in the view. Below the box: the canonical (printed) form, the
// result's kind and size, where it is, and, when it does not resolve, its issues with their hints. Each worked example is
// a story of its own, with a few queries to try; every string is put in as text, never markup.

const Mono = 'font: var(--msui-text-xs) / var(--msui-leading-xs) var(--msui-font-mono);';
const Muted = 'color: var(--msui-text-muted);';

/** What a story shows: the view its data comes from, and an example (the free story has none). */
interface QueryModel {
  readonly view: ExampleView;
  readonly example?: QueryExample;
}

const viewOf = (id: string) => FeatureViews.find((v) => v.id === id)!;

/** The queries the free story offers: the first of each example on p53's view. */
const FreeQueries = QueryExamples.filter((e) => e.view === 'query').flatMap((e) => e.queries);

/** How long typing rests before the query runs. */
const DebounceMs = 100;

class QueryRoot extends UiElement<QueryModel> {
  static readonly tag = 'seqx-query-root';

  protected render({ view, example }: QueryModel): HTMLElement {
    const root = document.createElement('div');
    root.style.cssText =
      'display: grid; gap: var(--msui-space-4); padding: var(--msui-space-5); min-height: 100vh; align-content: start; background: var(--msui-bg); color: var(--msui-text); font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);';
    const about = document.createElement('p');
    about.style.cssText = `margin: 0; max-width: 72ch; ${Muted}`;
    about.textContent = example?.about ?? view.summary;

    const queries = example?.queries ?? FreeQueries;
    let emphasis: Emphasis = example?.emphasis ?? 'select';
    const text = new State(queries[0].query);
    // forward reference: the runner runs what `runNow` below draws
    const runner = new QueryRunner(() => runNow(), DebounceMs);

    const box = document.createElement('div');
    box.style.cssText = 'display: flex; gap: var(--msui-space-2); align-items: center; flex-wrap: wrap;';
    const input = TextInput({
      value: text,
      placeholder: 'An address, such as /data/p53/(range :from 94 :to 292)',
      size: 'md',
      'aria-label': 'Query',
      onInput: (v) => {
        text.set(v);
        runner.soon();
      },
      // Enter applies the text at once, as a chip does
      onSubmit: (v) => {
        text.set(v);
        runner.now();
      },
    }) as HTMLElement;
    input.style.cssText += ' flex: 1 1 28em; min-width: 16em; font-family: var(--msui-font-mono);';
    const mode = Segmented<Emphasis>({
      options: [
        { value: 'select', label: 'Select' },
        { value: 'highlight', label: 'Highlight' },
      ],
      value: emphasis,
      size: 'sm',
      'aria-label': 'What the result does in the view',
      onChange: (v) => {
        emphasis = v;
        runner.now();
      },
    }) as HTMLElement;
    box.append(input, mode);

    const chips = document.createElement('div');
    chips.style.cssText = 'display: flex; gap: var(--msui-space-2); flex-wrap: wrap;';
    for (const q of queries) {
      chips.append(
        Button({
          size: 'sm',
          children: q.label,
          title: q.query,
          onClick: () => {
            text.set(q.query);
            runner.now();
          },
        }) as HTMLElement,
      );
    }

    const canvasHost = document.createElement('div');
    canvasHost.style.cssText = 'min-width: 0;';
    // told to a screen reader as one line, whole (`aria-atomic`); the report below it is read on demand, not announced
    const status = document.createElement('p');
    status.style.cssText = 'margin: 0; font-weight: 600;';
    status.setAttribute('role', 'status');
    status.setAttribute('aria-live', 'polite');
    status.setAttribute('aria-atomic', 'true');
    const out = document.createElement('div');
    out.style.cssText = 'display: grid; gap: var(--msui-space-1); max-width: 96ch;';
    root.append(about, box, chips, canvasHost, status, out);

    const viewer = composeViewer(view, { target: canvasHost });
    // no run after the story is gone (a debounce, a load still pending); the query's marks go before the viewer does
    let disposed = false;
    onCleanup(() => {
      disposed = true;
      runner.dispose();
      try {
        clearQuery(viewer);
      } finally {
        viewer.retire();
      }
    });

    const line = (label: string, value: string, mono = true) => {
      const row = document.createElement('div');
      row.style.cssText = 'display: flex; gap: var(--msui-space-3); align-items: baseline; flex-wrap: wrap;';
      const k = document.createElement('span');
      k.style.cssText = `min-width: 7em; ${Muted}`;
      k.textContent = label;
      const v = document.createElement(mono ? 'code' : 'span');
      if (mono) v.style.cssText = Mono;
      v.textContent = value;
      row.append(k, v);
      return row;
    };

    const draw = (r: QueryReport) => {
      status.textContent = statusLine(r);
      const rows: HTMLElement[] = [];
      if (r.status !== 'empty') {
        rows.push(line('Printed', r.printed ?? 'It does not parse'));
        for (const s of r.systems) rows.push(line('On', `${s.system}  ${readableRanges(s.ranges)}`));
        for (const x of r.rects)
          rows.push(
            line('Rect', `${x.a.system} ${readableRanges(x.a.ranges)}  ×  ${x.b.system} ${readableRanges(x.b.ranges)}`),
          );
        for (const i of r.issues) {
          const at = i.offset !== undefined ? ` (at character ${i.offset + 1})` : '';
          const row = line(i.severity, `${i.message}${at}`, false);
          if (i.severity === 'error') (row.lastElementChild as HTMLElement).style.color = 'var(--msui-danger)';
          rows.push(row);
          if (i.hint) rows.push(line('Hint', i.hint, false));
          rows.push(line('Code', i.code));
        }
      }
      out.replaceChildren(...rows);
    };

    let last: QueryReport | undefined;
    const runNow = () => {
      if (disposed) return;
      last = runQuery(viewer, text.value, emphasis);
      draw(last);
    };
    // a query waiting for its data runs again when the view rests with it arrived
    onCleanup(rerunWhileWaiting(viewer, () => last, runNow));
    viewer.load(view.build()).then(
      () => runner.now(),
      (e: unknown) => {
        if (disposed) return;
        draw({
          text: text.value,
          status: 'failed',
          systems: [],
          rects: [],
          issues: [
            { code: 'story.query.load', severity: 'error', message: e instanceof Error ? e.message : String(e) },
          ],
        });
      },
    );
    return root;
  }
}

defineElements([QueryRoot]);

const meta: Meta = { title: 'Features/Query' };
export default meta;

function story(model: QueryModel): StoryObj {
  return {
    render: () => {
      const el = document.createElement(QueryRoot.tag) as QueryRoot;
      el.model = model;
      return el;
    },
  };
}

const exampleStory = (id: string): StoryObj => {
  const example = QueryExamples.find((e) => e.id === id)!;
  return story({ view: viewOf(example.view), example });
};

/** A query box over p53: type an address; it is parsed, checked and resolved as typed, its result selected or highlighted. */
export const QueryBox: StoryObj = story({ view: viewOf('query') });

/** Items by category: the domains a field names, as elements drawn by their extents. */
export const ItemsByCategory: StoryObj = exampleStory('category');

/** Ranges and slices: labels 1-based, `index` 0-based, a single residue, a range on the gene. */
export const RangesAndSlices: StoryObj = exampleStory('ranges');

/** A result on the gene: the protein's selection projected through the codon map. */
export const ResultOnTheGene: StoryObj = exampleStory('on-gene');

/** Merge, intersect, except: combinators over positions, and the issue a mixed query gives. */
export const AndOrNot: StoryObj = exampleStory('combined');

/** An alignment's rows and columns: a row, columns, and a rect of both. */
export const AlignmentRowsAndColumns: StoryObj = exampleStory('alignment');
