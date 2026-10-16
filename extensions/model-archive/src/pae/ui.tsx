/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Select, Show } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { SceneNode } from '@molstar/prism-viewer/model/scene';
import { SceneRow } from '@molstar/prism-viewer/views/structure/scene-row';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { combineLatest, distinctUntilChanged, map, of, startWith } from 'rxjs';
import { type PaeMatrix, paePairTarget, paeRangeLabel, paeRangeTarget } from './data.js';
import type { PaeSize } from './operations.js';
import { PaePlot } from './plot.js';

// The PAE card (parity.md §6, "PAE plot"; v5's Predicted Aligned Error panel, anew): the predicted aligned error of a
// model as a heat map (plot.tsx), as AlphaFold DB draws it. The cell under the pointer is lit in the 3D view and the
// sequence, a click selects its two residues, a box marks its aligned and scored ranges in both until cleared (marks.md);
// all through the API as the person. The matrix comes from the `model-archive` area; the card opens from the toolbar
// while a loaded model has a PAE, or by `model-archive.show-pae`; the plot is also a row under such a structure in the
// Structure panel, which opens the card at it.

export const PaeCard = 'model-archive.overlay.pae';
/** The PAE's row in the Structure panel, under a structure that has one. */
export const PaeRow = 'model-archive.row.pae';

/** "149 residues · 2 chains": what the row says of a PAE (its size, not its matrix: the card reads that). */
const paeDescription = (s: PaeSize | undefined) =>
  s && `${s.residues} residues${s.chains > 1 ? ` · ${s.chains} chains` : ''}`;

const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';

interface Shown {
  readonly ref: string;
  readonly label: string;
}

export function paeUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api['model-archive'];
  if (!area) return {};
  // The structures that have a PAE (by the file's metric list: no matrix is read for this), and those AlphaFold DB
  // has one for (an entry without it yet): asked again when the state changes (a fetched PAE is a state cell).
  const structures = combineLatest([ctx.sources.structures, area.changes.pipe(startWith(undefined))]).pipe(
    map(([list]) => list),
  );
  const withPae = structures.pipe(
    map((list) => list.filter((s) => area.hasPae(s.ref)).map((s): Shown => ({ ref: s.ref, label: s.label }))),
  );
  const fetchable = structures.pipe(
    map((list) => list.filter((s) => area.canFetchPae(s.ref)).map((s): Shown => ({ ref: s.ref, label: s.label }))),
  );
  const fetchPae = (structure: string) => void ctx.api['model-archive']?.fetchPae({ structure });
  const FetchButton = (props: { shown: Shown }) => (
    <Button size="sm" onClick={() => fetchPae(props.shown.ref)}>
      Fetch PAE from AlphaFold DB{props.shown.label ? ` (${props.shown.label})` : ''}
    </Button>
  );
  // The structure the card shows (by ref; none: the first with a PAE): its picker's, or the Structure panel row's.
  const chosen = observableState<string | undefined>(undefined);
  return {
    overlays: [
      {
        id: PaeCard,
        title: 'Predicted aligned error',
        render: domView(() => {
          const shown = combineLatest([withPae, chosen.asReadonly]).pipe(
            map(([list, ref]) => list.find((s) => s.ref === ref) ?? list[0]),
            distinctUntilChanged((a, b) => a?.ref === b?.ref),
          );
          // The card's text column (its title's: the card's inset and a row's own), on both sides: the plot is no row.
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0; padding: 0 calc(var(--msp-panel-inset) + var(--msp-inset)) calc(var(--msp-panel-inset) + var(--msp-inset))">
              <Show when={withPae.pipe(map((l) => (l.length > 1 ? l : undefined)))}>
                {(list) => (
                  <Select
                    size="sm"
                    aria-label="Structure"
                    options={list.map((s) => [s.ref, s.label] as const)}
                    value={shown.pipe(map((s) => s?.ref ?? ''))}
                    onChange={(ref) => chosen.set(ref)}
                  />
                )}
              </Show>
              <Show
                when={shown}
                fallback={() => (
                  <div style="display: grid; gap: var(--msp-space-2); justify-items: start">
                    <p style={`margin: 0; ${Quiet}`}>
                      No loaded model has a PAE. Model Archive's ModelCIF files carry one; AlphaFold DB ships it apart.
                    </p>
                    <Show when={fetchable}>
                      {(list) => (
                        <>
                          {list.map((s) => (
                            <FetchButton shown={s} />
                          ))}
                        </>
                      )}
                    </Show>
                  </div>
                )}
              >
                {(s) => {
                  const m = area.paeMatrix(s.ref);
                  if (!m)
                    return (
                      <div style="display: grid; gap: var(--msp-space-2); justify-items: start">
                        <p style={`margin: 0; ${Quiet}`}>{s.label} has no PAE.</p>
                        <FetchButton shown={s} />
                      </div>
                    );
                  return <PaeView ctx={ctx} structure={s.ref} matrix={m} />;
                }}
              </Show>
            </div>
          );
        }),
      },
    ],
    // In the Structure panel too: a row under a structure that has a PAE; a click opens the card at it (a plot needs
    // more room than the inspector has).
    rows: [
      {
        id: PaeRow,
        when: (node) => node.kind === 'structure' && area.hasPae(node.id),
        render: domView(({ ctx: view }) => {
          const structure = (view.subject as SceneNode).id;
          return (
            <SceneRow
              label="Predicted aligned error"
              description={of(paeDescription(area.paeSize(structure)))}
              onOpen={() => {
                chosen.set(structure);
                void ctx.api.view.open({ card: PaeCard });
              }}
            />
          );
        }),
      },
    ],
    toolbar: [
      {
        id: 'model-archive.toolbar.pae',
        icon: 'gauge',
        title: 'Predicted aligned error',
        group: 'scene',
        overlay: PaeCard,
        visible: combineLatest([withPae, fetchable]).pipe(
          map(([a, b]) => a.length + b.length > 0),
          distinctUntilChanged(),
        ),
      },
    ],
  };
}

/**
 * The plot on a structure, acting through the API as the person: the cell under the pointer lit in the 3D view and the
 * sequence (`highlight.set`), a click selecting its two residues, a box marking its aligned range (highlight tone) and
 * scored range (accent tone) in both (`highlight.mark`, marks.md) until cleared.
 */
export function PaeView(props: { ctx: ViewerUIContext; structure: string; matrix: PaeMatrix }) {
  const { ctx, structure, matrix: m } = props;
  const api = ctx.api;
  return (
    <PaePlot
      matrix={m}
      colors={ctx.sources.marks.colors}
      onHover={(cell) =>
        void (cell
          ? api.highlight.set({ target: paePairTarget(m, cell.i, cell.j), structure, sequence: true })
          : api.highlight.clear())
      }
      onPick={(i, j) => void api.select.set({ target: paePairTarget(m, i, j), structure })}
      onBox={(b) =>
        void (b
          ? api.highlight.mark({
              owner: PaeMarksOwner,
              structure,
              // The rest greyed, as AlphaFold DB does: the ranges stand out of any coloring.
              dimRest: true,
              // The scored last: painted over the aligned where they overlap, as the plot draws them.
              marks: [
                {
                  id: 'aligned',
                  tone: 'highlight',
                  label: `Aligned ${paeRangeLabel(m, b.i0, b.i1)}`,
                  target: paeRangeTarget(m, b.i0, b.i1),
                },
                {
                  id: 'scored',
                  tone: 'accent',
                  label: `Scored ${paeRangeLabel(m, b.j0, b.j1)}`,
                  target: paeRangeTarget(m, b.j0, b.j1),
                },
              ],
            })
          : api.highlight.unmark({ owner: PaeMarksOwner }))
      }
    />
  );
}

/** The PAE plot's marks' owner (marks.md). */
export const PaeMarksOwner = 'model-archive.pae';
