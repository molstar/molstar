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
import { derive, toState } from '@molstar/prism-ui/reactive/state';
import { type SceneNode, structureOfGroup } from '@molstar/prism-viewer/model/scene';
import { SceneRow } from '@molstar/prism-viewer/views/structure/scene-row';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { type PaeMatrix, paePairTarget, paeRangeLabel, paeRangeTarget } from './data.js';
import type { PaeSize } from './operations.js';
import { type PaeBox, PaePlot } from './plot.js';

// The PAE card (parity.md §6, "PAE plot"; v5's Predicted Aligned Error panel, anew): the predicted aligned error of a
// model as a heat map (plot.tsx), as AlphaFold DB draws it. The cell under the pointer is lit in the 3D view and the
// sequence, a click selects its two residues, a box marks its aligned and scored ranges in both until cleared (marks.md);
// all through the API as the person. The matrix comes from the `model-archive` area; the card opens from the toolbar
// while a loaded model has a PAE, or by `model-archive.show-pae`. The PAE is also a row in the Structure panel under a
// structure that has one or whose PAE can be fetched (as volume streaming's row): its click opens the card at it,
// fetching the PAE first when it is not here yet.

export const PaeCard = 'model-archive.overlay.pae';
/** The PAE's row in the Structure panel, under a structure that has one or can fetch it. */
export const PaeRow = 'model-archive.row.pae';

/** "149 residues · 2 chains": what the row says of a PAE (its size, not its matrix: the card reads that). */
const paeDescription = (s: PaeSize | undefined) =>
  s && `${s.residues} residues${s.chains > 1 ? ` · ${s.chains} chains` : ''}`;

const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';

const sameShown = (a: readonly Shown[], b: readonly Shown[]) =>
  a.length === b.length && a.every((x, i) => x.ref === b[i].ref && x.label === b[i].label);

interface Shown {
  readonly ref: string;
  readonly label: string;
}

export function paeUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api['model-archive'];
  if (!area) return {};
  // The structures that have a PAE (by the file's metric list: no matrix is read for this), and those whose PAE can be
  // fetched (Model Archive's associated file, an AlphaFold DB entry): asked again when the state changes (a fetched PAE is a state cell). The
  // area tells of its changes by a stream, not a state: a state told at each (never the same), at this edge.
  const changed = toState<void>(area.changes, undefined, () => false);
  const shownOf = (has: (ref: string) => boolean) =>
    derive(
      [ctx.sources.structures, changed],
      (list) => list.filter((s) => has(s.ref)).map((s): Shown => ({ ref: s.ref, label: s.label })),
      sameShown,
    );
  const withPae = shownOf((ref) => area.hasPae(ref));
  const fetchable = shownOf((ref) => area.canFetchPae(ref));
  const fetchPae = (structure: string) => void ctx.api['model-archive']?.fetchPae({ structure });
  const FetchButton = (props: { shown: Shown }) => (
    <Button size="sm" onClick={() => fetchPae(props.shown.ref)}>
      Fetch PAE{props.shown.label ? ` (${props.shown.label})` : ''}
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
          const shown = derive(
            [withPae, chosen.asReadonly],
            (list, ref) => list.find((s) => s.ref === ref) ?? list[0],
            (a, b) => a?.ref === b?.ref,
          );
          // The card's text column (its title's: the card's inset and a row's own), on both sides: the plot is no row.
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0; padding: 0 calc(var(--msp-panel-inset) + var(--msp-inset)) calc(var(--msp-panel-inset) + var(--msp-inset))">
              <Show when={derive(withPae, (l) => (l.length > 1 ? l : undefined))}>
                {(list) => (
                  <Select
                    size="sm"
                    aria-label="Structure"
                    options={list.map((s) => [s.ref, s.label] as const)}
                    value={derive(shown, (s) => s?.ref ?? '')}
                    onChange={(ref) => chosen.set(ref)}
                  />
                )}
              </Show>
              <Show
                when={shown}
                fallback={() => (
                  <div style="display: grid; gap: var(--msp-space-2); justify-items: start">
                    <p style={`margin: 0; ${Quiet}`}>
                      No loaded model has a PAE. A ModelCIF file may carry one; Model Archive and AlphaFold DB often
                      ship it apart.
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
    // In the Structure panel too: a row in the Other group of a structure that has a PAE or can fetch one; a click opens
    // the card at it (a plot needs more room than the inspector has), fetching the PAE first when it is not here yet.
    rows: [
      {
        id: PaeRow,
        when: (node) => {
          const structure = node.kind === 'other' ? structureOfGroup(node) : undefined;
          return !!structure && (area.hasPae(structure) || area.canFetchPae(structure));
        },
        // A fetched PAE is a state cell, not a change of the scene: the tree asks `when` again on the area's changes.
        whenChanged: changed,
        render: domView(({ ctx: view }) => {
          const structure = structureOfGroup(view.subject as SceneNode)!;
          const fetching = observableState(false);
          const description = derive([changed, fetching.asReadonly], (_, busy) =>
            area.hasPae(structure)
              ? paeDescription(area.paeSize(structure))
              : busy
                ? 'Fetching…'
                : `Not loaded · ${area.paeService(structure) ?? 'elsewhere'}`,
          );
          const open = () => {
            chosen.set(structure);
            void ctx.api.view.open({ card: PaeCard });
          };
          return (
            <SceneRow
              label="PAE"
              description={description}
              onOpen={async () => {
                if (area.hasPae(structure)) return open();
                if (fetching.value) return;
                fetching.set(true);
                const r = await ctx.api['model-archive']?.fetchPae({ structure });
                fetching.set(false);
                if (r?.ok) open();
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
        visible: derive([withPae, fetchable], (a, b) => a.length + b.length > 0),
      },
    ],
  };
}

/**
 * The plot on a structure, acting through the API as the person: the cell under the pointer lit in the 3D view and the
 * sequence (`highlight.set`), a click selecting its two residues (again: unselecting them), a box marking its aligned range (highlight tone) and
 * scored range (accent tone) in both (`highlight.mark`, marks.md) until cleared.
 */
export function PaeView(props: { ctx: ViewerUIContext; structure: string; matrix: PaeMatrix }) {
  const { ctx, structure, matrix: m } = props;
  const api = ctx.api;
  const area = api['model-archive'];
  // A click selects the cell's two residues; on the pair already selected (and nothing else), it unselects them.
  // A pair or a box, one at a time: a box takes off the pair the plot selected last (the sequence is not held on a pair
  // the person moved on from), and while it is shown a click selects nothing (the plot). Only while the selection is
  // still the one the plot made (`selectionToken`): a selection made since is the person's, and left alone.
  let picked: { target: string; token: object | undefined } | undefined;
  // The calls one after another, in the order the plot asked for them: a box, a ✕ or the card closing while a click's or
  // another box's call is on its way never lands before it (the marks and the selection end as the plot shows them).
  let queue: Promise<unknown> = Promise.resolve();
  const inTurn = (fn: () => Promise<unknown>) => {
    queue = queue.then(fn, fn);
    return queue;
  };
  const togglePair = (target: string) =>
    inTurn(async () => {
      const r = await api.select.set({ target, structure });
      picked = { target, token: area?.selectionToken(structure) };
      if (r.ok && r.unchanged) {
        await api.select.set({ target, structure, mode: 'remove' });
        picked = undefined;
      }
    });
  const unpick = async () => {
    const was = picked;
    picked = undefined;
    if (!was || !area || area.selectionToken(structure) !== was.token) return;
    await api.select.set({ target: was.target, structure, mode: 'remove' });
  };
  const mark = (b: PaeBox) =>
    inTurn(async () => {
      await unpick();
      await api.highlight.mark({
        owner: PaeMarksOwner,
        structure,
        // The rest greyed, as AlphaFold DB does: the ranges stand out of any coloring.
        dimRest: true,
        // Both ranges in view in the sequence (it zooms out to them), as a selection's are.
        focusSequence: true,
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
      });
    });
  const unmark = () => inTurn(() => api.highlight.unmark({ owner: PaeMarksOwner }));

  return (
    <PaePlot
      matrix={m}
      colors={ctx.sources.marks.colors}
      onHover={(cell) =>
        void (cell
          ? api.highlight.set({ target: paePairTarget(m, cell.i, cell.j), structure, sequence: true })
          : api.highlight.clear())
      }
      onPick={(i, j) => void togglePair(paePairTarget(m, i, j))}
      onBox={(b) => void (b ? mark(b) : unmark())}
    />
  );
}

/** The PAE plot's marks' owner (marks.md). */
export const PaeMarksOwner = 'model-archive.pae';
