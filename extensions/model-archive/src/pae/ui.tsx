/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, Legend, onCleanup, Select, Show } from '@molstar/prism-ui';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import { type PaeMatrix, PaeColorStyles, paeCellAt, paeCellLabel, paeChains, paeImage, paePairTarget } from './data.js';

// The PAE card (parity.md §6, "PAE plot"; v5's Predicted Aligned Error panel, anew): the predicted aligned error of a
// model as a heat map, one pixel per residue pair, drawn once into a canvas. The pointer names the cell under it (read
// and written once per frame), a click selects its two residues (`select.set`, as the person). The matrix comes from
// the `model-archive` area; the card opens from the toolbar while a loaded model has a PAE, or by
// `model-archive.show-pae`.

export const PaeCard = 'model-archive.overlay.pae';

const Quiet = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';

interface Shown {
  readonly ref: string;
  readonly label: string;
}

export function paeUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api['model-archive'];
  if (!area) return {};
  // The structures that have a PAE (by the file's metric list: no matrix is read for this).
  const withPae = ctx.sources.structures.pipe(
    map((list) => list.filter((s) => area.hasPae(s.ref)).map((s): Shown => ({ ref: s.ref, label: s.label }))),
  );
  return {
    overlays: [
      {
        id: PaeCard,
        title: 'Predicted aligned error',
        render: domView(() => {
          const chosen = observableState<string | undefined>(undefined);
          const shown = combineLatest([withPae, chosen.asReadonly]).pipe(
            map(([list, ref]) => list.find((s) => s.ref === ref) ?? list[0]),
            distinctUntilChanged((a, b) => a?.ref === b?.ref),
          );
          return (
            <div style="display: grid; gap: var(--msp-space-3); min-width: 0; padding: 0 var(--msp-panel-inset) var(--msp-panel-inset)">
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
                  <p style={`margin: 0; ${Quiet}`}>
                    No loaded model has a PAE. Model Archive's ModelCIF files carry one; AlphaFold DB ships it apart.
                  </p>
                )}
              >
                {(s) => {
                  const m = area.paeMatrix(s.ref);
                  if (!m) return <p style={`margin: 0; ${Quiet}`}>{s.label} has no PAE.</p>;
                  return (
                    <PaePlot
                      matrix={m}
                      onPick={(i, j) => void ctx.api.select.set({ target: paePairTarget(m, i, j), structure: s.ref })}
                    />
                  );
                }}
              </Show>
            </div>
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
        visible: withPae.pipe(
          map((l) => l.length > 0),
          distinctUntilChanged(),
        ),
      },
    ],
  };
}

/** The heat map with the chains along it, the cell under the pointer, and the scale. */
export function PaePlot(props: { matrix: PaeMatrix; onPick(i: number, j: number): void }) {
  const m = props.matrix;
  const canvas = (
    <canvas
      width={m.size}
      height={m.size}
      role="img"
      aria-label={`${m.name}: ${m.size} by ${m.size} residues`}
      style="display: block; width: 100%; aspect-ratio: 1; image-rendering: pixelated; cursor: crosshair; border-radius: var(--msp-radius-sm); background: var(--msp-color-well)"
    />
  ) as HTMLCanvasElement;
  const g = canvas.getContext('2d');
  if (g) {
    const image = g.createImageData(m.size, m.size);
    paeImage(m, image.data);
    g.putImageData(image, 0, 0);
  }

  const readout = (<p style={`margin: 0; min-height: 1.5em; ${Quiet}`} aria-live="polite" />) as HTMLElement;
  const idle = 'Point at a cell; a click selects its two residues.';
  readout.textContent = idle;

  // Pointer-rate input, handled once per frame: the frame reads the plot's box, then writes the line (styling.md §5).
  let frame = 0;
  let x = 0;
  let y = 0;
  let inside = false;
  let cell = -1;
  const update = () => {
    frame = 0;
    const r = canvas.getBoundingClientRect();
    const at = inside && r.width > 0 ? paeCellAt(m, (x - r.left) / r.width, (y - r.top) / r.height) : undefined;
    const next = at ? at.i * m.size + at.j : -1;
    if (next === cell) return;
    cell = next;
    readout.textContent = at ? paeCellLabel(m, at.i, at.j) : idle;
  };
  const schedule = (e: PointerEvent, over: boolean) => {
    x = e.clientX;
    y = e.clientY;
    inside = over;
    if (!frame) frame = requestAnimationFrame(update);
  };
  canvas.addEventListener('pointermove', (e) => schedule(e, true));
  canvas.addEventListener('pointerleave', (e) => schedule(e, false));
  canvas.addEventListener('click', (e) => {
    const r = canvas.getBoundingClientRect();
    const at = paeCellAt(m, (e.clientX - r.left) / r.width, (e.clientY - r.top) / r.height);
    if (at) props.onPick(at.i, at.j);
  });
  onCleanup(() => {
    if (frame) cancelAnimationFrame(frame);
  });

  const chains = paeChains(m);
  return (
    <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
      {canvas}
      {chains.length > 1 ? (
        <div style={`display: flex; ${Quiet}`} aria-label="Chains">
          {chains.map((c) => (
            <span style={`flex: ${c.to - c.from} 1 0; min-width: 0; overflow: hidden; text-align: center`}>
              {c.label}
            </span>
          ))}
        </div>
      ) : null}
      {readout}
      <Legend
        ramps={[
          {
            stops: [
              { t: 0, color: PaeColorStyles.low },
              { t: 1, color: PaeColorStyles.high },
            ],
            lo: '0 Å',
            hi: `${m.max.toFixed(0)} Å`,
            caption: 'Expected position error',
          },
        ]}
      />
    </div>
  );
}
