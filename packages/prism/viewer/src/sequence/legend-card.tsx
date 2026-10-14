/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { For, Legend } from '@molstar/prism-ui';
import { map } from 'rxjs';
import type { SequenceViewer } from './model.js';

/** The legend of the sequence's colors, in Prism's `Legend`: the card the Sequence menu's "Legend…" opens. */
export function LegendCard(props: { viewer: SequenceViewer }) {
  const { legend } = props.viewer;
  return (
    <div style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-3)">
      <Legend
        title={legend.pipe(map((l) => l.title))}
        items={legend.pipe(map((l) => l.items))}
        ramps={legend.pipe(map((l) => l.ramps))}
      />
      <For each={legend.pipe(map((l) => l.sections ?? []))}>
        {(section) => <Legend title={section.title} items={section.items} />}
      </For>
    </div>
  );
}
