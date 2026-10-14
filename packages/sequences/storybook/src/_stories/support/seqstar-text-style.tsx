/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A case's text views in either built-in text preset, the reader's pick: `wrapped` (the default) or `compact`; and one
// of the select behavior's params overridden by the host (`shiftPositions`), as a viewer spec would.

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { useSeqViewer } from '@molstar/seq-viewer/react';
import { type ReactNode, useMemo, useState } from 'react';
import { Choices } from '../ui/index.js';
import { textStyled } from './seqstar-cases.js';

const STYLES = [
  {
    id: 'wrapped',
    label: 'Wrapped',
    title: 'The `wrapped` preset: spaced to read, the group numbers above each line',
  },
  {
    id: 'compact',
    label: 'Compact',
    title: 'The `compact` preset: half the group gap, lines and lanes touching, still groups of 10',
  },
];

/** With `on`, the document in the picked text preset and its `Choices`; else the document as it is and nothing. */
export function useTextStyle(doc: SeqViewSpecDoc, on: boolean): [SeqViewSpecDoc, ReactNode] {
  const [style, setStyle] = useState('wrapped');
  const shown = useMemo(() => (on && style !== 'wrapped' ? textStyled(doc, style) : doc), [doc, on, style]);
  if (!on) return [doc, null];
  return [shown, <Choices key="text-style" label="Text" choices={STYLES} active={style} onPick={setStyle} />];
}

const SHIFT = [
  { id: 'extend', label: 'Extends', title: 'The default: Shift-click extends a range of positions from its far end' },
  {
    id: 'toggle',
    label: 'Adds or removes',
    title: '`shiftPositions: "toggle"`: Shift-click adds or removes a position',
  },
];

/** Shift-click on positions as the host sets the select behavior's `shiftPositions` (inside the viewer's provider). */
export function ShiftClickChoice() {
  const { context } = useSeqViewer();
  const [rule, setRule] = useState('extend');
  const pick = (id: string) => {
    setRule(id);
    context.behaviors.get('select')?.setParams({ shiftPositions: id });
  };
  return <Choices label="Shift-click or Shift-drag on residues" choices={SHIFT} active={rule} onPick={pick} />;
}
