/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { button, type DomPart, type RowProps as DomRowProps, setAttr } from '@molstar/seq-viewer/dom';
import type { RowProps } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { useState } from 'react';
import { familyWithTree, syntheticFamily } from '../../../viewer/_test/families.js';
import { computedTreeDoc, realignedDoc } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard, type ViewerOptions } from './support/seqstar-viewer.js';
import s from './support/seqstar-viewer.module.css';
import { type Args, viewerControls } from './support/seqstar-viewer-args.js';
import { Tooltip } from './ui/index.js';

/**
 * Alignments: the rows' vertical navigation, a tree beside the rows, guide trees computed at load and a realigned family,
 * and a control beside each row.
 */
const meta = {
  title: 'Seq*/Views and layouts/Alignments',
  ...viewerControls,
  // Said here too: Storybook reads a story file's layout from its own meta, not through a spread.
  parameters: { ...viewerControls.parameters, layout: 'fullscreen' },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

const VERTICAL: readonly [string, ViewerOptions][] = [
  ['scrollbar (the default)', { verticalNav: 'scrollbar' }],
  ["scrollbar, verticalZoom: its thumb's ends zoom rows", { verticalNav: 'scrollbar', verticalZoom: true }],
  ['overview', { verticalNav: 'overview' }],
  ["both (the `msa` preset's default)", { verticalNav: 'both' }],
  ["both, verticalZoom: handles on the window's top and bottom", { verticalNav: 'both', verticalZoom: true }],
  ['none (Alt-wheel and the keyboard still zoom rows)', { verticalNav: 'none' }],
];
/** A synthetic SH2-like family (lib/seqstar/viewer/_test/families.ts: case 3's rows mutated; not real sequences). */
const FAMILY = syntheticFamily(160, 10);
/** The rows' vertical navigation choices (render props `verticalNav`, `verticalZoom`) on 160 rows, 10 in view. */
export const FamilyVertical: Story = {
  name: 'Vertical navigation (160 rows)',
  render: (a) => (
    <PrismPage
      title="Vertical navigation on a synthetic family of 160 rows"
      caption={
        <>
          Synthetic: case 3's SH2 rows with deterministic substitutions and gaps (not real sequences); a box section of
          10 rows. The `msa` preset's choice per card; Alt-wheel and the keyboard zoom rows in each.
        </>
      }
    >
      {VERTICAL.map(([label, o]) => (
        <ViewerCard
          key={label}
          doc={FAMILY}
          label={`verticalNav: ${label}`}
          options={{ ...o, navigator: a.navigator, ui: a.ui }}
        />
      ))}
    </PrismPage>
  ),
};

const TREE_FAMILY = familyWithTree(160, 16, { columns: true, sites: true });
/** The tree beside the rows: click a clade's branch to collapse it into a wedge, again to open. */
export const FamilyTree: Story = {
  name: 'A tree beside the rows',
  render: (a) => (
    <PrismPage
      title="A synthetic family of 160 rows with a synthetic tree"
      caption={
        <>
          Synthetic: case 3's SH2 rows mutated and a generated tree (a clade per seed, balanced; not a real phylogeny),
          synthetic data columns and two sites on rows (lanes). Rows in tree order; click a clade's branch or connector
          to collapse it, the wedge to open it; hover a node to see its rows; Alt-wheel zooms rows (below 4 px they draw
          as one band).
        </>
      }
    >
      <ViewerCard
        doc={TREE_FAMILY}
        label="A family of 160 rows with a tree"
        options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }}
      />
    </PrismPage>
  ),
};

/** Guide trees: `nj` and `upgma` node providers on a committed alignment, computed at load from its distances. */
export const ComputedTrees: Story = {
  name: 'Guide trees and realignment',
  render: (a) => (
    <PrismPage
      title="A guide tree computed from the alignment"
      caption={
        <>
          One synthetic alignment of 13 rows (evolved down a known tree with unequal rates, in a fixed shuffle, so the
          data order is not the tree's), no tree in the document: <code>tree {'{ provider: "nj" }'}</code> and{' '}
          <code>upgma</code> compute theirs from p-distances at load, and the rows follow their leaf order. Both find
          the distant outgroup O1 and the B clade; they differ on the fast-evolving A5. NJ is rooted at the midpoint of
          its longest path and keeps A5's long branch of its own (about 0.21 against 0.01 to 0.06 for the rest of its
          clade); UPGMA assumes a clock, so every leaf ends at one height and A5 is pulled to a sister of its clade.
          Click a branch to collapse a clade.
        </>
      }
    >
      {(['nj', 'upgma'] as const).map((provider) => (
        <ViewerCard
          key={provider}
          doc={computedTreeDoc(provider)}
          label={provider === 'nj' ? 'Neighbor joining' : 'UPGMA'}
          options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }}
        />
      ))}
      <div className={f.kicker}>
        Realigning: a second family, 12 rows from four clades with insertions and deletions (B lost 5 residues, C gained
        GSGS, D lost 4), its rows given as just their residues padded at the end (first card: the indels shift
        everything after them), and as <code>alignment {'{ provider: "realign", inputs: "given" }'}</code>
        (second): biomsa's progressive alignment of the rows' own residues at load, each new row bound to its input row,
        labels and order kept, an NJ tree computed over the result, rows in tree order. The columns line up again after
        each indel.
      </div>
      {[false, true].map((realigned) => (
        <ViewerCard
          key={String(realigned)}
          doc={realignedDoc(realigned)}
          label={realigned ? 'Realigned' : 'As given'}
          options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }}
        />
      ))}
    </PrismPage>
  ),
};

/** The host's control per row (the `RowControl` slot): a 3D toggle, as Stories would show a member in 3D. */
function RowToggle({ dataRow }: RowProps) {
  const [on, setOn] = useState(false);
  return (
    <Tooltip content={`Show row ${dataRow + 1} in 3D`}>
      <button
        type="button"
        className={s.rowToggle}
        aria-pressed={on}
        aria-label={`Show row ${dataRow + 1} in 3D`}
        onClick={() => setOn(!on)}
      >
        3D
      </button>
    </Tooltip>
  );
}
const ROW_UI = { components: { RowControl: RowToggle } };
/**
 * The same control as a `DomPart` for the `<seq-viewer>` variant. It renders in the element's shadow root, where the
 * page's CSS does not reach: it borrows the skin's own `seqv-choice` (pressed state included).
 */
const domRowToggle: DomPart<DomRowProps> = (host, _ctx, { dataRow }) => {
  const label = `Show row ${dataRow + 1} in 3D`;
  const b = button(
    { class: 'seqv-choice', 'aria-pressed': false, 'aria-label': label },
    () => {
      setAttr(b, 'aria-pressed', b.getAttribute('aria-pressed') !== 'true');
    },
    '3D',
  );
  host.append(b);
  return { dispose: () => b.remove() };
};
const ROW_DOM = { RowControl: domRowToggle };
/** A control beside each visible row, through the `RowControl` slot. */
export const RowControls: Story = {
  name: 'Row controls',
  render: (a) => (
    <PrismPage
      title="A control beside each row"
      caption={
        <>
          The host's component in the `RowControl` slot, in the header's last 40 px (render prop `controlWidth`); only
          the rows in view get one, and they follow the rows' scroll by a transform (no React render per scroll).
        </>
      }
    >
      <ViewerCard
        doc={TREE_FAMILY}
        label="A family of 160 rows, a control beside each"
        ui={ROW_UI}
        domSlots={ROW_DOM}
        options={{
          navigator: a.navigator,
          rowsColor: a.rowsColor,
          ui: a.ui,
          render: { controlWidth: 40, headerWidth: 230 },
        }}
      />
    </PrismPage>
  ),
};
