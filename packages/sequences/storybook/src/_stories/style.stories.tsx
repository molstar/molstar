/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useSeqState, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { readPromise } from './support/seqstar-read.js';
import { type CSSProperties, type ReactNode, Suspense } from 'react';
import { compiledCase } from './support/seqstar-cases.js';
import { BAD, FLOORS, HEADS, STYLE_FAMILY, styled } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { IssueList } from './support/seqstar-issues.js';
import { PrismPage } from './support/seqstar-prism.js';
import st from './support/seqstar-style.module.css';
import { ViewerCard } from './support/seqstar-viewer.js';

/**
 * A section's `style` as a preset name or an inline preset: `preset` names the base and a curated
 * list of keys overrides it for that section, mapped onto the viewer's render props. The
 * document's key beats the preset's value; a bad key is a warning in the issues, the rest applies.
 */
const meta = { title: 'Seq*/Views and layouts/Section style', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

const json = (v: unknown) => JSON.stringify(v).replaceAll(',', ', ').replaceAll(':', ': ');
const Label = ({ style }: { style: unknown }) => <span className={f.code}>"style": {json(style)}</span>;

function Frame({ title, kicker, children }: { title: string; kicker: string; children: ReactNode }) {
  return (
    <PrismPage title={title} caption={kicker}>
      <Suspense fallback={<div className={f.status}>Loading…</div>}>{children}</Suspense>
    </PrismPage>
  );
}

function HeadsCards() {
  const c = readPromise(compiledCase('case3'));
  return (
    <div className={f.stack}>
      {HEADS.map((style) => (
        <ViewerCard
          key={style.column_heads}
          doc={styled(c.doc, 'family.main', style)}
          label={<Label style={style} />}
        />
      ))}
    </div>
  );
}
/** Case 3 with the data column's head (Species) on the last pinned track's line, and on a row of its own. */
export const ColumnHeads: Story = {
  name: 'Column heads',
  render: () => (
    <Frame
      title="Case 3: the data columns' heads, shared or on their own row"
      kicker="column_heads: shared (the default) puts Species on the last pinned track's line (FLVR arginine); own_row gives the heads a row of their own between the pinned tracks and the rows, which start one row lower."
    >
      <HeadsCards />
    </Frame>
  ),
};

/** The family's four data columns at 640 px: each floor keeps a data area of at least that many px. */
export const DataFloor: Story = {
  name: 'Data floor at a narrow width',
  render: () => (
    <Frame
      title="The data area's floor: data columns yield to a narrow width"
      kicker="A synthetic family with four data columns (Len, Id., Species, Clade) in a 640 px card. data_floor 0: no column yields, and at this width they leave the letters nothing, so the view waits for a width (the overview still maps the rows); 64 (msa's default, left out): columns yield from the last while the data area would be narrower; 160: more yield. Widen the window to see them come back."
    >
      <div className={st.narrow} style={{ '--style-card-width': '640px' } as CSSProperties}>
        {FLOORS.map((style) => (
          <ViewerCard
            key={String(style.data_floor)}
            doc={styled(STYLE_FAMILY, 'fam.main', style)}
            label={<Label style={style} />}
          />
        ))}
      </div>
    </Frame>
  ),
};

/** The context's issues as the default UI's badge lists them (`issues.current`), drawn open. */
function ShownIssues() {
  const { context } = useSeqViewer();
  const issues = useSeqState(context.issues.current);
  return (
    <div className={f.stack}>
      <h3 className={f.subheading}>issues.current</h3>
      <IssueList issues={issues} />
    </div>
  );
}
function InvalidCard() {
  const c = readPromise(compiledCase('case3'));
  return (
    <ViewerCard doc={styled(c.doc, 'family.main', BAD)} label={<Label style={BAD} />} options={{ parts: ['Issues'] }}>
      <ShownIssues />
    </ViewerCard>
  );
}
/** A misspelled key and a bad value: two warnings naming what is allowed; the section takes `row_height` (26 px). */
export const InvalidKey: Story = {
  name: 'An invalid key',
  render: () => (
    <Frame
      title="Case 3 with a misspelled style key and a bad value"
      kicker="colum_heads is no key and data_floor -10 is below 0: each is a warning (svs.section.style) naming the keys or the values allowed, listed in the issues (the badge, and drawn open below). The document still loads and the section is drawn with the rest: msa with 26 px rows."
    >
      <InvalidCard />
    </Frame>
  ),
};
