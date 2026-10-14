/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec, type PatchEvent, type SeqViewerContext } from '@molstar/seq-viewer';
import { InteractionsViewerRegistry } from '@molstar/seq-viewer/interactions/viewer';
import { SeqViewer, useSeqState, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { readPromise } from './support/seqstar-read.js';
import { Suspense, useEffect, useMemo, useState } from 'react';
import { fixtureLoaders, MolstarFixture } from '../../../core/_test/format/load-fixtures.js';
import { compiledCase } from './support/seqstar-cases.js';
import { antibodyPlain, antibodyWithSnapshots } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import p from './support/seqstar-patches.module.css';
import { PrismCanvas, PrismPage } from './support/seqstar-prism.js';
import { STORY_CHROME, type UiChoice, ViewerCard } from './support/seqstar-viewer.js';
import { type Args, viewerControls } from './support/seqstar-viewer-args.js';
import { ScrollArea, TextAction } from './ui/index.js';

/**
 * Patches and snapshots: every change a person makes in the viewer is a small ref-addressed patch.
 */
const meta = {
  title: 'Seq*/Interaction/Patches and snapshots',
  ...viewerControls,
  // Said here too: Storybook reads a story file's layout from its own meta, not through a spread.
  parameters: { ...viewerControls.parameters, layout: 'fullscreen' },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

/** The last proposal: its label, how it was routed, its ops, and the time to the next frame. */
function PatchLog() {
  const { context: ctx } = useSeqViewer();
  const [last, setLast] = useState<{ e: PatchEvent; ms: number | null } | null>(null);
  const canUndo = useSeqState(ctx.history.canUndo);
  useEffect(() => {
    const sub = ctx.events.patch.subscribe((e) => {
      const t0 = performance.now();
      setLast({ e, ms: null });
      requestAnimationFrame(() => setLast((l) => (l?.e === e ? { e, ms: performance.now() - t0 } : l)));
    });
    return () => sub.unsubscribe();
  }, [ctx]);
  return (
    <div className={p.change}>
      <div className={f.row}>
        <span className={p.title}>Last change</span>
        <button
          type="button"
          className={p.chip}
          disabled={!canUndo}
          onClick={() => ctx.actions.run('history.undo', null)}
        >
          Undo
        </button>
      </div>
      {last ? (
        <>
          <span>
            {last.e.label} · {last.e.routed} · {last.e.patch.ops.length} ops · frame{' '}
            {last.ms === null ? '…' : `${last.ms.toFixed(1)} ms`}
          </span>
          <ScrollArea orientation="both">
            <pre className={p.ops}>{last.e.patch.ops.map((o) => JSON.stringify(o)).join('\n')}</pre>
          </ScrollArea>
        </>
      ) : (
        <span className={p.quiet}>Hover a track header for its menu, or use a switch or a snapshot.</span>
      )}
    </div>
  );
}

/** What the patches are made with, left out by `clean`: the views, the headers (their switches and hidden chips) and the track menu. */
const PATCH_PARTS = ['ViewTabs', 'SectionHeader', 'TrackActions'] as const;
function PatchesCase({ ui }: { ui: UiChoice }) {
  const c = readPromise(compiledCase('case2'));
  const doc = useMemo(() => antibodyWithSnapshots(c.doc), [c]);
  return (
    <ViewerCard doc={doc} label="Fv · trastuzumab" options={{ ui, parts: PATCH_PARTS }}>
      <PatchLog />
    </ViewerCard>
  );
}
/** Filter view, hide and remove, switches, snapshots and undo on the real viewer. */
export const ViewerPatches: Story = {
  name: 'Filter view, hide and remove, snapshots, undo',
  render: (a) => (
    <PrismPage
      title="Filter view, hide and remove, patch snapshots"
      caption={
        <>
          Case 2 with three snapshots added. Each action is a small ref-addressed patch proposed to the document port
          (here the context itself, with undo); the graph reruns only what the patch touched.
        </>
      }
    >
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <PatchesCase ui={a.ui} />
      </Suspense>
    </PrismPage>
  ),
};

// ---- The same patches through `SeqViewer`'s `onReady`, on case 2 ----

/** Undo, from the context the viewer made (a host's own control beside the UI). */
function Undo({ context: ctx }: { context: SeqViewerContext }) {
  const canUndo = useSeqState(ctx.history.canUndo);
  return (
    <TextAction tone="quiet" disabled={!canUndo} onClick={() => ctx.actions.run('history.undo', null)}>
      Undo
    </TextAction>
  );
}
function SeqViewerPatchesView({ chrome }: { chrome: 'full' | 'clean' }) {
  const doc = useMemo(antibodyPlain, []);
  const [ctx, setCtx] = useState<SeqViewerContext | null>(null);
  return (
    <PrismPage
      title="Filter view, hide and remove, patch snapshots"
      caption={
        <>
          Case 2 with three snapshots added. Each action is a small ref-addressed patch proposed to the document port
          (here the context itself, with undo).
        </>
      }
    >
      <section className={f.card}>
        <div className={f.row}>
          <h2 className={f.heading}>Fv · trastuzumab</h2>
          {ctx && <Undo context={ctx} />}
        </div>
        <PrismCanvas>
          <SeqViewer spec={spec} ports={ports} document={doc} ui={{ chrome }} onReady={setCtx} />
        </PrismCanvas>
      </section>
    </PrismPage>
  );
}
export const SeqViewerPatches: StoryObj<{ chrome: 'full' | 'clean' }> = {
  name: 'The same through SeqViewer',
  args: { chrome: STORY_CHROME },
  argTypes: {
    chrome: { control: 'inline-radio', options: ['full', 'clean'], description: "This UI's parts (a UI option)" },
  },
  render: (a) => <SeqViewerPatchesView chrome={a.chrome} />,
};

// The fixtures' structures and the interactions renderers; read at mount.
const base = DefaultSeqViewerSpec();
const spec = { ...base, registry: [...(base.registry ?? []), { svs: MolstarFixture }, InteractionsViewerRegistry] };
const ports = { loaders: fixtureLoaders() };
