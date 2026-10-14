/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { resolveView } from "@molstar/seq-core/format";
import type { Issue } from "@molstar/seq-core/model";
import { materialize, NightingaleSupported } from "@molstar/seq-nightingale";
import { mount } from "@molstar/seq-nightingale/mount";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { Suspense, use, useEffect, useMemo, useRef, useState } from "react";
import { type CompiledCase, compiledCase, reg } from "./support/seqstar-cases";
import f from "./support/seqstar-format.module.css";
import { IssueList } from "./support/seqstar-issues";
import { PrismPage } from "./support/seqstar-prism";
import { BuilderCode, type NavigatorChoice, type UiChoice, ViewerCard } from "./support/seqstar-viewer";
import s from "./support/seqstar-viewer.module.css";
import { type Args, viewerControls } from "./support/seqstar-viewer-args";
import { ScrollArea } from "./ui";

/**
 * The same case and view in our viewer and in Nightingale.
 */
const meta = {
  title: "Seq*/Cases/Default UIs",
  ...viewerControls,
  // Said here too: Storybook reads a story file's layout from its own meta, not through a spread.
  parameters: { ...viewerControls.parameters, layout: "fullscreen" },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

/** Nightingale on one view, as the Nightingale stories mount it (no event wiring here). */
function NightingaleView({ c, view }: { c: CompiledCase; view: string }) {
  const r = useMemo(() => {
    const { view: resolved, issues } = resolveView(c.doc, c.set, view, reg, NightingaleSupported);
    const m = materialize(resolved, c.set);
    const lost = [...issues.filter((i) => i.code.startsWith("svs.render.")), ...m.issues];
    return { view: resolved, root: m.root, lost };
  }, [c, view]);
  const host = useRef<HTMLDivElement>(null);
  const [failed, setFailed] = useState<Issue[]>([]);
  useEffect(() => {
    const el = host.current;
    if (!el) return;
    const m = mount(r.root, el, { view: r.view, set: c.set });
    let live = true;
    m.ready.then((issues) => live && setFailed(issues));
    return () => {
      live = false;
      m.dispose();
    };
  }, [r, c.set]);
  return (
    <section className={f.card}>
      <h2 className={f.heading}>
        Nightingale <span className={f.code}>{view}</span>
      </h2>
      <ScrollArea orientation="horizontal">
        <div ref={host} className={f.nightingale} />
      </ScrollArea>
      <h3 className={f.subheading}>What Nightingale cannot draw</h3>
      <IssueList issues={[...r.lost, ...failed]} />
    </section>
  );
}

function SideBySide({ name, navigator, ui }: { name: string; navigator: NavigatorChoice; ui: UiChoice }) {
  const c = use(compiledCase(name));
  const v = c.views[0];
  return (
    <div className={s.pair}>
      <ViewerCard doc={c.doc} label={`Our viewer: ${v.label}`} options={{ navigator, view: v.ref, ui }} />
      <NightingaleView c={c} view={v.ref} />
    </div>
  );
}

/** The same case and view in our viewer and in Nightingale (the format's second viewer), side by side. */
export const SideBySideNightingale: Story = {
  name: "Side by side with Nightingale",
  render: (a) => (
    <PrismPage
      title="Our viewer and Nightingale on case 7, Barnase–barstar"
      caption="The first view of the case, resolved for each viewer's supported set."
    >
      <BuilderCode name="case7" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <SideBySide name="case7" navigator={a.navigator} ui={a.ui} />
      </Suspense>
    </PrismPage>
  ),
};
