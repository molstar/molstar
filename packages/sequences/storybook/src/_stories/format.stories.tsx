/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ResolvedView, resolveView, supportedAll, toJSON } from "@molstar/seq-core/format";
import type { SeqSet } from "@molstar/seq-core/model";
import { NightingaleSupported } from "@molstar/seq-nightingale";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { Suspense, use, useMemo } from "react";
import { CASE_TITLES, type CompiledCase, compiledCase, reg } from "./support/seqstar-cases.js";
import s from "./support/seqstar-format.module.css";
import { IssueList } from "./support/seqstar-issues.js";
import { PrismPage } from "./support/seqstar-prism.js";
import { makeBundle, openFixture } from "./support/seqstar-svsx.js";
import { ScrollArea } from "./ui/index.js";

/**
 * The format: each acceptance case as its document, its
 * validation and compile issues on the test fixtures, and each view resolved, summarized: sections, blocks, and
 * layers with the type declared → the type drawn by a viewer that draws everything, and by Nightingale.
 */
const meta = { title: "Seq*/Cases/Format", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

const AXIS = ["blocks", "columns"];
const ARRANGEMENT = ["side", "stacked"];

function Summary({ view, set, nightingale }: { view: ResolvedView; set: SeqSet; nightingale: ResolvedView }) {
  const str = (id: number) => (id >= 0 ? set.strings.get(id) : "–");
  const S = view.sections;
  const L = view.layers;
  const drawn = (v: ResolvedView, l: number) => (v.layers.type[l] >= 0 ? str(v.layers.type[l]) : "omitted");
  return (
    <>
      {Array.from({ length: S.count }, (_, k) => {
        const blocks = Array.from({ length: S.blockCount[k] }, (_, i) => S.blockStart[k] + i);
        const layers = Array.from({ length: L.count }, (_, l) => l).filter(
          (l) => view.tracks.section[L.track[l]] === k,
        );
        return (
          <div key={S.ref[k]} className={s.card}>
            <h3 className={s.subheading}>
              Section <span className={s.mono}>{str(S.ref[k])}</span>: {AXIS[S.axis[k]]} axis,{" "}
              {ARRANGEMENT[S.arrangement[k]]}
              {S.hidden[k] ? ", hidden" : ""}
              {S.scrollGroup[k] >= 0 ? `, scroll group ${str(S.scrollGroup[k])}` : ""}
            </h3>
            <table className={s.table}>
              <thead>
                <tr>
                  <th>Block</th>
                  <th>Length</th>
                  <th>Shown</th>
                  <th>Numbering</th>
                </tr>
              </thead>
              <tbody>
                {blocks.map((b) => (
                  <tr key={b}>
                    <td>{str(view.blocks.label[b])}</td>
                    <td>{view.blocks.length[b]}</td>
                    <td>{view.blocks.positions[b] ? "a subset" : "all"}</td>
                    <td className={s.mono}>{str(view.blocks.numbering[b]?.name ?? -1)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
            <table className={s.table}>
              <thead>
                <tr>
                  <th>Track</th>
                  <th>Layer</th>
                  <th>Declared → drawn</th>
                  <th>Nightingale draws</th>
                  <th>Placements</th>
                </tr>
              </thead>
              <tbody>
                {layers.map((l) => (
                  <tr key={l}>
                    <td>{str(view.tracks.label[L.track[l]])}</td>
                    <td className={s.mono}>{str(L.ref[l])}</td>
                    <td className={s.mono}>
                      {str(L.declared[l])} → {drawn(view, l)}
                    </td>
                    <td className={s.mono}>{drawn(nightingale, l)}</td>
                    <td>{L.placementCount[l]}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        );
      })}
    </>
  );
}

function ViewSummary({ c, view }: { c: CompiledCase; view: { ref: string; label: string } }) {
  const r = useMemo(() => {
    const all = resolveView(c.doc, c.set, view.ref, reg, supportedAll(reg));
    const nightingale = resolveView(c.doc, c.set, view.ref, reg, NightingaleSupported).view;
    return { ...all, nightingale };
  }, [c, view.ref]);
  return (
    <section className={s.card}>
      <h2 className={s.heading}>
        View {view.label} <span className={s.code}>{view.ref}</span>
      </h2>
      <IssueList issues={r.issues} />
      <Summary view={r.view} set={c.set} nightingale={r.nightingale} />
    </section>
  );
}

function Document({ name }: { name: string }) {
  const c = use(compiledCase(name));
  const json = useMemo(() => toJSON(c.doc, { pretty: true }), [c]);
  return (
    <>
      <section className={s.card}>
        <h2 className={s.heading}>The document</h2>
        <ScrollArea className={s.box} orientation="both">
          <pre className={s.pre}>{json}</pre>
        </ScrollArea>
        <h3 className={s.subheading}>Validation</h3>
        <IssueList issues={c.validation} />
        <h3 className={s.subheading}>Compile, on the fixtures</h3>
        <IssueList issues={c.compile} />
      </section>
      {c.views.map((v) => (
        <ViewSummary key={v.ref} c={c} view={v} />
      ))}
    </>
  );
}

function Case({ name }: { name: string }) {
  return (
    <PrismPage title={CASE_TITLES[name]} caption="SeqViewSpec 1.0: the document, its issues, and each view resolved.">
      <Suspense fallback={<div className={s.status}>Compiling…</div>}>
        <Document name={name} />
      </Suspense>
    </PrismPage>
  );
}

export const Case1: Story = { name: "p53: text and tracks", render: () => <Case name="case1" /> };
export const Case2: Story = { name: "Antibody Fv, Kabat and IMGT", render: () => <Case name="case2" /> };
export const Case3: Story = { name: "SH2 family", render: () => <Case name="case3" /> };
export const Case4: Story = { name: "tRNA family, base pairs", render: () => <Case name="case4" /> };
export const Case5: Story = { name: "DNA construct", render: () => <Case name="case5" /> };
export const Case6: Story = { name: "Imatinib in ABL", render: () => <Case name="case6" /> };
export const Case7: Story = { name: "Barnase–barstar", render: () => <Case name="case7" /> };
export const Case8: Story = { name: "AlphaFold DB complex", render: () => <Case name="case8" /> };
export const Case9: Story = { name: "PF00042 by view", render: () => <Case name="case9" /> };
export const Case10: Story = { name: "The Stories drafts", render: () => <Case name="case10" /> };
export const Case12: Story = { name: "Designed dimer", render: () => <Case name="case12" /> };
export const Case13: Story = { name: "p53 features, UniProt-style", render: () => <Case name="case13" /> };
export const Case14: Story = { name: "Protein over its coding DNA", render: () => <Case name="case14" /> };
export const Case15: Story = { name: "Indels in protein designs", render: () => <Case name="case15" /> };
export const Case16: Story = { name: "A plasmid map", render: () => <Case name="case16" /> };
export const Case17: Story = { name: "A genotype matrix", render: () => <Case name="case17" /> };
export const Case18: Story = { name: "Finding things in p53", render: () => <Case name="case18" /> };
export const Case19: Story = { name: "A chromosome by region", render: () => <Case name="case19" /> };
export const Case20: Story = { name: "A protein on its genome", render: () => <Case name="case20" /> };
export const Case21: Story = { name: "A gene on its reverse strand", render: () => <Case name="case21" /> };
export const Case22: Story = { name: "A family triaged", render: () => <Case name="case22" /> };

const opened = openFixture();
const bundle = makeBundle();

function FileTable({ files }: { files: ReadonlyMap<string, Uint8Array> }) {
  return (
    <table className={s.table}>
      <thead>
        <tr>
          <th>File in the archive</th>
          <th>Bytes</th>
        </tr>
      </thead>
      <tbody>
        {[...files].map(([path, bytes]) => (
          <tr key={path}>
            <td className={s.mono}>{path}</td>
            <td>{bytes.length}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}

function Json({ text }: { text: string }) {
  return (
    <ScrollArea className={s.box} orientation="both">
      <pre className={s.pre}>{text}</pre>
    </ScrollArea>
  );
}

function OpenedArchive() {
  const o = use(opened);
  const json = useMemo(() => (o.svsx.doc ? toJSON(o.svsx.doc, { pretty: true }) : ""), [o]);
  return (
    <>
      <section className={s.card}>
        <h2 className={s.heading}>
          The archive <span className={s.code}>kinase.svsx</span>, {o.size} bytes
        </h2>
        <FileTable files={o.svsx.files} />
        <p className={s.status}>
          Main file <span className={s.mono}>{o.svsx.main}</span> (MessagePack, found before{" "}
          <span className={s.mono}>index.svsj</span>); relative URIs resolve against{" "}
          <span className={s.mono}>{o.svsx.base || "the archive's root"}</span>.
        </p>
        <h3 className={s.subheading}>The document, as the archive carries it</h3>
        <Json text={json} />
        <h3 className={s.subheading}>Validation, with the archive as the base</h3>
        <IssueList issues={o.svsx.issues} />
        <h3 className={s.subheading}>Compile over the archive's files: nothing else is asked</h3>
        <IssueList issues={o.compiled?.compile ?? []} />
        {o.compiled ? (
          <div className={s.status}>
            {o.compiled.set.molecules.length} molecule(s), {o.compiled.set.annotations.length} annotation(s) loaded
          </div>
        ) : null}
      </section>
      {o.compiled?.views.map((v) => (
        <ViewSummary key={v.ref} c={o.compiled as CompiledCase} view={v} />
      ))}
    </>
  );
}

function Bundled() {
  const b = use(bundle);
  return (
    <>
      <section className={s.card}>
        <h2 className={s.heading}>A document with absolute URIs</h2>
        <Json text={b.before} />
        <h3 className={s.subheading}>Fetched once each, through a fixture fetch (no network)</h3>
        <ul className={s.issues}>
          {b.requests.map((r) => (
            <li key={r}>
              <span className={s.mono}>{r}</span>, {b.sources[r]?.length ?? 0} bytes
            </li>
          ))}
        </ul>
      </section>
      <section className={s.card}>
        <h2 className={s.heading}>The bundle, {b.size} bytes: each URI rewritten to its path in the archive</h2>
        <FileTable files={b.read.files} />
        <Json text={b.after} />
        <h3 className={s.subheading}>Read back and validated</h3>
        <IssueList issues={b.read.issues} />
      </section>
    </>
  );
}

/** `.svsx`: one zip, a document and the files it reads, opened with no network. */
export const Svsx: Story = {
  name: "A document with its data (.svsx)",
  render: () => (
    <PrismPage
      title="A document with its data (.svsx)"
      caption={
        <>
          A zip with the document as index.svsm and the files it reads, opened from a committed fixture; then the bundle
          step that makes one from a document with absolute URIs.
        </>
      }
    >
      <Suspense fallback={<div className={s.status}>Opening…</div>}>
        <OpenedArchive />
        <Bundled />
      </Suspense>
    </PrismPage>
  ),
};
