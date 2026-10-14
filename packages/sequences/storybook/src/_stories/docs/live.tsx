/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the Seq* docs pages compute live, from the packages themselves, so the pages never drift from the code: the
// Introduction's example run three ways, the format's validation and generated reference, SeqQL's generated reference and
// skill text.

import {
  createRegistries,
  KernelRegistry,
  MvsRegistry,
  printReference,
  register,
  type SeqViewSpecDoc,
  validate,
} from "@molstar/seq-core/format";
import type { Issue, SeqSet } from "@molstar/seq-core/model";
import { checkText, createLanguage, printQuery, referenceText, renderIssue, skillText } from "@molstar/seq-core/query";
import { SeqViewer } from "@molstar/seq-viewer/react";
import { Markdown, Source, Unstyled } from "@storybook/addon-docs/blocks";
import { type ReactNode, useEffect, useMemo, useRef, useState } from "react";
import f from "../support/seqstar-format.module.css";
import { IssueList } from "../support/seqstar-issues";
import { PrismCanvas } from "../support/seqstar-prism";
import s from "./docs.module.css";
import { exampleLink, mountNightingaleSide } from "./linked";
import { doc, mountDomExample, svgExample } from "./minimal-example";

/** Whether the page is dark (the docs container sets `data-msp-theme` from the toolbar; blue is a dark theme). */
const readDark = () => {
  const theme = document.documentElement.dataset.mspTheme ?? document.documentElement.dataset.theme;
  return theme === "dark" || theme === "blue";
};
function useDark() {
  const [dark, setDark] = useState(readDark);
  useEffect(() => {
    const o = new MutationObserver(() => setDark(readDark()));
    o.observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme", "data-msp-theme"] });
    setDark(readDark()); // the container sets it in a layout effect, after this component first rendered
    return () => o.disconnect();
  }, []);
  return dark;
}

/** Storybook's `Source`, dark on a dark page. */
export function Code({ code, language = "tsx" }: { code: string; language?: "tsx" | "json" | "bash" }) {
  return <Source code={code} language={language} dark={useDark()} />;
}

/** A card: its label, then the example on a canvas in the page's theme. */
function Panel({ label, children }: { label: string; children: ReactNode }) {
  return (
    <section className={s.panel}>
      <div className={s.label}>{label}</div>
      <PrismCanvas>{children}</PrismCanvas>
    </section>
  );
}

function DomExample() {
  const host = useRef<HTMLDivElement>(null);
  useEffect(() => (host.current ? mountDomExample(host.current) : undefined), []);
  return <div ref={host} />;
}

function SvgExample() {
  const [svg, setSvg] = useState<string | null>(null);
  useEffect(() => {
    let live = true;
    void svgExample().then((text) => live && setSvg(text));
    return () => {
      live = false;
    };
  }, []);
  return svg ? (
    // biome-ignore lint/security/noDangerouslySetInnerHtml: renderSvg's string from our own document, not user input
    <div className={s.svg} dangerouslySetInnerHTML={{ __html: svg }} />
  ) : (
    <div className={s.note}>Drawing…</div>
  );
}

/**
 * The page's first sight: the minimal example's document in the React UI (the snippet's `ReactExample`, with `onReady`
 * to give its context to the link with Nightingale below).
 */
export function ExampleViewer() {
  useEffect(() => () => exampleLink.setViewer(null), []);
  return (
    <Unstyled>
      <div className={s.stack}>
        <Panel label="@molstar/seq-viewer/react: <SeqViewer document={doc} />">
          <SeqViewer document={doc} onReady={(ctx) => exampleLink.setViewer(ctx)} />
        </Panel>
      </div>
    </Unstyled>
  );
}

/**
 * The same document in Nightingale, the format's second viewer (`@molstar/seq-nightingale`): its view resolved for what
 * Nightingale draws, materialized as Nightingale's components and mounted (`mountNightingaleExample`'s steps, on the
 * set the viewer above shows, so the two are linked: `linked.ts`); then what it reports it cannot draw.
 */
export function ExampleNightingale() {
  const host = useRef<HTMLDivElement>(null);
  const [lost, setLost] = useState<Issue[] | null>(null);
  useEffect(() => {
    const el = host.current;
    if (!el) return;
    let live = true;
    const noop = () => {};
    let dispose = noop;
    let drawn = 0;
    let set: SeqSet | null = null;
    // Nightingale sizes its parts from its width when it mounts: mounted before the docs page settles its layout, its
    // letters track stays empty. So it mounts again when the width changes (or the viewer's set does).
    const draw = () => {
      if (!live || !set || el.clientWidth === 0 || (el.clientWidth === drawn && dispose !== noop)) return;
      drawn = el.clientWidth;
      dispose();
      const m = mountNightingaleSide(el, set);
      dispose = () => {
        exampleLink.setSide(null);
        m.dispose();
        dispose = noop;
      };
      void m.ready.then((failed) => {
        if (!live || dispose === noop) return;
        // Linked once drawn: Nightingale's first events are its drawing's.
        exampleLink.setSide(m.side);
        setLost(failed);
      });
    };
    const stopSet = exampleLink.onSet((next) => {
      if (next === set) return;
      set = next;
      dispose();
      draw();
    });
    let timer: ReturnType<typeof setTimeout> | undefined;
    const resized = new ResizeObserver(() => {
      clearTimeout(timer);
      timer = setTimeout(draw, 100);
    });
    resized.observe(el);
    return () => {
      live = false;
      clearTimeout(timer);
      resized.disconnect();
      stopSet();
      dispose();
    };
  }, []);
  return (
    <Unstyled>
      <div className={s.stack}>
        <section className={s.panel}>
          <div className={s.label}>@molstar/seq-nightingale: mountNightingaleExample(host)</div>
          <div ref={host} className={f.nightingale} />
          <div className={s.note}>
            {lost === null ? (
              "Drawing…"
            ) : lost.length === 0 ? (
              "Nightingale reports nothing it cannot draw."
            ) : (
              <>
                What Nightingale reports it cannot draw:
                <IssueList issues={lost} />
              </>
            )}
          </div>
        </section>
      </div>
    </Unstyled>
  );
}

/** The rest of the minimal example run: the DOM UI and the headless SVG of the same document. */
export function MinimalExample() {
  return (
    <Unstyled>
      <div className={s.stack}>
        <Panel label="@molstar/seq-viewer/dom: <seq-viewer>">
          <DomExample />
        </Panel>
        <Panel label="@molstar/seq-viewer/svg: renderSvg(doc, { width: 640, background: true })">
          <SvgExample />
        </Panel>
      </div>
    </Unstyled>
  );
}

function kernel() {
  const reg = createRegistries();
  register(reg, [KernelRegistry, MvsRegistry]);
  return reg;
}

/** JSON indented, with arrays of plain values kept on one line. */
const compactJson = (v: unknown) =>
  JSON.stringify(v, null, 2).replace(
    /\[\s+([^[\]{}]*?)\s+\]/g,
    (_, inner: string) => `[${inner.replace(/\s*\n\s*/g, " ")}]`,
  );

/** The example document as JSON, as the builder wrote it. */
export function ExampleJson() {
  return (
    <Unstyled>
      <pre className={s.code}>{compactJson(doc)}</pre>
    </Unstyled>
  );
}

/** The example with two mistakes, and what `validate` says about them. */
function broken(): SeqViewSpecDoc {
  const json = JSON.stringify(doc)
    .replace('"data":"domains"', '"data":"domain"')
    .replace('"labels":true', '"lables":true');
  return JSON.parse(json) as SeqViewSpecDoc;
}
export function ValidationExample() {
  const issues = useMemo(() => validate(broken(), kernel()), []);
  return (
    <Unstyled>
      <pre className={s.code}>{issues.map((i) => `${i.severity} ${i.code} ${i.path}\n  ${i.message}`).join("\n")}</pre>
    </Unstyled>
  );
}

/** The format's reference, generated from the parameter tables and the kernel's and MVS's registries. */
export function FormatReference() {
  const text = useMemo(() => printReference(kernel()), []);
  return (
    <details className={s.details}>
      <summary>The generated reference ({text.split("\n").length.toLocaleString("en-US")} lines)</summary>
      <Markdown options={{ disableParsingRawHTML: true }}>{text}</Markdown>
    </details>
  );
}

/** A query read, bound and checked against the base language: the JSON (full names), and every issue. */
const bound = (text: string) => checkText(text, createLanguage().language);

/** Queries as written and as the printer writes them back; with a JSON, the stored form of the first. */
export function QueryExamples({ texts, json = false }: { texts: readonly string[]; json?: boolean }) {
  const out = useMemo(
    () =>
      texts.map((t) => {
        const r = bound(t);
        return r.expr === null ? r.issues.map(renderIssue).join("\n") : printQuery(r.expr);
      }),
    [texts],
  );
  const first = useMemo(() => (json && texts[0] ? bound(texts[0]).expr : null), [json, texts]);
  return (
    <Unstyled>
      <pre className={s.code}>{out.join("\n")}</pre>
      {first !== null && <pre className={s.code}>{compactJson(first)}</pre>}
    </Unstyled>
  );
}

/** What the reader and binder say about queries with mistakes, rendered as a person or the agent sees them. */
export function QueryErrors({ texts }: { texts: readonly string[] }) {
  const out = useMemo(
    () => texts.map((t) => `${t}\n${bound(t).issues.map(renderIssue).join("\n")}`).join("\n\n"),
    [texts],
  );
  return (
    <Unstyled>
      <pre className={s.code}>{out}</pre>
    </Unstyled>
  );
}

/** SeqQL's base table, counted from the language the binder uses. */
export function QueryCounts() {
  const { symbols, aliases } = useMemo(() => createLanguage().language, []);
  return (
    <>
      {symbols.size} base symbols, {aliases.size} short names
    </>
  );
}

/** SeqQL's reference, generated from the symbol table: every symbol's full entry. */
export function QueryReference() {
  const text = useMemo(() => referenceText(createLanguage().language), []);
  return (
    <details className={s.details}>
      <summary>The generated reference ({text.split("\n").length.toLocaleString("en-US")} lines)</summary>
      <Markdown options={{ disableParsingRawHTML: true }}>{text}</Markdown>
    </details>
  );
}

/** The agent's skill as generated for the base language, with no host options. */
export function QuerySkill() {
  const text = useMemo(() => skillText(createLanguage().language), []);
  return (
    <Unstyled>
      <pre className={s.code}>{text}</pre>
    </Unstyled>
  );
}
