/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type Apply,
  bindQuery,
  checkText,
  createLanguage,
  type Expression,
  printQuery,
  type QueryIssue,
  readQuery,
  referenceText,
  renderIssue,
  skillText,
} from "@molstar/seq-core/query";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { useMemo, useState } from "react";
import s from "./support/seqstar-format.module.css";
import { PrismPage } from "./support/seqstar-prism.js";
import q from "./support/seqstar-query.module.css";
import { ScrollArea, TextAction, TextTabs } from "./ui/index.js";

/**
 * SeqQL's surface: an S-expression box and what the reader, the binder
 * and the checker make of it: the issues under the line and column they point at (inside a `sel` string too), the
 * expansion the runtime will compile beside the box, the canonical print a document stores, the JSON with full names,
 * and how each name resolved. A second story shows the generated skill and reference. Nothing runs here; the runtime has its own gallery.
 */
const meta = { title: "Seq*/Query and search/SeqQL", parameters: { layout: "fullscreen" } } satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

const LANGUAGE = createLanguage().language;
const EXAMPLES: readonly (readonly [string, string])[] = [
  [
    "Slots",
    "(pick (items :of domains :where (> (length) 30)) :where (and (matches (label) kinase) (>= (field identity :of stats) 0.9)))",
  ],
  ["Operand kind", "(and 1brs/A (not (partners 1brs/B :via contacts)) (letters C))"],
  ["onto", '(onto (sel "1brs/A[10-20]") msa :keep true)'],
  ["A mix", "(and 1brs/A (> (length) 3))"],
  ["Did you mean", "(items :in domains :where (letter))"],
  ["Brackets", "(range 1brs/A[10-20])"],
  ["sel", '(and (sel "1brs/A[10-20, 45], $selection@0:5") (not (sel "domains[#PF00069]")))'],
  ["sel errors", '(or (sel "1brs/A[10-20")\n    (sel "1brs/B@20:9"))'],
  ["Checker", "(pick (each (items :of domains) (expand (current) -1))\n  :where (> (count 123) (seq.position.index)))"],
];

const at = (e: Expression, path: string) =>
  path
    .split("/")
    .slice(1)
    .reduce<unknown>((x, k) => (x as Record<string, unknown>)[k], e) as Apply;

/** The text with each issue under its line: a caret at its column, then the issue. */
function inline(text: string, issues: readonly QueryIssue[]): string {
  const lines = text.split("\n");
  const out = lines.flatMap((line, k) => [
    line,
    ...issues
      .filter((i) => i.at?.line === k + 1)
      .map((i) => {
        const pad = " ".repeat((i.at?.column ?? 1) - 1);
        return `${pad}^ ${renderIssue(i).replace(/^\S+ /, "").replace("\n  ", `\n${pad}  `)}`;
      }),
  ]);
  const placeless = issues.filter((i) => !i.at).map(renderIssue);
  return [...out, ...placeless].join("\n");
}

function analyze(text: string) {
  const checked = checkText(text, LANGUAGE);
  const read = readQuery(text);
  if (!read.node) return { ...checked, names: [] };
  const b = bindQuery(read.node, LANGUAGE);
  const names = [...b.source].flatMap(([path, n]) =>
    n.kind === "list" && b.expr !== null
      ? [
          {
            path,
            written: String(n.items[0]?.kind === "atom" && n.items[0].value),
            place: n.at,
            full: at(b.expr, path).head,
          },
        ]
      : [],
  );
  return { ...checked, names };
}

function QueryLanguage() {
  const [text, setText] = useState(EXAMPLES[0]?.[1] ?? "");
  const r = useMemo(() => analyze(text), [text]);
  return (
    <PrismPage title="Query language">
      <div className={s.row}>
        {EXAMPLES.map(([label, e]) => (
          <TextAction key={label} onClick={() => setText(e)}>
            {label}
          </TextAction>
        ))}
      </div>
      <div className={q.pair}>
        <textarea
          className={q.box}
          value={text}
          onChange={(e) => setText(e.target.value)}
          spellCheck={false}
          rows={5}
          aria-label="SeqQL"
        />
        <section className={s.card}>
          <h2 className={s.heading}>Expansion, as the runtime compiles it</h2>
          {r.expanded === null ? (
            <div className={s.none}>none while there are errors</div>
          ) : (
            <pre className={q.text}>{printQuery(r.expanded)}</pre>
          )}
        </section>
      </div>
      <section className={s.card}>
        <h2 className={s.heading}>Issues</h2>
        {r.issues.length === 0 ? (
          <div className={s.none}>none</div>
        ) : (
          <pre className={q.text}>{inline(text, r.issues)}</pre>
        )}
      </section>
      {r.expr !== null && (
        <>
          <section className={s.card}>
            <h2 className={s.heading}>Canonical print</h2>
            <pre className={q.text}>{printQuery(r.expr)}</pre>
          </section>
          <section className={s.card}>
            <h2 className={s.heading}>Names</h2>
            <table className={s.table}>
              <thead>
                <tr>
                  <th>At</th>
                  <th>Written</th>
                  <th>Bound to</th>
                  <th>JSON path</th>
                </tr>
              </thead>
              <tbody className={s.mono}>
                {r.names.map((n) => (
                  <tr key={n.path}>
                    <td>{`${n.place.line}:${n.place.column}`}</td>
                    <td>{n.written}</td>
                    <td>{n.full}</td>
                    <td>{n.path || "/"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </section>
          <section className={s.card}>
            <h2 className={s.heading}>JSON, full names</h2>
            <ScrollArea className={s.box} orientation="both">
              <pre className={s.pre}>{JSON.stringify(r.expr, null, 2)}</pre>
            </ScrollArea>
          </section>
        </>
      )}
    </PrismPage>
  );
}

export const ReadBindPrint: Story = { name: "Read, bind, check, print", render: () => <QueryLanguage /> };

/** The host's lines of the skill: Stories' ref line and inputs. */
const OPTIONS = {
  refs: "Refs: 1brs/A (a chain's sequence), p53 (a molecule), msa (an alignment), msa/3 (its row 3), domains (an annotation), ctx.site.hinge (a Context entity).",
  inputs: [
    { name: "selection", says: "the person's selection" },
    { name: "focus", says: "the Focus" },
  ],
};

function SkillAndReference() {
  const [tab, setTab] = useState("skill");
  const texts = useMemo(() => ({ skill: skillText(LANGUAGE, OPTIONS), reference: referenceText(LANGUAGE) }), []);
  const shown = tab === "skill" ? texts.skill : texts.reference;
  const words = shown.split(/\s+/).filter(Boolean).length;
  return (
    <PrismPage
      title="Skill and reference"
      caption="Generated from the symbol table: the skill the agent reads, and every symbol's full entry, loaded on demand."
    >
      <TextTabs
        value={tab}
        onValueChange={setTab}
        tabs={[
          { value: "skill", label: "Skill" },
          { value: "reference", label: "Reference" },
        ]}
        end={<span className={s.kicker}>{`${words} words`}</span>}
      />
      <ScrollArea className={q.page} orientation="both">
        <pre className={s.pre}>{shown}</pre>
      </ScrollArea>
    </PrismPage>
  );
}

export const Skill: Story = { name: "Skill and reference", render: () => <SkillAndReference /> };
