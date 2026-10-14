/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `query` computed provider (seqspec:seqql#in-documents; seqspec:format#computed-nodes, QD1): a stored SeqQL query
// checked statically (the query language's codes as `svs.query.*`, `search` and the host's state refused, its refs the node's
// dependencies with `sel` expanded, a stale text a warning, the version rule) and run at load into a ranges annotation,
// after the queries it names.

import { type Annotation, Issue, type SeqSet } from "@molstar/seq-core/model";
import { bindQuery, createLanguage, type Expression, Q, readQuery, type SeqQLQuery } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";
import { compile, type SeqViewSpecDoc, validate } from "../../src/format/index.ts";
import { CASES } from "./cases.ts";
import { fixtureLoaders, registries } from "./load-fixtures.ts";

/** A query as stored, bound but not checked (so a test can store a wrong one). */
const stored = (text: string): SeqQLQuery => ({
  seqql: "1.0",
  expr: bindQuery(readQuery(text).node as never, createLanguage().language).expr as Expression,
  text,
});
/** Case 1 (p53 and its regions) with query nodes added to its data. */
function withQueries(queries: Record<string, SeqQLQuery>): SeqViewSpecDoc {
  const doc = CASES.case1();
  const data = (doc.root.children ?? []).find((c) => c.kind === "data") as { children: unknown[] };
  for (const [ref, query] of Object.entries(queries))
    data.children.push({ kind: "computed", ref, params: { provider: "query", query } });
  return doc;
}
const reg = registries();
const said = (doc: SeqViewSpecDoc) =>
  validate(doc, reg).map((i) => `${i.severity} ${i.code} ${i.path} ${i.message}`.replace(/\/root\/children\/0/, ""));
const annotationOf = (set: SeqSet, ref: string) =>
  set.annotations.find((a) => a && set.strings.get(a.key) === ref) as Annotation;

describe("the query provider's static checks", () => {
  it("accepts a query, and names what is wrong with one in svs.query.* codes", () => {
    expect(said(withQueries({ big: Q.read("(items :of regions :where (> (length) 50))") }))).toEqual([]);
    const issues = said(
      withQueries({
        a: stored("(search GXGXXG :in p53)"),
        b: stored("(and p53 (input selection))"),
        c: stored("(expand (current) 5)"),
        d: stored("(slice p53 9)"),
        e: { ...Q.read("(range p53 10-20)"), seqql: "2.0" },
        f: { ...Q.read("(range p53 10-20)"), text: "(range p53 10-20 )" },
      }),
    );
    expect(issues.map((i) => i.split(" ").slice(0, 3).join(" "))).toEqual([
      "error svs.query.search /children/9/params/query/expr",
      "error svs.query.input /children/10/params/query/expr/args/1/0",
      "error svs.query.input /children/11/params/query/expr/args/0",
      "error svs.query.invalid /children/12/params/query/expr/args/1/0",
      "error svs.query.version /children/13/params/query/seqql",
      "warning svs.query.stale /children/14/params/query/text",
    ]);
    expect(issues[0]).toMatch(/a document query cannot search; run the search, then save its hits/);
    expect(issues[1]).toMatch(/data only: it cannot read what the host passes/);
    expect(issues[4]).toMatch(/this query is SeqQL 2.0; this reader reads 1.x$/);
    expect(issues[5]).toMatch(/the printed expr, which wins: \(range p53 10-20\)$/);
  });

  it("gives the refs it names, sel's expanded, as the node's: unknown ones and cycles are the format's errors", () => {
    expect(
      said(withQueries({ a: Q.read('(or (sel "p53[1-5], nope[1-5]") regions)') })).map((i) => i.slice(0, 77)),
    ).toEqual(['error svs.ref.unresolved /children/9/params/query/expr/args/0 "nope" names no']);
    const cycle = said(withQueries({ a: Q.read("(and p53 b)"), b: Q.read("(and p53 a)") }));
    expect(cycle.map((i) => i.split(" ").slice(0, 2).join(" "))).toEqual(["error svs.ref.cycle"]);
  });
});

describe("the query provider at load", () => {
  it("runs each query after the ones it names: one ranges item per group, with its item's id and label", async () => {
    const doc = withQueries({
      outside: Q.read("(and p53 (not big))"),
      big: Q.read("(items :of regions :where (> (length) 50))"),
      gap: Q.read("(range p53 1-5 :numbering nope)"),
    });
    const { set, issues } = await compile(doc, reg, fixtureLoaders());
    const big = annotationOf(set, "big");
    expect([big.geometry, big.count, big.itemCount]).toEqual(["ranges", 1, 1]);
    expect([big.start[0], big.stop[0], set.strings.get(big.itemKey[0])]).toEqual([101, 292, "q0"]);
    expect(String(set.strings.get(big.rowLabel?.[0] ?? -1))).toMatch(/binding/i);
    const outside = annotationOf(set, "outside");
    expect([outside.count, Array.from(outside.start), Array.from(outside.stop)]).toEqual([2, [0, 292], [101, 393]]);
    // a query written wrong in a way only its run can see is invalid input, an error that says what to write; its
    // node is left out, nothing cascading (GV-09)
    expect(annotationOf(set, "gap")).toBeUndefined();
    expect(issues.map((i) => `${i.severity} ${i.code}`)).toEqual(["error svs.query.invalid"]);
    expect(issues[0].message).toMatch(/^seqql.numbering.unknown: p53\/unit has no numbering "nope"/);
    void Issue;
  });

  it("a run's issues by kind: a bad regex invalid, an empty result a note", async () => {
    const doc = withQueries({
      regex: Q.read('(items :of regions :where (matches (label) "("))'),
      none: Q.read('(items :of regions :where (= (label) "no such label"))'),
    });
    const { issues } = await compile(doc, reg, fixtureLoaders());
    expect(issues.map((i) => `${i.severity} ${i.code} ${i.message.slice(0, i.message.indexOf(":"))}`)).toEqual([
      "error svs.query.invalid seqql.type",
      "info svs.query.note seqql.empty",
    ]);
  });
});
