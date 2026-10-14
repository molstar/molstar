/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// Insertion codes in the model and SeqQL (seqspec:model#numberings; seqspec:seqql#the-sel-shorthand): a heavy chain in
// Kabat-style author numbers (52A, 82A-82C, 100A-100B) on its unit, `label` on its molecule. A label with a code names
// its residue, a range by plain numbers holds the inserted ones between them, `(label)` reads the codes, an absent code
// says what the labels are; `by_numbering` between two chains pairs codes with codes.

import { type SeqSet, SeqSets, type SetInput } from "@molstar/seq-core/model";
import { compileQuery, createQueryEnv, SeqSelection } from "@molstar/seq-core/model/query";
import { checkText } from "@molstar/seq-core/query";
import { describe, expect, it } from "vitest";

const VH =
  "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
/** Kabat-style labels over `length` residues, with the codes at the given numbers. */
function kabat(length: number, codes: Record<number, string[]>): (number | string)[] {
  const out: (number | string)[] = [];
  for (let n = 1; out.length < length; n++) out.push(n, ...(codes[n] ?? []).map((c) => `${n}${c}`));
  return out.slice(0, length);
}
const H1 = kabat(VH.length, { 52: ["A"], 82: ["A", "B", "C"], 100: ["A", "B"] });
const H2 = kabat(VH.length, { 52: ["A", "B"], 82: ["A", "B", "C"] });

function chains(): SetInput {
  const mol = (key: string) => ({ key, residues: VH, numberings: [{ name: "label", method: "offset", first: 1 }] });
  const auth = (labels: (number | string)[]) => [{ name: "auth", method: "labels", labels, isDefault: true }];
  return {
    molecules: [mol("e1"), mol("e2")],
    units: [
      { key: "H", molecule: "e1", numberings: auth(H1) },
      { key: "K", molecule: "e2", numberings: auth(H2) },
    ],
    mappings: [{ key: "kabat", from: "H", to: "K", method: "by_numbering", numbering: "auth" }],
  } as SetInput;
}

const env = createQueryEnv([]);
function run(set: SeqSet, text: string) {
  const c = checkText(text, env.language);
  expect(c.issues).toEqual([]);
  const r = compileQuery(c.expr as never, env).run({ set, inputs: {} });
  return { summary: SeqSelection.summary(r.selection, { maxParts: 8 }), said: r.issues.map((i) => i.message) };
}

describe("insertion codes", () => {
  const { set, issues } = SeqSets.build(chains());
  it("build cleanly", () => expect(issues).toEqual([]));

  it.each([
    ['(sel "H[52A]")', "1 group, 1 position: H 52A"],
    ['(sel "H[auth:52A]")', "1 group, 1 position: H 52A"],
    // A range by plain numbers holds the inserted residues between them; one from a code starts there.
    ['(sel "H[52-53]")', "1 group, 3 positions: H 52-53"],
    ['(sel "H[52A-53]")', "1 group, 2 positions: H 52A-53"],
    ['(sel "H[82-83]")', "1 group, 5 positions: H 82-83"],
    ['(sel "H[100A-101]")', "1 group, 3 positions: H 100A-101"],
    // `label` (label_seq_id) counts the inserted residues: its 53 is auth 52A.
    ['(sel "H[label:53]")', "1 group, 1 position: H 52A"],
    ['(positions :in H :where (matches (label) "[A-Z]$"))', "3 groups, 6 positions: H 52A, 82A-82C, 100A-100B"],
  ])("%s", (text, summary) => {
    expect(run(set, text)).toEqual({ summary, said: [] });
  });

  it("a code the chain does not have selects nothing and says what its labels are", () => {
    const r = run(set, '(sel "H[52B]")');
    expect(r.summary).toBe("selects nothing");
    expect(r.said[0]).toMatch(/^H has no labels 52B; its labels run 1-/);
  });

  it("by_numbering pairs codes with codes: H's 52A is K's 52A; K's 52B and H's 100A-100B pair nothing", () => {
    const m = set.mappings.find((x) => x && set.strings.get(x.key) === "kabat");
    if (!m) throw new Error("no mapping");
    const label = (labels: (number | string)[], i: number) => (i < 0 ? null : String(labels[i]));
    const at = (l: string) => H1.findIndex((x) => String(x) === l);
    expect(["52", "52A", "53", "82C", "100", "100A", "100B", "101"].map((l) => label(H2, m.forward[at(l)]))).toEqual([
      "52",
      "52A",
      "53",
      "82C",
      "100",
      null,
      null,
      "101",
    ]);
  });
});
