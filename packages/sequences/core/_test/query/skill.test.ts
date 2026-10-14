/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// QL2b: the generated skill and reference: every symbol in the skill, every example reads, binds
// and checks, every hint's code emitted, A.3's snapshot, an extension's symbols in its host's skill.
import { describe, expect, it } from 'vitest';
import {
  checkText,
  createLanguage,
  defineSymbol,
  type Language,
  printQuery,
  Q,
  referenceText,
  type SkillOptions,
  type SymbolDef,
  skillText,
} from '../../src/query/index.js';

/** Stories' options and its `chain` macro, as Appendix A.3 renders them. */
const chain = defineSymbol({
  name: 'stories.chain',
  alias: 'chain',
  kind: 'macro',
  args: [{ name: 'C', type: 'str', positional: true, doc: 'the chain id' }],
  result: 'query',
  signature: '(chain C)',
  doc: 'units of chain C',
  examples: [{ text: '(chain A)', says: 'chain A' }],
  hints: [],
  expand: (call) =>
    Q.seq.generator.systems({
      where: Q.core.rel.eq(Q.seq.system.field('auth_asym_id', { ns: 'mvs' }), call.args?.['0'] as string),
    }),
});
const STORIES: SkillOptions = {
  refs: "Refs: 1brs/A (a chain's sequence), p53 (a molecule), msa (an alignment), msa/3 (its row 3), domains (an annotation), ctx.site.hinge (a Context entity).",
  inputs: [
    { name: 'selection', says: "the person's selection" },
    { name: 'focus', says: 'the Focus' },
  ],
};
const STORIES_LANGUAGE = createLanguage([{ symbols: [chain] }]).language;

/**
 * A query that makes each static hint's code appear on its symbol (seqspec:seqql#errors); the others are the runtime's.
 */
const EMITTED: Record<string, string> = {
  'core.logic.and seqql.type': '(items :where (core.logic.and x (> (length) 3)))',
  'core.logic.or seqql.type': '(items :where (core.logic.or x (> (length) 3)))',
  'core.rel.lt seqql.type': '(items :where (< (label) 3))',
  'core.rel.lte seqql.type': '(items :where (<= (label) 3))',
  'core.rel.gr seqql.type': '(items :where (> (label) 3))',
  'core.rel.gre seqql.type': '(items :where (>= (label) 3))',
  'core.rel.in seqql.arg.missing': '(items :where (in (label)))',
  'core.str.matches seqql.arg.missing': '(items :where (matches (label)))',
  'seq.combinator.intersect seqql.type': '(and 1brs/A (> (length) 3))',
  'seq.modifier.complement seqql.arg.count': '(not 1brs/A 1brs/B)',
  'seq.macro.sel seqql.arg.missing': '(sel)',
  'seq.macro.sel seqql.type': '(sel 45)',
  'seq.macro.sel seqql.sel.syntax': '(sel "1brs/A[10-20]x")',
  'seq.macro.sel seqql.sel.mixed': '(sel "1brs/A[10,#x]")',
  'seq.generator.systems seqql.arg.unknown': '(systems :frobs x)',
  'seq.generator.items seqql.arg.unknown': '(items :frobs x)',
  'seq.generator.letters seqql.arg.missing': '(letters)',
  'seq.generator.search seqql.arg.missing': '(search)',
  'seq.generator.current seqql.current': '(current)',
  'seq.modifier.range seqql.arg.order': '(range 1brs/A :numbering auth 10-20)',
  'seq.modifier.slice seqql.arg.count': '(slice 1brs/A 9)',
  'seq.modifier.slice seqql.type': '(slice 1brs/A 45A 50)',
  'seq.modifier.expand seqql.arg.count': '(expand 1brs/A 5 6)',
  'seq.modifier.onto seqql.type': '(onto 1brs/A msa :kind codon)',
  'seq.modifier.group seqql.arg.missing': '(group 1brs/A)',
  'seq.filter.pick seqql.arg.missing': '(pick (all))',
  'seq.filter.ids seqql.arg.missing': '(ids domains)',
  'seq.system.field seqql.field.unknown': '(systems :where (= (field auth) A))',
};
const RUNTIME = new Set(
  'ref.unknown ref.kind numbering.unknown slice.range slice.clipped field.unknown input.unknown classes.unknown limit.positions limit.runs limit.groups limit.time label.unknown onto.ambiguous onto.unmapped lazy empty version'
    .split(' ')
    .map((c) => `seqql.${c}`),
);
const symbolsOf = (language: Language): SymbolDef[] => [...language.symbols.values()];

describe('the skill', () => {
  it("renders the v1 core with Stories' options as its snapshot (updated on purpose)", async () => {
    await expect(skillText(STORIES_LANGUAGE, STORIES)).toMatchFileSnapshot('snapshots/skill-stories.txt');
  });

  it('renders the v1 core alone as packages/sequences/_spec/reference/seqql-skill.md', async () => {
    const core = createLanguage().language;
    await expect(
      `# SeqQL skill\n\nGenerated from the symbol table; do not edit.\n\n\`\`\`text\n${skillText(core)}\`\`\`\n`,
    ).toMatchFileSnapshot('../../../_spec/reference/seqql-skill.md');
  });

  it('lists every symbol of the language, by its signature', () => {
    const skill = skillText(STORIES_LANGUAGE);
    const missing = symbolsOf(STORIES_LANGUAGE).filter((s) => !skill.includes(s.signature));
    expect(missing.map((s) => s.name)).toEqual([]);
    expect(skill).not.toMatch(/\{(refs|inputs|table)\}|\nRefs:|\nInputs:/);
  });

  it("puts an extension's symbols in its host's skill and reference, by kind and slot", () => {
    const evalue = defineSymbol({
      name: 'x.item.evalue',
      alias: 'evalue',
      kind: 'property',
      reads: 'item',
      args: [],
      result: 'num',
      signature: '(evalue)',
      doc: "a hit's e-value",
      examples: [{ text: '(items :of hits :where (< (evalue) 0.001))', says: 'significant hits', skill: 3 }],
      hints: [],
    });
    const near = defineSymbol({ ...chain, name: 'x.near', alias: null, signature: '(x.near Q n)', doc: 'near Q' });
    const { language, issues } = createLanguage([{ symbols: [evalue, near] }]);
    expect(issues).toEqual([]);
    const skill = skillText(language).split('\n');
    expect(skill.find((l) => l.startsWith('item '))).toBe(
      'item (items, partners :where): (id) (label) (category) (value) (type) (annotation) (length) (evalue)',
    );
    expect(skill.filter((l) => l.startsWith('(x.near') || l.startsWith('(current)'))).toEqual([
      '(current)  inside each, the group being run on',
      '(x.near Q n)  near Q',
    ]);
    const examples = skill.slice(skill.indexOf('EXAMPLES') + 2);
    expect(examples.slice(2, 4)).toEqual([
      '(pick (items :of domains) :where (matches (label) kinase))  kinase domains',
      '(items :of hits :where (< (evalue) 0.001))  significant hits',
    ]);
    expect(referenceText(language)).toContain('### evalue\n`x.item.evalue` · property · result: number · reads: item');
  });
});

describe("the table's metadata, as the skill and the checker read it", () => {
  it('every example reads, binds and checks, and prints as written', () => {
    const examples = symbolsOf(STORIES_LANGUAGE).flatMap((s) => s.examples.map((e) => e.text));
    expect(examples.length).toBeGreaterThan(60);
    for (const text of examples) {
      const r = checkText(text, STORIES_LANGUAGE);
      expect([text, r.issues]).toEqual([text, []]);
      expect(printQuery(r.expr as never, STORIES_LANGUAGE)).toBe(text);
    }
  });

  it("every hint's code is one the binder or checker emits on its symbol, or the runtime's", () => {
    const language = createLanguage().language;
    for (const s of symbolsOf(language)) {
      for (const h of s.hints) {
        const key = `${s.name} ${h.code}`;
        const text = EMITTED[key] ?? (h.code === 'seqql.slot' ? `(${s.name})` : undefined);
        if (text === undefined) {
          expect([key, RUNTIME.has(h.code)]).toEqual([key, true]);
          continue;
        }
        const issues = checkText(text, language).issues;
        const on = issues.filter((i) => i.code === h.code && i.signature?.includes(s.signature));
        expect([key, on.length > 0]).toEqual([key, true]);
      }
    }
  });
});

describe('the reference', () => {
  it("renders every symbol's full entry", async () => {
    const reference = referenceText(STORIES_LANGUAGE);
    for (const s of symbolsOf(STORIES_LANGUAGE)) expect(reference).toContain(`\`${s.name}\` · ${s.kind}`);
  });

  it('renders the v1 core as packages/sequences/_spec/reference/seqql.md', async () => {
    await expect(referenceText(createLanguage().language)).toMatchFileSnapshot('../../../_spec/reference/seqql.md');
  });
});
