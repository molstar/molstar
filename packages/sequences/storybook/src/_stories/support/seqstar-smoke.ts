/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What each Seq* story opens, for the smoke test (packages/sequences/storybook/_test/smoke.test.ts and smoke-more.test.ts,
// half each): per story, the documents it shows, how it loads them, the patches its buttons make and the warnings it
// shows on purpose. The test opens each on a headless viewer and fails on an error, on a warning not listed here, and
// on a visible section with nothing in it; and it fails when a story file has an export this list does not name, so a
// new story is never forgotten. No React here.
//
// A story is named by its file and export: `viewer.stories.tsx#Case1`. Add a new story's entry beside its neighbors:
//   - `cases`: the documents it shows (a function, so nothing is built when the story is not run);
//   - `same`: the key of the story whose documents it shows, which is smoked there;
//   - `skip`: why it has no document to open (a live source, a page of text).

import type { Loaders, Patch, SeqViewSpecDoc } from '@molstar/seq-core/format';
import { type SeqViewerContext, type SeqViewerRegistryEntry, sectionInfo } from '@molstar/seq-viewer';
import { renderSvg } from '@molstar/seq-viewer/svg';
import pf00042 from '../../../../core/_test/fixtures/pf00042.json?raw';
import { CASES } from '../../../../core/_test/format/cases.js';
import { genomeDoc, synthGenome } from '../../../../core/_test/format/genome-fixture.js';
import { familyWithTree, syntheticFamily } from '../../../../viewer/_test/families.js';
import { keyedFamilyDoc, MODES, modePatch } from '../../../../viewer/_test/keyed-family.js';
import { lazyStoryDoc, lazyStoryLoaders } from '../../../../viewer/_test/lazy-story-fixtures.js';
import { splitUnitsDoc } from '../../../../viewer/_test/split-units-fixtures.js';
import { described, hideEmpty, textStyled } from './seqstar-cases.js';
import { circularGenomeDoc, PLASMID_RING, rotateRing, zoomDetail } from './seqstar-docs-circular.js';
import {
  alignmentGallery,
  customDoc,
  HOST,
  manyDoc,
  NUCLEIC,
  PROTEIN,
  pickerDoc,
  rangeDoc,
  sharedDoc,
  valuesGallery,
} from './seqstar-docs-colors.js';
import { sashimiDoc, scalesDoc, sizesDoc } from './seqstar-docs-sizes.js';
import {
  atlasDoc,
  atlasLoaders,
  CODING_PROTEIN,
  codingLoaders,
  codingStoryDoc,
  dynamicDoc,
  dynamicKnobs,
  dynamicLoaders,
  hoverResidue,
  items,
  keyedDoc,
  keyedLoaders,
  PEAKS,
  STREAMED_AT,
  STREAMED_SECTION,
  streamedDoc,
  streamedKnobs,
  streamedLoaders,
  transcriptLoaders,
  transcriptStoryDoc,
} from './seqstar-docs-streamed.js';
import {
  antibodyPlain,
  antibodyWithSnapshots,
  asText,
  BAD,
  BOOKMARKS,
  codingPair,
  computedTreeDoc,
  describedCase10,
  exonsOnly,
  FLOORS,
  familyDoc,
  flipped,
  HEADS,
  LOCATION_SECTION,
  largeAlignment,
  letterMapStrands,
  linkStrands,
  locationSetsDoc,
  logoDoc,
  longSequence,
  pDemoRuns,
  plasmid,
  realignedDoc,
  STYLE_FAMILY,
  styled,
  tracksDoc,
} from './seqstar-docs-views.js';
import { calmodulinDotPlot } from './seqstar-dot-plots.js';
import { FILE_FORMATS_SMOKE } from './seqstar-file-formats-smoke.js';
import { GENOTYPE_DENSITIES, GENOTYPE_MODES, genotypeDensityPatch, genotypeModePatch } from './seqstar-genotypes.js';
import {
  IDEOGRAM_PLACES,
  IDEOGRAM_REGISTRY,
  IDEOGRAM_SECTION,
  ideogramDoc,
  ideogramLoaders,
} from './seqstar-ideogram.js';
import { INDEL_LIT, INDEL_SELECTIONS, litPerSection, selectCommand } from './seqstar-indels.js';
import { p53MolstarCase } from './seqstar-molstar-smoke.js';
import { storyPatch } from './seqstar-patch.js';
import { AskedLog, PREFETCH_SECTION, PREFETCH_SNAPSHOTS, prefetchDoc, prefetchLoaders } from './seqstar-prefetch.js';
import { TRIAGE_MODES, triagePatch } from './seqstar-triage.js';

/** A warning a story shows on purpose: its code, and text its message must hold when given. */
export interface AllowedWarning {
  readonly code: string;
  readonly includes?: string;
}

/**
 * A step after the document is open: what a button does. A patch goes through `storyPatch`, as the stories send it.
 * `run` may return `false` (a refused action, an expectation that did not hold) or issues; both are findings unless
 * allowed.
 */
export interface SmokeStep {
  readonly label: string;
  readonly run: (ctx: SeqViewerContext) => Promise<unknown> | unknown;
}

/** One document a story opens, with the way it opens it. */
export interface SmokeCase {
  readonly label: string;
  readonly doc: SeqViewSpecDoc;
  /** The loaders for the document's sources (default: the format's fixtures). */
  readonly loaders?: Loaders;
  /** Registry entries after the stories' own (a host's color scheme). */
  readonly registry?: readonly SeqViewerRegistryEntry[];
  /** The view the story shows first (default: the document's). */
  readonly view?: string;
  /** What the story's buttons do, each after the snapshots, in order. */
  readonly steps?: readonly SmokeStep[];
  /** Warnings the story shows on purpose. */
  readonly allow?: readonly AllowedWarning[];
  /** The story streams: after its steps the stream manager has asked for regions, and none failed. */
  readonly streams?: true;
  /** False: the story registers the structure extension itself (`@molstar/seq-molstar-poc`), so the fixture entries stay out. */
  readonly fixtureStructures?: false;
}

/** What the test supplies for documents that come from files. */
export interface SmokeEnv {
  /** A file of the repository, by its path from the repository root. */
  readFile(path: string): Uint8Array;
}

export interface SmokeStory {
  /** `file#Export`, the file under packages/sequences/storybook/src/_stories/. */
  readonly key: string;
  readonly cases?: (env: SmokeEnv) => readonly SmokeCase[] | Promise<readonly SmokeCase[]>;
  /** The key of the story with the same documents. */
  readonly same?: string;
  /** Why there is nothing to open. */
  readonly skip?: string;
}

/** A patch as a button sends it. */
const patchStep = (label: string, patch: Patch): SmokeStep => ({
  label,
  run: async (ctx) => (await storyPatch(ctx, patch, label)).issues,
});

/** A section header's switch as a click sends it: its action with the section and the value; then the param holds it. */
const sectionSwitch = (section: string, action: string, value: string): SmokeStep => ({
  label: `${action} ${value}`,
  run: async (ctx) => {
    await ctx.actions.run(action, { section, value }, 'host:smoke');
    await ctx.graph.update();
    return sectionInfo(ctx, section).switches.find((x) => x.action === action)?.value === value;
  },
});

// ---- The format's cases ----

/**
 * The warnings the acceptance cases show on purpose: each exercises a diagnostic. Notes (`info`: case 3's and case
 * 15's inserts) are no warnings.
 */
const CASE_WARNINGS: Readonly<Record<string, readonly AllowedWarning[]>> = {
  // A descending row header (NAME/137-27) gives no row numbering.
  case4: [{ code: 'svs.numbering.descending' }],
  // Case 7 names a provider no one registered ("hyd"), and a chain-each view whose data reaches no block.
  case7: [{ code: 'svs.provider.unavailable' }, { code: 'svs.layer.unplaced' }],
};
const caseOf = (name: string, make: () => SeqViewSpecDoc): SmokeCase => ({
  label: name,
  doc: described(name, make()),
  allow: CASE_WARNINGS[name],
});
const format =
  (name: string): SmokeStory['cases'] =>
  () => [caseOf(name, CASES[name])];
const FORMAT_CASES = ['case1', 'case2', 'case3', 'case4', 'case5', 'case6', 'case7', 'case8', 'case9'] as const;

// ---- Pieces used by several stories ----

const FAMILY_160 = syntheticFamily(160, 10);
const TREE_FAMILY = familyWithTree(160, 16, { columns: true, sites: true });
const KEYED_FAMILY = (JSON.parse(pf00042) as { rows: [string, number, string][] }).rows.map(([label, , gapped]) => ({
  label,
  gapped,
}));

export const SMOKE: readonly SmokeStory[] = [
  // ---- Seq*/Viewer ----
  ...FORMAT_CASES.map(
    (c, k): SmokeStory => ({
      key: `viewer.stories.tsx#Case${k + 1}`,
      // Case 1's text view in either text preset (its `Text` choice).
      cases:
        c === 'case1'
          ? () => [
              caseOf(c, CASES.case1),
              { ...caseOf(c, () => textStyled(CASES.case1(), 'compact')), label: 'case1, compact text' },
            ]
          : format(c),
    }),
  ),
  { key: 'viewer.stories.tsx#Case10', cases: format('case10') },
  {
    key: 'viewer.stories.tsx#Case12',
    cases: () => [
      {
        ...caseOf('case12', CASES.case12),
        // The section header's switches (the story places the header): the dimer as its two copies, and back; the heme
        // as its atoms, and back (the Ligand switch stays offered while the atoms are shown).
        steps: [
          sectionSwitch('design-view.main', 'section.copies', 'each'),
          sectionSwitch('design-view.main', 'section.copies', 'one'),
          sectionSwitch('design-view.main', 'section.ligandForm', 'atoms'),
          sectionSwitch('design-view.main', 'section.ligandForm', 'position'),
        ],
      },
    ],
  },
  {
    key: 'viewer.stories.tsx#Case13',
    cases: () => [
      caseOf('case13', CASES.case13),
      // Its rows choice: only the residues a variant has.
      {
        ...caseOf('case13', () => {
          const d = CASES.case13();
          return { ...d, root: hideEmpty(d.root) };
        }),
        label: 'case 13, split.empty hide',
      },
    ],
  },
  { key: 'viewer.stories.tsx#Case14', cases: format('case14') },
  {
    key: 'viewer.stories.tsx#Case15',
    cases: () => [
      {
        ...caseOf('case15', CASES.case15),
        // Its designs' inserts land on the wild type as boundaries, and D4's deleted H102 on the junction it leaves.
        allow: [...(CASE_WARNINGS.case15 ?? []), { code: 'svs.layer.insert' }],
        // In each view, a selection lights what the maps carry it to (the runs written out in seqstar-indels.ts).
        steps: Object.entries(INDEL_LIT).flatMap(([view, lit]) =>
          INDEL_SELECTIONS.map(
            (sel): SmokeStep => ({
              label: `${view}: select ${sel.label}`,
              run: async (ctx) => {
                ctx.views.setCurrent(view);
                await ctx.graph.update();
                const ok = await ctx.commands.run(selectCommand(sel), 'host:smoke');
                return ok && JSON.stringify(litPerSection(ctx)) === JSON.stringify(lit[sel.id]);
              },
            }),
          ),
        ),
      },
    ],
  },
  {
    key: 'viewer.stories.tsx#Case16',
    cases: () => [
      {
        ...caseOf('case16', CASES.case16),
        steps: [
          { label: 'rotate the origin', run: async (ctx) => (await storyPatch(ctx, rotateRing(450))).issues },
          { label: 'zoom the detail to the cloning site', run: zoomDetail },
          {
            // With its overlay: the ring's window, its `--seq-window` fill and its grips (the default skin's inks).
            label: 'export the ring as SVG, with the window',
            run: async (ctx) => {
              const svg = await renderSvg(ctx, {
                width: 640,
                background: true,
                overlay: true,
                sections: [PLASMID_RING],
              });
              return svg.includes('fill="#5b6370" fill-opacity="0.1"') && svg.includes('fill="#7c838e"');
            },
          },
        ],
      },
    ],
  },
  {
    key: 'viewer.stories.tsx#Case17',
    cases: () => [
      {
        ...caseOf('case17', CASES.case17),
        steps: [
          ...GENOTYPE_DENSITIES.map((d) => patchStep(d.label, genotypeDensityPatch(d.key) as Patch)),
          ...GENOTYPE_MODES.map((m) => patchStep(m.label, genotypeModePatch(m) as Patch)),
          {
            // The last mode keeps about 24 samples: few enough to draw as cells (a band a device pixel or more).
            label: 'the last mode draws cells',
            run: (ctx) => {
              const prepared = ctx.graph.assembled.value?.view?.prepared ?? [];
              const matrix = prepared.map((p) => p?.state as { stateOf?: Int32Array } | null).find((x) => x?.stateOf);
              const n = matrix?.stateOf?.length ?? 0;
              return n >= 20 && n <= 30;
            },
          },
        ],
      },
    ],
  },
  { key: 'viewer.stories.tsx#Case18', cases: format('case18') },
  {
    key: 'viewer.stories.tsx#Case19',
    cases: () => [
      {
        // Streamed by region under the ideogram navigator; its snapshots move through chr17's places.
        ...caseOf('case19', CASES.case19),
        registry: IDEOGRAM_REGISTRY,
        streams: true,
      },
    ],
  },
  {
    key: 'viewer.stories.tsx#Case20',
    cases: () => [
      {
        ...caseOf('case20', CASES.case20),
        streams: true,
        steps: [
          { label: 'residues 450–550', run: (ctx) => ctx.viewports.get(CODING_PROTEIN)?.zoomTo(0, 450, 550) },
          { label: "residue 500 hovered (its codon's letters asked)", run: (ctx) => hoverResidue(ctx, 500) },
          { label: "residue 500's codon, split by its intron", run: (ctx) => hoverResidue(ctx, 500) === 'T|CA' },
          {
            label: 'T1 is a coordinate system, shown by its section',
            run: (ctx) => ctx.graph.assembled.value?.set.byName('genes_src/T1', { kind: -1, id: -1, row: -1 }) === true,
          },
        ],
      },
    ],
  },
  { key: 'viewer.stories.tsx#Case21', cases: format('case21') },
  {
    key: 'viewer.stories.tsx#Case22',
    cases: () => [
      {
        ...caseOf('case22', CASES.case22),
        steps: TRIAGE_MODES.map((m) => patchStep(m.label, triagePatch(m) as Patch)),
      },
    ],
  },
  { key: 'alignments.stories.tsx#FamilyVertical', cases: () => [{ label: '160 rows', doc: FAMILY_160 }] },
  { key: 'alignments.stories.tsx#FamilyTree', cases: () => [{ label: '160 rows with a tree', doc: TREE_FAMILY }] },
  {
    key: 'alignments.stories.tsx#ComputedTrees',
    cases: () => [
      { label: 'NJ', doc: computedTreeDoc('nj') },
      { label: 'UPGMA', doc: computedTreeDoc('upgma') },
      { label: 'as given', doc: realignedDoc(false) },
      { label: 'realigned', doc: realignedDoc(true) },
    ],
  },
  { key: 'alignments.stories.tsx#RowControls', same: 'alignments.stories.tsx#FamilyTree' } as SmokeStory,
  {
    key: 'large-documents.stories.tsx#LargeAlignment',
    cases: () => [{ label: '10,000 x 300', doc: largeAlignment(10_000, 300) }],
  },
  {
    key: 'large-documents.stories.tsx#WideAlignment',
    cases: () => [{ label: '1,000 x 2,000', doc: largeAlignment(1_000, 2_000) }],
  },
  {
    key: 'large-documents.stories.tsx#LongSequence',
    cases: () => [{ label: '35,000 residues', doc: longSequence(35_000) }],
  },
  {
    key: 'large-documents.stories.tsx#GenomeScale',
    cases: () => [{ label: '3M nucleotides', doc: genomeDoc(synthGenome()) }],
  },
  {
    key: 'patches.stories.tsx#ViewerPatches',
    cases: () => [{ label: 'case 2 with snapshots', doc: antibodyWithSnapshots(CASES.case2()) }],
  },
  { key: 'side-by-side.stories.tsx#SideBySideNightingale', same: 'viewer.stories.tsx#Case7' },

  // ---- Seq*/Format: the same documents, shown as data ----
  ...[
    ...FORMAT_CASES,
    'case10',
    'case12',
    'case13',
    'case14',
    'case15',
    'case16',
    'case17',
    'case18',
    'case19',
    'case20',
    'case21',
    'case22',
  ].map(
    (c): SmokeStory => ({
      key: `format.stories.tsx#Case${c.slice(4)}`,
      same: `viewer.stories.tsx#Case${c.slice(4)}`,
    }),
  ),
  {
    key: 'format.stories.tsx#Svsx',
    skip: "No viewer: it shows the archive reader and the bundler, and the archive's relative paths do not validate in a viewer without its archive; the format's svsx tests cover both.",
  },

  // ---- Seq*/Nightingale: the same cases ----
  ...[
    ...FORMAT_CASES,
    'case10',
    'case12',
    'case13',
    'case14',
    'case15',
    'case16',
    'case17',
    'case18',
    'case19',
    'case20',
    'case21',
    'case22',
  ].map(
    (c): SmokeStory => ({
      key: `nightingale.stories.tsx#Case${c.slice(4)}`,
      same: `viewer.stories.tsx#Case${c.slice(4)}`,
    }),
  ),

  // ---- Seq*/React UI, Seq*/DOM UI, Seq*/Export, Seq*/UX drafts ----
  { key: 'react-ui.stories.tsx#Case10', same: 'viewer.stories.tsx#Case10' },
  { key: 'react-ui.stories.tsx#Composed', same: 'viewer.stories.tsx#Case10' },
  {
    key: 'patches.stories.tsx#SeqViewerPatches',
    cases: () => [{ label: 'case 2 with plain snapshots', doc: antibodyPlain() }],
  },
  {
    key: 'dom-ui.stories.tsx#DescriptionLinks',
    cases: () => [{ label: 'case 10 with markdown descriptions', doc: describedCase10() }],
  },
  {
    key: 'dom-ui.stories.tsx#DomComposed',
    cases: () => [{ label: 'case 10 with markdown descriptions', doc: describedCase10() }],
  },
  { key: 'export.stories.tsx#SvgAndPng', same: 'large-documents.stories.tsx#GenomeScale' },
  {
    key: 'ux-drafts.stories.tsx#Workspace',
    skip: "Hosted in apps/storybook (Stories' workspace screens): it opens viewer.stories.tsx#Case10's document.",
  },
  {
    key: 'ux-drafts.stories.tsx#Copies',
    skip: "Hosted in apps/storybook (Stories' workspace screens): it opens viewer.stories.tsx#Case7's document.",
  },
  {
    key: 'ux-drafts.stories.tsx#Reader',
    skip: "Hosted in apps/storybook (Stories' workspace screens): it opens viewer.stories.tsx#Case10's document.",
  },
  {
    key: 'ux-drafts.stories.tsx#Phone',
    skip: "Hosted in apps/storybook (Stories' workspace screens): it opens viewer.stories.tsx#Case10's document.",
  },
  {
    key: 'ux-drafts.stories.tsx#Navigation',
    skip: "Hosted in apps/storybook (Stories' workspace screens): it opens viewer.stories.tsx#Case6's document.",
  },

  // ---- Seq*/Cyclic ----
  {
    key: 'cyclic.stories.tsx#Linear',
    cases: () => [{ label: 'pDemo', doc: plasmid() }],
  },
  {
    key: 'cyclic.stories.tsx#Ring',
    cases: () => [
      {
        label: 'pDemo, ring',
        doc: plasmid(),
        view: 'ring',
        steps: [
          {
            // The marks the story sets once the ring is in: a selection and a highlight across the origin.
            label: 'selection and highlight across the origin',
            run: (ctx) => {
              ctx.selection.set(pDemoRuns(0, 225, 240, 0, 0, 12), 'api');
              ctx.highlight.set('story', pDemoRuns(0, 235, 240, 0, 0, 6));
            },
          },
          {
            // The figure the story draws beside the ring: `renderSvg` with its overlay.
            label: 'export the ring as SVG',
            run: async (ctx) =>
              (await renderSvg(ctx, { width: 480, overlay: true, hover: true, background: true })).includes('<svg'),
          },
        ],
      },
    ],
  },

  {
    key: 'circular-maps.stories.tsx#CircularGenome',
    cases: () => [{ label: 'a circular genome', doc: circularGenomeDoc() }],
  },

  // ---- Seq*/Hover ----
  {
    key: 'hover.stories.tsx#Emphasis',
    cases: () => [
      { label: 'every emphasized type', doc: tracksDoc('p53, every emphasized type') },
      { label: "a profile's cells", doc: logoDoc() },
    ],
  },
  {
    key: 'hover.stories.tsx#Off',
    cases: () => [{ label: 'emphasis off', doc: tracksDoc('p53, hover emphasis off', { hover_emphasis: false }) }],
  },
  {
    key: 'hover.stories.tsx#RowMarks',
    cases: () => [{ label: 'per-row marks', doc: splitUnitsDoc({ highlights: true }) }],
  },

  // ---- Seq*/Letter maps, Seq*/Links, Seq*/Pieces ----
  {
    key: 'letter-maps.stories.tsx#Maps',
    cases: () => [
      { label: 'strands', doc: letterMapStrands() },
      { label: 'exons only', doc: exonsOnly() },
      { label: 'as text', doc: asText() },
    ],
  },
  {
    key: 'links.stories.tsx#Follow',
    cases: () => [
      { label: 'a protein following its gene', doc: codingPair() },
      { label: 'the reverse strand following', doc: linkStrands() },
    ],
  },
  { key: 'pieces.stories.tsx#Direction', cases: () => [{ label: 'a construct, flipped', doc: flipped() }] },
  {
    key: 'dot-plots.stories.tsx#SelfComparison',
    cases: () => [{ label: 'calmodulin against itself', doc: calmodulinDotPlot() }],
  },

  // ---- Seq*/Style ----
  {
    key: 'style.stories.tsx#ColumnHeads',
    cases: () =>
      HEADS.map((style) => ({
        label: JSON.stringify(style),
        doc: styled(CASES.case3(), 'family.main', style),
        allow: CASE_WARNINGS.case3,
      })),
  },
  {
    key: 'style.stories.tsx#DataFloor',
    cases: () =>
      FLOORS.map((style) => ({ label: JSON.stringify(style), doc: styled(STYLE_FAMILY, 'fam.main', style) })),
  },
  {
    key: 'style.stories.tsx#InvalidKey',
    cases: () => [
      {
        label: 'a misspelled key and a bad value',
        doc: styled(CASES.case3(), 'family.main', BAD),
        // The story's point: two warnings that name what is allowed.
        allow: [...(CASE_WARNINGS.case3 ?? []), { code: 'svs.section.style' }],
      },
    ],
  },

  // ---- Seq*/Color schemes ----
  {
    key: 'color-schemes.stories.tsx#Gallery',
    cases: () => [
      {
        label: 'protein letters',
        doc: alignmentGallery('SH2 family (case 3)', 'PF00017.seed.sto', PROTEIN),
        allow: CASE_WARNINGS.case3,
      },
      {
        label: 'nucleotide letters',
        doc: alignmentGallery('tRNA family (case 4)', 'RF00005.seed.sto', NUCLEIC),
        allow: CASE_WARNINGS.case4,
      },
      { label: 'values and categories', doc: valuesGallery() },
    ],
  },
  {
    key: 'color-schemes.stories.tsx#Picker',
    cases: () => [
      {
        label: 'SH2 conservation and members',
        doc: pickerDoc(),
        registry: HOST,
        allow: CASE_WARNINGS.case3,
        steps: [
          ...[
            ['members', 'taylor'],
            ['members', 'story:pastel'],
            ['cons.color', 'magma'],
            ['members', null],
            ['cons.color', null],
          ].map(
            ([layer, scheme]): SmokeStep => ({
              label: `seqv:colorScheme ${layer} ${scheme}`,
              run: async (ctx) =>
                (await ctx.commands.run({ command: 'seqv:colorScheme', params: { layer, scheme } }, 'host:links')) ||
                false,
            }),
          ),
        ],
      },
    ],
  },
  { key: 'color-schemes.stories.tsx#ManyTracks', cases: () => [{ label: '25 tracks', doc: manyDoc() }] },
  {
    key: 'color-schemes.stories.tsx#CustomPalette',
    cases: () => [{ label: 'custom', doc: customDoc(), allow: CASE_WARNINGS.case3 }],
  },
  {
    key: 'color-schemes.stories.tsx#SharedScales',
    cases: () => (['own', 'group', 'explicit'] as const).map((mode) => ({ label: mode, doc: sharedDoc(mode) })),
  },
  { key: 'color-schemes.stories.tsx#ColorRange', cases: () => [{ label: 'range', doc: rangeDoc() }] },
  { key: 'size-themes.stories.tsx#BySize', cases: () => [{ label: 'sized by a field', doc: sizesDoc() }] },
  { key: 'size-themes.stories.tsx#ByScale', cases: () => [{ label: 'linear, sqrt, log', doc: scalesDoc() }] },
  { key: 'size-themes.stories.tsx#Sashimi', cases: () => [{ label: 'sashimi', doc: sashimiDoc() }] },

  // ---- Seq*/Keyed sets ----
  {
    key: 'keyed-sets.stories.tsx#OnAlignment',
    cases: () => [
      {
        label: 'PF00042',
        doc: keyedFamilyDoc(KEYED_FAMILY),
        steps: MODES.map((m) => patchStep(m.label, modePatch(m) as Patch)),
      },
    ],
  },

  // ---- Seq*/Dynamic nodes ----
  {
    key: 'dynamic-nodes.stories.tsx#FirstPart',
    cases: () => {
      dynamicKnobs.latency = 0;
      return [
        {
          label: 'a 2 Mb chromosome',
          streams: true,
          doc: dynamicDoc(),
          loaders: dynamicLoaders(),
          steps: [
            ...Object.entries(PEAKS).map(
              ([k, at]): SmokeStep => ({
                label: `zoom to peak ${k}`,
                run: (ctx) => ctx.viewports.get('dyn.main')?.zoomTo(0, at - 4000, at + 4000),
              }),
            ),
          ],
        },
      ];
    },
  },
  {
    key: 'dynamic-nodes.stories.tsx#SecondPart',
    cases: () => {
      dynamicKnobs.latency = 0;
      const asked: string[] = [];
      return [
        {
          label: 'templates at max_depth 4',
          doc: keyedDoc(4),
          loaders: keyedLoaders(asked),
          steps: [
            {
              // The templates reach the nested sources: their manifests are asked for (the genes', then theirs).
              label: 'nested sources were asked for',
              run: () => asked.length > 1,
            },
            {
              label: 'add a gene, as a patch would',
              run: (ctx) =>
                (ctx.streams as unknown as { nodes: { add(k: string, n: unknown[]): void } }).nodes.add('genes', [
                  items('X1', "Gene X1 (a patch's)", 50_000, 60_000),
                ]),
            },
          ],
        },
        {
          label: 'templates at max_depth 2',
          doc: keyedDoc(2),
          loaders: keyedLoaders([]),
          // Sources deeper than max_depth are reported, never asked: the story says so.
          allow: [{ code: 'svs.stream.depth' }],
        },
      ];
    },
  },

  // ---- Seq*/Streamed data ----
  {
    key: 'streamed-data.stories.tsx#Letters',
    cases: () => {
      streamedKnobs.latency = 0;
      streamedKnobs.failure = 0;
      return [
        {
          label: 'a 20 Mb chromosome',
          streams: true,
          doc: streamedDoc(),
          loaders: streamedLoaders(),
          steps: [1000, 100].map(
            (span): SmokeStep => ({
              label: `zoom to ${span} positions`,
              run: (ctx) => ctx.viewports.get(STREAMED_SECTION)?.zoomTo(0, STREAMED_AT, STREAMED_AT + span),
            }),
          ),
        },
      ];
    },
  },
  {
    key: 'streamed-data.stories.tsx#ThroughAPath',
    cases: () => {
      streamedKnobs.latency = 0;
      streamedKnobs.failure = 0;
      return [
        {
          label: "a protein over three exons of chr2, and the gene's section",
          doc: codingStoryDoc(),
          loaders: codingLoaders(),
          streams: true,
          steps: [
            { label: 'residues 450–550', run: (ctx) => ctx.viewports.get(CODING_PROTEIN)?.zoomTo(0, 450, 550) },
            { label: "residue 500 hovered (its codon's letters asked)", run: (ctx) => hoverResidue(ctx, 500) },
            { label: "residue 500's codon, split by its intron", run: (ctx) => hoverResidue(ctx, 500) === 'T|CA' },
          ],
        },
      ];
    },
  },
  {
    key: 'streamed-data.stories.tsx#OnItsTranscript',
    cases: () => {
      streamedKnobs.latency = 0;
      streamedKnobs.failure = 0;
      return [
        {
          label: "T1's genome and its own residues, P1 mapped onto it",
          doc: transcriptStoryDoc(),
          loaders: transcriptLoaders(),
          streams: true,
          steps: [
            {
              label: 'T1 is a coordinate system, shown by its section',
              run: (ctx) =>
                ctx.graph.assembled.value?.set.byName('genes_src/T1', { kind: -1, id: -1, row: -1 }) === true,
            },
          ],
        },
      ];
    },
  },
  {
    key: 'streamed-data.stories.tsx#LazyAlignment',
    cases: () => [{ label: '20,000 rows', doc: lazyStoryDoc(), loaders: lazyStoryLoaders(), streams: true }],
  },
  {
    key: 'prefetch.stories.tsx#Ahead',
    cases: () => [
      {
        label: 'a chromosome with three snapshots',
        doc: prefetchDoc(),
        loaders: prefetchLoaders(new AskedLog()),
        streams: true,
        steps: [
          {
            label: 'pan right, zoom in',
            run: (ctx) => {
              ctx.viewports.get(PREFETCH_SECTION)?.panBy(300, 0, 'api');
              ctx.viewports.get(PREFETCH_SECTION)?.panBy(300, 0, 'api');
            },
          },
          ...PREFETCH_SNAPSHOTS.map(
            ([title]): SmokeStep => ({ label: `snapshot ${title}`, run: (ctx) => void ctx.snapshots.apply(title) }),
          ),
        ],
      },
    ],
  },
  {
    key: 'ideogram.stories.tsx#Karyotype',
    cases: () => [
      {
        label: '24 chromosomes under the ideogram navigator',
        doc: ideogramDoc(),
        loaders: ideogramLoaders(),
        registry: IDEOGRAM_REGISTRY,
        streams: true,
        steps: IDEOGRAM_PLACES.map(
          ([label, block, a, e]): SmokeStep => ({
            label,
            run: (ctx) => ctx.viewports.get(IDEOGRAM_SECTION)?.zoomTo(block, a, e),
          }),
        ),
      },
    ],
  },
  {
    key: 'jbrowse.stories.tsx#Live',
    skip: 'Live: reads hg38 files over the network. Its document is opened offline by packages/sequences/viewer/_test/jbrowse-story.test.ts.',
  },

  // ---- Seq*/Performance ----
  {
    key: 'performance.stories.tsx#GlyphAtlas',
    cases: () => [{ label: '100,000 x 5,000, lazy', doc: atlasDoc(), loaders: atlasLoaders, streams: true }],
  },

  // ---- Seq*/Search, Seq*/Query, Seq*/Query language ----
  { key: 'search.stories.tsx#AppComposedSearch', same: 'viewer.stories.tsx#Case1' },
  {
    key: 'search.stories.tsx#LocationSets',
    cases: () => [
      {
        label: 'p53 with location sets',
        doc: locationSetsDoc(),
        steps: [
          patchStep('restyle the first entry', {
            format: 'seq_view_spec_patch',
            ops: [
              {
                op: 'set',
                ref: LOCATION_SECTION,
                params: {
                  highlights: [
                    { data: 'cysteines', style: 'band' },
                    { data: 'zinc', style: 'accent' },
                    { data: BOOKMARKS },
                  ],
                },
              },
            ],
          }),
          patchStep('remove an entry', {
            format: 'seq_view_spec_patch',
            ops: [{ op: 'set', ref: LOCATION_SECTION, params: { highlights: [{ data: 'zinc', style: 'accent' }] } }],
          }),
          patchStep('add the domains', {
            format: 'seq_view_spec_patch',
            ops: [
              {
                op: 'set',
                ref: LOCATION_SECTION,
                params: { highlights: [{ data: 'zinc', style: 'accent' }, { data: 'domains' }] },
              },
            ],
          }),
          {
            label: 'step the bookmarks',
            run: (ctx) =>
              ctx.actions.run('bookmark.step', { section: LOCATION_SECTION, dir: 1 }, 'host:story') || false,
          },
        ],
      },
    ],
  },
  { key: 'query.stories.tsx#Gallery', cases: () => [{ label: 'case 7 with a barnase family', doc: familyDoc() }] },
  { key: 'query-language.stories.tsx#ReadBindPrint', skip: 'A text box: no document.' },
  { key: 'query-language.stories.tsx#Skill', skip: 'A page of generated text: no document.' },

  // ---- Seq*/Mol* ----
  { key: 'molstar.stories.tsx#P53', cases: async (env) => [await p53MolstarCase(env)] },

  // ---- Seq*/File formats ----
  ...FILE_FORMATS_SMOKE,
];
