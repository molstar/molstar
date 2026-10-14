/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from "@molstar/seq-core/format";
import { nodePath, SeqCommands, type TrackActionList, trackDataRefs } from "@molstar/seq-viewer";
import { button, type DomPart, type RowProps as DomRowProps, setAttr } from "@molstar/seq-viewer/dom";
import { type RowProps, useLiveVersion, useSeqViewer } from "@molstar/seq-viewer/react";
import { renderSvg } from "@molstar/seq-viewer/svg";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { createContext, Suspense, use, useCallback, useContext, useMemo, useState } from "react";
import { CASE_TITLES, type CompiledCase, compiledCase, described, hideEmpty } from "./support/seqstar-cases";
import { PLASMID_RING, PSYN_LENGTH, rotateRing, zoomDetail } from "./support/seqstar-docs-circular";
import { CODING_PROTEIN } from "./support/seqstar-docs-streamed";
import f from "./support/seqstar-format.module.css";
import {
  GENOTYPE_DENSITIES,
  GENOTYPE_MODES,
  GENOTYPE_SAMPLES,
  GENOTYPE_SITES,
  type GenotypeMode,
  genotypeDensityPatch,
  genotypeModePatch,
} from "./support/seqstar-genotypes";
import { IDEOGRAM_REGISTRY } from "./support/seqstar-ideogram";
import { INDEL_SELECTIONS, INDEL_VIEWS, selectCommand } from "./support/seqstar-indels";
import { storyPatch } from "./support/seqstar-patch";
import { PrismPage, prism } from "./support/seqstar-prism";
import { Case18Viewer } from "./support/seqstar-search-part";
import { ShiftClickChoice, useTextStyle } from "./support/seqstar-text-style";
import { TRIAGE_MODES, type TriageMode, triagePatch } from "./support/seqstar-triage";
import {
  BuilderCode,
  type CleanPart,
  caseParts,
  STORY_CHROME,
  ViewerCard,
  type ViewerOptions,
} from "./support/seqstar-viewer";
import vs from "./support/seqstar-viewer.module.css";
import { type Args, viewerControls } from "./support/seqstar-viewer-args";
import { Choices, TextAction, Tooltip } from "./ui";

/**
 * The format's acceptance cases in our viewer: each case's document on the test fixtures, shown by the default React UI
 * (or the `<seq-viewer>` element) in Stories' clean chrome (the canvas, tooltip, snapshot player and descriptions with
 * command links; a case with several views gets its tabs, and the cases about a section's switches their section
 * headers: numbering, copies, the ligand's form). "The Stories drafts" has a
 * `chrome` control that adds the headers, track menu, legend and export.
 */
const meta = {
  title: "Seq*/Cases/Viewer",
  ...viewerControls,
  // Said here too: Storybook reads a story file's layout from its own meta, not through a spread.
  parameters: { ...viewerControls.parameters, layout: "fullscreen" },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

const docs = new Map<string, SeqViewSpecDoc>();
function caseDoc(name: string, c: CompiledCase) {
  const d = docs.get(name) ?? described(name, c.doc);
  docs.set(name, d);
  return d;
}

/** Parts `clean` leaves out that a case is about, besides the view tabs a multi-view case needs. */
const CASE_PARTS: Readonly<Record<string, readonly CleanPart[]>> = {
  case2: ["SectionHeader"], // the Kabat and IMGT numbering switch
  case6: ["SectionHeader"], // imatinib's ligand form (a position, or its atoms), and the numbering switch
  case7: ["SectionHeader"], // the numbering and copies switches
  case8: ["SectionHeader"], // the numbering switch
  case9: ["SectionHeader"], // the chains section's copies switch
  case12: ["SectionHeader"], // the dimer's copies (one chain, or its two copies as units) and the heme's ligand form
  case15: ["SectionHeader"], // the designs' own-axes sections: their WT numbering with insertion codes (93A…), and the numbering switch
  case21: ["SectionHeader"], // the sections' names: the gene on its own strand, the mRNA, the region as written
};
const casePartsFor = (name: string, c: CompiledCase) => caseParts(c, CASE_PARTS[name]);

/** Cases whose text view the reader can switch between the `wrapped` and `compact` presets. */
const TEXT_STYLE = new Set(["case1"]);

function CaseViewer({ name, options }: { name: string; options?: ViewerOptions }) {
  const c = use(compiledCase(name));
  const o = { ...options, parts: options?.parts ?? casePartsFor(name, c) };
  const [doc, textStyle] = useTextStyle(caseDoc(name, c), TEXT_STYLE.has(name));
  return (
    <>
      {textStyle}
      <ViewerCard doc={doc} label={c.views.map((v) => v.label).join(" · ")} options={o}>
        {TEXT_STYLE.has(name) ? <ShiftClickChoice /> : null}
      </ViewerCard>
    </>
  );
}

/** What to try in a case, under its kicker. */
const CASE_NOTES: Readonly<Record<string, string>> = {
  case1:
    "Select in the text as in an editor: click a residue (again, it clears), drag across lines (out of the canvas " +
    "too), Shift-click or Shift-drag to extend, click an underline for its region; Enter and Escape from the " +
    "keyboard. The Text choice switches the text's presets: Wrapped spaces it to read (a group gap of about a " +
    "cell, lines apart, the underlines a little under the letters); Compact packs it (half the group gap, lines and " +
    "lanes touching), still in groups of 10. Both are render props (`wrappedGroupGap`, `wrappedGap`, " +
    "`wrappedTrackGap`, `wrappedLaneHeight`) a host sets. The Shift-click choice sets the select behavior's " +
    "`shiftPositions`, as a host's viewer spec would.",
  case8:
    "Open the Aligned error (PAE) group: AlphaFold DB's PAE as a dot plot, the scored residues across and the " +
    "residues aligned on down, over the same window as the tracks (zoom or pan them, the plot follows), with the " +
    "interface's contacts as red dots. Hover a cell: both residues light, with a crosshair; drag to select a " +
    "rectangle, its two ranges.",
};

function Case({ name, options }: { name: string; options?: ViewerOptions }) {
  return (
    <PrismPage
      title={CASE_TITLES[name]}
      caption={
        <>
          @molstar/seq-viewer, its canvas engine and the default UI (
          {options?.ui === "DOM"
            ? "the <seq-viewer> element, @molstar/seq-viewer/dom"
            : "React, @molstar/seq-viewer/react"}
          ); fixtures, no network.
        </>
      }
    >
      {CASE_NOTES[name] ? <div className={f.kicker}>{CASE_NOTES[name]}</div> : null}
      <BuilderCode name={name} />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <CaseViewer name={name} options={options} />
      </Suspense>
    </PrismPage>
  );
}

export const Case1: Story = {
  name: "p53: text and tracks",
  render: (a) => <Case name="case1" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case2: Story = {
  name: "Antibody Fv, Kabat and IMGT",
  render: (a) => <Case name="case2" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case3: Story = {
  name: "SH2 family",
  render: (a) => <Case name="case3" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case4: Story = {
  name: "tRNA family, base pairs",
  render: (a) => <Case name="case4" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case5: Story = {
  name: "DNA construct",
  render: (a) => <Case name="case5" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case6: Story = {
  name: "Imatinib in ABL",
  render: (a) => <Case name="case6" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case7: Story = {
  name: "Barnase–barstar",
  render: (a) => <Case name="case7" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case8: Story = {
  name: "AlphaFold DB complex",
  render: (a) => <Case name="case8" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
export const Case9: Story = {
  name: "PF00042 by view",
  render: (a) => <Case name="case9" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
const coveredRefs = (doc: SeqViewSpecDoc | null, track: string) => {
  const n = nodePath(doc, track)?.at(-1);
  return n ? trackDataRefs(n) : [];
};
/** The track menu (`UI.TrackActions`): the built-ins, then a host's own action. `remove` is not in the default list. */
const CASE10_ACTIONS: TrackActionList = [
  "rows",
  "filter",
  "focus",
  "hide",
  "remove",
  {
    id: "select-covered",
    label: "Select what it covers",
    when: (ctx, _section, track) => coveredRefs(ctx.document.shown.value, track).length > 0,
    run: (ctx, _section, track) => {
      const refs = coveredRefs(ctx.document.shown.value, track);
      ctx.commands.dispatch(SeqCommands.Select, { target: { refs }, origin: "api" }).catch(() => {});
    },
  },
];
/** With `chrome`: `full`, or `clean` as Stories shows it: the canvas, hover, selection, tooltip, navigator. */
export const Case10: StoryObj<Args & { chrome: "full" | "clean" }> = {
  name: "The Stories drafts",
  parameters: {
    docs: {
      description: {
        story:
          'The track menu here is the host\'s own list (`UI.TrackActions`): the built-ins, "Remove track" among them, and "Select what it covers". "Remove track" is built in but not in the default menu; this story enables it.',
      },
    },
  },
  args: { chrome: STORY_CHROME },
  argTypes: {
    chrome: { control: "inline-radio", options: ["full", "clean"], description: "The viewer's `UI.Chrome`" },
  },
  render: (a) => (
    <Case
      name="case10"
      options={{
        navigator: a.navigator,
        rowsColor: a.rowsColor,
        ui: a.ui,
        chrome: a.chrome,
        parts: [], // the toggle shows the difference: clean is clean
        trackActions: CASE10_ACTIONS,
      }}
    />
  ),
};
export const Case12: Story = {
  name: "Designed dimer",
  render: (a) => <Case name="case12" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
const ROWS = [
  { id: "all", label: "All" },
  { id: "variants", label: "With variants" },
];
/** Case 13's viewer, its split layers' rows all, or only the residues a variant has (`split.empty: "hide"`). */
function Case13Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case13"));
  const [rows, setRows] = useState("all");
  const doc = useMemo(
    () => (rows === "all" ? caseDoc("case13", c) : { ...caseDoc("case13", c), root: hideEmpty(c.doc.root) }),
    [c, rows],
  );
  return (
    <ViewerCard
      key={rows}
      doc={doc}
      label={c.views.map((v) => v.label).join(" · ")}
      options={{ ...options, parts: casePartsFor("case13", c) }}
    >
      <div className={prism.controls}>
        <div className={prism.wide}>
          <Choices label="Variant rows (split.empty)" choices={ROWS} active={rows} onPick={setRows} />
        </div>
      </div>
    </ViewerCard>
  );
}
export const Case13: Story = {
  name: "p53 features, UniProt-style",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case13}
      caption={
        <>
          @molstar/seq-viewer, its canvas engine and the default UI; fixtures, no network. The variants grid has a row
          per residue a variant may change to; with variants, <code>split.empty: "hide"</code> keeps only the rows a
          variant has.
        </>
      }
    >
      <BuilderCode name="case13" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case13Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};
export const Case14: Story = {
  name: "Protein over its coding DNA",
  render: (a) => <Case name="case14" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />,
};
/** Case 15's views, and selections that show what each map carries across (best seen in Own axes). */
function Case15Controls() {
  const { context } = useSeqViewer();
  const [view, setView] = useState<string>(INDEL_VIEWS.wt);
  const [selected, setSelected] = useState("none");
  const views = [
    { id: INDEL_VIEWS.wt, label: "On the wild type" },
    { id: INDEL_VIEWS.aligned, label: "Aligned" },
    { id: INDEL_VIEWS.own, label: "Own axes" },
  ];
  const selections = [...INDEL_SELECTIONS.map((s) => ({ id: s.id, label: s.label })), { id: "none", label: "Nothing" }];
  const pickView = (id: string) => {
    setView(id);
    context.views.setCurrent(id);
  };
  const pickSelection = (id: string) => {
    setSelected(id);
    const s = INDEL_SELECTIONS.find((x) => x.id === id);
    if (s) void context.commands.run(selectCommand(s), "host:story");
    else context.selection.clear("host:story");
  };
  return (
    <div className={prism.controls}>
      <div className={prism.wide}>
        <Choices label="View" choices={views} active={view} onPick={pickView} />
      </div>
      <div className={prism.wide}>
        <Choices label="Select" choices={selections} active={selected} onPick={pickSelection} />
      </div>
    </div>
  );
}
function Case15Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case15"));
  return (
    <ViewerCard doc={caseDoc("case15", c)} label="WT barnase and four designs" options={{ ...options, maxHeight: 720 }}>
      <Case15Controls />
    </ViewerCard>
  );
}
/** Case 15: insertions and deletions with the format's own primitives (mappings, an alignment, boundaries). */
export const Case15: Story = {
  name: "Indels in protein designs",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case15}
      caption={
        <>
          Barnase&apos;s mature chain (UniProt P00648, numbered 1–110; its secondary structure, active site and the
          residues touching barstar in 1BRS) and four synthetic designs, each written as edits on the wild type and
          named as HGVS names protein changes: D1 substitutions only; D2 a loop shortened (P64_S67del); D3 two
          N-terminal residues removed and a loop grafted (D93_W94insGNGKDG); D4 an N-terminal Met, His102 and its
          neighbors deleted (catalytically dead) and a C-terminal His tag. Each design is a <code>sequence</code> with a
          gapped <code>pairwise</code> mapping onto the wild type and a <code>wt</code> numbering that names its
          residues in the wild type&apos;s terms, inserts with letters (93A–93F), as Kabat names an antibody&apos;s;
          nothing in the format is specific to indels (frameshifts belong to the DNA, not here). <b>On the wild type</b>
          : a row per design placed through its mapping (<code>letters</code> with <code>of</code>), substitutions on a
          blue cell, deletions as a red bar, insertions as a green mark on their boundary (hover it for the letters);
          below, the same edits as <code>variants</code>, a deletion one marker on a line over its residues, an
          insertion its length where it fits. <b>Aligned</b>: the inserts in columns of their own, lower case on a
          veiled cell and marked on the ruler; the wild type&apos;s numbering skips them. <b>Own axes</b>: each design
          on its own residues, numbered in the wild type&apos;s terms, with the wild type&apos;s features placed through
          its mapping: D4&apos;s H102 lands on the junction its deletion leaves, and the β-sheet spans the graft (a
          range takes its first to last paired residue); each design&apos;s deletions show as the junctions they leave.
          On the wild type, the designs&apos; inserted stretches, data on their own residues, land on the boundaries
          they sit at. Select H102: it lights in every design that keeps it, at its own index there, and in D4, which
          deleted it, as a caret (a bar between two cells) at the junction the deletion leaves. Select D3&apos;s graft:
          it lights in D3 (in Aligned, its insert columns), and as a caret where it sits on the wild type and the other
          designs; hover an insertion&apos;s mark and it lights a caret too. D4&apos;s tag, beyond what the mappings
          pair, lights only in D4.
        </>
      }
    >
      <BuilderCode name="case15" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case15Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};
/** Case 16's actions: rotate the ring's origin, zoom the detail to the cloning site, export the ring as SVG. */
function Case16Controls() {
  const { context } = useSeqViewer();
  const [origin, setOrigin] = useState(0);
  const [svg, setSvg] = useState<string | null>(null);
  const rotate = useCallback(() => {
    const next = (origin + 450) % PSYN_LENGTH;
    setOrigin(next);
    void storyPatch(context, rotateRing(next), "rotate the origin");
  }, [context, origin]);
  const exported = useCallback(() => {
    void renderSvg(context, { width: 640, background: true, overlay: true, sections: [PLASMID_RING] }).then(setSvg);
  }, [context]);
  return (
    <>
      <div className={f.row}>
        <TextAction tone="quiet" onClick={rotate}>
          Rotate the origin (now {origin + 1})
        </TextAction>
        <TextAction tone="quiet" onClick={() => zoomDetail(context)}>
          Zoom the detail to the cloning site
        </TextAction>
        <TextAction tone="quiet" onClick={exported}>
          Export the ring as SVG
        </TextAction>
      </div>
      {svg && (
        // biome-ignore lint/security/noDangerouslySetInnerHtml: renderSvg's string from our own document, not user input
        <div className={vs.figure} dangerouslySetInnerHTML={{ __html: svg }} />
      )}
    </>
  );
}
function Case16Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case16"));
  return (
    <ViewerCard doc={caseDoc("case16", c)} label="pSyn, 2,686 nt" options={options}>
      <Case16Controls />
    </ViewerCard>
  );
}
/** Case 16: a plasmid map through `b.preset.plasmidMap`, a ring with its linear detail and the window between them. */
export const Case16: Story = {
  name: "A plasmid map",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case16}
      caption={
        <>
          <code>b.preset.plasmidMap(view, pSyn, {"{ features, sites, pairs }"})</code> writes a circular section and a
          linear one following it, with generic params only. pSyn is synthetic (2,686 nt, a pUC19-like layout over
          random letters). Each feature has an arrowhead by its <code>strand</code> (<code>direction_field</code>), on
          the ring and in the detail; the f1 region crosses the origin. Labels sit outside the ring, past the
          ruler&apos;s numbers, each with a leader from its feature: sorted by angle and spread a line apart per side,
          the nine cut sites of the cloning site fanned out; hover a label and its feature is hovered, as on the ring
          (the tooltip, the detail), click it to select it. The cut sites are ticks on the features&apos; ring, in the
          skin&apos;s faint ink. Three repeat pairs are chords through the hub. The ring shows the detail&apos;s view as
          a window: a sector with a grip at each end on the ruler&apos;s baseline. Zoom the detail (or wheel in it, or
          focus the ring and press <kbd>+</kbd>), then drag the window round the ring to pan the detail, or a grip to
          zoom it, the other end held. Rotate the origin: the ring starts elsewhere, the labels follow, and a window
          across the new origin wraps. Export the ring as SVG: the figure below the viewer, with its window.
        </>
      }
    >
      <BuilderCode name="case16" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case16Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};
/** Case 17's choices: the matrix layer's order and filter over its sample table, and its density. */
/** The orders: each patches the matrix layer's `order` and `filter`, as the Stories agent would. */
function Case17Modes() {
  const { context } = useSeqViewer();
  const [mode, setMode] = useState("data");
  const pick = (key: string) => {
    const m = GENOTYPE_MODES.find((x) => x.key === key) as GenotypeMode;
    setMode(m.key);
    void storyPatch(context, genotypeModePatch(m) as never, m.label);
  };
  return (
    <div className={prism.wide}>
      <Choices
        label="Order and filter"
        choices={GENOTYPE_MODES.map((m) => ({ id: m.key, label: m.label }))}
        active={mode}
        onPick={pick}
      />
    </div>
  );
}

/** How a device row of many samples picks its color: patches the layer's `density`. */
function Case17Densities() {
  const { context } = useSeqViewer();
  const [density, setDensity] = useState("majority");
  const pick = (key: string) => {
    const d = GENOTYPE_DENSITIES.find((x) => x.key === key) as (typeof GENOTYPE_DENSITIES)[number];
    setDensity(d.key);
    void storyPatch(context, genotypeDensityPatch(d.key) as never, d.label);
  };
  return (
    <div className={prism.wide}>
      <Choices
        label="A pixel row's color (density)"
        choices={GENOTYPE_DENSITIES.map((d) => ({ id: d.key, label: d.label }))}
        active={density}
        onPick={pick}
      />
    </div>
  );
}

function Case17Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case17"));
  return (
    <ViewerCard
      doc={caseDoc("case17", c)}
      label={`${GENOTYPE_SAMPLES.toLocaleString("en-US")} samples × ${GENOTYPE_SITES} sites`}
      options={options}
    >
      <div className={prism.controls}>
        <Case17Modes />
        <Case17Densities />
      </div>
    </ViewerCard>
  );
}
/** Case 17: dense matrix states (KS2), a matrix's keyed states ordered and filtered by a key table. */
export const Case17: Story = {
  name: "A genotype matrix",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case17}
      caption={
        <>
          {GENOTYPE_SAMPLES.toLocaleString("en-US")} synthetic samples as the states of a <code>matrix</code> over{" "}
          {GENOTYPE_SITES} variant sites, heterozygous amber and homozygous for the alternative allele red (a{" "}
          <code>field_palette</code>; the reference genotype is not colored). The alternative allele&apos;s frequency
          varies by population and in blocks of 20 sites. The samples are thinner than a pixel, so each device row shows
          one color for its samples, by the layer&apos;s <code>density</code>: the most common color of its colored
          samples, faded by the share of its samples with the reference genotype (<code>majority</code>), or the
          strongest genotype among them, so a rare homozygote shows (<code>maximum</code>). Alt and the wheel zoom the
          samples vertically. A sample table keyed by <code>state</code> gives each its <code>population</code> and its{" "}
          <code>alt_alleles</code>; the choices patch the layer&apos;s <code>order</code>, <code>filter</code> and{" "}
          <code>density</code>. By population, the blocks line up; the last order keeps about 24 samples, which draw as
          cells (a band each; labels only from 9 px a band).
        </>
      }
    >
      <BuilderCode name="case17" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case17Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};
/** Case 18: finding things in p53: SeqQL queries run at load, the Search part, location sets and bookmarks. */
export const Case18: Story = {
  name: "Finding things in p53",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case18}
      caption={
        <>
          p53&apos;s Tracks view: the letters and three feature tracks (domains, zinc ligands, DNA-contact residues),
          each item with a name to search by. Two SeqQL queries run at load as <code>computed</code> nodes: the
          cysteines, and the domains holding four or more of them (a query naming the other). The section&apos;s{" "}
          <code>highlights</code> show them, as any annotation: the rich domain as a band behind the tracks, the
          cysteines lit, the zinc ligands accented, and the bookmarks (an inline <code>ranges</code> table the bookmark
          actions edit and step). The Search part is the viewer&apos;s; the search is the story&apos;s own: item names
          from two characters, anywhere in the name, then the letters from one, or a regular-expression motif such as{" "}
          <code>{"[KR]{2}.{10,12}[KR]{3,}"}</code> (the bipartite nuclear localization signal). Every result lights and
          ticks in the navigator; Enter, Shift+Enter and ‹ › step through them; Keep holds them as a location set the
          legend lists under Highlights. Place the part in the viewer&apos;s bar (the UI spec gives <code>search</code>)
          or above the viewer, as an app would (<code>{"<Search context={ctx} search={fn} />"}</code>).
        </>
      }
    >
      <BuilderCode name="case18" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case18Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};

/** What the stream manager has done: regions asked, in flight, held. */
function StreamStats() {
  const { context: ctx } = useSeqViewer();
  useLiveVersion(ctx.streams.stats);
  const st = ctx.streams.stats.current;
  return (
    <div className={f.code}>
      {st.requests} regions asked, {st.inFlight} in flight, {st.errors} failed; {(st.bytes / 1024).toFixed(0)} KiB held
    </div>
  );
}
/** Case 19's viewer: the ideogram navigator, a registered entry (the cytobands are the story's data). */
function Case19Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case19"));
  return (
    <ViewerCard
      doc={caseDoc("case19", c)}
      label="chr17 of GRCh38, 83 Mb, streamed by region"
      options={{ ...options, navigator: "ideogram", registry: IDEOGRAM_REGISTRY }}
    >
      <StreamStats />
    </ViewerCard>
  );
}
/** Case 19: a chromosome by region, streamed, with the ideogram as its navigator. */
export const Case19: Story = {
  name: "A chromosome by region",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case19}
      caption={
        <>
          chr17 of GRCh38 as one block of a streamed genome: a letters stream&apos;s manifest names the chromosomes and
          their lengths, and nothing more is asked until the view needs it. The navigator is the ideogram (a registered
          entry, from hg38&apos;s committed cytobands): the chromosome&apos;s bands, its centromere in red, the
          view&apos;s window over them; press to jump, drag to pan, drag an end to zoom. The view opens on 17p13.1; the
          snapshots go to the whole chromosome, 200 kb of 17p13.1 and 50 positions in a gene&apos;s exon. Zoomed out,
          genes and variants are density from each region&apos;s counts and the Sequence track is hatched; at 200 kb the
          genes arrive as items (exons joined, colored by category) and the variants by their alternative; at 50
          positions the letters, colored by nucleotide. phyloP is binned at every zoom. Every track is streamed by
          region from the fixtures (synthetic letters, genes, variants and scores; no network); the line under the
          viewer counts the regions asked.
        </>
      }
    >
      <BuilderCode name="case19" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case19Viewer options={{ rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};

/** Case 20's zooms on P1's axis. */
function Case20Zooms() {
  const { context: ctx } = useSeqViewer();
  const go = (a: number, e: number) => ctx.viewports.get(CODING_PROTEIN)?.zoomTo(0, a, e);
  return (
    <div className={f.row}>
      <TextAction tone="quiet" onClick={() => go(0, 1334)}>
        Whole protein
      </TextAction>
      <TextAction tone="quiet" onClick={() => go(450, 550)}>
        Residues 450–550 (a codon split at 500)
      </TextAction>
      <TextAction tone="quiet" onClick={() => go(880, 920)}>
        Residues 880–920
      </TextAction>
    </div>
  );
}
function Case20Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case20"));
  return (
    <ViewerCard
      doc={caseDoc("case20", c)}
      label="P1, 1,334 residues coded by three exons over 220 kb of chr2"
      options={{ ...options, maxHeight: 720 }}
    >
      <Case20Zooms />
      <StreamStats />
    </ViewerCard>
  );
}
/** Case 20: a protein on its genome, streamed through a codon map, and the transcript's item coordinate system. */
export const Case20: Story = {
  name: "A protein on its genome",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case20}
      caption={
        <>
          P1 (synthetic, 1,334 residues) is coded by three exons over 220 kb of chr2 of a streamed genome, joined to it
          by a <code>codon</code> mapping. Its section shows its domains and chr2&apos;s streamed genes, variants and
          phyloP on its own axis, through that mapping: only the regions its exons reach are asked, a variant lands on
          the residue its position codes, and phyloP over a residue is read over its codon. The second section is chr2
          at the gene: P1&apos;s residues over their codons (<code>letters</code> with <code>of</code>) and its domains
          along the exons, a domain across an intron in two parts. The genes stream also serves T1, the transcript that
          codes P1, with its map (a codon map&apos;s params): when its region arrives T1 becomes an item coordinate
          system, <code>genes_src/T1</code>, and the third section shows its own residues, the genome&apos;s genes,
          variants and phyloP placed through its map. Hover a residue of P1 or T1: its label shows the codon that codes
          it, read from the genome&apos;s streamed letters (asked when first hovered), with its positions on chr2; a
          codon an intron splits shows both parts (residue 500&apos;s T|CA, 900&apos;s A|TG). Hover a domain or a
          variant: it lights in every section that places it.
        </>
      }
    >
      <BuilderCode name="case20" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case20Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};

/** Case 21: a gene on the reverse strand, read in its own direction, with its splice junctions and its mRNA. */
export const Case21: Story = {
  name: "A gene on its reverse strand",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case21}
      caption={
        <>
          GENE2 is synthetic: four exons on the reverse strand of a 1,400-nucleotide region, its protein translated by a
          minus-strand <code>codon</code> mapping (the CDS&apos;s ranges read downward, complemented). The first section
          is the gene as a piece of the region in its own direction (<code>range</code> and <code>direction: -1</code>
          on its block): read 5′ → 3′, the ruler counting down, its letters complemented by a letter map (the data stays
          the + strand&apos;s). Its transcript is one item (four rows of one <code>id</code>, joined), arrowheads by its{" "}
          <code>strand</code>; the protein&apos;s residues over their codons; an RNA-seq-like coverage; and its splice
          junctions as a sashimi plot (<code>b.preset.sashimi</code>): arcs stroked and labeled by read count, the
          gene&apos;s (−) hanging below their line, an exon-skipping one thin, and two weak junctions of an antisense
          transcript (+) rising above theirs. The second section is the mRNA as its own track: the same piece showing
          only the exons (<code>b.preset.exonsOnly</code>, the introns breaks), its letters as RNA (the{" "}
          <code>complement_rna</code> letter map). The third is the region as written, the + strand: the gene&apos;s
          arrowheads point left. The sections follow the first; hover or select anywhere and the same nucleotides light
          in all three.
        </>
      }
    >
      <BuilderCode name="case21" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <CaseViewer name="case21" options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};

/** Case 22's kept rows: the host's own triage state, which the row control toggles. */
const KeptRows = createContext<{ kept: ReadonlySet<number>; toggle: (row: number) => void }>({
  kept: new Set(),
  toggle: () => {},
});
/** The host's control per row (the `RowControl` slot): keep the row, as a curator marks the members to take on. */
function KeepToggle({ dataRow }: RowProps) {
  const { kept, toggle } = useContext(KeptRows);
  const on = kept.has(dataRow);
  const label = on ? `Row ${dataRow + 1} kept: click to drop it` : `Keep row ${dataRow + 1}`;
  return (
    <Tooltip content={label}>
      <button
        type="button"
        className={vs.rowToggle}
        aria-pressed={on}
        aria-label={label}
        onClick={() => toggle(dataRow)}
      >
        Keep
      </button>
    </Tooltip>
  );
}
const KEEP_UI = { components: { RowControl: KeepToggle } };
/** The same control for the `<seq-viewer>` variant, in its shadow root (the skin's `seqv-choice`). */
const domKeepToggle: DomPart<DomRowProps> = (host, _ctx, { dataRow }) => {
  const label = `Keep row ${dataRow + 1}`;
  const b = button(
    { class: "seqv-choice", "aria-pressed": false, "aria-label": label },
    () => setAttr(b, "aria-pressed", b.getAttribute("aria-pressed") !== "true"),
    "Keep",
  );
  host.append(b);
  return { dispose: () => b.remove() };
};
const KEEP_DOM = { RowControl: domKeepToggle };

/** The triage choices: each patches the rows layer's `order` and `filter`, as the Stories agent would. */
function Case22Triage({ kept }: { kept: number }) {
  const { context } = useSeqViewer();
  const [mode, setMode] = useState("tree");
  const pick = (key: string) => {
    const m = TRIAGE_MODES.find((x) => x.key === key) as TriageMode;
    setMode(m.key);
    void storyPatch(context, triagePatch(m) as never, m.label);
  };
  return (
    <div className={prism.controls}>
      <div className={prism.wide}>
        <Choices
          label="Triage (order and filter)"
          choices={TRIAGE_MODES.map((m) => ({ id: m.key, label: m.label }))}
          active={mode}
          onPick={pick}
        />
      </div>
      <div className={f.code}>{kept} rows kept</div>
    </div>
  );
}
function Case22Viewer({ options }: { options?: ViewerOptions }) {
  const c = use(compiledCase("case22"));
  const [kept, setKept] = useState<ReadonlySet<number>>(new Set());
  const store = useMemo(
    () => ({
      kept,
      toggle: (row: number) =>
        setKept((k) => {
          const next = new Set(k);
          if (!next.delete(row)) next.add(row);
          return next;
        }),
    }),
    [kept],
  );
  return (
    <KeptRows.Provider value={store}>
      <ViewerCard
        doc={caseDoc("case22", c)}
        label="PF00042, 32 globins"
        ui={KEEP_UI}
        domSlots={KEEP_DOM}
        options={{ ...options, verticalNav: "both", maxHeight: 480, render: { controlWidth: 44, headerWidth: 380 } }}
      >
        <Case22Triage kept={kept.size} />
      </ViewerCard>
    </KeptRows.Provider>
  );
}
/** Case 22: a family triaged: a guide tree, key tables on its rows, the vertical overview and a control per row. */
export const Case22: Story = {
  name: "A family triaged",
  render: (a) => (
    <PrismPage
      title={CASE_TITLES.case22}
      caption={
        <>
          Pfam&apos;s PF00042 seed (32 globins, committed as a fixture) in an <code>msa</code> section of 16 rows. A
          guide tree is computed at load (<code>tree {'{ provider: "nj" }'}</code>, neighbor joining on p-distances) and
          the rows follow it; click a branch to fold a clade. Three key tables give the rows fields, shown as data
          columns: the review status (the document&apos;s own <code>rows</code> table: three Swiss-Prot entries,
          HBA_HUMAN, MYG_CHICK and YNX4_YEAST, the rest TrEMBL&apos;s, read from the entry names),{" "}
          <code>row_stats</code>&apos; identity to HBA_HUMAN and coverage, and <code>redundancy</code>&apos;s clusters
          at 50%. The triage choices patch the rows layer&apos;s <code>order</code> and <code>filter</code>: the
          reviewed first, the reviewed only, the rows covering 80% of the columns, one per cluster; the tree is drawn
          only in tree order. Beside the scrollbar, the vertical overview: every row in miniature, the rows in view a
          window to drag (Alt-wheel zooms rows). Each row has the host&apos;s control (the <code>RowControl</code>{" "}
          slot): Keep, the curator&apos;s mark, held by the story and counted under the viewer.
        </>
      }
    >
      <BuilderCode name="case22" />
      <Suspense fallback={<div className={f.status}>Compiling…</div>}>
        <Case22Viewer options={{ navigator: a.navigator, rowsColor: a.rowsColor, ui: a.ui }} />
      </Suspense>
    </PrismPage>
  ),
};
