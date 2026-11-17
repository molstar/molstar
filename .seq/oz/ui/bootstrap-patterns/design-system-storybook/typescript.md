# Design system and Storybook: TypeScript

The concept, invariants and conformance cases are in [README.md](README.md); this file is the TypeScript implementation.

## Layout

```
lib/
  design-tokens/            # tokens.ts: the single definition (compiled to CSS vars)
  ui/                       # shared presentational components: props in, events out, no models
    button/
      button.tsx
      button.module.css
      button.stories.tsx
packages/product-x/ui/      # connected views + their stories; composes @ui, no private styling primitives
prototypes/
  checkout-redesign-reference/   # temporary design-time mocks, a package outside the build (see below)
tooling/
  storybook/                # the one Storybook of the repo: aggregates every package by section
    .storybook/main.ts
    docs/                   # narrative pages that have no other home
```

## One Storybook, sections declared by packages

```ts
// tooling/storybook/.storybook/main.ts: one specifier per package that declares a title; folder paths never show
import { readFileSync } from "node:fs";
import { globSync } from "node:fs";                          // Node 22+
const root = "../../..";
type Manifest = { name: string; storybook?: { title: string } };

const packages = globSync(`${root}/{lib,packages}/**/package.json`, { exclude: (p) => p.includes("node_modules") })
  .map((file) => ({ dir: file.replace(/\/package\.json$/, ""), m: JSON.parse(readFileSync(file, "utf8")) as Manifest }))
  .filter(({ m }) => m.storybook?.title)
  .sort((a, b) => a.m.storybook!.title.localeCompare(b.m.storybook!.title));

const sections = [
  { directory: "../docs", files: "**/*.mdx", titlePrefix: "Docs" },
  ...globSync(`${root}/prototypes/*-reference`).map((dir) => ({ directory: `${dir}/src`, files: "**/*.stories.tsx",
    titlePrefix: `Reference/${dir.split("/").pop()!.replace(/-reference$/, "")}` })),
  // e.g. lib/ui → "Design system/Components", packages/product-x/ui → "App/Screens"; src/ and stories/ never show
  ...packages.map(({ dir, m }) => ({ directory: `${dir}/src`, files: "**/*.stories.tsx", titlePrefix: m.storybook!.title })),
];
export default { framework: "@storybook/react-vite", stories: sections, addons: ["@storybook/addon-docs", "@storybook/addon-a11y"] };
```

A story file's auto title is its path below `src/` without folders a reader doesn't need: name files for what they
show (`button.stories.tsx` → "Design system/Components/Button") and set a short `title` only to regroup.

## Icons: Lucide plus custom, one format (§one-icon-set)

```ts
// lib/icons/src/index.ts: icons as data, framework-free; Lucide's node format for both sources
import { icons as lucide } from "lucide";                       // IconNode = [tag, attrs][] on a 24x24 grid
export type IconNode = readonly (readonly [tag: string, attrs: Record<string, string | number>])[];

const custom = {
  Branch2: [["circle", { cx: 6, cy: 6, r: 2 }], ["circle", { cx: 18, cy: 18, r: 2 }], ["path", { d: "M6 8v4a4 4 0 0 0 4 4h6" }]],
} satisfies Record<string, IconNode>;

export const icons = { ...lucide, ...custom } as const;
export type IconName = keyof typeof icons;
// the format every node is drawn in; the Icon component applies it, custom nodes may not override it
export const iconFrame = { viewBox: "0 0 24 24", fill: "none", stroke: "currentColor", "stroke-width": 2,
  "stroke-linecap": "round", "stroke-linejoin": "round" } as const;
```

```tsx
// lib/ui/icon/icon.tsx: the only renderer; size from tokens, color from the surrounding text
export function Icon({ name, size = "md", label }: { name: IconName; size?: "sm" | "md" | "lg"; label?: string }) {
  return (
    <svg {...iconFrame} className={styles[size]} role={label ? "img" : undefined} aria-label={label} aria-hidden={!label}>
      {icons[name].map(([tag, attrs], i) => createElement(tag, { key: i, ...attrs }))}
    </svg>
  );
}
```

The gate check (`tooling/scripts/icons-check.ts`) walks `custom`: only `path`, `circle`, `rect`, `line`, `polyline`,
`polygon`, `ellipse`; no `fill`, `stroke`, `stroke-width`, `transform` or `style` attributes; coordinates inside 0-24.
A story renders the whole registry as a grid (Lucide and custom side by side) so a custom icon that looks off shows; it
lives beside `Icon` in `lib/ui` (the icons package may not depend on the renderer). The full registry is large (all of
Lucide), so an app bundles only the names its UI draws: it aliases the icons package to a module of named imports, and
a quick test fails when a UI source names an icon missing there.

## Tokens defined once

```ts
// lib/design-tokens/tokens.ts
export const tokens = {
  color: { bg: "#0b0d10", fg: "#e8ebef", accent: "#4f8cff", danger: "#e5484d" },
  space: { 1: 4, 2: 8, 3: 12, 4: 16, 6: 24 },
  radius: { sm: 4, md: 8, xl: 16, pill: 999 },
  size: { icon: 16, control: 32, "h1": 28, "h2": 22, "h3": 18 },   // sizes too, so components carry no raw lengths
  font: { body: "14px/1.45 system-ui", mono: "13px/1.4 ui-monospace" },
  motion: { fast: "120ms", base: "200ms" },
} as const;

// takes any token set, so dark/compact overrides use the same function
export const toCssVars = (t: Record<string, Record<string, string | number>>): string =>
  Object.entries(t).flatMap(([group, vals]) =>
    Object.entries(vals).map(([k, v]) => `--${group}-${k}:${typeof v === "number" ? `${v}px` : v};`)).join("");
```

Emission on `:root` with `[data-theme=dark]` and `[data-density=compact]` overrides, and breakpoints as constants: see [README.md](README.md#shape).

## Styling: CSS Modules

```css
/* lib/ui/button/button.module.css: locally scoped class names; only tokens, no raw values */
.root {
  padding: var(--space-2) var(--space-4);
  border-radius: var(--radius-sm);
  background: var(--color-accent);
  color: var(--color-fg);
  font: var(--font-body);
  transition: background var(--motion-fast);
}
.root:disabled { opacity: 0.6; }
.root[data-variant="danger"] { background: var(--color-danger); }
```

```tsx
// lib/ui/button/button.tsx: presentational; plain props, no models
import styles from "./button.module.css";

export interface ButtonProps { label: string; busy?: boolean; variant?: "primary" | "danger"; onPress: () => void }

export function Button({ label, busy = false, variant = "primary", onPress }: ButtonProps) {
  return (
    <button className={styles.root} data-variant={variant} disabled={busy} aria-busy={busy} onClick={onPress}>
      {label}
    </button>
  );
}

// dynamic values go through CSS variables, not generated styles
<div className={styles.bar} style={{ "--progress": `${p * 100}%` } as CSSProperties} />
```

Invariants: see [README.md](README.md#invariants).

- Lint: reject raw colors and spacing/font literals in `*.module.css` (stylelint), and runtime CSS-in-JS imports in `lib/ui`.

## Two kinds of components: presentational and connected

**Presentational** components (the design system, like `Button` above) are props in, events out. They know nothing about models, and may keep *transient local UI state* (hover, focus, open/closed, a text draft) in the view.

```tsx
// transient UI state stays local: nobody else needs it, and it is never saved, shared, replayed or inspected
export function Disclosure({ summary, children }: { summary: string; children: ReactNode }) {
  const [open, setOpen] = useState(false);
  return <section><button aria-expanded={open} onClick={() => setOpen(!open)}>{summary}</button>{open && children}</section>;
}
```

**Connected** views live in the product's UI package (`packages/<product>/ui`) at feature boundaries: they bind a model to presentational components ([§reactive-models](../reactive-models/README.md)) and pass plain values and callbacks down.

```tsx
// packages/product-x/ui/save-button.tsx: the only place that touches the model
interface DocModel { readonly state: ReadonlyBehavior<{ readonly saving: boolean }>; readonly save: () => void }

export function SaveButton({ model }: { model: DocModel }) {
  const saving = useObservableProp(model.state, (s) => s.saving);     // the slice it uses, not the whole state
  return <Button label="Save" busy={saving} onPress={model.save} />;
}
```

Where state lives: see [reactive-models](../reactive-models/README.md#shape) (§thin-view-layer).

## Stories

Presentational components use plain args; connected views use a model on the `local()` composition ([§ports-and-adapters](../../../architecture/bootstrap-patterns/ports-and-adapters/README.md)).

```tsx
// lib/ui/button/button.stories.tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { expect, fn, userEvent, within } from "storybook/test";

const meta = { title: "ui/Button", component: Button, tags: ["autodocs"], args: { label: "Save", onPress: fn() } } satisfies Meta<typeof Button>;
export default meta;
type Story = StoryObj<typeof meta>;

export const Default: Story = {};
export const Busy: Story = { args: { busy: true } };
export const Danger: Story = { args: { variant: "danger" } };
export const LongLabel: Story = { args: { label: "A very long label ".repeat(6) } };

// story as interaction test (medium tier)
export const Press: Story = {
  play: async ({ canvasElement, args }) => {
    await userEvent.click(within(canvasElement).getByRole("button"));
    await expect(args.onPress).toHaveBeenCalledOnce();
  },
};
```

```tsx
// packages/product-x/ui/save-button.stories.tsx: connected view on the in-memory composition, no backend
export default { title: "product-x/SaveButton" };

export const SaveFlow: StoryObj = {
  loaders: [async () => {
    const scope = new Scope();
    const app = createApp(local(scope), scope);
    const { id } = (await app.invoke("doc.create", { title: "Hello" })) as { id: string };   // seed through a command (§test-tiers)
    return { model: openDoc(app, id) };                                                       // helper assumed: builds the DocModel
  }],
  render: (_args, { loaded }) => <SaveButton model={loaded.model as DocModel} />,
};
```

To cover every node kind, generate one story per `vocabularySamples()` entry ([§conformance](../../../testing/bootstrap-patterns/conformance/README.md)).

## Accessibility and visual checks (gate)

Through the Storybook test integration of the version in use (today the vitest addon in browser mode; earlier the
standalone test runner shown here); the checks are the same: interaction, a11y, reference stories excluded.

```ts
// tooling/storybook/.storybook/test-runner.ts
import type { TestRunnerConfig } from "@storybook/test-runner";
import { injectAxe, checkA11y } from "axe-playwright";

const config: TestRunnerConfig = {
  tags: { exclude: ["reference"] },                      // reference mocks are never tested or baselined (§reference-mocks)
  async preVisit(page) { await injectAxe(page); },
  async postVisit(page) { await checkA11y(page, "#storybook-root", { detailedReport: true }); },
};
export default config;
```

```bash
npm run storybook:build          # must succeed at the gate
npx test-storybook --ci          # interaction + a11y for all real stories
# visual regression baselines (Chromatic, Playwright snapshots, ...): per story, reviewed on change
```

## Documentation hub, including non-UI modules

Docs are rendered from the durable sources and the registries, not copied (§rendered-docs).

```mdx
{/* tooling/storybook/docs/modules/storage.mdx */}
import { Meta, Markdown } from "@storybook/addon-docs/blocks";
import spec from "../../../../lib/storage/_spec/module.md?raw";
import { CommandCatalog, ConformanceStatus } from "../components";

<Meta title="modules/storage" />

<Markdown>{spec}</Markdown>

## Conformance
<ConformanceStatus suite="DocumentStore" />

## Commands
<CommandCatalog prefix="doc." />
```

```tsx
// tooling/storybook/docs/components.tsx: generated from the command registry (§commands), so docs cannot drift
const registry = registryFor(local());
export const CommandCatalog = ({ prefix }: { prefix: string }) =>
  <Table rows={registry.catalog().filter((c) => c.name.startsWith(prefix))} cols={["name", "effect", "description"]} />;
```

The rules render the same way: one `import x from ".../rules/<file>.md?raw"` and `<Markdown>{x}</Markdown>` per rules file. Error code catalogs ([§error-codes](../../../architecture/bootstrap-patterns/error-codes/README.md)) render as a table per module.

## Reference/mock section (temporary, design-time)

```tsx
// prototypes/checkout-redesign-reference/summary.stories.tsx: throwaway; hardcoded data, ad-hoc CSS, no models/ports
export default {
  tags: ["reference"],
  parameters: { reference: { design: ".design/checkout/redesign.md" },   // required: owning design
                chromatic: { disableSnapshot: true } },
};
export const Mock = { render: () => <div style={{ padding: 16, background: "#111", color: "#fff" }}>Total <b>$42.00</b></div> };
```

```ts
// tooling/scripts/reference-gate.ts (merge gate): entries name an existing design; nothing imports them
import { execSync } from "node:child_process";
import { existsSync, readFileSync } from "node:fs";
const sh = (c: string) => execSync(c, { encoding: "utf8" }).trim();
const problems: string[] = [];

for (const f of sh("git ls-files 'prototypes/*-reference/**.stories.tsx'").split("\n").filter(Boolean)) {
  const design = /design:\s*"([^"]+)"/.exec(readFileSync(f, "utf8"))?.[1];
  if (!design || !existsSync(design)) problems.push(`${f}: names no existing .design/ file`);
}
const importers = sh(`git grep -l -E "from ['\\"](@scope/[a-z0-9-]+-reference|[^'\\"]*-reference/)" -- lib packages apps services || true`);
if (importers) problems.push(`imports from a reference package:\n${importers}`);
process.exit(problems.length ? (console.error(problems.join("\n")), 1) : 0);
```

Folding in: build the real component and its stories from the reference, link the design to the real story ids, and delete the reference entry (a checklist item of the work that implements the design, §design-plan-implement).
