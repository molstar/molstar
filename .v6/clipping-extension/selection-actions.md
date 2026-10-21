# Selection actions and the tree filtered by the selection (working notes)

Order: after Prism particles, before the clipping extension, which adds a row to the menu (David, 2026-10-21).

## The bar

`views/canvas/selection-bar.tsx` today is a summary plus Focus, Make a component, More… and Clear. It becomes two parts:

- **The summary is the menu button** ("9 residues ▾"). It opens a menu like the canvas Expand/Fullscreen popover, each row a title and one line of description:
  - **Focus**: bring the selection into view.
  - **Show in structure tree**: the tree filtered by the selection (below).
  - **Make a component**: a copy of the selection as its own component.
  - **More in the Selection panel…**: color, measure, superpose (and clipping).
- **Clear** (×) stays a button.

**Extensions add rows.** `ViewerExtension.selectionActions` is headless, as `nodeMenu` is (clipping plan, D4):
- `{ id, title, description, icon, order, when?(selection), run(ctx) }`;
- the clipping extension adds "Clip…", which opens its section scoped to the selection.

The actions stay commands of the `selection` group (the command menu offers them), and Focus keeps its shortcut.

## The tree filtered by the selection

The scene tree has a text filter (`views/structure/scene-tree.tsx`, `filter`):
- **Hits** are the matching rows and their ancestors.
- **The search reaches into Contents** (entities, chains).

**The filter gets a non-text mode.**
- Today `SceneTree` takes `filter: ReadonlyState<string>`. It becomes `ReadonlyState<SceneFilter>`, where `SceneFilter = { text: string; selection?: boolean }`.
- The search box shows the non-text part as a chip before the text. Backspace in an empty box, or the chip's ×, removes it.
- Text and the chip combine: a row is a hit when it passes both. The mode alone, with empty text, is a filter too.
- The tree's search state, now text only, keeps both, so the chip survives a tab switch as the text does.

A second kind of filter: **"In the selection"**.
- **What counts as a hit:** the rows whose loci overlap the selection. That is the structures, their components and representations that cover it, and in Contents the entities and chains that hold it.
- **Rows unfold** as the search's hits do. There is no scrolling through thousands of rows, since only what holds the selection is there.
- **The filter shows as a chip** in the tree's search field ("In the selection ×"). Typed text narrows it further.
- **It follows the selection while on:** a new selection, a new filter. It goes off with ×, or when the selection is cleared. Off, the tree is as it was.
- **"Show in structure tree"** opens the Structure panel with it on, and picks the first component hit.

Where the hits come from: the scene source computes them, since it owns the plugin, as `scene.selectionHits: ReadonlyState<ReadonlySet<SceneNode> | undefined>`:
- components and representations from `resolveCovered`'s rule (loci remapped onto the component);
- entities and chains from the contents' loci.

It is recomputed on a selection change only while the filter is on. Idle is idle (Prism AGENTS.md).

## Tests and stories

- **Tests:**
  - the hits for a ligand selection on 1cbs: its structure, the ligand component, its representations, Contents → its entity and chain;
  - the filter follows a new selection and goes off with the clear;
  - the bar's menu runs each action;
  - an extension's row shows.
- **Stories:** the bar with its menu open, and the tree filtered by a selection on 1cbs, wide and narrow, light and dark.
