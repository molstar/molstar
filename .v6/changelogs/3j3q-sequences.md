## Seq\*: many copies of one chain

- A merged block of many copies no longer searches a path from every copy to every other: a unit reaches only itself
  among its copies unless a mapping ends on a unit. Prism's sequence over 3J3Q (1,356 copies): a 5 s freeze on load →
  none (the document's resolve, headless: 6.5 s → 0.2 s).
- The `secondary` layer draws an opaque feature equal to one already drawn at its start once: a merged block's copies
  alike cost one shape (3J3Q's 18,965 rows, unmerged: 20,345 → 50 shapes a frame).
- `molstar:secondary_structure` takes `copies: "one"`: one row set per molecule merged from its chains (per residue the
  class most chains have). Prism asks for it while copies are merged, so 3J3Q's track shows 14 elements, not every
  chain's stacked; showing the chains apart gives each its own again.
