# adopt-principles: UX foundations, structure maps, the map model (2026-11-17)

## Asked

- UX foundations like `@molstar/ui`'s, aimed at sequence data ("when to squircles"), a nicer UI, and good default
  behavior with "consistent logic for selections, toggling, etc."
- Keep a Nightingale reference rendering of a subset, as a check that the renderer is general.
- Maps onto structures: "basically SIFTS, but in abstract form".
- A correction: "computed" is how a map is sourced, not its kind (a computed codon map is a codon map), with runtimes
  bound at viewer composition; and compressed encodings, "not always we need full expanded table".

## Decided

- Maps have a kind, a sourcing and an encoding, answered separately; one interface; structure systems and the
  `structure` map kind ([design/maps.md](../../design/maps.md)).
- UX foundations and one interaction grammar, provisional until the foundation stories and the grammar's tests settle
  the numbers ([design/ux-foundations.md](../../design/ux-foundations.md)).
- The Nightingale reference adapter, a later private package ([design/rendering.md](../../design/rendering.md)).

## Delivered

Uncommitted on `packages-seq`: the designs above; the building blocks and the core and viewer specs updated; the
Storybook pages Building blocks and UX foundations; TODO and status. The three scaffold checks and `pnpm format:check`
pass.

## Learned

The first reading of "computed" folded sourcing into the kind; the author's correction separated three axes, which
also made the structure map an ordinary kind rather than a special case.

## Left open

The foundations' numbers and the grammar's details (provisional); the edge-case outcomes awaiting confirmation; the
root nouns.
