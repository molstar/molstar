# Drift: by region (lane I)

1. **Node refs**: protocol.md has nodes referring to each other by `(node :key K)` against the answering source; the
   language's `ref` takes absolute addresses only, so runtimes write absolute ones. Options: keep (recommended), or a
   relative `ref` form (a `packages/core` repr change).
2. **`for_each.node_kind`**, not `kind` (a reserved param key in the description language). Record in the design.
3. **An item system's length** is `floor(coding / 3)` (a stop counts as a residue). Options: subtract a stop when the
   letters show one, or an explicit length param. Recommend leaving it until letters arrive.
4. **Summaries per attached source**, not per data entry (data-store.md): a coarse tile per entry would need unbounded
   level-1 tiles. Recommend per source, bounded by `maxLength`; update the design.
5. **Paths between spec systems** always take the static model's path, even where a dynamic map would be shorter, so
   the static tiers never move. Recommend keeping it.
6. **Template re-expansion**: every template re-expands on each overlay version (the design: only the arrival's own).
   Cheap today; a followup.
7. **Instance labels**: the template's label plus its top ancestor's key; nodes have no label param. Option: add one
   to the description.

## From the review (strong), decided with defaults in the fix round

- Nodes are also held while shown (a placed instance, a node whose item system a shown layer uses, a node a held node
  names): `nodes.holdShown`, default true (dynamic-nodes.md's lifetime rule).
- Nested sources are wanted by region for the shown instances, not only the template's root source.
- `maxTemplateTracks` (default 200) caps expansion, with an issue.
- Summaries over `maxLength`: binned on the fly from a bounded level-1 read (or not served); the pyramid counted in the
  store's budget.
- Followups: per-template re-expansion (EC-178); the overlay keyed on what it reads of the model (a view-only edit
  rebuilds it today).
