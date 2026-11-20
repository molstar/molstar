# Drift: a section's `via` counts only its own layers' paths

**What**: lane MX found that a section's `via` (maps.md, structure-syncing.md: the maps a section's projections use)
is counted only against the paths its own layers take; a `via` naming a map the section's hover projection would use
toward another view, but no layer of its own, is `seq.resolve.via.unused`. MX added P04637's domains to the 2OCJ
chains' section so its `via` is used.

**Options**:
1. As built: `via` scopes the section's layers only; a projection toward other views takes the default path.
2. `via` scopes every path from the section: its layers and its projections (hover, selection, focus, follow), and is
   unused only when none of them could take it.

**Recommendation**: 2. structure-syncing.md names `via` for hover and selection projection; a structure section with
no annotation layer is the common case.
