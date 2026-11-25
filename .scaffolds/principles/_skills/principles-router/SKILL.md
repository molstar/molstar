---
name: principles-router
description: "Load before code or design work in an area the project's manifest covers: it says which module (the package's or the instance's) to read."
---

# Principles router

Pick the modules for the work from the project's manifest ("Read before" column): `.<instance>/principles/README.md`,
the path the instance skill names. If no instance skill is loaded, load it first. The package's own module list is
the catalogue in [README.md](../../README.md#modules-the-catalogue), the starting point a manifest is made from, not
what the router reads. A module of the instance named like a package module extends it: read the package's first.
Then for each:

1. Its `README.md`: the header (scope, status, what it depends on).
2. Its `rules.md`, if it has one: "In brief" first, then the full section before changing code it governs.
3. Its code patterns, if it has them: `bootstrap-patterns/` in the package (take in where it fits); in an instance,
   `stable-patterns/` (follow as written) and `sketch-patterns/` (harden, report drift). A bootstrap pattern the
   manifest's adaptations say is superseded is read from its instance copy. What that means is in
   [_core/system.md](../../_core/system.md#code-patterns).

How work is run is the `principles` skill ([../principles/SKILL.md](../principles/SKILL.md)).
