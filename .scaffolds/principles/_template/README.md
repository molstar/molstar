# Templates

Skeletons only, one per kind of module (see [../_core/system.md](../_core/system.md)). To add a module, copy the
skeleton for its kind to the instance, `.<instance>/principles/<area>/`, fill in the header, write the content, and
add a row to the instance's manifest, `.<instance>/principles/README.md` (with its "Read before" line, which the
`principles-router` skill reads). Starter content lives in the modules themselves, not here.

| Skeleton | Kind |
|---|---|
| [rule-set/](rule-set/README.md) | rule set: `rules.md` with sections of rules, and its code patterns in `stable-patterns/` and `sketch-patterns/` (a project's own module, in its instance; the package's modules use `bootstrap-patterns/`; a folder per pattern: language-neutral README plus one file per language) |
| [pattern-collection/](pattern-collection/README.md) | pattern collection: patterns that are the module itself (look and feel), by source folder |
| [layout/](layout/README.md) | layout: how a codebase is laid out |
| [instance-skill/](instance-skill/SKILL.md) | the instance skill: one per project, named after the repository (adopting, "Link the skills") |
| [instance-manifest/](instance-manifest/README.md) | the instance manifest: one per project, `.<instance>/principles/README.md` (adopting, step 1) |
