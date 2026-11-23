---
name: <instance>
description: "Load first in any session in this repository: what <repo> is, where its instance of the principles keeps things, and its own modules. Then load principles and principles-router."
---

<!-- The instance skill: one per project, named after the repository, created by the bootstrap (adopting pattern,
"Link the skills"). It lives in the instance (`.<instance>/_skills/<instance>/SKILL.md`), is linked from
`.agents/skills/<instance>`, and links and copies nothing beyond pointers. Replace the <...> placeholders and remove
this comment. -->

# <instance>

<One or two sentences: what this repository is and what it builds.>

It is an instance of the principles package: the package is `<path to the package>` (copied from <home> at <commit or
date>); this repository's own modules, designs, sources, records and work folders are in `.<instance>/`.

1. **Read the manifest**, `.<instance>/principles/README.md`: the folder map (where the package's `.design/`,
   `.records/`, `.execution/` ... are here), the modules (the package's and this project's own) and the adaptations.
2. **Load `principles`** (how work is run), then **`principles-router`**, which picks modules from that manifest.
3. <Anything a session here must know before it starts: the one command that runs the checks, hot files, who the
   author is.>
