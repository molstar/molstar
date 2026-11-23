# Drift: a right-click on a ring's ribbon

**Where**: the genome ring (`ext-rondo`), CL2. **Found by**: CL2's review, 2026-11-24.

**The design says**: answer B (`drafts/novice-questions.md`): a right-click on a stem opens the stem menu; over the
ring's empty space, the browser's menu; the track's actions (Move up, Reverse, Hide track) stay in its header menu.
It says nothing of a ribbon.

**The code does**: a right-click on a ribbon opens the track's menu (Hide track, Reverse, Collapse). On the genome
ring, which has no track header in view, that is the only pointer route to the track's actions.

**Question for the author**: keep it (a ribbon's right-click is the track's menu), or give a ribbon its own menu (its
pair: show the pair, select its ends) and move the track's actions elsewhere?

**Recommended**: keep it for version 1, and say so in the ring's design and spec. Why: the ring has no header to hold
the track's actions, and a ribbon menu is a new design no finding asked for.
