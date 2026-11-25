# The novice pass's findings, fourth slice (2026-11-21)

A light agent played a first-time user over the Storybook of the main checkout at the stage 4 gate (light, and dark for
Track headers and Legends), one task script a story (`briefs/novice-pass.md`). Services (SV) are not built: skipped. The
Secondary click on the "2 hidden" chip and the Domains marks could not be sent by the browser tool, so the chip's list
is unseen (a gap, not a finding). The pane repaints only after a mouse move, so a screenshot taken straight after an
action can be one step stale; every finding below was checked against the DOM or a second look.

29 findings. Counted by kind: **bug 5** (F5, F7, F13, F15, F23), **missing design 11** (F1, F3, F8, F9, F11, F14, F21,
F22, F25, F28, F29), **look 13** (F2, F4, F6, F10, F12, F16, F17, F18, F19, F20, F24, F26, F27). What worked: hover both
ways (with and without latency, the stale messages did not take the hover), selection both ways, Esc, double-click
focus, color by this track, Host data's drop and "Read the file instead", the host widget following Move down, the
extensibility toggle, dragging the coverage strip's window.

Friction: 1 = none, 2 = slowed, 3 = hunted for a minute, 4 = gave up on the first idea, 5 = blocked. Nothing blocked.

| # | Story | Tried | What happened | Friction | Kind | Placed |
| --- | --- | --- | --- | --- | --- | --- |
| F1 | Structure link › Hover both ways | hover a residue of P04637, then a cell of a chain | both lit on the other side, but nothing names the residue or its number; the only words are the fake's message log below the fold | 2 | missing design | the hover card (features.md: planned): asked again, with the third slice's list |
| F2 | Structure link › Hover with latency | sweep across the fake and back with a 600 ms delay | works; nothing shows that the messages are late (the log has no times) | 1 | look | TODO: time-stamp the fake's log |
| F3 | Structure link › Selection both ways | "Select chain", then "Within 5 Å" | "Select chain" selects all four copies and replaces the Seq* selection (776); "Within 5 Å" changed nothing (still 776) and says of what | 3 | missing design | design question: what does the fake's "Within 5 Å" do on a selection that is already the chain; NP: label the buttons by effect |
| F4 | Structure link › Selection both ways | drag 100-135 on P04637 | the tint crosses the ruler and hides the "100" tick label | 2 | look | TODO (the ring's tint over numbers, third slice 9, is the same) |
| F5 | Structure link › Focus and back | double-click a residue, then Alt+Left | the Seq* view goes back; the fake's camera stays at 146 and the log shows no second message | 3 | bug | fixes.md: Back sends the structure the previous focus (or the story says the camera is not followed) |
| F6 | Structure link › Focus and back | read the zoomed view | the ruler under the zoom is a stray bracket with no number | 2 | look | NP: ruler at high zoom (a tick at the focus) |
| F7 | Structure link › Chain maps | scroll to the fake | its box is painted empty (the rows are in the DOM, 8 rows at the right place, hover events arrive); every other story of the folder paints them | 3 | bug | fixes.md: the fake's rows blank in Chain maps (cause not found: not the layout, no opacity; the second `Seq*` section may paint over it) |
| F8 | Structure link › Chain maps | hover residue ~160 | the light lands at x for ~107 in both sections: the "shifted by one" of the text cannot be seen at this scale, and the toggle's effect cannot be read | 2 | missing design | NP: the story hovers a zoomed window so one residue wide shifts show |
| F9 | Structure link, Track headers › Color by this track | click the box icon on Domains | the fake is colored (one run, "Colored 1 runs"); the icon has no pressed state, nothing says which track colors; the box icon moves with the label's width | 2 | missing design | design question: a marked state for the coloring track (viewer-models) |
| F10 | the two Color by this track stories | open both | the same picture twice | 1 | look | TODO: one of the two goes |
| F11 | Host data › Served by the fake | "Drop chain" | chain D's tracks empty, nothing in the view says why; the issue is below the fold, as `seq.resolve.system.unnamed` | 3 | missing design | design question: an empty section says its issue on the view (the hidden-track chip's place) |
| F12 | Host data | read the B-factor tracks | no value axis, no unit, no legend; the tracks end at 195 of 219 | 2 | look | NP: the story's legend; the end gap is the data (check) |
| F13 | Track headers › Built-ins, Built-ins as DOM parts | point at Domains | the hide slot and menu cut the label to "Domair" (canvas) or overprint it (DOM parts) | 2 | bug | fixes.md: the label column's width under the slots (the host-entry story, with a wider column, is fine) |
| F14 | Track headers › Built-ins | hover and click the "1 filtered" mark | the card says "1 filtered", not what; the click clears the filter at once | 2 | missing design | design question: the card names the filter; a click that clears wants the card to say so |
| F15 | Track headers › Built-ins, Hidden tracks chip | click a mark, then Cmd+Z and Ctrl+Z | nothing comes back, the toolbar's Back is disabled; the text says the key brings it back | 3 | bug | fixes.md: the header's click leaves the focus where the undo key does not reach (or the key is the canvas's only: then say so) |
| F16 | Track headers › Short tracks | read the 14 px Domains | colors are lost (a grey-blue heat with one dark patch, no legend); the marks are 5 px squares; the menu is secondary-click only | 3 | look | TODO: short tracks keep the layer's colors; a mark's hit target |
| F17 | Track headers › Ring through its menu | "Menu of the forward strand" | opens near the top left, does not name its track or point at it | 2 | look | NP: the menu's title is the track's label |
| F18 | Track headers › Host widget from the geometry | read, click badge 2 | works; the text names `viewer.headers`; the menu's "Collapse or expand hidden stretches" is opaque; the menu covers its track | 2 | look | NP: the story's words; TODO: the menu's item |
| F19 | Track headers › Host entry | click the copy slot | the copy happened (log), the icon shows nothing | 1 | look | TODO |
| F20 | Extensibility › Toggled at run time | register the extension | the sites become lollipops and lose the labels (S6, S15…) the fallback had; the Domains legend lists `dna-binding`, not drawn | 2 | look | NP: the lollipop keeps its label; the legend lists what is in view |
| F21 | Extensibility | Alt+double-click on a lollipop | it zooms; nothing shows a site (the text says it "shows the sites") | 2 | missing design | NP: the text says "zooms to a site"; or a readout |
| F22 | Legends › Legend placements | read the band under the canvas | compact swatches carry no names: which pale square is "special"? | 3 | missing design | design question: a compact legend's card on hover (the hover card again) |
| F23 | Legends | read Variants | the band shows "Variants" twice, the second with no swatches; the side panel lists 0.000001, 0.05, 0.1 with no circles; no unit | 3 | bug | fixes.md: the size legend's key (circles) in the compact and the side forms |
| F24 | Legends | compare light and dark | `helix +6` in light, `+5` in dark; the first hydrophobicity swatch is white on a light ground | 2 | look | TODO: the chip count by width; a swatch outline |
| F25 | Navigators › Coverage strip | drag the strip at fit, scroll, Ctrl+wheel | nothing (the window is the whole strip); the strip is the picture of the Conservation track below it; the zoom took the keys panel to find | 3 | missing design | design question: a window at fit is drawn as one; NP: the story starts zoomed |
| F26 | Navigators (all), the keys panel | open it from the keys button | it is clipped at the viewport's bottom (the zoom rows are cut) and stays over the view | 2 | look | TODO: the panel scrolls or opens upward |
| F27 | Navigators › Host drawn | read the strip | the labels NLS and TET overprint | 3 | look | NP: the host story's labels (the model gives the room; the host's DOM does not use it) |
| F28 | Navigators › Host drawn | hover DBD | the readout says the window (`p53:1–393 of p53`), not what is hovered; the "legend" highlight stays after the pointer left | 2 | missing design | NP: the story's readout names the hovered range |
| F29 | Archives › Plasmid saved, Refused | read | static listings: no button saves or opens; Refused has no text above its four cases | 2 | missing design | NP: a save and an open button; one line of text |

## Placed

| Where | Findings |
| --- | --- |
| fixes.md (bugs) | F5, F7, F13, F15, F23 |
| design questions for the author | F1, F22 (the hover card, again), F3, F9, F11, F14, F25 |
| NP (the stories) | F8, F12, F17, F18, F20, F21, F27, F28, F29, F3 |
| TODO | F2, F4, F6, F10, F16, F19, F24, F26 |

The three that would help a first user most: F22 and F1 (the hover card: nothing on the page names what a swatch,
a residue or a mark is), F15 (the undo key that does nothing after a header click), F11 (an empty section that does
not say why).
