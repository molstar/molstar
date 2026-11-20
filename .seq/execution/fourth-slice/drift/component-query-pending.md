# Drift: `(elements :component X)` on a source-backed sequence waits for its letters

**What**: lane MX: a query by `component` on a sequence whose letters come from a source stays pending until the
letters are read (the component is the letter's). A bookmark over a long remote sequence shows nothing until then.

**Options**:
1. As built: pending, staged as pending (EC-169), and resolved once the letters arrive; the query asks for them.
2. Read the component from the source's index without the letters (where a format holds it).

**Recommendation**: 1, as built (it follows EC-169); a format with a component index may add 2 later (TODO).
