# Drift: a layout allocation regression after stage 1

The allocation tier's "p53 on its gene, whole protein" pan, bisected over stage 1's merges (`frame.quick.ts`):

| Merge | Pan | Zoom |
|---|---|---|
| AX `1e7e87870` | 42.3 KB | 42.6 KB |
| SL `2457f3ab9` | 42.6 KB | 42.9 KB |
| TH `1772436b9` | 72.4 KB | 68.6 KB |

TH added ~30 KB a frame to the layout (TH did not run the tier); its bound is 64 KB. TX's interface adds ~3 KB more.
Being fixed by lane TA (`briefs/ta-theme-allocation.md`); the bound is not moved.
