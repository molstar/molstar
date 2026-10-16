# Accessibility and touch

Required for every view: phones (compact mode), keyboard, screen readers, colorblind-safe and high-contrast colors. Best
practice (WCAG 2.2 AA, the WAI-ARIA patterns), checked per view in tests and Storybook. The code is the truth for how.

- **Contrast from the tokens**: the accent and status colors are text too, at least 4.5:1 on a panel; a high-contrast
  theme (black and white, a yellow accent, borders) and `prefers-contrast` (`styles/theme.css`). No meaning by color
  alone.
- **Targets**: at least 24px; 44px in compact mode (`--msp-touch-target`, `views/shell/shell.module.css`); the canvas's
  own buttons 38px there (`views/shell/canvas.module.css`).
- **Visible focus**, keyboard only (`:focus-visible`).
- **Keyboard**: every view usable without a pointer. The command menu is a modal combobox (Tab stays in, focus returns
  to the opener); trees are ARIA trees (`dom/tree-keys.ts`: Tab reaches a tree once, arrows move, open, close; Enter
  does what a click does); cards and dialogs trap and restore focus (`dom/focus.ts`); a closed panel gives focus to its
  tab; menu buttons follow the menu button pattern; the rail's tabs are a tab list. Keys can be remapped from the
  keyboard in Help (a recorder takes the next chord, Escape cancels; the result is said politely:
  [commands.md](commands.md) §6).
- **Every field has a name**: a params row's control is labelled by its label's id
  (`components/params/params-form.tsx`).
- **The canvas without sight or pointer**: picking, measuring and focusing by target (select, measure, focus take
  targets); the canvas has a text alternative (what is loaded and selected, in words; `viewer/src/ui/shell.tsx`).
- **Screen readers**: one polite live region for results (`ApiResult.did`), the latest message only
  (`managers/notifications.ts`, the canvas chrome's status line); logs are not live; the agent's reply is said once when
  finished, not as it streams; the sequence's keyboard cursor is said.
- **Touch**: a long press gives what hover gives a mouse (the hover label, the sequence's tooltip;
  `viewer/src/sources/canvas-host.ts`); ordered picks take one tap at a time, or targets.
- **Reduced motion**: `prefers-reduced-motion` turns the shell's transitions and the camera's animated moves into cuts
  (`viewer/src/api/camera.ts`, `viewer/src/sources/compose-plugin.ts`).
- **Data colors**: a colorblind-safe palette (Okabe–Ito) beside Prism's and Mol\*'s classic one, chosen in Settings (the
  viewer's `_spec/viewer-spec.md` §2.8).
