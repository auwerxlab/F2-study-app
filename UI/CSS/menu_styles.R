# F2 title-bar gradient (palette + motion from UI/CSS/navbar_background.R) used as:
#   - the accordion sub-titles' background (Analysis Parameters, Plot Parameters, ...)
#   - the LEFT MENU's selected-tab FONT colour ONLY (no sidebar background at all).
tags$style(HTML("
/* Shared colour-shifting motion (matches the title bar's 15s cadence) */
@keyframes menu-gradient {
  0%   { background-position: 0% 50%; }
  50%  { background-position: 100% 50%; }
  100% { background-position: 0% 50%; }
}

/* ---- Accordion sub-titles: animated gradient background, slightly transparent ---- */
.accordion-button {
  background-color: transparent;
  background-image: linear-gradient(-45deg,
    rgba(102,126,234,0.26), rgba(118,75,162,0.26),
    rgba(240,147,251,0.26), rgba(245,87,108,0.26));
  background-size: 400% 400%;
  animation: menu-gradient 15s ease infinite;
  color: #33294a;
  font-weight: 500;
}
.accordion-button:hover {
  background-image: linear-gradient(-45deg,
    rgba(102,126,234,0.34), rgba(118,75,162,0.34),
    rgba(240,147,251,0.34), rgba(245,87,108,0.34));
}
.accordion-button:not(.collapsed) {
  color: #33294a;
  background-image: linear-gradient(-45deg,
    rgba(102,126,234,0.40), rgba(118,75,162,0.40),
    rgba(240,147,251,0.40), rgba(245,87,108,0.40));
  box-shadow: none;                      /* drop the default divider line */
}
.accordion-button:focus {
  box-shadow: none;
  border-color: rgba(118,75,162,0.45);
}
[data-bs-theme=\"dark\"] .accordion-button,
[data-bs-theme=\"dark\"] .accordion-button:not(.collapsed) {
  color: #f0e8f7;
}

/* ---- Methods: plain grey \">\" that rotates to \"v\" on open (native <details>) ---- */
.method-details > summary {
  list-style: none;                 /* drop the default disclosure triangle (Firefox) */
  cursor: pointer;
  color: #6c757d;
  font-weight: 500;
  padding: 4px 0;
  user-select: none;
}
.method-details > summary::-webkit-details-marker { display: none; }  /* Chrome/Safari */
.method-details > summary::before {
  content: \">\";
  display: inline-block;
  width: 1em;
  color: #6c757d;
  transition: transform 0.15s ease;
}
.method-details[open] > summary::before { transform: rotate(90deg); }  /* > becomes v */
.method-details .method-body { padding: 6px 0 2px 1.1em; }
[data-bs-theme=\"dark\"] .method-details > summary,
[data-bs-theme=\"dark\"] .method-details > summary::before { color: #adb5bd; }

/* ---- LEFT MENU: selected tab = coloured FONT only, no background ---- */
/* Solid F2 purple fallback (also the icon + label). */
.sidebar .btn.btn-outline-secondary.active-tab,
.sidebar .btn.action-button.active-tab,
.sidebar .btn.active-tab .fa,
.sidebar .btn.active-tab span,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary.active-tab,
[data-bs-theme=\"dark\"] .sidebar .btn.active-tab .fa,
[data-bs-theme=\"dark\"] .sidebar .btn.active-tab span {
  background: none !important;
  background-color: transparent !important;
  border-color: transparent !important;
  box-shadow: none !important;
  font-weight: 600 !important;
  color: #764ba2 !important;
}
/* Where supported, paint the font with the animated F2 gradient instead. */
@supports ((-webkit-background-clip: text) or (background-clip: text)) {
  .sidebar .btn.btn-outline-secondary.active-tab,
  .sidebar .btn.action-button.active-tab,
  [data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary.active-tab {
    background-image: linear-gradient(-45deg, #667eea, #764ba2, #f093fb, #f5576c) !important;
    background-size: 400% 400% !important;
    -webkit-background-clip: text !important;
    background-clip: text !important;
    animation: menu-gradient 15s ease infinite;
  }
  .sidebar .btn.active-tab,
  .sidebar .btn.active-tab .fa,
  .sidebar .btn.active-tab span,
  [data-bs-theme=\"dark\"] .sidebar .btn.active-tab .fa,
  [data-bs-theme=\"dark\"] .sidebar .btn.active-tab span {
    -webkit-text-fill-color: transparent !important;
  }
}
"))
