# Ensure sidebar menu button text is white in dark mode and readable on hover
tags$style(HTML("
/* Base state */
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary {
  color: #ffffff !important;
  border-color: rgba(255, 255, 255, 0.35) !important;
  background-color: transparent !important;
}
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary .fa,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary span {
  color: #ffffff !important;
} /* Hover, focus, active */
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:hover,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:focus,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:active,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary.active,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:focus-visible {
  color: #ffffff !important;
  background-color: rgba(255, 255, 255, 0.12) !important;
  border-color: rgba(255, 255, 255, 0.5) !important;
  box-shadow: none !important;
}
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:hover .fa,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:hover span,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:focus .fa,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:focus span,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:active .fa,
[data-bs-theme=\"dark\"] .sidebar .btn.btn-outline-secondary:active span {
  color: #ffffff !important;
}
"))