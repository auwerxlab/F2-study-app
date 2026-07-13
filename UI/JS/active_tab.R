# Toggle the .active-tab class on the sidebar menu button that matches the
# current tab. Driven from the server via session$sendCustomMessage("setActiveTab",
# paste0("tab_", current_tab())), so it also follows programmatic navigation
# (e.g. the version-footer link jumping to Info).
tags$script(HTML("
Shiny.addCustomMessageHandler('setActiveTab', function(tabId) {
  document.querySelectorAll('.active-tab').forEach(function(el) {
    el.classList.remove('active-tab');
  });
  if (tabId) {
    var el = document.getElementById(tabId);
    if (el) el.classList.add('active-tab');
  }
});
"))
