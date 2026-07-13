tags$script(HTML("
  $(document).ready(function() {
    // Wait for bs_themer to load, then collapse it
    setTimeout(function() {
      $('#bsthemerAccordion').removeClass('show');
    }, 500);
  });
"))