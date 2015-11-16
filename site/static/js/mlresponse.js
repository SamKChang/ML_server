/* inline figure wrapper */
/*
$(document).ready(function() {
  $(".search_icon").each(function(){
    $(this).add($(this).prev())
      .wrapAll('<div class="imageInputWrapper"></div>');
  }); 
});
*/
/* show/hide focus function */
$(document).ready(function () {
  $("#ml_single").show();
});
$(document).ready(function() {
  $(".container").click(function(e) {
    if($(e.target).is('h2')){
      $(this).siblings('div').children("form").hide(200);
      $(this).children("form").show(200);
    }
  });
});

$(document).ready(function() {
  $("#crystal").datepicker();
});
