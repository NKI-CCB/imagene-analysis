$(document).ready(function(){
  $('.datatable>table').DataTable({
    order: []
  });
  $('pre code.sourceCode').each(function(i, block) {
    hljs.highlightBlock(block);
  });

  $('section.collapsed h2 ~ ').hide();
  $('section h2').click(function(e){
    $(this).parent().toggleClass('collapsed');
    $(this).siblings().toggle();
  });
});
