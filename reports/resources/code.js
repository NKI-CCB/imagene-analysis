$(document).ready(function(){
    $('.datatable>table').DataTable({
        columnDefs: [
            {   targets: 0,
                render: ""
            },
            {   targets: '_all',
                render: $.fn.dataTable.render.number( ',', '.', 3)
            }
        ]
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
