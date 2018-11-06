$(document).ready(function () {
    $('#reference-lib-block').hide();
});

$(document).ready(function () {
    $(document).on('change', '#use_reference_lib', function () {
        if ($(this).is(':checked')) {
            $('#reference-lib-block').show();
        }
        else {
            $('#reference-lib-block').hide();
        }
    });
});

$(document).ready(function () {
    $('#results-table').DataTable();
});
