setHook(
    packageEvent("languageserver", "onLoad"),
    function(...) {
        options(languageserver.default_linters = lintr::with_defaults(assignment_linter = NULL,commas_linter=NULL,infix_spaces_linter=NULL,spaces_left_parentheses_linter=NULL))
    }
)
