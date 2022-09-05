#!/usr/bin/awk -f

#
# taken from: https://fortran-lang.org/learn/building_programs/project_make
#

BEGIN {
    # Fortran is case insensitive, disable case sensitivity for matching
    IGNORECASE = 1
}

# Match a module statement
# - the first argument ($1) should be the whole word module
# - the second argument ($2) should be a valid module name
$1 ~ /^module$/ &&
$2 ~ /^[a-zA-Z][a-zA-Z0-9_]*$/ {
    # count module names per file to avoid having modules twice in our list
    if (modc[FILENAME,$2]++ == 0) {
        # add to the module list, the generated module name is expected
        # to be lowercase, the FILENAME is the current source file
        mod[++im] = sprintf("%s.mod = $(%s)", tolower($2), FILENAME)
    }
}

# Match a use statement
# - the first argument ($1) should be the whole word use
# - the second argument ($2) should be a valid module name
$1 ~ /^(!@[a-zA-Z0-9_]+[\s]+)?use$/ &&
$2 ~ /^[a-zA-Z][a-zA-Z0-9_]*,?$/ {
    # Remove a trailing comma from an optional only statement
    gsub(/,/, "", $2)
    # count used module names per file to avoid using modules twice in our list
    if (usec[FILENAME,$2]++ == 0) {
        # add to the used modules, the generated module name is expected
        # to be lowercase, the FILENAME is the current source file
        use[++iu] = sprintf("$(%s) += $(%s.mod)", FILENAME, tolower($2))
    }
}

# Match an include statement
# - the first argument ($1) should be the whole word include
# - the second argument ($2) can be everything, as long as delimited by quotes
$1 ~ /^(#:?)?include$/ &&
$2 ~ /^["'].+["']$/ {
    # Remove quotes from the included file name
    gsub(/'|"/, "", $2)
    # count included files per file to avoid having duplicates in our list
    if (incc[FILENAME,$2]++ == 0) {
        # Add the included file to our list, this might be case-sensitive
        inc[++ii] = sprintf("$(%s) += %s", FILENAME, $2)
    }
}

# Finally, produce the output for make, loop over all modules, use statements
# and include statements, empty lists are ignored in awk
END {
    for (i in mod) print mod[i]
    for (i in use) print use[i]
    for (i in inc) print inc[i]
}
