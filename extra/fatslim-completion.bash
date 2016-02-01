#!/binbash

_fatslim_help_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\naggregates\napl\nbenchmark\nhelp\nmembranes\nself-test\nthickness\nversion' -- $c)); return 0; fi
}

_fatslim_version_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose' -- $c)); return 0; fi
}

_fatslim_selftest_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose' -- $c)); return 0; fi
}

_fatslim_benchmark_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\n-c\n--conf\n--trajectory\n-t\n--index\n-n\n--hg-group\n--nthreads\n--begin\n-b\n--end\n-e' -- $c)); return 0; fi
case "$p" in
--conf|-c) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.gro' -f -- $c ; compgen -S '/' -d $c)) ;;
--trajectory|-t) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*@(.xtc|.trr)' -f -- $c ; compgen -S '/' -d $c)) ;;
--index|-n) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
esac
}

_fatslim_aggregates_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\n-c\n--conf\n--trajectory\n-t\n--index\n-n\n--hg-group\n--nthreads\n--begin\n-b\n--end\n-e\n--cutoff\n--output\n-o\n--output-index\n--output-index-hg' -- $c)); return 0; fi
case "$p" in
--conf|-c) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.gro' -f -- $c ; compgen -S '/' -d $c)) ;;
--trajectory|-t) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*@(.xtc|.trr)' -f -- $c ; compgen -S '/' -d $c)) ;;
--index|-n) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--output|-o) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.xvg' -f -- $c ; compgen -S '/' -d $c)) ;;
--output-index) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--output-index-hg) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
esac
}

_fatslim_membranes_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\n-c\n--conf\n--trajectory\n-t\n--index\n-n\n--hg-group\n--nthreads\n--begin\n-b\n--end\n-e\n--cutoff\n--output\n-o\n--output-index\n--output-index-hg' -- $c)); return 0; fi
case "$p" in
--conf|-c) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.gro' -f -- $c ; compgen -S '/' -d $c)) ;;
--trajectory|-t) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*@(.xtc|.trr)' -f -- $c ; compgen -S '/' -d $c)) ;;
--index|-n) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--output|-o) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.xvg' -f -- $c ; compgen -S '/' -d $c)) ;;
--output-index) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--output-index-hg) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
esac
}

_fatslim_thickness_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\n-c\n--conf\n--trajectory\n-t\n--index\n-n\n--hg-group\n--nthreads\n--begin\n-b\n--end\n-e\n--cutoff\n--idfreq\n--thickness-cutoff\n--plot-thickness\n--export-thickness-raw\n' -- $c)); return 0; fi
case "$p" in
--conf|-c) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.gro' -f -- $c ; compgen -S '/' -d $c)) ;;
--trajectory|-t) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*@(.xtc|.trr)' -f -- $c ; compgen -S '/' -d $c)) ;;
--index|-n) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--plot-thickness) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.xvg' -f -- $c ; compgen -S '/' -d $c)) ;;
--export-thickness-raw) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.csv' -f -- $c ; compgen -S '/' -d $c)) ;;
esac
}

_fatslim_apl_compl() {
local IFS=$'\n'
local c=${COMP_WORDS[COMP_CWORD]}
local n
for ((n=1;n<COMP_CWORD;++n)) ; do [[ "${COMP_WORDS[COMP_CWORD-n]}" == -* ]] && break ; done
local p=${COMP_WORDS[COMP_CWORD-n]}
COMPREPLY=()
if (( $COMP_CWORD <= 1 )) || [[ $c == -* ]]; then COMPREPLY=( $(compgen -S ' '  -W $'-h\n--help\n--debug\n-v\n--verbose\n-c\n--conf\n--trajectory\n-t\n--index\n-n\n--hg-group\n--nthreads\n--begin\n-b\n--end\n-e\n--cutoff\n--idfreq\n--plot-apl\n--export-apl-raw\n--apl-by-type\n--apl-limit\n--apl-cutoff\n--plot-area' -- $c)); return 0; fi
case "$p" in
--conf|-c) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.gro' -f -- $c ; compgen -S '/' -d $c)) ;;
--trajectory|-t) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*@(.xtc|.trr)' -f -- $c ; compgen -S '/' -d $c)) ;;
--index|-n) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.ndx' -f -- $c ; compgen -S '/' -d $c)) ;;
--plot-apl) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.xvg' -f -- $c ; compgen -S '/' -d $c)) ;;
--export-apl-raw) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.csv' -f -- $c ; compgen -S '/' -d $c)) ;;
--plot-area) (( $n <= 1 )) && COMPREPLY=( $(compgen -S ' ' -X '!*.xvg' -f -- $c ; compgen -S '/' -d $c)) ;;
esac
}


_fatslim_compl() {
local i c m
local IFS=$'\n'
COMPREPLY=()
unset COMP_WORDS[0]
for (( i=1; i<COMP_CWORD; ++i )) ; do
[[ "${COMP_WORDS[i]}" != -* ]] && break
unset COMP_WORDS[i]
done
if (( i == COMP_CWORD )); then
c=${COMP_WORDS[COMP_CWORD]}
COMPREPLY=( $(compgen -S ' ' -W $'-h\n--help\n--version\nhelp\nversion\nself-test\naggregates\nmembranes\nthickness\napl\nbenchmark' -- $c) )
return 0
fi
m=${COMP_WORDS[i]}
COMP_WORDS=( "${COMP_WORDS[@]}" )
COMP_CWORD=$((COMP_CWORD-i))
case "$m" in
    help) _fatslim_help_compl ;;
    version) _fatslim_version_compl ;;
    self-test) _fatslim_selftest_compl ;;
    benchmark) _fatslim_benchmark_compl ;;
    aggregates) _fatslim_aggregates_compl ;;
    membranes) _fatslim_membranes_compl ;;
    thickness) _fatslim_thickness_compl ;;
    apl) _fatslim_apl_compl;;
esac
}

complete -o nospace -F _fatslim_compl fatslim
