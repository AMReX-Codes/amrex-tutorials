#!/usr/bin/env bash
#
# Check the layout of the tutorials and reject legacy AMReX patterns.

set -u
cd "$(git rev-parse --show-toplevel)"

TOPICS="Base AMR LinearSolvers EB Particles FFT TimeIntegration Coupling Fortran Python ML Applications BuildSystems"
status=0
fail () { echo "ERROR: $*"; status=1; }

for topic in ${TOPICS}; do
    [ -d "${topic}" ] || continue
    for dir in "${topic}"/*/; do
        dir=${dir%/}
        [ -f "${dir}/README.md" ] || fail "${dir}: missing README.md"
        case "${topic}" in
            Python|BuildSystems|Applications) ;;
            *)
                for f in CMakeLists.txt GNUmakefile Make.package; do
                    [ -f "${dir}/${f}" ] || fail "${dir}: missing ${f}"
                done ;;
        esac
    done
done

# Every GNUmakefile must default AMREX_HOME to an amrex checkout next to this repository.
while IFS= read -r mk; do
    depth=$(dirname "${mk}" | tr -cd '/' | wc -c)
    expected="AMREX_HOME ?= $(printf '../%.0s' $(seq 0 $((depth + 1))))amrex"
    grep -qxF "${expected}" "${mk}" || fail "${mk}: expected '${expected}'"
done < <(git ls-files ${TOPICS} | grep '/GNUmakefile$')

# Legacy patterns (BL_PROFILE is still the profiling macro and is allowed)
pattern='Gpu::LaunchRegion|AMREX_LAUNCH_|AMREX_FOR_[1-4]D|AMREX_PARALLEL_FOR_|AMREX_HOST_DEVICE_FOR_|AMREX_HOST_DEVICE_PARALLEL_FOR_|\bBL_(ASSERT|TO_FORTRAN|FORT|SPACEDIM|USE_|LANG|REAL|ARLIM)|AMREX_ARLIM|\.copy *\('
if git ls-files ${TOPICS} | grep -E '\.(cpp|H|h)$' | xargs -r grep -nE "${pattern}"; then
    fail "legacy AMReX pattern(s) found above"
fi
if git ls-files ${TOPICS} | grep -v '^Fortran/' | grep -E '\.(f|F|f90|F90)$'; then
    fail "Fortran file(s) found above"
fi

exit ${status}
